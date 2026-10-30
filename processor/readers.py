"""
Platform readers. Each reads one sample into a raw AnnData (counts as CSR, every
feature including controls, platform ids as obs_names, coordinates in obs["x"] /
obs["y"]) plus a `meta` dict for the contract (platform, kind, unit, grid, spot
diameter, files read, warnings). processor/contract.py turns it into the contract
shape; nothing here filters or normalises.
"""
from __future__ import annotations

import json
import logging
import math
import re
import struct
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from processor.inputs import Candidate, sanitise

log = logging.getLogger("processor-st-prep")

CHUNK_ROWS = 20_000
CONTROL_PATTERNS = re.compile(
    r"^(NegControlProbe|NegControlCodeword|BLANK_|Blank-|UnassignedCodeword|DeprecatedCodeword|"
    r"NegPrb|SystemControl|FalseCode|Negative\d)", re.IGNORECASE)
POSITION_COLUMNS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


@dataclass
class Read:
    adata: ad.AnnData
    meta: dict
    sub_id: str | None = None        # for files holding several samples
    warnings: list[str] = field(default_factory=list)


def member(values, pool) -> np.ndarray:
    """Boolean mask: which of `values` are in `pool` (hash lookup; np.isin on strings is quadratic)."""
    return np.asarray(pd.Index(values).isin(pd.Index(pool)), dtype=bool)


def _rel(p: Path, root: Path) -> str:
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)


# --- 10x HDF5 ---------------------------------------------------------------------------------

def read_10x_h5(path: Path) -> tuple[sp.csr_matrix, np.ndarray, pd.DataFrame]:
    """
    A 10x Genomics HDF5 feature-barcode matrix (Cell Ranger / Space Ranger / Xenium, format
    version 2+) as (counts CSR obs x features, barcodes, var with gene_ids / feature_types / genome).
    The file stores features x barcodes as CSC; its transpose is the CSR we want, without a copy.
    """
    with h5py.File(path, "r") as h:
        if "matrix" not in h:
            raise ValueError(f"{path.name}: not a 10x HDF5 matrix (no 'matrix' group)")
        m = h["matrix"]
        n_features, n_barcodes = (int(v) for v in m["shape"][:])
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        barcodes = m["barcodes"][:].astype(str)
        ft = m["features"]
        names = ft["name"][:].astype(str)
        ids = ft["id"][:].astype(str) if "id" in ft else names
        types = ft["feature_type"][:].astype(str) if "feature_type" in ft else np.full(len(names), "Gene Expression")
        var = pd.DataFrame({"gene_ids": ids, "feature_types": types}, index=pd.Index(names, name=None))
        if "genome" in ft:
            var["genome"] = ft["genome"][:].astype(str)
    csc = sp.csc_matrix((data, indices, indptr), shape=(n_features, n_barcodes))
    return csc.T.tocsr(copy=False), barcodes, var


def read_10x_mtx(folder: Path) -> tuple[sp.csr_matrix, np.ndarray, pd.DataFrame]:
    """A 10x MTX folder (matrix.mtx[.gz], features.tsv[.gz] or genes.tsv, barcodes.tsv[.gz]) as read_10x_h5 returns."""
    import gzip
    from scipy.io import mmread

    def pick(*names: str) -> Path:
        for n in names:
            if (folder / n).is_file():
                return folder / n
        raise ValueError(f"{folder.name}: no {names[0]}")

    mtx = pick("matrix.mtx.gz", "matrix.mtx")
    opener = gzip.open if mtx.suffix == ".gz" else open
    with opener(mtx, "rb") as fh:
        M = mmread(fh)                       # features x barcodes, COO
    X = sp.csr_matrix(M.T)
    feats = pd.read_csv(pick("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"), sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(pick("barcodes.tsv.gz", "barcodes.tsv"), sep="\t", header=None, dtype=str)[0].to_numpy()
    names = feats[1].to_numpy() if feats.shape[1] > 1 else feats[0].to_numpy()
    var = pd.DataFrame({"gene_ids": feats[0].to_numpy(),
                        "feature_types": feats[2].to_numpy() if feats.shape[1] > 2 else np.full(len(names), "Gene Expression")},
                       index=pd.Index(names))
    if X.shape != (len(barcodes), len(names)):
        raise ValueError(f"{folder.name}: matrix is {X.shape[::-1]} but {len(names)} features x {len(barcodes)} barcodes listed")
    return X, barcodes, var


def read_10x_matrix(path: Path) -> tuple[sp.csr_matrix, np.ndarray, pd.DataFrame]:
    return read_10x_mtx(path) if path.is_dir() else read_10x_h5(path)


# --- Visium / Visium HD -----------------------------------------------------------------------

def read_positions(path: Path) -> pd.DataFrame:
    """tissue_positions(_list).csv (Space Ranger 1.x without header, 2.x with) or tissue_positions.parquet."""
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        with open(path, encoding="utf-8-sig") as fh:
            first = fh.readline()
        if first.lower().startswith("barcode"):
            df = pd.read_csv(path, header=0, encoding="utf-8-sig")
        else:
            df = pd.read_csv(path, header=None, names=POSITION_COLUMNS, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in POSITION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: positions file lacks {', '.join(missing)}")
    df["barcode"] = df["barcode"].astype(str)
    return df.set_index("barcode")


def load_image(path: Path) -> np.ndarray:
    from PIL import Image

    with Image.open(path) as im:
        if im.mode not in ("RGB", "RGBA", "L"):
            im = im.convert("RGB")
        return np.asarray(im, dtype=np.uint8)


def _visium_like(matrix: Path, positions: Path, files: dict[str, Path], root: Path, platform: str,
                 bin_size_um: int | None, library_name: str) -> Read:
    warnings: list[str] = []
    X, barcodes, var = read_10x_matrix(matrix)
    pos = read_positions(positions)
    have = member(barcodes, pos.index)
    if not have.all():
        n = int((~have).sum())
        warnings.append(f"{n} barcodes of the matrix have no position in {positions.name} and were dropped")
        X = X[np.flatnonzero(have)]
        barcodes = barcodes[have]
    pos = pos.loc[barcodes]
    obs = pd.DataFrame(index=pd.Index(barcodes.astype(str)))
    obs["x"] = pos["pxl_col_in_fullres"].to_numpy(dtype=np.float64)
    obs["y"] = pos["pxl_row_in_fullres"].to_numpy(dtype=np.float64)
    obs["in_tissue"] = pos["in_tissue"].to_numpy(dtype=np.int32)
    obs["array_row"] = pos["array_row"].to_numpy(dtype=np.int32)
    obs["array_col"] = pos["array_col"].to_numpy(dtype=np.int32)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    scalefactors: dict = {}
    if "scalefactors" in files:
        with open(files["scalefactors"]) as fh:
            scalefactors = json.load(fh)
    else:
        warnings.append("no scalefactors_json.json: spot diameter estimated from the spot spacing")
    images = {}
    for role in ("hires", "lowres"):
        if role in files:
            try:
                images[role] = load_image(files[role])
            except Exception as exc:  # noqa: BLE001 - an unreadable image is not fatal
                warnings.append(f"tissue_{role}_image.png unreadable: {exc}")
    if not images:
        warnings.append("no tissue_hires_image.png / tissue_lowres_image.png")
    entry: dict = {"scalefactors": scalefactors}
    if images:
        entry["images"] = images
    md: dict = {}
    if matrix.is_file():
        with h5py.File(matrix, "r") as h:
            md = {k: (v.decode() if isinstance(v, bytes) else str(v)) for k, v in h.attrs.items()
                  if k in ("chemistry_description", "software_version")}
    if md:
        entry["metadata"] = md
    adata.uns["spatial"] = {library_name: entry}

    spot = scalefactors.get("spot_diameter_fullres")
    if not spot:
        spot = estimate_diameter(obs.loc[obs["in_tissue"] == 1, ["x", "y"]].to_numpy()
                                 if (obs["in_tissue"] == 1).sum() > 2 else obs[["x", "y"]].to_numpy(),
                                 spacing_fraction=0.55 if platform == "visium" else 1.0)
    meta = {
        "platform": platform, "kind": "spot" if platform == "visium" else "bin", "coordinate_unit": "pixel",
        "spot_diameter": float(spot), "grid": "hex" if platform == "visium" else "square",
        "files": {k: _rel(v, root) for k, v in files.items() if v.exists()},
        "n_positions": int(len(pos)), "n_positions_in_tissue": int((pos["in_tissue"] == 1).sum()),
        "image": {r: (r in images) for r in ("hires", "lowres")}, "scalefactors": scalefactors,
        "software_version": md.get("software_version"),
    }
    if bin_size_um is not None:
        meta["bin_size_um"] = int(bin_size_um)
        if scalefactors.get("microns_per_pixel"):
            meta["microns_per_pixel"] = float(scalefactors["microns_per_pixel"])
    return Read(adata, meta, warnings=warnings)


def read_visium(c: Candidate, root: Path, params: dict) -> list[Read]:
    files = dict(c.files)
    return [_visium_like(files.pop("matrix"), files.pop("positions"), files | {"matrix": c.files["matrix"], "positions": c.files["positions"]},
                         root, "visium", None, c.sample_id)]


def read_visium_hd(c: Candidate, root: Path, params: dict) -> list[Read]:
    want = f"square_{int(params['BIN_SIZE_UM']):03d}um"
    bins = sorted(k for k in c.files if k.startswith("square_"))
    if want not in c.files:
        raise ValueError(f"no {want} under binned_outputs (available: {', '.join(bins)}); set BIN_SIZE_UM")
    d = c.files[want]
    files: dict[str, Path] = {}
    for n in ("filtered_feature_bc_matrix.h5", "filtered_feature_bc_matrix", "raw_feature_bc_matrix.h5", "raw_feature_bc_matrix"):
        if (d / n).is_file() or ((d / n).is_dir() and any((d / n / m).is_file() for m in ("matrix.mtx.gz", "matrix.mtx"))):
            files["matrix"] = d / n
            break
    if "matrix" not in files:
        raise ValueError(f"{want}: no filtered_feature_bc_matrix.h5")
    for n in ("tissue_positions.parquet", "tissue_positions.csv"):
        if (d / "spatial" / n).is_file():
            files["positions"] = d / "spatial" / n
            break
    if "positions" not in files:
        raise ValueError(f"{want}: no spatial/tissue_positions.parquet")
    for role, name in (("scalefactors", "scalefactors_json.json"), ("hires", "tissue_hires_image.png"),
                       ("lowres", "tissue_lowres_image.png")):
        for sd in (d / "spatial", c.source / "spatial"):
            if (sd / name).is_file():
                files[role] = sd / name
                break
    return [_visium_like(files["matrix"], files["positions"], files, root, "visium-hd", int(params["BIN_SIZE_UM"]), c.sample_id)]


# --- imaging platforms ------------------------------------------------------------------------

def estimate_diameter(xy: np.ndarray, spacing_fraction: float = 1.0) -> float:
    """
    A typical object diameter from the coordinates: the median nearest-neighbour distance
    (times `spacing_fraction`: Visium spots are 55 µm on a 100 µm pitch, hence 0.55).
    """
    from scipy.spatial import cKDTree

    xy = np.asarray(xy, dtype=np.float64)
    if len(xy) < 2:
        return 1.0
    if len(xy) > 50_000:
        xy = xy[np.random.default_rng(0).choice(len(xy), 50_000, replace=False)]
    d, _ = cKDTree(xy).query(xy, k=2)
    nn = d[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    return float(np.median(nn) * spacing_fraction) if len(nn) else 1.0


def diameter_from_area(area: np.ndarray) -> float | None:
    a = np.asarray(area, dtype=np.float64)
    a = a[np.isfinite(a) & (a > 0)]
    return float(2.0 * math.sqrt(np.median(a) / math.pi)) if len(a) else None


def read_xenium(c: Candidate, root: Path, params: dict) -> list[Read]:
    warnings: list[str] = []
    X, barcodes, var = read_10x_h5(c.files["matrix"])
    cells_path = c.files["cells"]
    if cells_path.suffix == ".parquet":
        cells = pd.read_parquet(cells_path)
    else:
        cells = pd.read_csv(cells_path)
    if "cell_id" not in cells.columns:
        raise ValueError(f"{cells_path.name}: no cell_id column")
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    have = member(barcodes, cells.index)
    if not have.all():
        warnings.append(f"{int((~have).sum())} cells of the matrix are not in {cells_path.name} and were dropped")
        X = X[np.flatnonzero(have)]
        barcodes = barcodes[have]
    cells = cells.loc[barcodes]
    obs = pd.DataFrame(index=pd.Index(barcodes.astype(str)))
    for col in ("x_centroid", "y_centroid"):
        if col not in cells.columns:
            raise ValueError(f"{cells_path.name}: no {col} column")
    obs["x"] = cells["x_centroid"].to_numpy(dtype=np.float64)
    obs["y"] = cells["y_centroid"].to_numpy(dtype=np.float64)
    obs["in_tissue"] = np.ones(len(obs), dtype=np.int32)
    for col in ("cell_area", "nucleus_area"):
        if col in cells.columns:
            obs[col] = cells[col].to_numpy(dtype=np.float64)
    for col in ("segmentation_method", "fov_name", "region"):
        if col in cells.columns:
            obs[col] = cells[col].astype(str).to_numpy()
    adata = ad.AnnData(X=X, obs=obs, var=var)
    spot = diameter_from_area(obs["cell_area"]) if "cell_area" in obs else None
    meta = {"platform": "xenium", "kind": "cell", "coordinate_unit": "micron",
            "spot_diameter": spot if spot else estimate_diameter(obs[["x", "y"]].to_numpy()), "grid": "none",
            "files": {k: _rel(v, root) for k, v in c.files.items()}}
    return [Read(adata, meta, warnings=warnings)]


def _read_dense_counts_csv(path: Path, id_col: int | str | None = 0, drop_cols: tuple[str, ...] = ()
                           ) -> tuple[sp.csr_matrix, list[str], list[str], pd.DataFrame]:
    """
    A dense cell x gene CSV (MERSCOPE cell_by_gene.csv, CosMx *_exprMat_file.csv) read in
    chunks, each turned into CSR, so no dense copy of the whole table exists.
    Returns (CSR, row ids, gene columns, the dropped id columns as a DataFrame).
    """
    blocks, ids, extras = [], [], []
    genes: list[str] | None = None
    reader = pd.read_csv(path, index_col=id_col, chunksize=CHUNK_ROWS, low_memory=False)
    for chunk in reader:
        if genes is None:
            genes = [g for g in chunk.columns if g not in drop_cols]
        extras.append(chunk[list(drop_cols)].copy() if drop_cols else None)
        ids.extend(chunk.index.astype(str).tolist())
        vals = chunk[genes].to_numpy()
        if vals.dtype.kind == "f" and np.all(np.mod(vals[np.isfinite(vals)], 1) == 0):
            vals = np.nan_to_num(vals).astype(np.int32)
        blocks.append(sp.csr_matrix(vals))
    if genes is None:
        raise ValueError(f"{path.name}: empty table")
    X = sp.vstack(blocks, format="csr") if len(blocks) > 1 else blocks[0]
    extra = pd.concat([e for e in extras if e is not None]) if drop_cols else pd.DataFrame(index=ids)
    return X, ids, genes, extra


def read_merscope(c: Candidate, root: Path, params: dict) -> list[Read]:
    warnings: list[str] = []
    X, ids, genes, _ = _read_dense_counts_csv(c.files["counts"], id_col=0)
    meta_df = pd.read_csv(c.files["metadata"], index_col=0, low_memory=False)
    meta_df.index = meta_df.index.astype(str)
    have = member(np.asarray(ids), meta_df.index)
    if not have.all():
        warnings.append(f"{int((~have).sum())} cells of cell_by_gene.csv are not in cell_metadata.csv and were dropped")
        X = X[np.flatnonzero(have)]
        ids = [i for i, h in zip(ids, have) if h]
    meta_df = meta_df.loc[ids]
    for col in ("center_x", "center_y"):
        if col not in meta_df.columns:
            raise ValueError(f"{c.files['metadata'].name}: no {col} column")
    obs = pd.DataFrame(index=pd.Index(ids))
    obs["x"] = meta_df["center_x"].to_numpy(dtype=np.float64)
    obs["y"] = meta_df["center_y"].to_numpy(dtype=np.float64)
    obs["in_tissue"] = np.ones(len(obs), dtype=np.int32)
    if "volume" in meta_df.columns:
        obs["cell_volume"] = meta_df["volume"].to_numpy(dtype=np.float64)
    if "fov" in meta_df.columns:
        obs["fov"] = meta_df["fov"].astype(str).to_numpy()
    types = np.where([g.lower().startswith("blank") for g in genes], "Blank Codeword", "Gene Expression")
    var = pd.DataFrame({"gene_ids": genes, "feature_types": types}, index=pd.Index(genes))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    spot = None
    if all(col in meta_df.columns for col in ("min_x", "max_x", "min_y", "max_y")):
        ext = ((meta_df["max_x"] - meta_df["min_x"]) + (meta_df["max_y"] - meta_df["min_y"])).to_numpy(dtype=np.float64) / 2
        ext = ext[np.isfinite(ext) & (ext > 0)]
        spot = float(np.median(ext)) if len(ext) else None
    if spot is None and "volume" in meta_df.columns:
        v = meta_df["volume"].to_numpy(dtype=np.float64)
        v = v[np.isfinite(v) & (v > 0)]
        spot = float(2 * (3 * np.median(v) / (4 * math.pi)) ** (1 / 3)) if len(v) else None
    meta = {"platform": "merscope", "kind": "cell", "coordinate_unit": "micron",
            "spot_diameter": spot if spot else estimate_diameter(obs[["x", "y"]].to_numpy()), "grid": "none",
            "files": {k: _rel(v, root) for k, v in c.files.items()}}
    return [Read(adata, meta, warnings=warnings)]


COSMX_META_COLS = ("fov", "cell_ID", "cell", "cell_id", "slide_ID", "slide_ID_numeric")


def _cosmx_ids(df: pd.DataFrame) -> np.ndarray:
    if "cell" in df.columns:
        return df["cell"].astype(str).to_numpy()
    if "cell_id" in df.columns:
        return df["cell_id"].astype(str).to_numpy()
    if "cell_ID" in df.columns and "fov" in df.columns:
        return (df["cell_ID"].astype(str) + "_" + df["fov"].astype(str)).to_numpy()
    if "cell_ID" in df.columns:
        return df["cell_ID"].astype(str).to_numpy()
    raise ValueError("no cell / cell_ID column")


def read_cosmx(c: Candidate, root: Path, params: dict) -> list[Read]:
    warnings: list[str] = []
    px = float(params["COSMX_PIXEL_UM"])
    head = pd.read_csv(c.files["counts"], nrows=0)
    id_cols = tuple(col for col in COSMX_META_COLS if col in head.columns)
    X, _, genes, id_df = _read_dense_counts_csv(c.files["counts"], id_col=None, drop_cols=id_cols)
    ids = _cosmx_ids(id_df)
    meta_df = pd.read_csv(c.files["metadata"], low_memory=False)
    meta_df.index = pd.Index(_cosmx_ids(meta_df))
    have = member(ids, meta_df.index)
    if not have.all():
        warnings.append(f"{int((~have).sum())} cells of the expression matrix are not in the metadata file and were dropped")
        X = X[np.flatnonzero(have)]
        ids = ids[have]
    meta_df = meta_df.loc[ids]
    for col in ("CenterX_global_px", "CenterY_global_px"):
        if col not in meta_df.columns:
            raise ValueError(f"{c.files['metadata'].name}: no {col} column")
    obs = pd.DataFrame(index=pd.Index(ids))
    obs["x"] = meta_df["CenterX_global_px"].to_numpy(dtype=np.float64) * px
    obs["y"] = meta_df["CenterY_global_px"].to_numpy(dtype=np.float64) * px
    obs["in_tissue"] = np.ones(len(obs), dtype=np.int32)
    if "Area" in meta_df.columns:
        obs["cell_area"] = meta_df["Area"].to_numpy(dtype=np.float64) * px * px
    if "fov" in meta_df.columns:
        obs["fov"] = meta_df["fov"].astype(str).to_numpy()
    types = []
    for g in genes:
        gl = g.lower()
        if gl.startswith("negprb") or gl.startswith("negative"):
            types.append("Negative Control Probe")
        elif gl.startswith("systemcontrol"):
            types.append("System Control")
        elif gl.startswith("falsecode"):
            types.append("False Code")
        else:
            types.append("Gene Expression")
    var = pd.DataFrame({"gene_ids": genes, "feature_types": types}, index=pd.Index(genes))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    spot = diameter_from_area(obs["cell_area"]) if "cell_area" in obs else None
    meta = {"platform": "cosmx", "kind": "cell", "coordinate_unit": "micron",
            "spot_diameter": spot if spot else estimate_diameter(obs[["x", "y"]].to_numpy()), "grid": "none",
            "pixel_um": px, "files": {k: _rel(v, root) for k, v in c.files.items()}}
    return [Read(adata, meta, warnings=warnings)]


# --- AnnData / SpatialData --------------------------------------------------------------------

ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{5,}")
SYMBOL_COLUMNS = ("feature_name", "gene_symbols", "gene_symbol", "symbol", "gene_name", "gene_names", "SYMBOL")


def _counts_layer(adata: ad.AnnData, warnings: list[str]) -> sp.csr_matrix:
    """Raw counts of an AnnData: layers["counts"] when present, else X (or raw.X) if integer-valued."""
    def integer_valued(M) -> bool:
        M = sp.csr_matrix(M) if not sp.issparse(M) else M.tocsr()
        d = M.data[:200_000] if M.nnz > 200_000 else M.data
        return bool(np.all(np.mod(d, 1) == 0)) and (d.size == 0 or float(d.max()) < 2**31)

    if "counts" in adata.layers:
        return sp.csr_matrix(adata.layers["counts"])
    X = adata.X
    if X is None:
        raise ValueError("AnnData has no X")
    if integer_valued(X):
        return sp.csr_matrix(X)
    if adata.raw is not None and adata.raw.n_vars == adata.n_vars and integer_valued(adata.raw.X):
        warnings.append("X is not integer-valued; raw counts taken from .raw.X")
        return sp.csr_matrix(adata.raw.X)
    warnings.append("X is not integer-valued and there is no layers['counts']: the values are used as they are")
    return sp.csr_matrix(X)


def _symbols(var: pd.DataFrame) -> tuple[pd.Index, pd.DataFrame, str | None]:
    """Gene symbols as var_names: when var_names are Ensembl ids and a symbol column exists, swap them."""
    names = var.index.astype(str)
    probe = names[: min(len(names), 200)]
    n_ens = sum(1 for n in probe if ENSEMBL_RE.match(n))
    if n_ens >= 0.8 * max(1, len(probe)):
        for col in SYMBOL_COLUMNS:
            if col in var.columns:
                sym = var[col].astype(str).to_numpy()
                var = var.copy()
                if "gene_ids" not in var.columns:
                    var["gene_ids"] = names.to_numpy()
                return pd.Index(sym), var, col
    if "gene_ids" not in var.columns:
        var = var.copy()
        var["gene_ids"] = names.to_numpy()
    return pd.Index(names), var, None


def _adata_to_read(adata: ad.AnnData, source: Path, root: Path, platform: str, sub_id: str | None,
                   coords: np.ndarray, warnings: list[str], meta_extra: dict | None = None) -> Read:
    X = _counts_layer(adata, warnings)
    names, var, sym_col = _symbols(adata.var)
    if sym_col:
        warnings.append(f"var_names were Ensembl ids; gene symbols taken from var['{sym_col}']")
    var = var.copy()
    var.index = names
    if "feature_types" not in var.columns:
        var["feature_types"] = np.where([bool(CONTROL_PATTERNS.match(n)) for n in names], "Control", "Gene Expression")
    obs = pd.DataFrame(index=adata.obs_names.astype(str))
    obs["x"] = coords[:, 0].astype(np.float64)
    obs["y"] = coords[:, 1].astype(np.float64)
    src = adata.obs
    obs["in_tissue"] = (src["in_tissue"].to_numpy(dtype=np.int32) if "in_tissue" in src.columns
                        else np.ones(len(obs), dtype=np.int32))
    for col in ("array_row", "array_col"):
        if col in src.columns:
            obs[col] = src[col].to_numpy(dtype=np.int32)
    for col in ("cell_area", "nucleus_area", "cell_volume"):
        if col in src.columns:
            obs[col] = src[col].to_numpy(dtype=np.float64)
    if "annotation" in src.columns:
        obs["annotation"] = src["annotation"].astype(str).to_numpy()
    out = ad.AnnData(X=X, obs=obs, var=var)
    prev = adata.uns.get("pennsieve") if isinstance(adata.uns, dict) or hasattr(adata.uns, "get") else None
    meta: dict = {"platform": platform, "files": {"source": _rel(source, root)}}
    if isinstance(prev, dict) and prev.get("platform"):
        # a file this app (or a step after it) wrote: keep what it knows
        for k in ("platform", "kind", "coordinate_unit", "spot_diameter", "grid", "bin_size_um"):
            if k in prev:
                meta[k] = prev[k]
        meta["source_platform"] = prev.get("platform")
    spatial_uns = adata.uns.get("spatial") if hasattr(adata.uns, "get") else None
    if isinstance(spatial_uns, dict) and spatial_uns:
        key = next(iter(spatial_uns))
        entry = spatial_uns[key]
        if isinstance(entry, dict) and ("scalefactors" in entry or "images" in entry):
            out.uns["spatial"] = {key: entry}
            sf = entry.get("scalefactors", {}) or {}
            meta.setdefault("kind", "spot")
            meta.setdefault("coordinate_unit", "pixel")
            meta.setdefault("grid", "hex")
            if sf.get("spot_diameter_fullres"):
                meta.setdefault("spot_diameter", float(sf["spot_diameter_fullres"]))
            meta["image"] = {r: r in (entry.get("images") or {}) for r in ("hires", "lowres")}
    if "kind" not in meta:
        meta["kind"] = "cell"
        meta["coordinate_unit"] = "micron"
        meta["grid"] = "none"
        warnings.append("no platform information in the file: kind=cell, coordinate_unit=micron assumed")
    if "spot_diameter" not in meta or not meta["spot_diameter"]:
        d = diameter_from_area(obs["cell_area"]) if "cell_area" in obs else None
        meta["spot_diameter"] = d if d else estimate_diameter(obs[["x", "y"]].to_numpy(),
                                                                0.55 if meta.get("grid") == "hex" else 1.0)
        warnings.append("spot_diameter estimated from the coordinates")
    if meta_extra:
        meta.update(meta_extra)
    return Read(out, meta, sub_id=sub_id, warnings=warnings)


def _split_by_sample(adata: ad.AnnData) -> list[tuple[str | None, ad.AnnData]]:
    if "sample_id" in adata.obs.columns:
        vals = adata.obs["sample_id"].astype(str)
        uniq = list(dict.fromkeys(vals.tolist()))
        if len(uniq) > 1:
            return [(sanitise(u), adata[(vals == u).to_numpy()]) for u in uniq]
    return [(None, adata)]


def read_h5ad(c: Candidate, root: Path, params: dict) -> list[Read]:
    adata = ad.read_h5ad(c.files["h5ad"])
    if "spatial" not in adata.obsm:
        raise ValueError("no obsm['spatial']")
    parts = _split_by_sample(adata)
    out = []
    for sub_id, part in parts:
        warnings: list[str] = []
        coords = np.asarray(part.obsm["spatial"])
        if coords.ndim != 2 or coords.shape[1] < 2:
            raise ValueError(f"obsm['spatial'] has shape {coords.shape}; expected (n, 2)")
        if coords.shape[1] > 2:
            warnings.append(f"obsm['spatial'] has {coords.shape[1]} columns; the first two are used")
            coords = coords[:, :2]
        out.append(_adata_to_read(part, c.source, root, "h5ad", sub_id, coords, warnings))
    if len(parts) > 1:
        for r in out:
            r.warnings.append(f"file holds {len(parts)} samples (obs['sample_id']); split")
    return out


# WKB geometry types
_WKB_POINT, _WKB_POLYGON, _WKB_MULTIPOLYGON = 1, 3, 6


def wkb_centroid(b: bytes) -> tuple[float, float]:
    """Centroid of a WKB Point / Polygon (exterior ring, area-weighted) / MultiPolygon. Raises ValueError otherwise."""
    def parse(off: int) -> tuple[tuple[float, float, float], int]:
        bo = "<" if b[off] == 1 else ">"
        gtype = struct.unpack_from(bo + "I", b, off + 1)[0] & 0xFF
        off += 5
        if gtype == _WKB_POINT:
            x, y = struct.unpack_from(bo + "dd", b, off)
            return (x, y, 1.0), off + 16
        if gtype == _WKB_POLYGON:
            n_rings = struct.unpack_from(bo + "I", b, off)[0]
            off += 4
            cx = cy = area = 0.0
            for r in range(n_rings):
                n = struct.unpack_from(bo + "I", b, off)[0]
                off += 4
                pts = np.frombuffer(b, dtype=bo + "f8", count=2 * n, offset=off).reshape(n, 2)
                off += 16 * n
                if r == 0:
                    x, y = pts[:, 0], pts[:, 1]
                    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
                    a = cross.sum() / 2.0
                    if abs(a) > 1e-12:
                        cx = ((x[:-1] + x[1:]) * cross).sum() / (6 * a)
                        cy = ((y[:-1] + y[1:]) * cross).sum() / (6 * a)
                        area = abs(a)
                    else:
                        cx, cy, area = float(x.mean()), float(y.mean()), 0.0
            return (cx, cy, area), off
        if gtype == _WKB_MULTIPOLYGON:
            n = struct.unpack_from(bo + "I", b, off)[0]
            off += 4
            sx = sy = sa = 0.0
            first = None
            for _ in range(n):
                (cx, cy, a), off = parse(off)
                first = first or (cx, cy)
                sx, sy, sa = sx + cx * a, sy + cy * a, sa + a
            if sa > 0:
                return (sx / sa, sy / sa, sa), off
            return (first[0], first[1], 0.0), off
        raise ValueError(f"unsupported WKB geometry type {gtype}")

    (x, y, _), _ = parse(0)
    return x, y


def _shapes_centroids(element_dir: Path) -> tuple[pd.DataFrame, float | None]:
    """Centroids (x, y) of a SpatialData shapes element (GeoParquet, WKB geometry), indexed by instance id,
    and 2 x median radius when the shapes are circles."""
    import pyarrow.parquet as pq

    pqf = element_dir / "shapes.parquet"
    if not pqf.is_file():
        raise ValueError(f"{element_dir.name}: no shapes.parquet")
    tbl = pq.read_table(pqf)
    df = tbl.to_pandas()
    if "geometry" not in df.columns:
        raise ValueError(f"{element_dir.name}: no geometry column")
    cent = np.array([wkb_centroid(g) for g in df["geometry"].to_numpy()], dtype=np.float64).reshape(-1, 2)
    out = pd.DataFrame({"x": cent[:, 0], "y": cent[:, 1]}, index=df.index.astype(str))
    diameter = None
    if "radius" in df.columns:
        r = df["radius"].to_numpy(dtype=np.float64)
        r = r[np.isfinite(r) & (r > 0)]
        diameter = float(2 * np.median(r)) if len(r) else None
    return out, diameter


def _points_coords(element_dir: Path) -> pd.DataFrame:
    import pyarrow.parquet as pq

    pqf = element_dir / "points.parquet"
    if not pqf.is_file():
        raise ValueError(f"{element_dir.name}: no points.parquet")
    df = pq.read_table(pqf, columns=["x", "y"]).to_pandas()
    df.index = df.index.astype(str)
    return df


def read_spatialdata(c: Candidate, root: Path, params: dict) -> list[Read]:
    """
    A SpatialData zarr store without the spatialdata package: the first table (anndata.read_zarr)
    and its coordinates from obsm["spatial"], else the centroids of the shapes (or points) element
    the table annotates (uns["spatialdata_attrs"]: region / region_key / instance_key).
    """
    tables_dir = c.files["tables"]
    names = sorted(d.name for d in tables_dir.iterdir() if d.is_dir() and not d.name.startswith("."))
    if not names:
        raise ValueError("no table in the store")
    name = names[0]
    warnings: list[str] = []
    if len(names) > 1:
        warnings.append(f"{len(names)} tables ({', '.join(names)}); read {name}")
    adata = ad.read_zarr(tables_dir / name)
    diameter = None
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"])[:, :2]
        region_note = "obsm['spatial']"
    else:
        attrs = adata.uns.get("spatialdata_attrs", {}) if hasattr(adata.uns, "get") else {}
        region = attrs.get("region")
        region_key, instance_key = attrs.get("region_key"), attrs.get("instance_key")
        if region is None or instance_key is None:
            raise ValueError("table has neither obsm['spatial'] nor spatialdata_attrs (region / instance_key)")
        regions = list(region) if isinstance(region, (list, tuple, np.ndarray)) else [region]
        obs_region = (adata.obs[region_key].astype(str).to_numpy() if region_key in adata.obs.columns
                      else np.full(adata.n_obs, regions[0]))
        inst = adata.obs[instance_key].astype(str).to_numpy()
        coords = np.full((adata.n_obs, 2), np.nan)
        for reg in regions:
            rows = np.flatnonzero(obs_region == str(reg))
            if not len(rows):
                continue
            for kind in ("shapes", "points"):
                ed = c.files["zarr"] / kind / str(reg)
                if ed.is_dir():
                    if kind == "shapes":
                        cdf, diameter = _shapes_centroids(ed)
                    else:
                        cdf = _points_coords(ed)
                    hit = cdf.reindex(inst[rows])
                    coords[rows] = hit[["x", "y"]].to_numpy()
                    break
            else:
                raise ValueError(f"region {reg!r} is neither a shapes nor a points element of the store")
        bad = ~np.isfinite(coords).all(axis=1)
        if bad.any():
            warnings.append(f"{int(bad.sum())} observations have no geometry in the shapes element; dropped")
            adata = adata[~bad].copy()
            coords = coords[~bad]
        region_note = f"centroids of {', '.join(map(str, regions))}"
    out = []
    for sub_id, part in _split_by_sample(adata):
        w = list(warnings)
        sub_coords = coords if sub_id is None else coords[member(adata.obs_names, part.obs_names)]
        extra = {"table": name, "coordinates": region_note}
        if diameter:
            extra["spot_diameter"] = diameter
        out.append(_adata_to_read(part, c.source, root, "spatialdata", sub_id, sub_coords, w, extra))
    return out


READERS = {
    "visium": read_visium, "visium-hd": read_visium_hd, "xenium": read_xenium, "merscope": read_merscope,
    "cosmx": read_cosmx, "h5ad": read_h5ad, "spatialdata": read_spatialdata,
}
