"""
A tiny synthetic input tree with one sample per supported platform, small enough for unit tests:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both]

writes data/fixture/linked/spatial/ (what a linked `spatial` input looks like) and/or
data/fixture/flat/ (the same files at the root of INPUT_DIR) with:

  sr1/                 Space Ranger 1.x layout: filtered_feature_bc_matrix.h5 holding every spot (in and out of
                       tissue), spatial/tissue_positions_list.csv (no header), scalefactors, hires/lowres PNGs
  sr2/outs/            Space Ranger 2.x layout: tissue_positions.csv (header), filtered + raw .h5, extra files
  sr3/                 Space Ranger 2.x with MTX folders (filtered_feature_bc_matrix/, raw_feature_bc_matrix/)
  xen/                 Xenium: cell_feature_matrix.h5 (genes + control probes/codewords) + cells.parquet
  mer/                 MERSCOPE: cell_by_gene.csv (Blank-* columns) + cell_metadata.csv
  cos/                 CosMx: S1_exprMat_file.csv (NegPrb*, SystemControl*) + S1_metadata_file.csv
  tab.h5ad             AnnData with Ensembl var_names, var["feature_name"] symbols, obsm["spatial"]
  sd.zarr/             SpatialData-like store: tables/table + shapes/cells/shapes.parquet (WKB points)
  nospatial.h5ad       an AnnData without obsm["spatial"] (ignored)
  notes.txt            a stray file (ignored)
  samples.csv          sample_id, path, condition, batch (one row, `ghost`, names a path that does not exist)
  annotations.csv      labels for some sr1 spots and xen cells
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import struct
from pathlib import Path

import numpy as np

GENES = ["ACTB", "CD3E", "CD8A", "MS4A1", "COL1A1", "EPCAM", "MT-CO1", "MT-ND1"]
ENSEMBL = [f"ENSG{i:011d}" for i in range(1, len(GENES) + 1)]
SR_SCALEFACTORS = {"spot_diameter_fullres": 60.0, "tissue_hires_scalef": 0.1, "fiducial_diameter_fullres": 95.0,
                   "tissue_lowres_scalef": 0.03}


def counts(n: int, g: int, rng: np.random.Generator, lam: float = 3.0) -> np.ndarray:
    return rng.poisson(lam, size=(n, g)).astype(np.int32)


def write_10x_h5(path: Path, X: np.ndarray, barcodes: list[str], names: list[str], ids: list[str],
                 types: list[str], software: str = "spaceranger-2.1.0") -> None:
    """A 10x HDF5 matrix (format version 2): features x barcodes as CSC."""
    import h5py
    import scipy.sparse as sp

    M = sp.csc_matrix(X.T)   # features x barcodes
    with h5py.File(path, "w") as h:
        h.attrs["filetype"] = "matrix"
        h.attrs["version"] = 2
        h.attrs["software_version"] = software
        h.attrs["chemistry_description"] = "Spatial 3' v1"
        m = h.create_group("matrix")
        m.create_dataset("data", data=M.data.astype(np.int32))
        m.create_dataset("indices", data=M.indices.astype(np.int64))
        m.create_dataset("indptr", data=M.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array(M.shape, dtype=np.int32))
        m.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        f = m.create_group("features")
        f.create_dataset("name", data=np.array(names, dtype="S"))
        f.create_dataset("id", data=np.array(ids, dtype="S"))
        f.create_dataset("feature_type", data=np.array(types, dtype="S"))
        f.create_dataset("genome", data=np.array(["GRCh38"] * len(names), dtype="S"))
        f.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


def write_10x_mtx(folder: Path, X: np.ndarray, barcodes: list[str], names: list[str], ids: list[str],
                  types: list[str]) -> None:
    import scipy.sparse as sp
    from scipy.io import mmwrite

    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "matrix.mtx.gz", "wb") as fh:
        mmwrite(fh, sp.csc_matrix(X.T), field="integer")
    with gzip.open(folder / "features.tsv.gz", "wt") as fh:
        for i, n, t in zip(ids, names, types):
            fh.write(f"{i}\t{n}\t{t}\n")
    with gzip.open(folder / "barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(f"{b}\n" for b in barcodes))


def write_png(path: Path, size: tuple[int, int], seed: int) -> None:
    from PIL import Image

    rng = np.random.default_rng(seed)
    Image.fromarray(rng.integers(100, 230, size=(size[1], size[0], 3), dtype=np.uint8), "RGB").save(path)


def visium_spots(n_rows: int, n_cols: int, in_tissue_fraction: float, rng: np.random.Generator):
    """Barcodes and a hex-like grid of spots: (barcode, in_tissue, array_row, array_col, pxl_row, pxl_col)."""
    rows = []
    k = 0
    for r in range(n_rows):
        for c in range(n_cols):
            k += 1
            bc = f"{''.join(rng.choice(list('ACGT'), 16))}-1"
            in_t = 1 if rng.random() < in_tissue_fraction else 0
            rows.append((bc, in_t, r, 2 * c + (r % 2), 1000 + r * 100, 1000 + (2 * c + (r % 2)) * 50))
    return rows


def make_visium(folder: Path, layout: str, rng: np.random.Generator) -> list[tuple]:
    """layout: 'sr1' (1.x, matrix holds every spot), 'sr2' (2.x, filtered + raw .h5), 'sr3' (2.x, MTX folders)."""
    spots = visium_spots(5, 6, 0.8, rng)
    spots[0] = (spots[0][0], 1, *spots[0][2:])     # at least one in tissue
    folder.mkdir(parents=True, exist_ok=True)
    sp_dir = folder / "spatial"
    sp_dir.mkdir()
    X_all = counts(len(spots), len(GENES), rng)
    all_bc = [s[0] for s in spots]
    in_bc = [s[0] for s in spots if s[1] == 1]
    X_in = X_all[[i for i, s in enumerate(spots) if s[1] == 1]]
    types = ["Gene Expression"] * len(GENES)
    if layout == "sr1":
        write_10x_h5(folder / "filtered_feature_bc_matrix.h5", X_all, all_bc, GENES, ENSEMBL, types, "spaceranger-1.1.0")
        with open(sp_dir / "tissue_positions_list.csv", "w", newline="") as fh:
            csv.writer(fh).writerows(spots)
    elif layout == "sr2":
        write_10x_h5(folder / "filtered_feature_bc_matrix.h5", X_in, in_bc, GENES, ENSEMBL, types)
        write_10x_h5(folder / "raw_feature_bc_matrix.h5", X_all, all_bc, GENES, ENSEMBL, types)
        with open(sp_dir / "tissue_positions.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
            w.writerows(spots)
        (sp_dir / "spatial_enrichment.csv").write_text("feature_id,feature_name\n")
        (sp_dir / "detected_tissue_image.jpg").write_bytes(b"\xff\xd8\xff\xd9")
    else:
        write_10x_mtx(folder / "filtered_feature_bc_matrix", X_in, in_bc, GENES, ENSEMBL, types)
        write_10x_mtx(folder / "raw_feature_bc_matrix", X_all, all_bc, GENES, ENSEMBL, types)
        with open(sp_dir / "tissue_positions.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
            w.writerows(spots)
    with open(sp_dir / "scalefactors_json.json", "w") as fh:
        json.dump(SR_SCALEFACTORS, fh)
    write_png(sp_dir / "tissue_hires_image.png", (40, 38), 1)
    write_png(sp_dir / "tissue_lowres_image.png", (12, 12), 2)
    return spots


XEN_CONTROLS = [("NegControlProbe_00001", "Negative Control Probe"), ("NegControlProbe_00002", "Negative Control Probe"),
                ("NegControlCodeword_0500", "Negative Control Codeword"), ("BLANK_0001", "Blank Codeword"),
                ("BLANK_0002", "Blank Codeword"), ("UnassignedCodeword_0001", "Unassigned Codeword")]


def make_xenium(folder: Path, rng: np.random.Generator, n: int = 40) -> dict:
    import pandas as pd

    folder.mkdir(parents=True, exist_ok=True)
    genes = GENES[:6]
    names = genes + [c for c, _ in XEN_CONTROLS]
    ids = ENSEMBL[:6] + [c for c, _ in XEN_CONTROLS]
    types = ["Gene Expression"] * 6 + [t for _, t in XEN_CONTROLS]
    X = np.hstack([counts(n, 6, rng, 4.0), counts(n, len(XEN_CONTROLS), rng, 0.3)])
    cell_ids = [str(i + 1) for i in range(n)]
    write_10x_h5(folder / "cell_feature_matrix.h5", X, cell_ids, names, ids, types, "xenium-1.0.1")
    area = rng.uniform(40, 120, n)
    cells = pd.DataFrame({
        "cell_id": np.array(cell_ids, dtype=object).astype(str), "x_centroid": rng.uniform(0, 500, n),
        "y_centroid": rng.uniform(0, 400, n), "transcript_counts": X[:, :6].sum(1),
        "control_probe_counts": X[:, 6:8].sum(1), "control_codeword_counts": X[:, 8:].sum(1),
        "total_counts": X.sum(1), "cell_area": area, "nucleus_area": area * 0.4,
    })
    cells["cell_id"] = cells["cell_id"].astype(np.int32)
    cells.to_parquet(folder / "cells.parquet", index=False)
    return {"X": X, "cell_ids": cell_ids, "control_counts": X[:, 6:].sum(1), "cells": cells}


def make_merscope(folder: Path, rng: np.random.Generator, n: int = 30) -> dict:
    import pandas as pd

    folder.mkdir(parents=True, exist_ok=True)
    genes = GENES[:5] + ["Blank-0", "Blank-1"]
    X = np.hstack([counts(n, 5, rng, 4.0), counts(n, 2, rng, 0.2)])
    ids = [str(100000 + i) for i in range(n)]
    df = pd.DataFrame(X, columns=genes, index=ids)
    df.index.name = "cell"
    df.to_csv(folder / "cell_by_gene.csv")
    cx, cy = rng.uniform(0, 1000, n), rng.uniform(0, 800, n)
    meta = pd.DataFrame({"fov": rng.integers(0, 3, n), "volume": rng.uniform(500, 2000, n), "center_x": cx, "center_y": cy,
                         "min_x": cx - 5, "max_x": cx + 5, "min_y": cy - 4, "max_y": cy + 4}, index=ids)
    meta.index.name = ""
    meta.to_csv(folder / "cell_metadata.csv")
    return {"X": X, "ids": ids, "control_counts": X[:, 5:].sum(1), "meta": meta}


def make_cosmx(folder: Path, rng: np.random.Generator, n: int = 30) -> dict:
    import pandas as pd

    folder.mkdir(parents=True, exist_ok=True)
    genes = GENES[:5] + ["NegPrb1", "NegPrb2", "SystemControl1"]
    X = np.hstack([counts(n, 5, rng, 4.0), counts(n, 3, rng, 0.2)])
    fov = np.repeat([1, 2], n // 2)
    cell_id = np.concatenate([np.arange(1, n // 2 + 1), np.arange(1, n // 2 + 1)])
    expr = pd.DataFrame(X, columns=genes)
    expr.insert(0, "cell_ID", cell_id)
    expr.insert(0, "fov", fov)
    expr.to_csv(folder / "S1_exprMat_file.csv", index=False)
    meta = pd.DataFrame({"fov": fov, "cell_ID": cell_id, "Area": rng.uniform(1000, 4000, n), "AspectRatio": 1.1,
                         "CenterX_local_px": rng.uniform(0, 4000, n), "CenterY_local_px": rng.uniform(0, 4000, n),
                         "CenterX_global_px": rng.uniform(0, 20000, n), "CenterY_global_px": rng.uniform(0, 20000, n),
                         "Width": 40, "Height": 40})
    meta.to_csv(folder / "S1_metadata_file.csv", index=False)
    ids = [f"{c}_{f}" for c, f in zip(cell_id, fov)]
    return {"X": X, "ids": ids, "control_counts": X[:, 5:].sum(1), "meta": meta}


def make_h5ad(path: Path, rng: np.random.Generator, n: int = 25, spatial: bool = True) -> dict:
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    X = counts(n, len(GENES), rng).astype(np.float32)
    obs = pd.DataFrame({"donor": ["d1"] * n}, index=[f"cell{i}" for i in range(n)])
    var = pd.DataFrame({"feature_name": GENES}, index=ENSEMBL)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    if spatial:
        a.obsm["spatial"] = np.column_stack([rng.uniform(0, 300, n), rng.uniform(0, 300, n)])
    a.write_h5ad(path)
    return {"X": X, "obs_ids": list(obs.index)}


def wkb_point(x: float, y: float) -> bytes:
    return struct.pack("<BIdd", 1, 1, x, y)


def wkb_polygon(pts: list[tuple[float, float]]) -> bytes:
    ring = pts + [pts[0]]
    return struct.pack("<BII", 1, 3, 1) + struct.pack("<I", len(ring)) + b"".join(struct.pack("<dd", *p) for p in ring)


def make_spatialdata(folder: Path, rng: np.random.Generator, n: int = 20) -> dict:
    """A SpatialData-shaped zarr store: tables/table (AnnData) + shapes/cells/shapes.parquet with WKB geometry."""
    import anndata as ad
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
    import scipy.sparse as sp

    folder.mkdir(parents=True, exist_ok=True)
    X = counts(n, 5, rng)
    inst = np.arange(n)
    obs = pd.DataFrame({"region": pd.Categorical(["cells"] * n), "cell_id": inst}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=GENES[:5]))
    a.uns["spatialdata_attrs"] = {"region": "cells", "region_key": "region", "instance_key": "cell_id"}
    (folder / "tables").mkdir()
    a.write_zarr(folder / "tables" / "table")
    cx, cy = rng.uniform(0, 100, n), rng.uniform(0, 100, n)
    # half points, half squares centred on (cx, cy): the centroid is the same either way
    geoms = [wkb_point(cx[i], cy[i]) if i % 2 == 0 else
             wkb_polygon([(cx[i] - 2, cy[i] - 2), (cx[i] + 2, cy[i] - 2), (cx[i] + 2, cy[i] + 2), (cx[i] - 2, cy[i] + 2)])
             for i in range(n)]
    (folder / "shapes" / "cells").mkdir(parents=True)
    tbl = pa.table({"radius": pa.array([3.0] * n), "geometry": pa.array(geoms, pa.binary()),
                    "__index_level_0__": pa.array(inst)})
    pq.write_table(tbl, folder / "shapes" / "cells" / "shapes.parquet")
    (folder / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')
    return {"X": X, "x": cx, "y": cy, "obs_ids": list(obs.index)}


def make_tree(root: Path, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    root.mkdir(parents=True, exist_ok=True)
    planted = {}
    planted["sr1"] = make_visium(root / "sr1", "sr1", rng)
    planted["sr2"] = make_visium(root / "sr2" / "outs", "sr2", rng)
    planted["sr3"] = make_visium(root / "sr3", "sr3", rng)
    planted["xen"] = make_xenium(root / "xen", rng)
    planted["mer"] = make_merscope(root / "mer", rng)
    planted["cos"] = make_cosmx(root / "cos", rng)
    planted["tab"] = make_h5ad(root / "tab.h5ad", rng)
    make_h5ad(root / "nospatial.h5ad", rng, n=5, spatial=False)
    planted["sd"] = make_spatialdata(root / "sd.zarr", rng)
    (root / "notes.txt").write_text("not a sample\n")
    with open(root / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "path", "condition", "batch"])
        w.writerows([["visium1", "sr1", "ctrl", "b1"], ["visium2", "sr2", "case", "b1"], ["visium3", "sr3", "case", "b2"],
                     ["xen1", "xen", "case", "b2"], ["mer1", "mer", "ctrl", "b2"], ["cos1", "cos", "ctrl", "b3"],
                     ["tab1", "tab.h5ad", "case", "b3"], ["sd1", "sd.zarr", "ctrl", "b3"],
                     ["ghost", "ghost", "case", "b9"]])
    in_spots = [s for s in planted["sr1"] if s[1] == 1]
    ann = [("visium1", s[0], "A" if i % 2 else "B") for i, s in enumerate(in_spots[:5])]
    ann += [("visium1", "NOTABARCODE-1", "A")]
    ann += [("xen1", cid, "T cell" if i % 3 else "Tumour") for i, cid in enumerate(planted["xen"]["cell_ids"][:10])]
    with open(root / "annotations.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "obs_id", "annotation"])
        w.writerows(ann)
    planted["annotations"] = ann
    return planted


def make_fixture(root: Path, layout: str = "both") -> dict[str, Path]:
    root = Path(root)
    out: dict[str, Path] = {}
    if layout in ("linked", "both"):
        d = root / "linked"
        shutil.rmtree(d, ignore_errors=True)
        make_tree(d / "spatial")
        out["linked"] = d
    if layout in ("flat", "both"):
        d = root / "flat"
        shutil.rmtree(d, ignore_errors=True)
        make_tree(d)
        out["flat"] = d
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    args = ap.parse_args()
    for k, v in make_fixture(Path(args.root), args.layout).items():
        print(f"{k}: {v}")
