"""
From a reader's raw AnnData to the contract (docs/conventions.md, "The AnnData contract"):
control features summed into obs["control_counts"] and dropped, gene symbols as unique
var_names, out-of-tissue spots dropped, the seeded subsample, QC metrics, obs_names
prefixed with the sample id, sample sheet columns and annotations on obs, obsm["spatial"],
uns["pennsieve"].
"""
from __future__ import annotations

import zlib
from collections import OrderedDict

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from processor.readers import CONTROL_PATTERNS, Read

# 10x / platform feature types that are controls: summed into obs["control_counts"], then dropped.
CONTROL_TYPES = {
    "Negative Control Probe", "Negative Control Codeword", "Blank Codeword", "Unassigned Codeword",
    "Deprecated Codeword", "Genomic Control", "System Control", "False Code", "Control",
}
IMAGING = {"xenium", "merscope", "cosmx"}
OBS_ORDER = ["sample_id", "x", "y", "in_tissue", "array_row", "array_col", "total_counts", "n_genes_by_counts",
             "pct_counts_mt", "cell_area", "nucleus_area", "cell_volume", "control_counts", "annotation"]


def as_counts(X, warnings: list[str]) -> sp.csr_matrix:
    """CSR with int32 data when the values are integers (float32 otherwise), int32 indices when they fit."""
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    X.sum_duplicates()
    if X.dtype.kind == "f":
        d = X.data
        if np.all(np.isfinite(d)) and np.all(np.mod(d, 1) == 0) and (d.size == 0 or d.max() < 2**31):
            X = X.astype(np.int32)
        else:
            warnings.append("counts are not integer-valued; kept as float32")
            X = X.astype(np.float32)
    elif X.dtype != np.int32:
        if X.nnz and X.data.max() >= 2**31:
            warnings.append("counts above 2^31; kept as float32")
            X = X.astype(np.float32)
        else:
            X = X.astype(np.int32)
    if X.nnz < 2**31 - 1 and X.shape[1] < 2**31 - 1:
        X.indices = X.indices.astype(np.int32)
        X.indptr = X.indptr.astype(np.int32)
    return X


def split_features(var: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(gene mask, control mask, other mask) over var."""
    types = var["feature_types"].astype(str).to_numpy()
    names = var.index.astype(str).to_numpy()
    gene = types == "Gene Expression"
    control = np.isin(types, list(CONTROL_TYPES)) | (gene & np.array([bool(CONTROL_PATTERNS.match(n)) for n in names]))
    gene = gene & ~control
    other = ~gene & ~control
    return gene, control, other


def subsample_index(n: int, max_obs: int, seed: int, sample_id: str) -> np.ndarray | None:
    """A seeded subset of row indices in source order; None when nothing is to be dropped."""
    if max_obs <= 0 or n <= max_obs:
        return None
    rng = np.random.default_rng([seed, zlib.crc32(sample_id.encode())])
    return np.sort(rng.choice(n, size=max_obs, replace=False))


def finalise(read: Read, sample_id: str, source_rel: str, params: dict, extra: dict[str, str],
             annotations: dict[str, str] | None, prep_version: str) -> tuple[ad.AnnData, dict]:
    """The contract AnnData and the prep.json record for one sample."""
    a = read.adata
    meta = read.meta
    warnings = list(read.warnings)
    info: dict = {"platform": meta["platform"], "kind": meta["kind"], "coordinate_unit": meta["coordinate_unit"],
                  "files": meta.get("files", {})}

    X = as_counts(a.X, warnings)
    var = a.var.copy()
    var.index = var.index.astype(str)
    n_obs_source, n_feat_source = X.shape

    # 1. features: genes kept, controls summed then dropped, anything else dropped
    gene, control, other = split_features(var)
    dropped: dict[str, dict] = {}
    control_counts = None
    if control.any():
        control_counts = np.asarray(X[:, np.flatnonzero(control)].sum(axis=1)).ravel()
        for t in sorted(set(var["feature_types"].astype(str).to_numpy()[control])):
            m = control & (var["feature_types"].astype(str).to_numpy() == t)
            dropped[t] = {"n_features": int(m.sum()), "total_counts": int(X[:, np.flatnonzero(m)].sum()),
                          "summed_into_control_counts": True}
    if other.any():
        for t in sorted(set(var["feature_types"].astype(str).to_numpy()[other])):
            m = other & (var["feature_types"].astype(str).to_numpy() == t)
            dropped[t] = {"n_features": int(m.sum()), "total_counts": int(X[:, np.flatnonzero(m)].sum()),
                          "summed_into_control_counts": False}
        warnings.append("non-gene features dropped: " + ", ".join(f"{t} ({d['n_features']})" for t, d in dropped.items()
                                                                  if not d["summed_into_control_counts"]))
    if not gene.any():
        raise ValueError("no Gene Expression features")
    X = X[:, np.flatnonzero(gene)]
    var = var.iloc[np.flatnonzero(gene)].copy()
    info["features"] = {"n_source": int(n_feat_source), "n_genes": int(gene.sum()), "dropped": dropped}

    # 2. observations: out-of-tissue, then the subsample
    obs = a.obs.copy()
    obs.index = obs.index.astype(str)
    if control_counts is not None or meta["platform"] in IMAGING:
        obs["control_counts"] = (control_counts if control_counts is not None else np.zeros(len(obs))).astype(np.int64)
    keep = np.ones(len(obs), dtype=bool)
    n_out = int((obs["in_tissue"].to_numpy() == 0).sum()) if "in_tissue" in obs else 0
    if n_out and not params["KEEP_OUT_OF_TISSUE"]:
        keep &= obs["in_tissue"].to_numpy() != 0
    n_after_tissue = int(keep.sum())
    idx = np.flatnonzero(keep)
    sub = subsample_index(n_after_tissue, params["MAX_OBS_PER_SAMPLE"], params["SEED"], sample_id)
    if sub is not None:
        idx = idx[sub]
    if len(idx) == 0:
        raise ValueError("no observations left" + (" (all out of tissue)" if n_out == n_obs_source else ""))
    if len(idx) < n_obs_source:
        X = X[idx]
        obs = obs.iloc[idx].copy()
    info["observations"] = {
        "n_source": int(n_obs_source), "n_out_of_tissue": n_out,
        "n_out_of_tissue_dropped": int(n_obs_source - n_after_tissue), "n_out_of_tissue_kept": int(n_out) if params["KEEP_OUT_OF_TISSUE"] else 0,
        "n_subsampled_to": int(len(idx)) if sub is not None else None, "n_written": int(len(idx)),
    }

    # 3. var: symbols made unique, mt flag
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names = adata.var_names.astype(str)
    n_dup = int(adata.var_names.duplicated().sum())
    if n_dup:
        adata.var_names_make_unique()
        warnings.append(f"{n_dup} duplicated gene symbols made unique")
    adata.var["mt"] = np.asarray(adata.var_names.str.startswith(("MT-", "mt-")), dtype=bool)
    for col in ("gene_ids", "feature_types"):
        adata.var[col] = adata.var[col].astype(str)

    # 4. QC metrics (scanpy)
    import scanpy as sc

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    # per-gene QC columns would go stale once st-qc drops observations; only the per-observation metrics stay
    adata.var = adata.var.drop(columns=["n_cells_by_counts", "mean_counts", "pct_dropout_by_counts", "total_counts"],
                               errors="ignore")
    adata.obs["total_counts"] = adata.obs["total_counts"].astype(np.float64)
    adata.obs["n_genes_by_counts"] = adata.obs["n_genes_by_counts"].astype(np.int64)
    adata.obs["pct_counts_mt"] = adata.obs["pct_counts_mt"].astype(np.float64)

    # 5. identities, sheet columns, annotations
    # an AnnData this app wrote before already carries "<sample_id>:" on its names; strip it, so a
    # contract file read again keeps its ids (and annotations still join on the barcode)
    pre = f"{sample_id}:"
    obs_ids = np.array([str(o)[len(pre):] if str(o).startswith(pre) else str(o) for o in adata.obs_names], dtype=object)
    adata.obs_names = pd.Index([f"{sample_id}:{o}" for o in obs_ids])
    adata.obs["sample_id"] = pd.Categorical([sample_id] * adata.n_obs)
    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype(np.int32)
    for col in ("array_row", "array_col"):
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].astype(np.int32)
    ann_info = None
    if annotations is not None:
        labels = pd.Series(obs_ids).map(annotations)
        n_hit = int(labels.notna().sum())
        if "annotation" in adata.obs:
            warnings.append("obs['annotation'] of the file replaced by annotations.csv")
        adata.obs["annotation"] = pd.Categorical(labels.to_numpy())
        ann_info = {"n_annotated": n_hit, "n_in_file": len(annotations),
                    "n_not_matched": int(len(set(annotations) - set(obs_ids)))}
        if n_hit == 0:
            warnings.append("annotations.csv has rows for this sample but no obs_id matched")
    elif "annotation" in adata.obs:
        adata.obs["annotation"] = pd.Categorical(adata.obs["annotation"].astype(str).to_numpy())
    info["annotation"] = ann_info
    for col, val in extra.items():
        if col in OBS_ORDER or col in ("obs_id",):
            warnings.append(f"sample sheet column {col!r} clashes with a contract column; not copied to obs")
            continue
        adata.obs[col] = pd.Categorical([val] * adata.n_obs)
    for col in list(adata.obs.columns):
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype) in ("str", "string"):
            adata.obs[col] = pd.Categorical(adata.obs[col].astype(str).to_numpy())
    order = [c for c in OBS_ORDER if c in adata.obs.columns] + [c for c in adata.obs.columns if c not in OBS_ORDER]
    adata.obs = adata.obs[order]

    # 6. spatial
    adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy(dtype=np.float64)
    if "spatial" in a.uns:
        entry = next(iter(a.uns["spatial"].values()))
        adata.uns["spatial"] = {sample_id: entry}
    pen = OrderedDict()
    pen["platform"] = str(meta["platform"])
    pen["kind"] = str(meta["kind"])
    pen["coordinate_unit"] = str(meta["coordinate_unit"])
    pen["spot_diameter"] = float(meta["spot_diameter"])
    pen["grid"] = str(meta.get("grid", "none"))
    if meta.get("bin_size_um") is not None:
        pen["bin_size_um"] = int(meta["bin_size_um"])
    pen["source"] = source_rel
    pen["prep_version"] = prep_version
    pen["sample_id"] = sample_id
    pen["n_obs_source"] = int(n_obs_source)
    if sub is not None:
        pen["subsampled_from"] = int(n_after_tissue)
    adata.uns["pennsieve"] = dict(pen)

    info.update({
        "spot_diameter": pen["spot_diameter"], "grid": pen["grid"],
        "image": meta.get("image", {"hires": False, "lowres": False}),
        "n_obs": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "n_obs_in_tissue": int((adata.obs["in_tissue"] == 1).sum()),
        "warnings": warnings,
    })
    for k in ("bin_size_um", "microns_per_pixel", "pixel_um", "software_version", "table", "coordinates",
              "n_positions", "n_positions_in_tissue", "source_platform"):
        if meta.get(k) is not None:
            info[k] = meta[k]
    return adata, info
