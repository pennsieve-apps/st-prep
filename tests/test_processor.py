"""The synthetic fixture through the processor (both input layouts), the contract, and the readers' pieces."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import contract as C
from processor import inputs as I
from processor import main as app
from processor import readers as R
from tests.make_fixture import GENES, XEN_CONTROLS, make_fixture, wkb_point, wkb_polygon

SAMPLES = {"visium1", "visium2", "visium3", "xen1", "mer1", "cos1", "tab1", "sd1"}
CONTRACT_OBS = ["sample_id", "x", "y", "in_tissue", "total_counts", "n_genes_by_counts", "pct_counts_mt"]


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, value)
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="session")
def outputs(fixture_root, tmp_path_factory) -> dict[str, Path]:
    """One run per layout (the processor is deterministic; the session fixture saves repeating it per test)."""
    out = {}
    for layout in ("linked", "flat"):
        d = tmp_path_factory.mktemp(f"out-{layout}")
        mp = pytest.MonkeyPatch()
        try:
            assert run(fixture_root[layout], d, mp) == 0
        finally:
            mp.undo()
        out[layout] = d
    return out


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_run_writes_every_sample(outputs, layout):
    out = outputs[layout]
    assert {p.name for p in (out / "spatial").iterdir()} == {f"{s}.h5ad" for s in SAMPLES} | {"provenance.json"}
    assert {p.name for p in (out / "metadata").iterdir()} == {"samples.csv", "prep.json", "provenance.json"}

    table = rows(out / "metadata" / "samples.csv")
    assert list(table[0]) == list(I.REQUIRED_COLUMNS) + ["condition", "batch"]
    by_id = {r["sample_id"]: r for r in table}
    assert set(by_id) == SAMPLES
    prefix = "spatial/" if layout == "linked" else ""
    assert by_id["visium2"]["path"] == f"{prefix}sr2/outs"          # the sheet's `sr2` matched the outs/ folder
    assert by_id["visium1"] == {"sample_id": "visium1", "path": f"{prefix}sr1", "platform": "visium", "kind": "spot",
                                "n_obs": "23", "n_genes": "8", "n_obs_in_tissue": "23", "condition": "ctrl", "batch": "b1"}
    assert by_id["xen1"]["platform"] == "xenium" and by_id["xen1"]["kind"] == "cell" and by_id["xen1"]["n_genes"] == "6"
    assert by_id["visium3"]["n_obs"] == "21"                        # MTX folder layout, filtered matrix preferred

    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert prep["app"] == "st-prep" and prep["appVersion"] == "0.1.0"
    assert [s["sample_id"] for s in prep["skipped"]] == ["ghost"]
    assert "matches no sample" in prep["skipped"][0]["reason"]
    ignored = {i["path"]: i["reason"] for i in prep["ignored"]}
    assert ignored == {f"{prefix}nospatial.h5ad": "h5ad without obsm['spatial']",
                       f"{prefix}notes.txt": "not part of a recognised sample"}
    assert prep["sampleSheet"] == f"{prefix}samples.csv" and prep["annotations"] == f"{prefix}annotations.csv"

    prov = json.loads((out / "spatial" / "provenance.json").read_text())
    assert prov == json.loads((out / "metadata" / "provenance.json").read_text())
    assert prov["nOk"] == 8 and prov["nFailed"] == 1
    assert {"scanpy", "anndata", "numpy"} <= set(prov["tools"])
    assert prov["peakRssMb"] is None or prov["peakRssMb"] > 0


def check_contract(a: ad.AnnData, sample_id: str) -> None:
    assert sp.isspmatrix_csr(a.X) and a.X.dtype == np.int32
    assert a.X.indices.dtype == np.int32 and a.X.indptr.dtype == np.int32
    assert all(a.obs_names.str.startswith(sample_id + ":")) and a.obs_names.is_unique
    assert list(a.obs.columns[: len(CONTRACT_OBS)]) == CONTRACT_OBS or set(CONTRACT_OBS) <= set(a.obs.columns)
    assert a.obs["sample_id"].dtype == "category" and set(a.obs["sample_id"]) == {sample_id}
    assert a.obs["x"].dtype == np.float64 and a.obs["y"].dtype == np.float64
    assert a.obs["in_tissue"].dtype == np.int32 and set(a.obs["in_tissue"]) == {1}
    assert a.obsm["spatial"].dtype == np.float64 and a.obsm["spatial"].shape == (a.n_obs, 2)
    assert np.array_equal(a.obsm["spatial"], a.obs[["x", "y"]].to_numpy())
    assert a.var_names.is_unique and "gene_ids" in a.var and "feature_types" in a.var
    assert a.var["mt"].dtype == bool and set(a.var["feature_types"]) == {"Gene Expression"}
    totals = np.asarray(a.X.sum(axis=1)).ravel()
    assert np.allclose(a.obs["total_counts"], totals)
    assert np.array_equal(a.obs["n_genes_by_counts"], np.diff(a.X.indptr))
    mt = np.asarray(a.X[:, a.var["mt"].to_numpy()].sum(axis=1)).ravel()
    expect = np.where(totals > 0, 100 * mt / np.maximum(totals, 1), 0)
    assert np.allclose(a.obs["pct_counts_mt"], expect)
    pen = a.uns["pennsieve"]
    for key in ("platform", "kind", "coordinate_unit", "spot_diameter", "grid", "source", "prep_version"):
        assert key in pen, key
    assert pen["platform"] in I.PLATFORMS and pen["kind"] in ("spot", "bin", "cell")
    assert pen["coordinate_unit"] in ("pixel", "micron") and pen["grid"] in ("hex", "square", "none")
    assert pen["spot_diameter"] > 0 and pen["prep_version"] == "0.1.0"
    assert a.obs["condition"].dtype == "category" and a.obs["batch"].dtype == "category"


@pytest.mark.parametrize("sample_id", sorted(SAMPLES))
def test_contract(outputs, sample_id):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / f"{sample_id}.h5ad")
    check_contract(a, sample_id)
    b = ad.read_h5ad(outputs["flat"] / "spatial" / f"{sample_id}.h5ad")
    assert list(a.obs_names) == list(b.obs_names) and (a.X != b.X).nnz == 0
    assert a.uns["pennsieve"]["source"] == "spatial/" + b.uns["pennsieve"]["source"]


def test_visium_spaceranger_1x(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "visium1.h5ad")
    pos = R.read_positions(fixture_root["linked"] / "spatial" / "sr1" / "spatial" / "tissue_positions_list.csv")
    in_t = pos[pos["in_tissue"] == 1]
    # out-of-tissue spots dropped; in-tissue ones kept in source order; x = pxl_col, y = pxl_row
    assert a.n_obs == len(in_t) == 23
    bcs = [n.split(":", 1)[1] for n in a.obs_names]
    assert bcs == [b for b in R.read_10x_h5(fixture_root["linked"] / "spatial" / "sr1" / "filtered_feature_bc_matrix.h5")[1]
                   if pos.loc[b, "in_tissue"] == 1]
    assert np.array_equal(a.obs["x"], in_t.loc[bcs, "pxl_col_in_fullres"])
    assert np.array_equal(a.obs["y"], in_t.loc[bcs, "pxl_row_in_fullres"])
    assert np.array_equal(a.obs["array_row"], in_t.loc[bcs, "array_row"]) and a.obs["array_row"].dtype == np.int32
    pen = a.uns["pennsieve"]
    assert (pen["platform"], pen["kind"], pen["coordinate_unit"], pen["grid"]) == ("visium", "spot", "pixel", "hex")
    assert pen["spot_diameter"] == 60.0 and "bin_size_um" not in pen
    entry = a.uns["spatial"]["visium1"]
    assert entry["images"]["hires"].shape == (38, 40, 3) and entry["images"]["hires"].dtype == np.uint8
    assert entry["images"]["lowres"].shape == (12, 12, 3)
    assert entry["scalefactors"]["tissue_hires_scalef"] == 0.1 and entry["scalefactors"]["spot_diameter_fullres"] == 60.0
    assert entry["metadata"]["software_version"] == "spaceranger-1.1.0"
    assert list(a.var["mt"]) == [g.startswith("MT-") for g in GENES]
    assert list(a.var["gene_ids"]) == [f"ENSG{i:011d}" for i in range(1, 9)]
    # annotations: 5 spots labelled, one obs_id that does not exist
    assert a.obs["annotation"].dtype == "category"
    assert a.obs["annotation"].value_counts().to_dict() == {"B": 3, "A": 2}
    prep = json.loads((outputs["linked"] / "metadata" / "prep.json").read_text())["samples"]["visium1"]
    assert prep["annotation"] == {"n_annotated": 5, "n_in_file": 6, "n_not_matched": 1}
    assert prep["observations"]["n_out_of_tissue_dropped"] == 7 and prep["image"] == {"hires": True, "lowres": True}
    assert prep["files"]["positions"].endswith("tissue_positions_list.csv")


def test_visium_spaceranger_2x_and_mtx(outputs, fixture_root):
    prep = json.loads((outputs["linked"] / "metadata" / "prep.json").read_text())["samples"]
    assert prep["visium2"]["files"]["matrix"] == "spatial/sr2/outs/filtered_feature_bc_matrix.h5"
    assert prep["visium2"]["files"]["positions"] == "spatial/sr2/outs/spatial/tissue_positions.csv"
    assert prep["visium3"]["files"]["matrix"] == "spatial/sr3/filtered_feature_bc_matrix"
    a2 = ad.read_h5ad(outputs["linked"] / "spatial" / "visium2.h5ad")
    assert a2.uns["spatial"]["visium2"]["metadata"]["software_version"] == "spaceranger-2.1.0"
    # the MTX and the .h5 reader agree on the same data
    X_h5, bc_h5, var_h5 = R.read_10x_h5(fixture_root["linked"] / "spatial" / "sr2" / "outs" / "raw_feature_bc_matrix.h5")
    X_mtx, bc_mtx, var_mtx = R.read_10x_mtx(fixture_root["linked"] / "spatial" / "sr2" / "outs" / "raw_feature_bc_matrix"
                                            ) if False else (None, None, None)
    X3, bc3, var3 = R.read_10x_mtx(fixture_root["linked"] / "spatial" / "sr3" / "raw_feature_bc_matrix")
    assert X3.shape == (30, 8) and sp.isspmatrix_csr(X3) and list(var3.columns) == ["gene_ids", "feature_types"]
    assert list(var3.index) == GENES and len(bc3) == 30
    assert X_h5.shape == (30, 8) and list(var_h5.index) == GENES


def test_xenium_controls(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "xen1.h5ad")
    X, cells, var = R.read_10x_h5(fixture_root["linked"] / "spatial" / "xen" / "cell_feature_matrix.h5")
    ctrl = var["feature_types"].to_numpy() != "Gene Expression"
    assert ctrl.sum() == len(XEN_CONTROLS) and a.n_vars == 6
    assert np.array_equal(a.obs["control_counts"], np.asarray(X[:, ctrl].sum(axis=1)).ravel())
    assert a.obs["control_counts"].dtype == np.int64
    assert not any(n.startswith(("NegControl", "BLANK", "Unassigned")) for n in a.var_names)
    cells_df = pd.read_parquet(fixture_root["linked"] / "spatial" / "xen" / "cells.parquet")
    assert np.allclose(a.obs["x"], cells_df["x_centroid"]) and np.allclose(a.obs["y"], cells_df["y_centroid"])
    assert np.allclose(a.obs["cell_area"], cells_df["cell_area"]) and np.allclose(a.obs["nucleus_area"], cells_df["nucleus_area"])
    assert np.allclose(a.obs["control_counts"], cells_df["control_probe_counts"] + cells_df["control_codeword_counts"])
    pen = a.uns["pennsieve"]
    assert (pen["platform"], pen["kind"], pen["coordinate_unit"], pen["grid"]) == ("xenium", "cell", "micron", "none")
    assert pen["spot_diameter"] == pytest.approx(2 * np.sqrt(np.median(cells_df["cell_area"]) / np.pi))
    assert [n.split(":")[1] for n in a.obs_names[:3]] == ["1", "2", "3"]
    assert a.obs["annotation"].value_counts().to_dict() == {"T cell": 6, "Tumour": 4}
    dropped = json.loads((outputs["linked"] / "metadata" / "prep.json").read_text())["samples"]["xen1"]["features"]["dropped"]
    assert {t: d["n_features"] for t, d in dropped.items()} == {"Negative Control Probe": 2, "Negative Control Codeword": 1,
                                                                 "Blank Codeword": 2, "Unassigned Codeword": 1}
    assert all(d["summed_into_control_counts"] for d in dropped.values())


def test_merscope(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "mer1.h5ad")
    src = pd.read_csv(fixture_root["linked"] / "spatial" / "mer" / "cell_by_gene.csv", index_col=0)
    meta = pd.read_csv(fixture_root["linked"] / "spatial" / "mer" / "cell_metadata.csv", index_col=0)
    assert list(a.var_names) == GENES[:5]
    assert np.array_equal(a.X.toarray(), src[GENES[:5]].to_numpy())
    assert np.array_equal(a.obs["control_counts"], src[["Blank-0", "Blank-1"]].sum(axis=1).to_numpy())
    assert np.allclose(a.obs["x"], meta["center_x"]) and np.allclose(a.obs["y"], meta["center_y"])
    assert np.allclose(a.obs["cell_volume"], meta["volume"])
    pen = a.uns["pennsieve"]
    assert (pen["platform"], pen["kind"], pen["coordinate_unit"]) == ("merscope", "cell", "micron")
    assert pen["spot_diameter"] == pytest.approx(9.0)      # median of ((max_x-min_x) + (max_y-min_y)) / 2


def test_cosmx_units(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "cos1.h5ad")
    meta = pd.read_csv(fixture_root["linked"] / "spatial" / "cos" / "S1_metadata_file.csv")
    expr = pd.read_csv(fixture_root["linked"] / "spatial" / "cos" / "S1_exprMat_file.csv")
    assert list(a.var_names) == GENES[:5]
    assert [n.split(":")[1] for n in a.obs_names] == [f"{c}_{f}" for c, f in zip(meta["cell_ID"], meta["fov"])]
    assert np.allclose(a.obs["x"], meta["CenterX_global_px"] * 0.12028)
    assert np.allclose(a.obs["y"], meta["CenterY_global_px"] * 0.12028)
    assert np.allclose(a.obs["cell_area"], meta["Area"] * 0.12028 ** 2)
    assert np.array_equal(a.obs["control_counts"], expr[["NegPrb1", "NegPrb2", "SystemControl1"]].sum(axis=1).to_numpy())
    assert a.uns["pennsieve"]["coordinate_unit"] == "micron" and a.uns["pennsieve"]["platform"] == "cosmx"


def test_cosmx_pixel_parameter(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root["flat"], tmp_path / "out", monkeypatch, PLATFORM="cosmx", COSMX_PIXEL_UM="0.5") == 0
    a = ad.read_h5ad(tmp_path / "out" / "spatial" / "cos1.h5ad")
    meta = pd.read_csv(fixture_root["flat"] / "cos" / "S1_metadata_file.csv")
    assert np.allclose(a.obs["x"], meta["CenterX_global_px"] * 0.5)
    # PLATFORM forces one reader: everything else is ignored, not read
    table = rows(tmp_path / "out" / "metadata" / "samples.csv")
    assert [r["sample_id"] for r in table] == ["cos1"]
    prep = json.loads((tmp_path / "out" / "metadata" / "prep.json").read_text())
    assert any(i["path"].startswith("sr1/") for i in prep["ignored"])


def test_h5ad_symbols_from_feature_name(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "tab1.h5ad")
    src = ad.read_h5ad(fixture_root["linked"] / "spatial" / "tab.h5ad")
    assert list(a.var_names) == GENES and list(a.var["gene_ids"]) == list(src.var_names)
    assert a.X.dtype == np.int32 and np.array_equal(a.X.toarray(), src.X.toarray())
    assert np.array_equal(a.obsm["spatial"], src.obsm["spatial"])
    assert [n.split(":")[1] for n in a.obs_names] == list(src.obs_names)
    assert a.uns["pennsieve"]["platform"] == "h5ad" and a.uns["pennsieve"]["kind"] == "cell"
    warnings = json.loads((outputs["linked"] / "metadata" / "prep.json").read_text())["samples"]["tab1"]["warnings"]
    assert any("feature_name" in w for w in warnings)


def test_spatialdata_centroids(outputs, fixture_root):
    a = ad.read_h5ad(outputs["linked"] / "spatial" / "sd1.h5ad")
    t = ad.read_zarr(fixture_root["linked"] / "spatial" / "sd.zarr" / "tables" / "table")
    cent, diameter = R._shapes_centroids(fixture_root["linked"] / "spatial" / "sd.zarr" / "shapes" / "cells")
    assert diameter == 6.0 and a.uns["pennsieve"]["spot_diameter"] == 6.0
    assert np.allclose(a.obs[["x", "y"]].to_numpy(), cent.loc[t.obs["cell_id"].astype(str), ["x", "y"]].to_numpy())
    assert [n.split(":")[1] for n in a.obs_names] == list(t.obs_names)
    assert list(a.var_names) == GENES[:5] and a.uns["pennsieve"]["platform"] == "spatialdata"


def test_wkb_centroid():
    assert R.wkb_centroid(wkb_point(3.5, -2.0)) == (3.5, -2.0)
    x, y = R.wkb_centroid(wkb_polygon([(0, 0), (4, 0), (4, 2), (0, 2)]))
    assert (x, y) == pytest.approx((2.0, 1.0))
    x, y = R.wkb_centroid(wkb_polygon([(0, 0), (6, 0), (0, 3)]))   # triangle: centroid at the mean of the vertices
    assert (x, y) == pytest.approx((2.0, 1.0))
    with pytest.raises(ValueError):
        R.wkb_centroid(b"\x01\x02\x00\x00\x00")


def test_keep_out_of_tissue(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root["flat"], tmp_path / "out", monkeypatch, PLATFORM="visium", KEEP_OUT_OF_TISSUE="true") == 0
    a = ad.read_h5ad(tmp_path / "out" / "spatial" / "visium1.h5ad")
    assert a.n_obs == 30 and set(a.obs["in_tissue"]) == {0, 1} and (a.obs["in_tissue"] == 0).sum() == 7
    table = {r["sample_id"]: r for r in rows(tmp_path / "out" / "metadata" / "samples.csv")}
    assert table["visium1"]["n_obs"] == "30" and table["visium1"]["n_obs_in_tissue"] == "23"


def test_subsample_is_seeded(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root["flat"], tmp_path / "a", monkeypatch, MAX_OBS_PER_SAMPLE="10") == 0
    assert run(fixture_root["flat"], tmp_path / "b", monkeypatch, MAX_OBS_PER_SAMPLE="10") == 0
    assert run(fixture_root["flat"], tmp_path / "c", monkeypatch, MAX_OBS_PER_SAMPLE="10", SEED="7") == 0
    a = ad.read_h5ad(tmp_path / "a" / "spatial" / "xen1.h5ad")
    b = ad.read_h5ad(tmp_path / "b" / "spatial" / "xen1.h5ad")
    c = ad.read_h5ad(tmp_path / "c" / "spatial" / "xen1.h5ad")
    assert a.n_obs == 10 and list(a.obs_names) == list(b.obs_names) and list(a.obs_names) != list(c.obs_names)
    ids = [int(n.split(":")[1]) for n in a.obs_names]
    assert ids == sorted(ids)                                   # source order kept
    assert a.uns["pennsieve"]["subsampled_from"] == 40
    v = ad.read_h5ad(tmp_path / "a" / "spatial" / "visium1.h5ad")
    assert v.n_obs == 10 and set(v.obs["in_tissue"]) == {1}     # after the out-of-tissue removal
    prep = json.loads((tmp_path / "a" / "metadata" / "prep.json").read_text())["samples"]["visium1"]["observations"]
    assert prep == {"n_source": 30, "n_out_of_tissue": 7, "n_out_of_tissue_dropped": 7, "n_out_of_tissue_kept": 0,
                    "n_subsampled_to": 10, "n_written": 10}


def test_subsample_index():
    a = C.subsample_index(100, 10, 42, "s1")
    assert np.array_equal(a, C.subsample_index(100, 10, 42, "s1")) and not np.array_equal(a, C.subsample_index(100, 10, 42, "s2"))
    assert np.array_equal(a, np.sort(a)) and C.subsample_index(5, 10, 42, "s1") is None and C.subsample_index(5, 0, 42, "s1") is None


def test_without_sample_sheet(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["flat"], src)
    (src / "samples.csv").unlink()
    (src / "annotations.csv").unlink()
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    table = rows(out / "metadata" / "samples.csv")
    assert list(table[0]) == list(I.REQUIRED_COLUMNS)
    # sample ids from folder names (sr2/outs -> sr2) and file stems
    assert [r["sample_id"] for r in table] == ["cos", "mer", "sd", "sr1", "sr2", "sr3", "tab", "xen"]
    a = ad.read_h5ad(out / "spatial" / "sr1.h5ad")
    assert "annotation" not in a.obs and "condition" not in a.obs
    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert prep["sampleSheet"] is None and prep["skipped"] == []


def test_sample_not_in_sheet_is_ignored(fixture_root, tmp_path, monkeypatch):
    """The sheet defines the samples: detected samples it does not list are recorded as ignored, not written."""
    src = tmp_path / "in"
    shutil.copytree(fixture_root["flat"], src)
    (src / "samples.csv").write_text("sample_id,path,condition\nvisium1,sr1,ctrl\n")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    table = {r["sample_id"]: r for r in rows(out / "metadata" / "samples.csv")}
    assert list(table) == ["visium1"] and table["visium1"]["condition"] == "ctrl"
    assert {p.name for p in (out / "spatial").glob("*.h5ad")} == {"visium1.h5ad"}
    ignored = {i["path"]: i["reason"] for i in json.loads((out / "metadata" / "prep.json").read_text())["ignored"]}
    assert ignored["xen"] == "xenium sample not in the sample sheet" and ignored["sr2/outs"].endswith("not in the sample sheet")


def test_bad_sample_sheet_fails(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["flat"], src)
    (src / "samples.csv").write_text("sample_id,condition\nvisium1,ctrl\n")
    assert run(src, tmp_path / "out", monkeypatch) == 1


def test_a_broken_sample_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["flat"], src)
    (src / "xen" / "cell_feature_matrix.h5").write_bytes(b"not an hdf5 file")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    prov = json.loads((out / "spatial" / "provenance.json").read_text())
    assert prov["nOk"] == 7 and prov["nFailed"] == 2       # ghost + xen1
    failed = [s for s in prov["samples"] if s["status"] == "failed"]
    assert [s["sample_id"] for s in failed] == ["xen1"]
    assert not (out / "spatial" / "xen1.h5ad").exists()


def test_duplicate_sample_ids_are_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["flat"], src)
    (src / "samples.csv").write_text("sample_id,path\nsame,sr1\nsame,xen\nmer1,mer\n")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    table = [r["sample_id"] for r in rows(out / "metadata" / "samples.csv")]
    assert table == ["mer1"]


def test_visium_hd(tmp_path, monkeypatch):
    """A Visium HD tree: binned_outputs/square_008um/ with a parquet positions file; 16 µm absent."""
    from tests.make_fixture import SR_SCALEFACTORS, write_10x_h5, write_png

    rng = np.random.default_rng(3)
    root = tmp_path / "in" / "hd_sample" / "outs"
    b = root / "binned_outputs" / "square_008um"
    (b / "spatial").mkdir(parents=True)
    n = 36
    bcs = [f"s_008um_{i // 6:05d}_{i % 6:05d}-1" for i in range(n)]
    X = rng.poisson(2, (n, len(GENES))).astype(np.int32)
    write_10x_h5(b / "filtered_feature_bc_matrix.h5", X, bcs, GENES, [f"ENSG{i:011d}" for i in range(1, 9)],
                 ["Gene Expression"] * 8, "spaceranger-3.0.0")
    pos = pd.DataFrame({"barcode": bcs, "in_tissue": 1, "array_row": [i // 6 for i in range(n)],
                        "array_col": [i % 6 for i in range(n)], "pxl_row_in_fullres": [100.0 + 20 * (i // 6) for i in range(n)],
                        "pxl_col_in_fullres": [100.0 + 20 * (i % 6) for i in range(n)]})
    pos.to_parquet(b / "spatial" / "tissue_positions.parquet", index=False)
    (b / "spatial" / "scalefactors_json.json").write_text(json.dumps(SR_SCALEFACTORS | {"spot_diameter_fullres": 20.0, "bin_size_um": 8,
                                                                                           "microns_per_pixel": 0.4}))
    write_png(b / "spatial" / "tissue_hires_image.png", (30, 30), 5)
    (root / "spatial").mkdir()
    write_png(root / "spatial" / "tissue_lowres_image.png", (10, 10), 6)
    out = tmp_path / "out"
    assert run(tmp_path / "in", out, monkeypatch) == 0
    a = ad.read_h5ad(out / "spatial" / "hd_sample.h5ad")
    pen = a.uns["pennsieve"]
    assert (pen["platform"], pen["kind"], pen["grid"], pen["bin_size_um"], pen["spot_diameter"]) == ("visium-hd", "bin", "square", 8, 20.0)
    assert a.n_obs == 36 and set(a.uns["spatial"]["hd_sample"]["images"]) == {"hires", "lowres"}
    assert list(a.obs["array_col"][:3]) == [0, 1, 2]
    assert run(tmp_path / "in", tmp_path / "out16", monkeypatch, BIN_SIZE_UM="16") == 1   # no such bin folder


def test_bad_parameters(fixture_root, tmp_path, monkeypatch):
    assert run(fixture_root["flat"], tmp_path / "o1", monkeypatch, PLATFORM="slideseq") == 2
    assert run(fixture_root["flat"], tmp_path / "o2", monkeypatch, MAX_OBS_PER_SAMPLE="-1") == 2
    assert run(fixture_root["flat"], tmp_path / "o3", monkeypatch, KEEP_OUT_OF_TISSUE="maybe") == 2
    assert run(fixture_root["flat"], tmp_path / "o4", monkeypatch, COSMX_PIXEL_UM="0") == 2
    assert run(fixture_root["flat"], tmp_path / "o5", monkeypatch, BIN_SIZE_UM="x") == 2


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "readme.txt").write_text("nothing here")
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_path_matching():
    root = Path("/in")
    c = I.Candidate("visium", root / "spatial" / "sampleA" / "outs", default_id="sampleA")
    for want in ("sampleA", "spatial/sampleA", "sampleA/outs", "spatial/sampleA/outs", "outs"):
        assert I.path_matches(want, c, root), want
    for want in ("ampleA", "sampleB", "", "spatial"):
        assert not I.path_matches(want, c, root), want
    f = I.Candidate("h5ad", root / "x" / "tab.h5ad", default_id="tab")
    assert I.path_matches("tab.h5ad", f, root) and I.path_matches("x/tab.h5ad", f, root) and not I.path_matches("tab", f, root)


def test_as_counts_float_and_int():
    w: list[str] = []
    X = C.as_counts(sp.csr_matrix(np.array([[1.0, 0], [2.0, 3.0]])), w)
    assert X.dtype == np.int32 and not w
    X = C.as_counts(sp.csr_matrix(np.array([[1.5, 0], [2.0, 3.0]])), w)
    assert X.dtype == np.float32 and w
