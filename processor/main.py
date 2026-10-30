"""
st-prep: the first step of the spatial transcriptomics workflows.

Finds every spatial transcriptomics sample under INPUT_DIR (Visium, Visium HD, Xenium,
MERSCOPE, CosMx, .h5ad, SpatialData .zarr), reads the optional sample sheet and
annotations, and writes

  OUTPUT_DIR/spatial/<sample_id>.h5ad     raw counts, coordinates, QC metrics, images (Visium), uns["pennsieve"]
  OUTPUT_DIR/spatial/provenance.json
  OUTPUT_DIR/metadata/samples.csv         sample_id, path, platform, kind, n_obs, n_genes, n_obs_in_tissue, sheet columns
  OUTPUT_DIR/metadata/prep.json           per sample: platform, files read, features and observations dropped, warnings
  OUTPUT_DIR/metadata/provenance.json

One sample at a time. Nothing is filtered or normalised except out-of-tissue Visium spots
and control features (summed into obs["control_counts"] first).
"""
from __future__ import annotations

from processor import threads  # noqa: F401  (sets the thread environment before numpy is imported)

import csv  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import os  # noqa: E402
import platform as _platform  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from collections import Counter  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from importlib import metadata as importlib_metadata  # noqa: E402
from pathlib import Path  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-st-prep")
log.setLevel(logging.INFO)
for noisy in ("numba", "matplotlib", "PIL", "h5py", "anndata", "scanpy"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-prep"
OUTPUTS = ("spatial", "metadata")
REPO_ROOT = Path(__file__).resolve().parent.parent
PLATFORM_CHOICES = ("auto", "visium", "visium-hd", "xenium", "merscope", "cosmx", "h5ad", "spatialdata")

DEFAULTS = {
    "SAMPLE_SHEET": "samples.csv",
    "ANNOTATIONS": "annotations.csv",
    "PLATFORM": "auto",
    "BIN_SIZE_UM": "8",
    "KEEP_OUT_OF_TISSUE": "false",
    "MAX_OBS_PER_SAMPLE": "0",
    "COSMX_PIXEL_UM": "0.12028",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def _bool(s: str, name: str) -> bool:
    v = s.strip().lower()
    if v in ("true", "1", "yes", "y"):
        return True
    if v in ("false", "0", "no", "n", ""):
        return False
    raise ParameterError(f"{name} must be true or false, got {s!r}")


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    p["SAMPLE_SHEET"] = raw["SAMPLE_SHEET"].strip()
    p["ANNOTATIONS"] = raw["ANNOTATIONS"].strip()
    p["PLATFORM"] = raw["PLATFORM"].strip().lower().replace("_", "-")
    if p["PLATFORM"] not in PLATFORM_CHOICES:
        raise ParameterError(f"PLATFORM must be one of {', '.join(PLATFORM_CHOICES)}, got {raw['PLATFORM']!r}")
    p["KEEP_OUT_OF_TISSUE"] = _bool(raw["KEEP_OUT_OF_TISSUE"], "KEEP_OUT_OF_TISSUE")
    try:
        p["BIN_SIZE_UM"] = int(raw["BIN_SIZE_UM"])
        p["MAX_OBS_PER_SAMPLE"] = int(raw["MAX_OBS_PER_SAMPLE"])
        p["SEED"] = int(raw["SEED"])
        p["THREADS"] = int(raw["THREADS"])
        p["COSMX_PIXEL_UM"] = float(raw["COSMX_PIXEL_UM"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["BIN_SIZE_UM"] <= 0 or p["MAX_OBS_PER_SAMPLE"] < 0 or p["THREADS"] < 0 or p["COSMX_PIXEL_UM"] <= 0:
        raise ParameterError("BIN_SIZE_UM and COSMX_PIXEL_UM must be > 0; MAX_OBS_PER_SAMPLE and THREADS >= 0")
    return p


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": _platform.python_version()}
    for pkg in ("scanpy", "anndata", "numpy", "scipy", "pandas", "h5py", "pyarrow", "pillow", "zarr"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)   # bytes on macOS, KiB on Linux
    except (ImportError, OSError):
        return None


def write_csv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=_json_default)
        fh.write("\n")


def _json_default(o):
    import numpy as np

    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    n_threads = threads.configure(str(params["THREADS"]))
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")

    from processor import contract as C
    from processor import inputs as I
    from processor import readers as R

    app_version = version_of(REPO_ROOT)
    rel = lambda p: I._rel(p, input_dir)  # noqa: E731

    # 1. what is there
    found = I.walk(input_dir, params["SAMPLE_SHEET"], params["ANNOTATIONS"], params["PLATFORM"])
    log.info("input: %d files, %d sample(s) detected (%s), sample sheet: %s, annotations: %s, threads %d",
             found.n_files, len(found.candidates),
             ", ".join(f"{k} x{v}" for k, v in sorted(Counter(c.platform for c in found.candidates).items())) or "none",
             rel(found.sheet) if found.sheet else "none", rel(found.annotations) if found.annotations else "none",
             n_threads)
    if not found.candidates:
        log.error("no spatial transcriptomics sample found under INPUT_DIR (%d files; %d ignored, e.g. %s)",
                  found.n_files, found.n_ignored,
                  "; ".join(f"{i['path']}: {i['reason']}" for i in found.ignored[:10]) or "nothing")
        return 1
    try:
        sheet = I.read_sheet(found.sheet)
    except (ValueError, OSError, UnicodeDecodeError) as exc:
        log.error("cannot use the sample sheet %s: %s", found.sheet, exc)
        return 1
    try:
        annotations = I.read_annotations(found.annotations)
    except (ValueError, OSError, UnicodeDecodeError) as exc:
        log.error("cannot use the annotations file %s: %s", found.annotations, exc)
        return 1
    sheet_errors = I.assign_identities(found.candidates, sheet, input_dir)
    for e in sheet_errors:
        log.error("sample sheet row %s (path %s): %s", e["sample_id"] or "?", e["path"], e["reason"])
    if sheet.path:
        # The sample sheet defines the samples: anything detected but not listed is left out (e.g. a
        # folder of preprocessed .h5ad copies next to the Space Ranger outputs), and recorded as ignored.
        for c in found.candidates:
            if c.sheet_row is None:
                log.warning("%s: %s sample not in %s; ignored", c.rel(input_dir), c.platform, sheet.path.name)
                found.ignored.append({"path": c.rel(input_dir), "reason": f"{c.platform} sample not in the sample sheet"})
        found.ignored.sort(key=lambda i: i["path"])
        found.n_ignored = len(found.ignored)
        found.candidates = [c for c in found.candidates if c.sheet_row is not None]
    ids = Counter(c.sample_id for c in found.candidates)
    dups = {s for s, k in ids.items() if k > 1}
    skipped: list[dict] = [{"sample_id": e["sample_id"], "path": e["path"], "reason": e["reason"]} for e in sheet_errors]
    if dups:
        for c in found.candidates:
            if c.sample_id in dups:
                log.error("%s: sample_id %s is used by %d samples; skipped (give each its own sample_id)",
                          c.rel(input_dir), c.sample_id, ids[c.sample_id])
                skipped.append({"sample_id": c.sample_id, "path": c.rel(input_dir), "reason": "duplicate sample_id"})
        found.candidates = [c for c in found.candidates if c.sample_id not in dups]

    spatial_dir, meta_dir = output_dir / "spatial", output_dir / "metadata"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    # 2. one sample at a time
    sample_rows: list[dict] = []
    prep_samples: dict[str, dict] = {}
    prov_samples: list[dict] = []
    written_ids: set[str] = set()
    for c in found.candidates:
        t0 = time.monotonic()
        src = c.rel(input_dir)
        for note in c.notes:
            log.warning("%s: %s", src, note)
        try:
            reads = R.READERS[c.platform](c, input_dir, params)
        except Exception as exc:  # noqa: BLE001 - one sample's failure must not sink the run
            log.error("%s (%s): FAILED: %s\n%s", c.sample_id, src, exc, traceback.format_exc())
            skipped.append({"sample_id": c.sample_id, "path": src, "reason": f"{type(exc).__name__}: {exc}"})
            prov_samples.append({"sample_id": c.sample_id, "input": src, "platform": c.platform, "status": "failed",
                                 "error": str(exc), "seconds": round(time.monotonic() - t0, 2)})
            continue
        if len(reads) > 1 and c.sheet_row is not None and c.sheet_row.get("sample_id"):
            log.warning("%s: holds %d samples; the sheet's sample_id %s is not used, obs['sample_id'] values are",
                        src, len(reads), c.sample_id)
        for r in reads:
            t1 = time.monotonic()
            sid = c.sample_id if r.sub_id is None else r.sub_id
            if sid in written_ids:
                log.error("%s: sample_id %s already written; skipped", src, sid)
                skipped.append({"sample_id": sid, "path": src, "reason": "duplicate sample_id"})
                continue
            try:
                adata, info = C.finalise(r, sid, src, params, c.extra, annotations.get(sid), app_version)
                out = spatial_dir / f"{sid}.h5ad"
                adata.write_h5ad(out, compression="gzip")
            except Exception as exc:  # noqa: BLE001
                log.error("%s (%s): FAILED: %s\n%s", sid, src, exc, traceback.format_exc())
                skipped.append({"sample_id": sid, "path": src, "reason": f"{type(exc).__name__}: {exc}"})
                prov_samples.append({"sample_id": sid, "input": src, "platform": c.platform, "status": "failed",
                                     "error": str(exc), "seconds": round(time.monotonic() - t1, 2)})
                continue
            written_ids.add(sid)
            seconds = round(time.monotonic() - t1, 2)
            info["source"] = src
            info["in_sample_sheet"] = (c.sheet_row is not None) if sheet.path else None
            info["notes"] = list(c.notes)
            info["seconds"] = seconds
            info["output"] = f"spatial/{sid}.h5ad"
            prep_samples[sid] = info
            for w in info["warnings"]:
                log.warning("%s: %s", sid, w)
            sample_rows.append({"sample_id": sid, "path": src, "platform": info["platform"], "kind": info["kind"],
                                "n_obs": info["n_obs"], "n_genes": info["n_genes"],
                                "n_obs_in_tissue": info["n_obs_in_tissue"], **c.extra})
            prov_samples.append({"sample_id": sid, "input": src, "platform": info["platform"], "status": "ok",
                                 "output": info["output"], "n_obs": info["n_obs"], "n_genes": info["n_genes"],
                                 "seconds": seconds})
            o = info["observations"]
            ctrl = {t: d["n_features"] for t, d in info["features"]["dropped"].items()}
            log.info("%s: %s %s -> %d obs x %d genes (%d source obs, %d out of tissue dropped%s; %d of %d features "
                     "dropped%s), %s%.1fs",
                     sid, info["platform"], src, info["n_obs"], info["n_genes"], o["n_source"],
                     o["n_out_of_tissue_dropped"],
                     f", subsampled to {o['n_subsampled_to']}" if o["n_subsampled_to"] else "",
                     info["features"]["n_source"] - info["n_genes"], info["features"]["n_source"],
                     (": " + ", ".join(f"{t} {n}" for t, n in ctrl.items())) if ctrl else "",
                     f"annotated {info['annotation']['n_annotated']}, " if info["annotation"] else "", seconds)
            del adata

    # 3. metadata
    write_csv(meta_dir / "samples.csv", list(I.REQUIRED_COLUMNS) + sheet.columns, sample_rows)
    finished = datetime.now(timezone.utc)
    prep_doc = {
        "app": APP_NAME, "appVersion": app_version, "generatedAt": finished.isoformat(timespec="seconds"),
        "parameters": params,
        "sampleSheet": rel(found.sheet) if found.sheet else None,
        "sampleSheetColumns": sheet.columns,
        "annotations": rel(found.annotations) if found.annotations else None,
        "samples": prep_samples, "skipped": skipped,
        "ignored": found.ignored, "nIgnored": found.n_ignored,
    }
    write_json(meta_dir / "prep.json", prep_doc)
    provenance = {
        "app": APP_NAME, "appVersion": app_version, "tools": tool_versions(),
        "startedAt": started.isoformat(timespec="seconds"), "finishedAt": finished.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2), "threads": n_threads,
        "parameters": params,
        "inputsFound": {"samples": [{"path": c.rel(input_dir), "platform": c.platform} for c in found.candidates],
                        "sampleSheet": rel(found.sheet) if found.sheet else None,
                        "annotations": rel(found.annotations) if found.annotations else None,
                        "nFiles": found.n_files, "nIgnored": found.n_ignored},
        "samples": prov_samples, "skipped": skipped, "nOk": len(sample_rows), "nFailed": len(skipped),
        "peakRssMb": peak_rss_mb(),
    }
    for d in (spatial_dir, meta_dir):
        write_json(d / "provenance.json", provenance)
    log.info("%d sample(s) written, %d skipped, peak RSS %s MB, %.1fs; wrote %s and %s", len(sample_rows), len(skipped),
             provenance["peakRssMb"], time.monotonic() - t_start, spatial_dir, meta_dir)
    if not sample_rows:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
