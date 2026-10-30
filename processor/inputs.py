"""
Finding the samples under INPUT_DIR: platform detection by the files present in
each folder (or per file for .h5ad / .zarr), the optional sample sheet and the
optional annotations file.
"""
from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("processor-st-prep")

PLATFORMS = ("visium", "visium-hd", "xenium", "merscope", "cosmx", "h5ad", "spatialdata")
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
REQUIRED_COLUMNS = ("sample_id", "path", "platform", "kind", "n_obs", "n_genes", "n_obs_in_tissue")
# Folder names Space Ranger / Xenium Ranger use for the output folder itself; a sample's default id
# comes from the folder above them.
OUTS_NAMES = {"outs", "binned_outputs", "output", "outputs"}
BIN_RE = re.compile(r"^square_(\d{3})um$")
MAX_IGNORED = 500


def sanitise(stem: str) -> str:
    """A sample_id from a folder or file name: anything outside [A-Za-z0-9_-] becomes '_'."""
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")
    s = re.sub(r"_+", "_", s)
    return s or "sample"


@dataclass
class Candidate:
    """One detected sample: a folder (or file) and the files the reader needs, by role."""
    platform: str
    source: Path                       # the sample's folder, or its file
    files: dict[str, Path] = field(default_factory=dict)
    default_id: str = ""
    notes: list[str] = field(default_factory=list)
    # set when identities are assigned
    sample_id: str = ""
    sheet_row: dict[str, str] | None = None
    extra: dict[str, str] = field(default_factory=dict)

    def rel(self, root: Path) -> str:
        return _rel(self.source, root)


@dataclass
class Found:
    candidates: list[Candidate] = field(default_factory=list)
    ignored: list[dict] = field(default_factory=list)      # {"path": ..., "reason": ...}
    n_ignored: int = 0
    sheet: Path | None = None
    annotations: Path | None = None
    n_files: int = 0


def _rel(path: Path, root: Path) -> str:
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def default_id_for(folder: Path, root: Path) -> str:
    """The sample id a folder gets without a sample sheet: its name, or its parent's for outs/ folders."""
    p = folder
    while p.name in OUTS_NAMES and p != root and p.parent != p:
        p = p.parent
    if p == root or not p.name:
        return "sample"
    return sanitise(p.name)


def _first(files: dict[str, Path], *names: str) -> Path | None:
    for n in names:
        if n in files:
            return files[n]
    return None


def _suffix_matches(files: dict[str, Path], suffix: str) -> dict[str, Path]:
    """Files whose name ends with suffix, keyed by the prefix before it."""
    return {name[: -len(suffix)]: p for name, p in files.items() if name.endswith(suffix)}


def _pairs(files: dict[str, Path], a_suffixes: tuple[str, ...], b_suffixes: tuple[str, ...]) -> list[tuple[str, Path, Path]]:
    """(prefix, a-file, b-file) pairs of two file kinds that share a prefix (MERSCOPE, CosMx exports)."""
    a: dict[str, Path] = {}
    for s in a_suffixes:
        a.update(_suffix_matches(files, s))
    b: dict[str, Path] = {}
    for s in b_suffixes:
        b.update(_suffix_matches(files, s))
    out = []
    for prefix in sorted(a):
        if prefix in b:
            out.append((prefix, a[prefix], b[prefix]))
        elif len(a) == 1 and len(b) == 1:   # one of each with different prefixes
            out.append((prefix, a[prefix], next(iter(b.values()))))
    return out


def _is_mtx_dir(d: Path) -> bool:
    return d.is_dir() and any((d / n).is_file() for n in ("matrix.mtx.gz", "matrix.mtx"))


def _has_obsm_spatial(path: Path) -> bool:
    import h5py

    try:
        with h5py.File(path, "r") as h:
            return "obsm" in h and "spatial" in h["obsm"]
    except OSError:
        return False


def _zarr_tables(path: Path) -> Path | None:
    for name in ("tables", "table"):
        d = path / name
        if d.is_dir() and any(c.is_dir() and not c.name.startswith(".") for c in d.iterdir()):
            return d
    return None


def detect_folder(folder: Path, dirnames: list[str], filenames: list[str], root: Path,
                  platform: str = "auto") -> tuple[list[Candidate], bool, list[dict]]:
    """
    Candidates in one folder, whether the folder is claimed as a whole (no descent), and files
    ignored with a reason. `platform` other than auto runs only that platform's detector.
    """
    files = {n: folder / n for n in filenames}
    want = (lambda p: platform in ("auto", p))
    cands: list[Candidate] = []
    ignored: list[dict] = []
    did = default_id_for(folder, root)

    # Visium HD: binned_outputs/square_<NNN>um/{filtered_feature_bc_matrix.h5, spatial/tissue_positions.parquet}
    if want("visium-hd") and "binned_outputs" in dirnames:
        binned = folder / "binned_outputs"
        bins = sorted(d.name for d in binned.iterdir() if d.is_dir() and BIN_RE.match(d.name))
        if bins:
            c = Candidate("visium-hd", folder, {"binned_outputs": binned}, did)
            c.files.update({b: binned / b for b in bins})
            return [c], True, ignored

    # Visium: *filtered_feature_bc_matrix.h5 + spatial/tissue_positions(_list).csv
    if want("visium"):
        mats = sorted(n for n in filenames if n.endswith("filtered_feature_bc_matrix.h5"))
        raws = sorted(n for n in filenames if n.endswith("raw_feature_bc_matrix.h5"))
        # Space Ranger's MTX folders (filtered_feature_bc_matrix/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz})
        mtx = {d: folder / d for d in dirnames if d.endswith("feature_bc_matrix") and _is_mtx_dir(folder / d)}
        files.update(mtx)
        mats += sorted(d for d in mtx if d.endswith("filtered_feature_bc_matrix"))
        raws += sorted(d for d in mtx if d.endswith("raw_feature_bc_matrix"))
        pos = None
        for d in (folder / "spatial", folder):
            for n in ("tissue_positions.csv", "tissue_positions_list.csv", "tissue_positions.parquet"):
                if (d / n).is_file():
                    pos = d / n
                    break
            if pos:
                break
        if (mats or raws) and pos:
            c = Candidate("visium", folder, {"positions": pos}, did)
            if mats:
                c.files["matrix"] = files[mats[0]]
            else:
                c.files["matrix"] = files[raws[0]]
                c.notes.append(f"no filtered_feature_bc_matrix; read {raws[0]} (all barcodes, in-tissue kept)")
            sp = pos.parent
            for role, name in (("scalefactors", "scalefactors_json.json"), ("hires", "tissue_hires_image.png"),
                               ("lowres", "tissue_lowres_image.png")):
                for d in (sp, folder / "spatial", folder):
                    if (d / name).is_file():
                        c.files[role] = d / name
                        break
            return [c], True, ignored
        if mats or raws:
            ignored.append({"path": _rel(files[(mats or raws)[0]], root),
                            "reason": "10x matrix without a spatial/tissue_positions file"})
            dirnames[:] = [d for d in dirnames if d not in mtx]

    # Xenium: cell_feature_matrix.h5 + cells.parquet | cells.csv.gz | cells.csv
    if want("xenium") and "cell_feature_matrix.h5" in files:
        cells = _first(files, "cells.parquet", "cells.csv.gz", "cells.csv")
        if cells is not None:
            return [Candidate("xenium", folder, {"matrix": files["cell_feature_matrix.h5"], "cells": cells}, did)], True, ignored
        ignored.append({"path": _rel(files["cell_feature_matrix.h5"], root), "reason": "Xenium matrix without cells.parquet / cells.csv.gz"})

    # MERSCOPE: *cell_by_gene.csv + *cell_metadata.csv
    if want("merscope"):
        pairs = _pairs(files, ("cell_by_gene.csv", "cell_by_gene.csv.gz"), ("cell_metadata.csv", "cell_metadata.csv.gz"))
        if pairs:
            out = []
            for prefix, a, b in pairs:
                sid = did if len(pairs) == 1 else sanitise(f"{did}_{prefix.rstrip('_')}")
                out.append(Candidate("merscope", folder if len(pairs) == 1 else a, {"counts": a, "metadata": b}, sid))
            return out, True, ignored

    # CosMx: *_exprMat_file.csv + *_metadata_file.csv
    if want("cosmx"):
        pairs = _pairs(files, ("_exprMat_file.csv", "_exprMat_file.csv.gz"), ("_metadata_file.csv", "_metadata_file.csv.gz"))
        if pairs:
            out = []
            for prefix, a, b in pairs:
                sid = did if len(pairs) == 1 else sanitise(f"{did}_{prefix}")
                out.append(Candidate("cosmx", folder if len(pairs) == 1 else a, {"counts": a, "metadata": b}, sid))
            return out, True, ignored

    # AnnData: *.h5ad with obsm["spatial"], one sample per file
    if want("h5ad"):
        for n in sorted(filenames):
            if n.lower().endswith(".h5ad"):
                if _has_obsm_spatial(files[n]):
                    cands.append(Candidate("h5ad", files[n], {"h5ad": files[n]}, sanitise(n[:-5])))
                else:
                    ignored.append({"path": _rel(files[n], root), "reason": "h5ad without obsm['spatial']"})
    return cands, False, ignored


def walk(root: Path, sheet_name: str, annotations_name: str, platform: str = "auto") -> Found:
    """Every sample under root (symlinks followed), the sample sheet and annotations file (shallowest wins)."""
    found = Found()
    claimed: list[Path] = []
    seen: set[str] = set()
    sheets: list[Path] = []
    annots: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        folder = Path(dirpath)
        real = os.path.realpath(dirpath)
        if real in seen:          # symlink loop
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        filenames = sorted(f for f in filenames if not f.startswith(".") and (folder / f).is_file())
        found.n_files += len(filenames)
        # the sample sheet and annotations are looked for anywhere (but not inside a sample folder)
        for f in filenames:
            if sheet_name and f == sheet_name:
                sheets.append(folder / f)
            if annotations_name and f == annotations_name:
                annots.append(folder / f)
        # SpatialData stores: a *.zarr folder with a tables element, one sample each; never descended
        zarrs = [d for d in dirnames if d.lower().endswith(".zarr")]
        for d in zarrs:
            z = folder / d
            tables = _zarr_tables(z)
            if platform in ("auto", "spatialdata") and tables is not None:
                found.candidates.append(Candidate("spatialdata", z, {"zarr": z, "tables": tables}, sanitise(d[:-5])))
            else:
                found.ignored.append({"path": _rel(z, root) + "/", "reason": "zarr store without a tables element"
                                      if tables is None else f"PLATFORM={platform}"})
        dirnames[:] = [d for d in dirnames if d not in zarrs]
        cands, whole, ignored = detect_folder(folder, dirnames, filenames, root, platform)
        found.candidates.extend(cands)
        found.ignored.extend(ignored)
        if whole:
            claimed.append(folder)
            dirnames[:] = []
            continue
        claimed_files = {p for c in cands for p in c.files.values()}
        noted = {Path(root, i["path"]) for i in ignored}
        for f in filenames:
            p = folder / f
            if p in claimed_files or p in noted or f == sheet_name or f == annotations_name:
                continue
            found.ignored.append({"path": _rel(p, root), "reason": "not part of a recognised sample"})
    for lst, attr in ((sheets, "sheet"), (annots, "annotations")):
        lst[:] = [p for p in lst if not any(p.is_relative_to(c) for c in claimed)] or lst
        if lst:
            lst.sort(key=lambda p: (len(p.relative_to(root).parts), str(p)))
            setattr(found, attr, lst[0])
            if len(lst) > 1:
                log.warning("%d files named %s; using %s", len(lst), lst[0].name, _rel(lst[0], root))
    found.candidates.sort(key=lambda c: c.rel(root))
    found.ignored.sort(key=lambda i: i["path"])
    found.n_ignored = len(found.ignored)
    if len(found.ignored) > MAX_IGNORED:
        found.ignored = found.ignored[:MAX_IGNORED]
    return found


# --- sample sheet ---------------------------------------------------------------------------

@dataclass
class SampleSheet:
    path: Path | None = None
    columns: list[str] = field(default_factory=list)     # extra columns, carried to samples.csv and obs
    rows: list[dict[str, str]] = field(default_factory=list)


def read_sheet(path: Path | None) -> SampleSheet:
    """The sample sheet: needs sample_id and path columns. Raises ValueError when it cannot be used."""
    if path is None:
        return SampleSheet()
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        header = [h.strip() for h in (reader.fieldnames or [])]
        for col in ("sample_id", "path"):
            if col not in header:
                raise ValueError(f"{path.name} has no '{col}' column (columns: {', '.join(header)})")
        rows = []
        for raw in reader:
            row = {(k or "").strip(): (v or "").strip() for k, v in raw.items() if k is not None}
            if not row.get("path") and not row.get("sample_id"):
                continue
            rows.append(row)
    extra = [h for h in header if h not in REQUIRED_COLUMNS and h]
    return SampleSheet(path=path, columns=extra, rows=rows)


def _norm(p: str) -> str:
    return p.replace("\\", "/").strip().strip("/")


def path_matches(sheet_path: str, candidate: Candidate, root: Path) -> bool:
    """A sheet `path` names a sample when it equals, or ends, the sample's path relative to the input,
    with or without a trailing outs/ folder, or its folder/file name."""
    want = _norm(sheet_path)
    if not want:
        return False
    rel = _norm(candidate.rel(root))
    forms = {rel}
    p = Path(rel)
    while p.name in OUTS_NAMES and p.parent != p and str(p.parent) != ".":
        p = p.parent
        forms.add(_norm(str(p)))
    for f in list(forms):
        forms.add(Path(f).name)
    return any(f == want or f.endswith("/" + want) for f in forms)


def assign_identities(cands: list[Candidate], sheet: SampleSheet, root: Path) -> list[dict]:
    """
    Give each candidate its sample_id and sheet columns. Returns the sheet rows that matched
    nothing (or several samples) as error records.
    """
    errors: list[dict] = []
    for c in cands:
        c.sample_id = c.default_id
        c.extra = {col: "" for col in sheet.columns}
    for row in sheet.rows:
        hits = [c for c in cands if path_matches(row.get("path", ""), c, root)]
        if not hits:
            errors.append({"sample_id": row.get("sample_id", ""), "path": row.get("path", ""),
                           "reason": f"path {row.get('path', '')!r} matches no sample under the input"})
            continue
        if len(hits) > 1:
            errors.append({"sample_id": row.get("sample_id", ""), "path": row.get("path", ""),
                           "reason": f"path {row.get('path', '')!r} matches {len(hits)} samples: "
                                     + ", ".join(h.rel(root) for h in hits)})
            continue
        c = hits[0]
        if c.sheet_row is not None:
            errors.append({"sample_id": row.get("sample_id", ""), "path": row.get("path", ""),
                           "reason": f"path {row.get('path', '')!r} matches the same sample as an earlier row"})
            continue
        c.sheet_row = row
        sid = row.get("sample_id") or c.default_id
        if not SAMPLE_ID_RE.match(sid):
            fixed = sanitise(sid)
            c.notes.append(f"sample_id {sid!r} sanitised to {fixed!r}")
            sid = fixed
        c.sample_id = sid
        c.extra = {col: row.get(col, "") for col in sheet.columns}
    return errors


# --- annotations ------------------------------------------------------------------------------

def read_annotations(path: Path | None) -> dict[str, dict[str, str]]:
    """annotations.csv (sample_id, obs_id, annotation) as {sample_id: {obs_id: annotation}}."""
    if path is None:
        return {}
    import pandas as pd

    df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    for col in ("sample_id", "obs_id", "annotation"):
        if col not in df.columns:
            raise ValueError(f"{path.name} has no '{col}' column (columns: {', '.join(df.columns)})")
    out: dict[str, dict[str, str]] = {}
    for sid, grp in df.groupby("sample_id", sort=False):
        out[str(sid)] = dict(zip(grp["obs_id"].astype(str).str.strip(), grp["annotation"].astype(str).str.strip()))
    return out
