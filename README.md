# st-prep

The first step of the three spatial transcriptomics workflows (Spatial Domain Mapping &
Condition Comparison; Cell-type Deconvolution & Niches; Single-cell Spatial Profiling). It
takes spatial transcriptomics data as the instrument software writes it, plus an optional
sample sheet and per-spot annotations, and turns every sample into one AnnData file of the
same shape, so `st-qc` and every step after it can read Visium spots, Visium HD bins and
Xenium / MERSCOPE / CosMx cells without knowing where they came from.

**For a biologist:** point it at the folder with your Space Ranger outputs (or Xenium,
MERSCOPE or CosMx exports) and a `samples.csv` saying which folder is which sample and
which condition it belongs to. It finds each sample, keeps the raw counts, the position of
every spot or cell, the H&E image for Visium, and adds a few quality numbers per spot/cell
(total counts, genes detected, % mitochondrial reads, control-probe counts). It removes only
what is not biology: Visium spots outside the tissue and the control probes of imaging panels
(their counts are kept per cell as `control_counts`, so `st-qc` can flag noisy cells). It does
not filter low-quality spots or normalise; that is `st-qc`'s job.

**For a computational analyst:** one `<sample_id>.h5ad` per sample following the AnnData
contract of the workflows (`docs/conventions.md` in `workflow-spatial-transcriptomics`): raw
integer counts as CSR int32 in `X`; `obs_names` `<sample_id>:<barcode|cell_id>`;
`obsm["spatial"]` in full-resolution pixels (Visium, Visium HD; scanpy/squidpy convention so
`uns["spatial"]` scale factors apply) or microns (imaging platforms); `uns["pennsieve"]` with
`platform`, `kind` (`spot`/`bin`/`cell`), `coordinate_unit`, `spot_diameter`, `grid`
(`hex`/`square`/`none`); QC metrics from `scanpy.pp.calculate_qc_metrics`. Samples are read and
written one at a time; matrices are never densified (dense CSV inputs are converted in chunks
of 20 000 rows).

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `spatial` | input | root | The data folder tree (searched recursively, symlinks followed), with the optional `samples.csv` and `annotations.csv` anywhere in it. |
| `spatial` | output | `OUTPUT_DIR/spatial/` | `<sample_id>.h5ad` per sample, `provenance.json`. |
| `metadata` | output | `OUTPUT_DIR/metadata/` | `samples.csv`, `prep.json`, `provenance.json`. |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | File name of the sample sheet, looked for anywhere in the input (the shallowest wins). Needs `sample_id` and `path` (a folder or file name, or a path relative to the input; it matches when it equals or ends the sample's path, with or without a trailing `outs/`). Every other column (`condition`, `batch`, `donor`, `tissue`, …) goes to `samples.csv` and onto every observation as a categorical `obs` column. **With a sheet, the sheet defines the samples**: anything detected but not listed is not processed and is listed in `prep.json → ignored` (e.g. a folder of preprocessed `.h5ad` copies next to the Space Ranger outputs, as in the TLS demo). A row whose path matches nothing, or several samples, is logged and listed in `prep.json → skipped`. Without a sheet, every detected sample is used, `sample_id` = the sanitised folder name (the folder above `outs/`) or file stem. |
| `ANNOTATIONS` | `annotations.csv` | File name of optional labels: `sample_id`, `obs_id` (barcode or cell id as in the platform's files), `annotation`. Joined onto `obs["annotation"]` (categorical; missing = NaN). The join counts are in `prep.json → samples.<id>.annotation`. |
| `PLATFORM` | `auto` | `auto` detects each sample by its files; `visium`, `visium-hd`, `xenium`, `merscope`, `cosmx`, `h5ad` or `spatialdata` only looks for that one. |
| `BIN_SIZE_UM` | `8` | Visium HD: which `binned_outputs/square_<NNN>um/` to read. |
| `KEEP_OUT_OF_TISSUE` | `false` | Visium: keep spots with `in_tissue = 0` (only present when the matrix holds them: a raw matrix, or an `.h5ad` that kept them). |
| `MAX_OBS_PER_SAMPLE` | `0` | Keep at most this many observations per sample: a seeded random subset, in source order. `0` keeps all. For quick test runs. |
| `COSMX_PIXEL_UM` | `0.12028` | Microns per pixel of CosMx global coordinates and `Area` (AtoMx export). |
| `SEED` | `42` | Seed of that subsample, combined with the `sample_id` (each sample gets its own reproducible subset). |
| `THREADS` | `0` | Threads for BLAS / numba; `0` = from the container's CPU quota. |

Bad parameters exit with code 2; no usable input, or every sample failing, exits with 1. A
sample that fails is logged and skipped.

## Platforms

Detection is by the files in each folder (one sample per folder, or per file for `.h5ad` and
`.zarr`). Everything not recognised is listed in `prep.json → ignored`.

| Platform | Detected by | Read | `kind` / unit / `grid` |
|---|---|---|---|
| Visium (Space Ranger 1.x–3.x) | `filtered_feature_bc_matrix.h5` **or** the MTX folder `filtered_feature_bc_matrix/{barcodes,features,matrix}` (`.tsv.gz`/`.mtx.gz`), plus `spatial/tissue_positions.csv` (2.x, header) or `tissue_positions_list.csv` (1.x, no header) | Matrix (10x HDF5 with h5py, or MTX with SciPy), positions, `scalefactors_json.json`, `tissue_hires_image.png`, `tissue_lowres_image.png`. `raw_feature_bc_matrix*` is used only when there is no filtered matrix; `spatial_enrichment.csv`, `detected_tissue_image.jpg`, `aligned_fiducials.jpg` and other images are not read. | `spot` / pixel / `hex` |
| Visium HD | `binned_outputs/square_<NNN>um/` with `filtered_feature_bc_matrix.h5` and `spatial/tissue_positions.parquet` | as Visium, for the `BIN_SIZE_UM` binning | `bin` / pixel / `square` |
| Xenium | `cell_feature_matrix.h5` + `cells.parquet` (or `cells.csv.gz`), or the same with a per-sample prefix (`<sample>_cell_feature_matrix.h5` + `<sample>_cells.parquet`, as 10x's downloads name them) | Matrix, `x_centroid`/`y_centroid`, `cell_area`, `nucleus_area` (µm, µm²) | `cell` / micron / `none` |
| MERSCOPE | `*cell_by_gene.csv` + `*cell_metadata.csv` | Counts (dense CSV → CSR in chunks), `center_x`/`center_y`, `volume` → `obs["cell_volume"]`, `fov` | `cell` / micron / `none` |
| CosMx | `*_exprMat_file.csv` + `*_metadata_file.csv` | Counts, `CenterX_global_px`/`CenterY_global_px` × `COSMX_PIXEL_UM`, `Area` × `COSMX_PIXEL_UM`² → `cell_area`; cell ids from a `cell` column, else `<cell_ID>_<fov>` | `cell` / micron / `none` |
| AnnData | `*.h5ad` with `obsm["spatial"]` | `layers["counts"]`, else integer-valued `X`, else integer `raw.X`; Ensembl `var_names` swapped for a symbol column (`feature_name`, `gene_symbols`, …); `uns["spatial"]` and an earlier `uns["pennsieve"]` kept; a file with several `obs["sample_id"]` values is split | from the file, else `cell` / micron |
| SpatialData | `*.zarr` with a `tables` element | The first table; coordinates from its `obsm["spatial"]`, else the centroids of the shapes (or points) element it annotates (`spatialdata_attrs`) | as AnnData |

**Control features**: 10x feature types other than `Gene Expression` (Negative Control Probe,
Negative Control Codeword, Blank / Unassigned / Deprecated Codeword, Genomic Control),
MERSCOPE `Blank-*`, CosMx `NegPrb*`, `SystemControl*`, `FalseCode*`, and names such as
`NegControlProbe_*` / `BLANK_*` are summed per observation into `obs["control_counts"]`
(int64; always present, possibly 0, for imaging platforms) and then dropped from `var`.
Other non-gene features (e.g. Antibody Capture) are dropped without being summed, with a
warning. Counts and feature numbers per type are in `prep.json → features.dropped`.

**SpatialData**: the `spatialdata` package is not bundled (its dependency tree, with dask,
xarray, geopandas and the OME-Zarr stack, would more than double the image). The reader opens
the store's table with `anndata.read_zarr` and computes shape centroids from the GeoParquet
WKB geometries itself (Point, Polygon, MultiPolygon). Images, labels and coordinate
transformations of the store are not read: coordinates are taken as they are in the table
or the shapes element (the element's intrinsic system).

## Outputs

**`spatial/<sample_id>.h5ad`**

| Slot | Content |
|---|---|
| `X` | Raw counts, CSR int32 (float32 with a warning if the source is not integer-valued). |
| `obs_names` | `<sample_id>:<barcode or cell id>` |
| `obs` | `sample_id`, `x`, `y`, `in_tissue` (int32), `array_row`, `array_col` (Visium/HD), `total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `total_counts_mt`, `cell_area`, `nucleus_area` (Xenium, CosMx), `cell_volume` (MERSCOPE), `control_counts` (imaging), `annotation` (when given), platform extras (`fov`, `segmentation_method`), then every sample sheet column. Order as in the source. |
| `var` | Gene symbols (`var_names_make_unique`), `gene_ids`, `feature_types` (`Gene Expression`), `genome` (10x), `mt` (bool, `MT-`/`mt-`). |
| `obsm["spatial"]` | `(n, 2)` float64 `[x, y]` = `obs[["x", "y"]]`. |
| `uns["spatial"][sample_id]` | Visium/HD: `images` (`hires`, `lowres`, uint8), `scalefactors` (as in `scalefactors_json.json`), `metadata` (`software_version`, `chemistry_description` from the h5). |
| `uns["pennsieve"]` | `platform`, `kind`, `coordinate_unit`, `spot_diameter` (Visium: `spot_diameter_fullres`; cells: diameter of a circle of the median cell area, else from MERSCOPE bounds or volume, else the median nearest-neighbour distance), `grid`, `bin_size_um` (HD), `source` (input path), `prep_version`, `sample_id`, `n_obs_source`, `subsampled_from` (when subsampled). |

**`metadata/samples.csv`**: `sample_id, path, platform, kind, n_obs, n_genes, n_obs_in_tissue`,
then every extra sample sheet column, unchanged.

**`metadata/prep.json`**: parameters, the sample sheet and annotations files used, and per
sample: platform, files read (relative paths), `features` (`n_source`, `n_genes`, `dropped` by
type with `n_features`, `total_counts`, `summed_into_control_counts`), `observations`
(`n_source`, `n_out_of_tissue`, `n_out_of_tissue_dropped`, `n_subsampled_to`, `n_written`),
coordinate unit, spot diameter, grid, image presence, `annotation` (`n_annotated`, `n_in_file`,
`n_not_matched`), Space Ranger version, warnings, seconds; `skipped` (sheet rows matching
nothing, duplicate ids, failed samples with the error); `ignored` (files and samples not used,
with the reason; at most 500 listed, `nIgnored` counts all).

**`provenance.json`** (both folders): app and version, tool versions, parameters, inputs
found, timings, per-sample status, peak RSS.

## Runtime

`runtime: cpu 2048, memory 8192 MiB, standard`. Measured with the linux/amd64 image under
amd64 emulation (Docker Desktop on an M-series Mac), `--memory 3g --cpus 2 --network none`,
full data (no subsampling); peak = cgroup `memory.peak` (includes page cache) and the
process's peak RSS:

| Demo | Samples | Observations | Genes | Wall | cgroup peak | peak RSS |
|---|---|---|---|---|---|---|
| Xenium FFPE breast Rep1 + Rep2 | 2 | 167 780 + 118 752 cells | 313 (228 control features dropped) | 27 s | 0.71 GB | 0.73 GB |
| TLS Visium (Space Ranger 2.1.0, MTX) | 8 | 27 715 spots | 17 845 | 60 s | 1.55 GB | 1.20 GB |
| Human Lymph Node (Space Ranger 1.1.0, h5) | 1 | 4 035 spots | 36 601 | 22 s | 1.07 GB | 1.05 GB |

Memory scales with the largest single sample (samples are processed one at a time), about
2–3× the size of its count matrix in memory, plus ~0.4 GB for the interpreter and libraries
and, for Visium, the decoded images. A ~500 k-cell Xenium slide with a 5 k-gene panel
(~150 M non-zeros) or a Visium HD 8 µm sample (~500 k bins × 18 k genes) needs roughly 4–6 GB;
8 GiB is chosen to leave headroom. 2 µm Visium HD bins (millions of bins) are out of scope.
Natively on the Mac the same runs take 5–11 s.

## Limitations

- It does not judge quality. Low-count spots, empty cells and doublets are all kept; `st-qc`
  filters them.
- Out-of-tissue spots are only dropped when the matrix contains them; Space Ranger's
  filtered matrix already holds only in-tissue spots, so `n_out_of_tissue` is usually 0.
- Gene symbols that occur twice (e.g. two Ensembl ids for one symbol) are made unique with a
  `-1` suffix; `var["gene_ids"]` keeps the original ids.
- Coordinates are taken as the platform reports them. Nothing is registered or aligned
  between samples; Visium pixel coordinates of different slides are not comparable.
- Xenium/MERSCOPE/CosMx images, transcripts and cell boundaries are not read (the cell table
  and matrix are enough for the workflows). Only one Visium HD bin size per run.
- CosMx flat files from instruments other than AtoMx may use a different pixel size; set
  `COSMX_PIXEL_UM`.
- `.h5ad` input without integer counts is used as it is, with a warning: downstream steps
  that need counts will then be working on whatever the file holds.

## Development

```sh
/opt/homebrew/bin/python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python      # pytest on the synthetic fixture, no Docker
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/
```

The fixture (`tests/make_fixture.py`) has one tiny sample per supported layout: Space Ranger
1.x `.h5` with `tissue_positions_list.csv` and out-of-tissue spots, 2.x `outs/` with
`tissue_positions.csv`, 2.x MTX folders, Visium HD (in the tests), Xenium with control
probes and codewords, MERSCOPE with `Blank-*`, CosMx with `NegPrb*`/`SystemControl*`, an
`.h5ad` with Ensembl ids, a SpatialData-like `.zarr` with WKB shapes, plus a sample sheet
(with a row that matches nothing) and annotations. See `sample-input/README.md` for real data.

## Licences

Apache-2.0 for this app's code (`LICENSE`). Bundled components are BSD/MIT/Apache/MPL
licensed; nothing GPL is in the image (see `THIRD_PARTY.md`).

## Citations

Verified on Crossref / DataCite and PubMed (title, full author list, year); also in
`CITATION.cff` with full author lists.

- **scanpy** (QC metrics): Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene
  expression data analysis. *Genome Biol* 2018;19:15. doi:10.1186/s13059-017-1382-0, PMID 29409532.
- **anndata** (file format): Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata:
  Access and store annotated data matrices. *J Open Source Softw* 2024;9(101):4371.
  doi:10.21105/joss.04371 (not in PubMed).
- **SpatialData** (store format read): Marconato L, Palla G, Yamauchi KA, et al. SpatialData:
  an open and universal data framework for spatial omics. *Nat Methods* 2025;22:58–62.
  doi:10.1038/s41592-024-02212-x, PMID 38509327.
- **Spatial transcriptomics / Visium**: Ståhl PL, Salmén F, Vickovic S, et al. Visualization
  and analysis of gene expression in tissue sections by spatial transcriptomics. *Science*
  2016;353:78–82. doi:10.1126/science.aaf2403, PMID 27365449.
- **Xenium demo data**: Janesick A, Shelansky R, Gottscho AD, et al. High resolution mapping
  of the tumor microenvironment using integrated single-cell, spatial and in situ analysis.
  *Nat Commun* 2023;14:8353. doi:10.1038/s41467-023-43458-x, PMID 38114474.
- **TLS Visium demo data**: Dawo S, Nonchev K, Silina K. 10x Visium Spatial Transcriptomics
  Dataset: Kidney (3) and Lung (5) Cancer with Tertiary Lymphoid Structures. Zenodo 2025.
  doi:10.5281/zenodo.14620362 (CC BY 4.0).
- **Human Lymph Node Visium demo data**: 10x Genomics, *Human Lymph Node*, Space Ranger 1.1.0
  (CC BY 4.0), files at https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/
  (the 10x Genomics dataset page "Human Lymph Node").
- **Xenium demo data files**: 10x Genomics, *FFPE Human Breast Cancer* (Xenium Onboard Analysis
  1.0.1, CC BY 4.0), files at https://cf.10xgenomics.com/samples/xenium/1.0.1/.

## Files side by side (Pennsieve data sources)

A Pennsieve data source stages the selected files side by side in one folder, without the dataset's
folders, so several samples' files must have distinct names. Xenium, MERSCOPE and CosMx samples can be
given with a per-sample file-name prefix (`breast_rep1_cell_feature_matrix.h5` +
`breast_rep1_cells.parquet`, `run1_cell_by_gene.csv` + `run1_cell_metadata.csv`,
`slide1_exprMat_file.csv` + `slide1_metadata_file.csv`); each prefix is one sample, and a sample
sheet's `path` can name it (`breast_rep1`). Without a sheet the prefix is the sample id. Visium and
Visium HD samples need their `spatial/` folder, so give them as AnnData `.h5ad` files (one per sample)
or one sample per run when staged this way.

