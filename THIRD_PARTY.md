# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following
third-party components, installed from PyPI at the versions pinned in `requirements.txt`.
All are permissively licensed (BSD, MIT, Apache-2.0, MPL-2.0, PSF); nothing GPL is in the
image (scanpy's optional igraph / leidenalg are not installed: this step does not cluster).

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| scanpy | 1.12.5 | QC metrics (`pp.calculate_qc_metrics`) | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata | 0.13.4 | The `.h5ad` files read and written; SpatialData tables (`read_zarr`) | BSD-3-Clause | https://github.com/scverse/anndata |
| h5py | 3.16.0 | Reads 10x HDF5 matrices (Space Ranger, Xenium) | BSD-3-Clause | https://www.h5py.org |
| NumPy | 2.5.4 | Arrays | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| SciPy | 1.18.1 | Sparse matrices (CSR), MTX reading, nearest-neighbour spacing | BSD-3-Clause | https://scipy.org |
| pandas | 3.1.0 | Positions, cell tables, sample sheet, annotations | BSD-3-Clause | https://pandas.pydata.org |
| pyarrow | 26.0.0 | Parquet (Xenium `cells.parquet`, Visium HD positions, SpatialData shapes) | Apache-2.0 | https://arrow.apache.org |
| Pillow | 13.0.0 | Reads the Visium tissue images | MIT-CMU | https://python-pillow.org |
| zarr, numcodecs, donfig, google-crc32c | 3.4.2, 0.17.1, 0.8.1.post1, 1.9.0 | SpatialData `.zarr` stores (via anndata) | MIT; MIT; MIT; Apache-2.0 | https://zarr.dev |
| scikit-learn, joblib, threadpoolctl | 1.9.1, 1.6.0, 3.7.0 | scanpy dependencies | BSD-3-Clause | https://scikit-learn.org |
| numba, llvmlite | 0.68.1, 0.50.0 | scanpy dependencies (JIT) | BSD-2-Clause; BSD-2-Clause AND Apache-2.0 WITH LLVM-exception | https://numba.pydata.org |
| umap-learn, pynndescent | 0.5.12, 0.7.0 | scanpy dependencies (unused here) | BSD-3-Clause; BSD-2-Clause | https://github.com/lmcinnes/umap |
| statsmodels, patsy, formulaic, interface_meta, wrapt | 0.15.0, 1.0.3, 1.2.2, 2.0.1, 2.5.2 | scanpy dependencies (unused here) | BSD-3-Clause; BSD-2-Clause; MIT; MIT; BSD-2-Clause | PyPI |
| matplotlib, seaborn, contourpy, cycler, kiwisolver, fonttools, pyparsing | 3.11.2, 0.13.2, 1.4.0, 0.12.1, 1.5.1, 4.68.0, 3.3.3 | scanpy dependencies (no plots are made) | Matplotlib License (PSF-based); BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; BSD-3-Clause; MIT; MIT | https://matplotlib.org |
| networkx, natsort, tqdm, session-info2, legacy-api-wrap, scverse-misc, fast-array-utils | 3.7, 8.4.0, 4.70.1, 0.4.2, 1.5, 0.1.8, 1.5.1 | scanpy / anndata dependencies | BSD-3-Clause; MIT; MPL-2.0 AND MIT; MPL-2.0; MPL-2.0; BSD-3-Clause; MPL-2.0 | PyPI |
| array-api-compat, narwhals, cloudpickle, msgspec | 1.16.0, 2.28.0, 3.1.2, 0.22.0 | anndata dependencies | MIT; MIT; BSD-3-Clause; BSD-3-Clause | PyPI |
| pydantic, pydantic-core, pydantic-settings, annotated-types, typing-inspection, python-dotenv, PyYAML | 2.14.2, 2.50.2, 2.16.2, 0.8.0, 0.4.4, 1.2.4, 6.0.3 | scverse-misc dependencies | MIT (python-dotenv BSD-3-Clause) | PyPI |
| packaging, python-dateutil, six, tzdata, typing_extensions, certifi | 26.3, 3.0.0, 1.17.0, 2026.6, 4.16.0, 2026.10.23 | dependencies | Apache-2.0 OR BSD-2-Clause; Apache-2.0/BSD dual; MIT; Apache-2.0; PSF-2.0; MPL-2.0 | PyPI |
| python:3.12.16-slim-bookworm | — | Base image (Debian 12 + CPython 3.12, pip) | PSF-2.0 (Python); Debian packages under their own free licences; pip MIT | https://hub.docker.com/_/python |

The SpatialData package (`spatialdata`) is **not** bundled: its dependency tree (dask,
xarray, geopandas/shapely, spatial-image, multiscale-spatial-image, ome-zarr, ...) would more
than double the image for a reader that only needs the table. The store is read with
anndata and pyarrow instead (see README, *SpatialData*).

Demo data used for testing (not in the image or the repository), all CC BY 4.0:
10x Genomics *Human Lymph Node* (Visium, Space Ranger 1.1.0); 10x Genomics *FFPE Human Breast
Cancer* Xenium replicates 1 and 2 (Janesick et al. 2023, doi:10.1038/s41467-023-43458-x);
Dawo S, Nonchev K, Silina K. *10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung
(5) Cancer with Tertiary Lymphoid Structures*, Zenodo 2025, doi:10.5281/zenodo.14620362.
