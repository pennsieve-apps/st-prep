# sample-input

st-prep reads raw spatial transcriptomics outputs. The three demo datasets of the workflows
are downloaded and checksummed by `scripts/fetch-sample-data.sh` in
`workflow-spatial-transcriptomics` (all CC BY 4.0):

```sh
cd ../workflow-spatial-transcriptomics
scripts/fetch-sample-data.sh all          # or tls | lymph | xenium
```

| Folder | What | Layout st-prep sees |
|---|---|---|
| `data/tls-visium/` | Dawo, Nonchev, Silina, Zenodo 2025, doi:10.5281/zenodo.14620362: 8 FFPE tumour sections (kidney ×3, lung ×5) with tertiary lymphoid structures, Space Ranger 2.1.0 | `TLS_VISIUM_USZ/10x_Visium/<id>/filtered_feature_bc_matrix/` (MTX only, no `.h5`), `spatial/tissue_positions.csv`; `samples.csv` (`sample_id, path, tissue, condition`) and `annotations.csv` (24 246 spots labelled TLS / Immune / Tumor / Normal / Unassigned). The `h5ad_preprocessed/` copies are not in the sheet and are ignored. |
| `data/lymph-node/` | 10x Genomics *Human Lymph Node*, Visium, Space Ranger 1.1.0 | `V1_Human_Lymph_Node/filtered_feature_bc_matrix.h5`, `spatial/tissue_positions_list.csv` (no header); `samples.csv`. Link only `V1_Human_Lymph_Node/` and `samples.csv`: `reference/` is the scRNA-seq reference for cell2location, not a spatial sample. |
| `data/xenium-breast/` | 10x Genomics *FFPE Human Breast Cancer*, Xenium Onboard Analysis 1.0.1, Rep1 and Rep2 (Janesick et al. 2023) | `Rep1/`, `Rep2/` with `cell_feature_matrix.h5`, `cells.parquet`; `samples.csv` (`sample_id, path, condition, replicate`) |

Run the image on one of them (read-only input, outputs in `out/spatial` and `out/metadata`):

```sh
D=../workflow-spatial-transcriptomics/data
docker run --rm --platform linux/amd64 --memory 3g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v $D/xenium-breast:/data/input:ro -v $PWD/out:/data/output st-prep:dev
```

Add `-e MAX_OBS_PER_SAMPLE=20000` for a quick subsampled run. The outputs of these runs feed
`st-qc` (all three workflows), then `st-domains`, `cell2location` or `st-celltype`; the
workflow builders keep them under `data/work/{tls,lymph,xenium,xenium-20k}/prep/`.

## Your own data

Put each sample's instrument output in its own folder (Space Ranger `outs/`, a Xenium output
folder, a MERSCOPE region folder with `cell_by_gene.csv` and `cell_metadata.csv`, a CosMx
flat-file export), add `samples.csv`:

```csv
sample_id,path,condition,donor
ctrl_1,run1/outs,control,D1
treated_1,run2/outs,treated,D1
```

and, optionally, `annotations.csv` (`sample_id, obs_id, annotation`, where `obs_id` is the
barcode or cell id as in the platform's files, e.g. `AAACAAGTATCTCCCA-1`).

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python    # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python        # Docker over data/fixture/linked/ -> data/output/
```
