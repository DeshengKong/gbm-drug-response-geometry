# Drug-response geometry along the glioblastoma trajectory

Analysis code for the Cell Reports manuscript **"Drug-response geometry along the glioblastoma trajectory"** (Kong, D.).

## Summary

This repository contains the code and result tables behind a multi-dataset computational reanalysis of glioblastoma (GBM) single-cell and spatial data. The central finding is *drug-response geometry*: the GBM cellular trajectory is a real population-level coordinate; drug-response **direction** is a function of trajectory position, while response **amplitude** is not predicted by baseline tumor or myeloid features and is only stratifiable at the tissue level; and **VEGFA** is the single ligand that tracks trajectory position across five independent datasets.

All analyses are reanalyses of public data. No new experiments were performed, and no survival or outcome labels exist in any dataset used, so nothing here can be read as efficacy.

## Data

Large raw datasets are **not** included in this repository. Download them from the accessions below before running the code.

| Dataset | Accession / location | Used for |
|---|---|---|
| de Jong et al., Nature 2026 — GBM snRNA-seq (12 IDH-wt primary GBM) and Table S10 meta-module gene lists | E-MTAB-17183 (h5ad) + Table S10 (DOI 10.1038/s41586-026-11234-8) | Trajectory axis definition (A-side) |
| Lavrenteva & Kim preprint — ex vivo drug perturbation of glioma slice cultures (Levitin et al. data, 391,444 cells, 43 treatment–vehicle comparisons, 19 TissueID, 9 drugs) | GEO GSE148842 (Levitin loom, DOI 10.64898/2026.11.03.763852) | Points 1–4: direction, amplitude, myeloid, stratification |
| Meylan et al., Cell 2026 — oncolytic HSV (oHSV) Xenium spatial, 8 patients paired pre/post | GEO GSE296577 (Seurat object `GSE296577_GBM_oHSV_seurat_Full_20260107.rds`, 3.8 GB); CODEX via EBI BioImage Archive S-BIAD1921 (DOI 10.1016/j.cell.2025.12.055, PMID 41679299) | Spatial validation: T-cell influx, VEGFA–hypoxia, hypoxic T-cell exclusion |
| Delta-24-RGD virotherapy cohort (MD Anderson, PI Fueyo) — Xenium + TRB, 15 human Xenium sections + 16 TRB samples | GEO GSE330794 (suppl. is a single 213 GB RAW tar; extract by byte range, see `03-delta24-replication/extract_tar.py`) | Independent replication: VEGFA–hypoxia, TCR clonal expansion |
| GBmap — GBM single-cell atlas (Ruiz-Moreno et al., Neuro-Oncology 2025) | CELLxGENE collection `999f2a15-3d7e-440b-96ae-2c806799c08c` (core 16 datasets, 110 patients, 338,564 cells; DOI 10.1093/neuonc/noaf113) | Atlas-scale validation of trajectory and ligand zoning |
| Wang et al., Nature Cancer 2022 — 29 primary–recurrent GBM pairs | GEO GSE174554 (DOI 10.1038/s43018-022-00475-x) | Recurrent-tumor test: mesenchymal shift and AP-1 |

Notes on excluded large files: the intermediate per-cell tables that exceed the size limit are regenerable from the accessions above — `colab_outputs/myeloid_panel_counts.npz` (27 MB), `meylan_analysis/mm2_tcell_tumor_cells.csv.gz` (59 MB), `meylan_analysis/mm2_genes_long.csv.gz` (7.3 MB), `gbmap/gbmap_per_cell.npz` (24 MB), and the raw data files themselves (GSE296577 RDS 3.8 GB, GSE330794 RAW tar 213 GB, GSE174554 RAW tar 1.65 GB, GBmap core h5ad).

## Repository layout

```
01-trajectory-drug-geometry/   Points 1–4 on the B-side data (Figures 1–4)
02-meylan-spatial/             Meylan oHSV Xenium spatial analysis (Figure 5)
03-delta24-replication/        Delta-24-RGD replication (Figure 6)
04-gbmap-wang/                 GBmap atlas validation + Wang recurrent test (Figure 7)
docs/                          Research memo and the v3 stratified-module spec
```

### `01-trajectory-drug-geometry/` — Figures 1–4

- `build_notebook.py` … `build_notebook5.py` — notebook generators. Each builds a Jupyter notebook (`gbm-point*.ipynb`) that runs the corresponding analysis on Google Colab; the executed notebooks (`*_output.ipynb`) are included.
- `gbm-point1-colab*.ipynb` — Point 1: TRAJ definition (late modules MM2/MM3/MM5 minus early modules MM1/MM4/MM6/MM7/MM8/MM10/MM11, scanpy `score_genes`), direction C vs trajectory position, amplitude vs trajectory, within-vs-composition decomposition.
- `gbm-point2-drug-module*.ipynb` — Point 2 v1: drug-module prototype (LOPO/LODO).
- `gbm-point2v2-point4*.ipynb` — Point 2 v2 (amplitude prediction, failed) and Point 4 (late-module score as tissue stratifier).
- `gbm-point2v3-stratified*.ipynb` — Point 2 v3: stratified (tertile) module, the final supported form.
- `gbm-point3-myeloid*.ipynb` — Point 3: myeloid composition/shift and ligand–trajectory coupling.
- `transition_matrix.py` + `transition_matrix_results_2026-11-08.md` — TRAJ-quintile compositional flow (partial reversion to bin 1, not to the most early state).
- `inspect_h5ad.py`, `inspect_loom.py`, `smoke_test.py` — data inspection and local smoke tests.
- `colab_outputs/`, `colab_outputs2v2/`, `colab_outputs2v3/` — result CSVs and figures.
- `MANIFEST.md` — data inventory (accessions, MD5s, file sizes).

### `02-meylan-spatial/` — Figure 5

- `extract_meylan.R`, `export_genes_csv.R`, `extract_mm2_genes.R` — R + SeuratObject extraction of metadata, centroids, and the 22-gene panel from the GSE296577 RDS.
- `build_notebook6.py`, `build_notebook7.py`, `build_notebook8.py` — notebook generators for the Q1/Q2/Q4 analysis (the executed notebooks are not included; regenerate with these scripts).
- `fixup_paired.py` — pairing fix (strips `_pre`/`_post` suffixes from `orig.ident` before pivoting).
- `find_tcr.R`, `find_tcr2.R`, `find_tcr3.R` — TCR clonotype search (negative result: no per-cell clonotype mapping exists in the RDS).
- `parallel_range_download.py`, `resume_rds_curl.py` — byte-range download of the 3.8 GB RDS.
- `meylan_mm2_tcell.py`, `meylan_mm2_tcell_followup.py`, `run_mm2_extract.py` — hypoxic T-cell exclusion analysis (tumor-cell to nearest-T distance, author mes2 vs tier2 discrete states).
- `inspect_rds.py` — RDS structure inspection.
- `analysis_design.md`, `meylan_results_2026-11-08.md`, `meylan_mm2_tcell_2026-11-08.md` — pre-registered designs and result memos.
- `q1_composition.csv`, `q1_distances.csv`, `q2a_ligand_by_tier2.csv`, `q2c_ligand_vs_mes2.csv`, `meylan_summary.json`, `meylan_analysis/mm2_tcell_per_section_*.csv` — result tables.

### `03-delta24-replication/` — Figure 6

- `extract_tar.py`, `resume_curl.py`, `download_list.txt` — byte-range extraction of the 61 needed files (15 sections × cell_feature_matrix.h5 + transcripts.parquet.gz + cells.parquet.gz, plus 16 TRB files) from the 213 GB `GSE330794_RAW.tar`.
- `trb_analysis.py` + `trb_analysis/` — TRB diversity/overlap (paired n=5).
- `xenium_analysis.py`, `xenium_fixup.py` + `xenium_analysis/` — Xenium marker-score cell typing, VEGFA–MM2 replication (14/14 sections), T-cell and myeloid fractions.
- `replication_design.md`, `replication_results_2026-11-08.md` — pre-registered design and results.

### `04-gbmap-wang/` — Figure 7

- `gbmap_analysis.py`, `gbmap_results_2026-11-08.md`, `analysis_design.md` — GBmap atlas validation (TRAJ vs PC1, rho=0.979; per-donor ligand–trajectory correlations). GSE148842 is excluded from the atlas (only one patient) for independence.
- `wang_analysis.py`, `wang_ap1.py`, `wang_results_2026-11-08.md`, `wang_ap1_results_2026-11-08.md`, `analysis_design.md` — Wang 29-pair recurrent test (S10 modules + direct AP-1 core-set test). Note: pair #33 is dropped (deposited MTX file `GSM5319529` is truncated).

### `docs/`

- `gbm-trajectory-drug-geometry-link-20261105.md` — research memo covering Points 1–4 and the systematic interpretation.
- `v3_module_spec.json` — machine-readable specification of the stratified drug-response module (late-score tertile cutoffs −0.0134 / 0.0619; stratum expected m12 0.038 / 0.071 / 0.093).

## Requirements

- Python 3.10+ with scanpy 1.12.5, anndata 0.13.4, scipy, scikit-learn, pandas, numpy (see the `*_output.ipynb` headers for the exact environment used).
- R 4.6.1 + SeuratObject 5.4.0 (only for the Meylan RDS extraction in `02-meylan-spatial/`).
- The analyses were run on Google Colab because of data size (the B-side loom is 7.4 GB; the Meylan RDS is 3.8 GB; this VM has 7 GB RAM). The `build_notebook*.py` generators build the notebooks; execute them on Colab (or any machine with enough memory). Local `smoke_test.py` reproduces the Colab numbers on small subsets.
- A local patch is needed for `google-colab-cli` 0.7.5: `_post_assignment` must send `data=b"{}"` with `Content-Type: application/json` (works around a Google frontend HTTP 411).

## Reproducing

Run in order; `01` first, then `02`, `03`, `04` independently (they share no intermediate files).

1. Download the datasets from the accessions in the Data section above.
2. `01-trajectory-drug-geometry/`: run `build_notebook.py` … `build_notebook5.py` to generate the notebooks, execute them on Colab, and compare the outputs in `colab_outputs/`, `colab_outputs2v2/`, `colab_outputs2v3/` with the included results. Then run `transition_matrix.py` for the compositional-flow refinement.
3. `02-meylan-spatial/`: download the RDS with `parallel_range_download.py` / `resume_rds_curl.py`, extract with the R scripts, then run the notebook generators and the `meylan_mm2_tcell*.py` scripts.
4. `03-delta24-replication/`: extract the 61 files from the RAW tar with `extract_tar.py` (byte ranges; see `download_list.txt`), then run `trb_analysis.py` and `xenium_analysis.py` (+ `xenium_fixup.py`).
5. `04-gbmap-wang/`: download the GBmap collection and GSE174554, then run `gbmap_analysis.py`, `wang_analysis.py`, and `wang_ap1.py`.

Pre-registered analysis designs (`analysis_design.md`, `replication_design.md`) were written before each analysis was run; the memos record where pre-set bars were not met (notably Point 2 v1/v2 amplitude prediction and the comparison-level v3 test).

## License

MIT — see `LICENSE`. Copyright (c) 2026 Desheng Kong.
