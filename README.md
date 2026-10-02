# Ferro-aging index in colorectal cancer — analysis pipeline

Analysis code for the manuscript:

> **A ferro-aging index in colorectal cancer is prognostic via its senescence/SASP–stromal arm, with TIMP1 as stromal hub**
> Submitted to *Cancer Science*.

This repository contains the **end-to-end, fully public-data pipeline** that regenerates
every figure and table reported in the manuscript from raw downloads to final outputs.
All results are **computational** and derived exclusively from public cohorts; there are no
patient-level private data and no wet-lab experiments in this pipeline.

## Scientific summary

Ferro-aging couples ferroptosis and cellular senescence, but in colorectal cancer (CRC) it
is untested which arm of such a composite carries the prognostic signal. From a curated
95-gene ferro-aging panel we derived a five-gene index (CDKN2A, CHEK1, MMP1, MMP3, TIMP1)
that stratifies overall survival in TCGA-COADREAD and validates in GSE39582. Compositional
ablation shows the prognostic weight is carried by the **senescence/SASP–stromal arm**, not
the ferroptosis arm, and localises to a **TIMP1⁺ CAF-enriched stroma**. The evidence is
computational and TIMP1 did not replicate externally, so the index warrants experimental and
prospective validation.

## Data sources

All inputs are downloaded automatically by `R_scripts/01_data_download.R`.

| Resource | Accession / portal | Role |
|---|---|---|
| TCGA-COADREAD | GDC (https://portal.gdc.cancer.gov) | Training cohort |
| GSE39582 | GEO (GPL570) | External survival validation |
| GSE132465 | GEO | Single-cell localisation |
| GSE281413 | GEO | GeoMx digital spatial profiling |
| MSI / TMB / consensus subtype | cBioPortal Public API | Molecular-subtype association |
| FerrDb V2 / CellAge | curated gene lists | Ferro-aging panel (fixed 95-gene core, defined in `03_FARG_input.R`) |
| DGIdb 4.0 | GraphQL API | Druggability annotation |

GDC and cBioPortal are portal resources and have **no GEO-style accession numbers**; the GEO
series above are accessioned.

## Requirements

- R ≥ 4.6 with the `BiocManager` infrastructure.
- Key packages: `survival`, `survminer`, `glmnet`, `timeROC`, `limma`, `clusterProfiler`,
  `GSVA`, `Seurat`, `scTenifoldKnk`, `rms`, `ComplexHeatmap`, `patchwork`, `ggplot2`,
  `dplyr`, `tidyr` (see `R_scripts/00_config.R` for the core set that is loaded automatically).
- Internet access for GDC / GEO / cBioPortal / DGIdb downloads.
- A fixed random seed (`SEED = 20260919`) makes every analysis reproducible.

## Running the pipeline

From the project root:

```bash
# 1. Download all raw public data (TCGA + GEO) and build the 95-gene panel cache
Rscript R_scripts/01_data_download.R
Rscript R_scripts/03_FARG_input.R

# 2. Run the whole staged analysis
STAGES=all Rscript R_scripts/99_run_all.R
```

`99_run_all.R` executes the pipeline in dependency-ordered stages:

| Stage | Scripts | Purpose |
|---|---|---|
| S0 | 01_data_download, 03_FARG_input | Data preparation |
| S1 | 04_DEG_intersection, 05_GO_KEGG | Differential expression + enrichment |
| S2 | 07_LASSO_multivariate, 08_ExternalValidation | Signature build + external validation |
| S3 | 09_immune_infiltration, 10_Clinical_Table1, 11_Nomogram_DCA, 12_Druggability_DGIdb, 13_hub_gene_TIMP1 | Microenvironment + clinical model + hub |
| S4 | 14_scRNA_celltype, 16_scTenifoldKnk | Single-cell + virtual knockout |
| S5 | 17_UnionVsArms, 18_ExternalAblation, 19_MSI_TMB_CMS | Compositional ablation + subtype |
| S6 | 20_spatial_DSP | Spatial co-localisation |
| S7 | 21_study_design | Study-design figure (Fig 1) |

Run a subset with `STAGES=S1,S2 Rscript R_scripts/99_run_all.R`.
Outputs are written to `figures/`, `tables/` and `results/`; these directories are created
automatically and are git-ignored.

## Repository layout

```
R_scripts/
  00_config.R            # paths, seed, parameters, shared helpers
  01_data_download.R     # GDC + GEO raw downloads
  03_FARG_input.R        # 95-gene ferro-aging panel (defined in-code)
  04_ … 21_              # staged analyses, one script per figure/table group
  99_run_all.R           # dependency-ordered orchestrator
```

## Notes and limitations

- Virtual-knockout (scTenifoldKnk) rankings are **hypothesis-generating**; only TIMP1 survives
  FDR correction, so the enrichment describes the manifold neighbourhood of the knocked-out gene.
- The consensus molecular subtype is reported as the PanCanAtlas subtype, not a formal
  Guinney CMS1–4 call (the reference classifier was not applied here).
- Cell-line IC50 deconvolution was not performed (the GDSC matrix was unreachable).

## Citation

If you use this pipeline, please cite the manuscript (see the manuscript's own reference list).

## License

This project is released under the MIT License — see [LICENSE](LICENSE).
