# Flavonoid intake and incident psoriasis: analysis code

This repository contains the R code used for the manuscript examining dietary flavonoid intake and risk of incident psoriasis using UK Biobank data.

## Repository structure

```text
.
├── _targets.R                         # main reproducible pipeline
├── run_pipeline.R                     # one-command runner
├── renv.lock                          # R package environment record/template
├── R/
│   ├── 00_setup.R
│   ├── 01_data_preparation.R
│   ├── 02_flavonoid_intake_construction.R
│   ├── 03_survival_and_covariates.R
│   ├── 04_food_source_and_spline_analysis.R
│   └── 05_cox_models_tables_interactions.R
├── data/
│   ├── README_data.md
│   └── raw/                           # not tracked; UK Biobank files go here locally
├── results/
│   ├── tables/
│   └── figures/
└── manuscript/
    ├── Supplementary_Methods_Code_Reproducibility.md
    └── Code_Availability_Statement.md
```

## How to run

1. Clone or download this repository.
2. Place the required input files in `data/raw/` as described in `data/README_data.md`.
3. Open R in the project root and run:

```r
source("run_pipeline.R")
```

Alternatively:

```r
install.packages(c("renv", "targets", "tarchetypes"))
renv::restore()
targets::tar_make()
```

## Data availability

The individual-level data used in this study are available from UK Biobank upon successful application. They are not publicly shared here due to data access restrictions and participant confidentiality requirements.

## Code availability

The analysis code is provided for transparency and reproducibility. After creating a public GitHub repository, archive the release using Zenodo to obtain a permanent DOI.

## Main R packages

- `dplyr`, `tidyr`, `purrr`, `stringr`
- `survival`, `rms`, `mediation`
- `tableone`, `broom`, `knitr`
- `ggplot2`, `ggsurvfit`, `forestmodel`
- `targets`, `tarchetypes`, `renv`

## Citation

Please cite the associated manuscript when using this code.
