# Phosphorus Microcosm Metatranscriptomics

![Version](https://img.shields.io/badge/version-v0.1.0-blue)
![R](https://img.shields.io/badge/R-metatranscriptomics-276DC3)
![Bioconductor](https://img.shields.io/badge/Bioconductor-edgeR-87B13F)
![Status](https://img.shields.io/badge/status-refactoring-yellow)

Reconstruction and refactoring of an R-based metatranscriptomic analysis workflow focused on microbial responses to phosphorus availability.

This repository is being developed from the original analysis scripts used during a marine microbial ecology research project. The goal is to transform the original research code into a structured, documented, and reproducible analytical workflow.

> **Status:** Work in progress.
> The original analysis code is currently being reorganized and validated before the workflow is considered fully reproducible.

## Associated publication

This repository reconstructs and refactors analyses underlying:

**Delgadillo-Nuño E. et al. (2026).**  
*Dealing with phosphorus deficiency: contrasting strategies in marine phytoplankton and bacteria.*  
ISME Communications 6(1), ycag035.

[Read the publication](https://doi.org/10.1093/ismeco/ycag035)

## Overview

The workflow analyzes prokaryotic and eukaryotic metatranscriptomic datasets and integrates functional annotation, gene-expression analysis, community-level statistics, and scientific visualization.

The analysis includes:

* Processing of prokaryotic and eukaryotic transcript count tables
* Functional and taxonomic annotation
* TPM-based transcript abundance calculations
* Differential gene expression analysis with **edgeR**
* Analysis of phosphorus-related metabolic functions
* Bray-Curtis dissimilarities
* Non-metric multidimensional scaling (**NMDS**)
* Permutational multivariate analysis of variance (**PERMANOVA**)
* Statistical summaries and publication-oriented figures

## Project goals

The original analysis was developed iteratively during the research process and contains exploratory code, repeated operations, historical alternatives, and dataset-specific assumptions.

This repository aims to:

1. Preserve the original analysis scripts.
2. Reproduce the original analytical results before modifying the workflow.
3. Separate data processing, statistical analysis, and visualization.
4. Replace duplicated code with reusable R functions.
5. Remove hard-coded paths and dataset-specific column indices where possible.
6. Document software dependencies and computational requirements.
7. Provide a transparent workflow that can be inspected and reproduced.

## Repository structure

```text
marine-p-deficiency-metat/
├── R/                    # Reusable R functions
├── scripts/              # Data processing and statistical analyses
├── figures/              # Scripts used to generate figures
├── data/
│   ├── raw/              # Original input data
│   ├── metadata/         # Sample and experimental metadata
│   ├── intermediate/     # Intermediate analysis objects
│   └── processed/        # Analysis-ready datasets
├── results/
│   ├── figures/          # Generated figures
│   └── tables/           # Generated statistical results
├── reports/              # Analysis reports and documentation
├── legacy/               # Original R / R Markdown scripts
└── README.md
```

## Planned workflow

```text
Raw data
   │
   ├── Prokaryotic transcript counts
   │
   └── Eukaryotic transcript counts
   │
   ▼
Data preparation and annotation
   │
   ├── Count matrices
   ├── TPM tables
   └── Functional / taxonomic annotations
   │
   ▼
Statistical analysis
   │
   ├── Differential gene expression
   ├── Phosphorus metabolism
   ├── Community dissimilarity
   └── NMDS / PERMANOVA
   │
   ▼
Processed results
   │
   ▼
Publication-oriented figures and tables
```

## Main tools

The workflow is primarily implemented in **R** and currently relies on packages from CRAN, Bioconductor, and external repositories.

Core tools include:

* `edgeR`
* `vegan`
* `ggplot2`
* `dplyr`
* `tidyr`
* `arrow`

Additional dependencies will be documented as the original scripts are refactored.

## Reproducibility strategy

Refactoring will be performed incrementally.

The original scripts will first be preserved under `legacy/`. Each analytical component will then be migrated into the new structure while checking that the resulting tables, statistical tests, and figures remain consistent with the original analysis.

Only after this validation step will legacy code and dependencies be modernized.

## Development status

* [x] Repository structure initialized
* [x] Original R environment reconstructed
* [x] Legacy package dependencies resolved
* [x] Original scripts archived
* [x] Prokaryotic preprocessing refactored
* [x] Eukaryotic preprocessing refactored
* [x] Differential expression workflow refactored
* [x] Multivariate analyses refactored
* [x] Figure-generation workflow separated
* [ ] Dependency versions frozen
* [ ] Full reproducibility validation

## Data availability

Input datasets and intermediate files are not yet included in the public repository structure.

Data provenance, availability, file sizes, and any redistribution restrictions will be documented during the reproducibility phase.

## Author

**Erick Delgadillo-Nuño**

Marine microbial ecology · Bioinformatics · Metatranscriptomics
