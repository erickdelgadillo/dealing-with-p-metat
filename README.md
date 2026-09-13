# Dealing with phosphorus deficiency

![R](https://img.shields.io/badge/R-paper%20figures-276DC3?logo=r&logoColor=white)
![Scope](https://img.shields.io/badge/scope-analysis--ready%20data-2E8B57)

Lightweight, reproducible figure workflow for the study **“Dealing with
phosphorus deficiency: contrasting strategies in marine phytoplankton and
bacteria.”**

This repository deliberately starts from small, analysis-ready datasets. Raw
count-table assembly, taxonomic and functional annotation, TPM calculation,
filtering, and differential-expression modelling are no longer performed here.
That upstream processing belongs in a separate workflow repository.

## Publication

Delgadillo-Nuño E., Teira E., Fernández E., Justel-Díez M., Di Leo D., Lundin
D., Pinhassi J. & Martínez-García S. (2026). *Dealing with phosphorus
deficiency: contrasting strategies in marine phytoplankton and bacteria.*
**ISME Communications**, 6(1), ycag035.

DOI: <https://doi.org/10.1093/ismeco/ycag035>

Raw sequencing data: ENA study `PRJEB94162`, samples `ERS25306826`–`ERS25306860`,
reads `ERR15316847`–`ERR15316881`.

## Reproduce the figures

From the repository root:

```bash
Rscript --vanilla analysis/paper_figures.R
```

The script reads only `data/derived/` and writes Figures 1–6 plus Supplementary
Figures 1–2 to `results/figures/`. It does not need `.Renviron`, external disks,
raw counts, annotation tables, or a saved R workspace.

Required R packages are `arrow`, `cowplot`, `dplyr`, `ggh4x`, `ggplot2`,
`magick`, `patchwork`, `readr`, and `scales`.

## Repository structure

```text
dealing-with-p-metat/
├── analysis/
│   └── paper_figures.R
├── data/
│   ├── README.md
│   └── derived/
│       ├── figure_*.csv
│       ├── figure_*.parquet
│       ├── supplementary_*.csv
│       ├── supplementary_*.parquet
│       └── SHA256SUMS
├── results/
│   └── figures/
└── README.md
```

The individual tables, their source versions, and verification status are
documented in `data/README.md`. `data/derived/SHA256SUMS` provides checksums for
the committed data snapshot.

## Provenance and scope change

The derived snapshot was assembled from the small plotting objects and final
local publication materials available before this refactor. The pre-refactor
state remains recoverable from commit `965c166` on `main`, and an independently
verified local backup retains the former `.RData` and `.Rhistory` files.

Removed from the active workflow:

- raw count, taxonomy, and annotation imports;
- TPM calculation and large table joins;
- prokaryotic and eukaryotic `edgeR` modelling;
- exploratory R Markdown notebooks and saved-session dependencies;
- non-paper diagnostic and bubble plots.

Retained in compact form:

- environmental summaries and the site map;
- paper-associated taxonomic composition and nMDS/PERMANOVA outputs;
- phosphorus-gene totals, taxonomic contributions, and proportions;
- complete DGE results needed for Supplementary Figure 2;
- selected DGE rows and annotations needed for Figures 4 and 5;
- final correlation inputs used by Figure 6.

## Known provenance notes

- The prokaryotic composition in Figure 2 uses the paper-associated
  `INTERES_Prok_WP2_TPMs_mean_Annotated_1KEGG_ko.parquet` source. It contains 45
  more rows than the later project copy and produces the published composition.
- No separate paper-final eukaryotic mean table was found locally. Its Figure 2
  composition was reconstructed from the current annotated eukaryotic TPM table
  and checked visually against the publication figure.
- The historical December 2024 correlation workbook and the operational
  workbook agree cell-for-cell on their shared sheets. The operational workbook
  additionally contains the small tables used by all six Figure 6 panels.
- Figures 4 and 5 reproduce the paper's displayed labels. The underlying DGE
  contrasts are stored as `R vs C`, `R+P vs C`, and `R+P vs R`; the publication
  labels the first two in the reverse textual order without reversing `logFC`.
  Both fields remain in the Parquet files so this distinction is explicit.
- An exploratory overwrite in the former integrated notebook caused its local
  Figure 2A to show only Rhodobacterales. The new script uses the verified full
  taxonomic-composition table and restores the all-taxa panel.
- The preliminary project map was replaced with the publication panel from the
  final Figure 1 source, removing its obsolete informal working title.

## Data integrity

Files under `data/derived/` are publication-sized outputs, not raw sequencing
data. New raw or intermediate data should be placed under `data/raw/`,
`data/intermediate/`, or `data/processed/`; those locations are ignored by Git.
