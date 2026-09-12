# Dealing with phosphorus deficiency

![R](https://img.shields.io/badge/R-metatranscriptomics-276DC3?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-edgeR-87B13F)
![Nextflow](https://img.shields.io/badge/Upstream-nf--core%2Fmetatdenovo-23aa62)
![Status](https://img.shields.io/badge/status-reproducibility%20refactoring-yellow)

> Reproducible R workflows for investigating contrasting phosphorus-acquisition strategies in marine phytoplankton and bacteria using metatranscriptomics.

This repository contains the reconstructed and refactored R workflows associated with the study **“Dealing with phosphorus deficiency: contrasting strategies in marine phytoplankton and bacteria.”**

The project investigates how natural marine microbial communities respond functionally to changes in phosphorus availability caused by high-N:P riverine inputs. A spring mesocosm experiment in the Ría de Vigo compared control seawater with P-deplete and P-replete riverine additions, allowing the transcriptional responses of prokaryotic and eukaryotic microbial communities to be examined under contrasting nutrient conditions.

The repository preserves the downstream analytical workflows used to investigate taxonomic composition, functional profiles, phosphorus-metabolism genes, differential gene expression, multivariate community patterns, and publication figures. Historical research code is being reorganized into a clearer and more reproducible structure while retaining the analyses associated with the published study.

## Scientific scope

The experiment was conducted during spring 2019 at the Toralla Marine Science Station (ECIMAT), Ría de Vigo, NW Spain, using three 500-L mesocosms maintained under in situ conditions.

Three experimental conditions were compared:

| Treatment | Description | Initial DIN | Initial DIP | Phosphorus status |
| --- | --- | ---: | ---: | --- |
| `C` | Control coastal seawater | ~0.3 µM | ~0.07 µM | P-deplete |
| `R` | Seawater + 10% river water | ~3 µM | ~0.07 µM | P-deplete / high N:P |
| `R+P` | Seawater + 10% river water + P addition | ~3 µM | ~1.6 µM | P-replete |

Metatranscriptomic samples were analyzed at **0 h and 72 h**.

The analyses represented in this repository include:

- Prokaryotic and eukaryotic metatranscriptomic datasets
- Taxonomic composition of expressed transcripts
- Functional annotation and transcript-abundance processing
- TPM-based expression analyses
- Phosphorus-acquisition and phosphorus-metabolism genes
- Differential gene expression with `edgeR`
- Bray-Curtis dissimilarities
- Non-metric multidimensional scaling (NMDS)
- PERMANOVA and community-level comparisons
- Integration with nutrient, chlorophyll-a, bacterial production, and primary-production measurements
- Publication-oriented figures and statistical summaries

The study revealed contrasting phosphorus-acquisition strategies between heterotrophic bacteria and eukaryotic phytoplankton under phosphorus-deficient conditions, highlighting different functional responses of the two microbial compartments to changes in nutrient stoichiometry.

## Associated publication

**Delgadillo-Nuño E., Teira E., Fernández E., Justel-Díez M., Di Leo D., Lundin D., Pinhassi J. & Martínez-García S. (2026).**  
*Dealing with phosphorus deficiency: contrasting strategies in marine phytoplankton and bacteria.*  
**ISME Communications**, 6(1), ycag035.

DOI: <https://doi.org/10.1093/ismeco/ycag035>

## Published metatranscriptomic workflow

The raw metatranscriptomic datasets were generated from microbial communities sampled from the mesocosms at 0 h and 72 h.

For the prokaryotic dataset, ribosomal RNA was depleted and the non-ribosomal RNA fraction was analyzed. For the eukaryotic dataset, poly(A)-selected RNA was used. Libraries were sequenced using an Illumina HiSeq 2500 platform.

Raw sequencing data were processed in the published study using the **nf-core/metatdenovo** workflow (`v.dev,1.0`, commit `8e0d117`), which included:

```text
Raw reads
   |
   +-- FastQC / MultiQC
   |
   +-- Cutadapt / Trim Galore
   |
   +-- de novo assembly with MEGAHIT
   |
   +-- ORF prediction
   |      |
   |      +-- Prokka          # prokaryotes
   |      +-- TransDecoder    # eukaryotes
   |
   +-- Read quantification
   |      |
   |      +-- BBMap
   |      +-- featureCounts
   |
   +-- Functional annotation
   |      |
   |      +-- eggNOG-mapper
   |
   +-- Taxonomic assignment
          |
          +-- prokaryotic reference resources
          +-- MAR-MEtASP / eukaryotic reference resources
```

This repository does **not** currently contain a reconstruction of the upstream `nf-core/metatdenovo` processing workflow. Its primary scope is the downstream R-based analysis of the resulting annotated metatranscriptomic datasets.

**PENDIENTE:** document whether the exact upstream pipeline configuration, parameter files, and intermediate outputs used for the publication can be made available or reconstructed independently.

## Published R workflow

The refactored analysis is organized around separate prokaryotic, eukaryotic, and integrated analysis notebooks.

| Notebook | Scope |
| --- | --- |
| `analysis/01_prokaryotes.Rmd` | Prokaryotic metatranscriptomic analysis for C, R, and R+P treatments at 0 h and 72 h; taxonomic and functional annotation, TPM processing, environmental metadata, and differential expression |
| `analysis/02_eukaryotes.Rmd` | Eukaryotic poly(A)-selected metatranscriptomic analysis; taxonomic and functional annotation, TPM processing, environmental metadata, and differential expression |
| `analysis/03_integrated_analysis_and_figures.Rmd` | Integrated historical analysis combining microbial compartments, multivariate analyses, phosphorus-metabolism summaries, and publication-oriented figures |
| `analysis/03_integrated_analysis_and_figures_clean.Rmd` | Cleaned/refactored version of the integrated analysis and figure-generation workflow |

The eukaryotic `R+P` treatment at 0 h contains two samples rather than three because one original sample was excluded from the published analysis due to insufficient sequencing depth.

The current refactoring strategy preserves the analytical logic of the publication while progressively separating preprocessing, statistical analysis, and visualization from the original exploratory research code.

## Running the published workflow

The repository is structured as an R project:

```text
marine-p-deficiency-metat.rproj
```

The main analysis notebooks are under:

```text
analysis/
```

A typical reconstruction workflow is:

```text
01_prokaryotes.Rmd
        |
        v
Prokaryotic processed analyses

02_eukaryotes.Rmd
        |
        v
Eukaryotic processed analyses

        \ /
         v

03_integrated_analysis_and_figures_clean.Rmd
        |
        v
Integrated statistics and publication figures
```

The notebooks currently rely on local metadata together with metatranscriptomic input tables that are not fully distributed through the public repository.

**PENDIENTE:** provide a tested, minimal command sequence for reproducing the complete analysis from a clean R environment once all required public input tables and dependency versions have been finalized.

For example, the final reproducibility workflow may be exposed through commands such as:

```r
rmarkdown::render("analysis/01_prokaryotes.Rmd")
rmarkdown::render("analysis/02_eukaryotes.Rmd")
rmarkdown::render("analysis/03_integrated_analysis_and_figures_clean.Rmd")
```

These commands should be considered provisional until full end-to-end validation is completed.

## Repository structure

```text
marine-p-deficiency-metat/
├── analysis/
│   ├── 01_prokaryotes.Rmd
│   ├── 02_eukaryotes.Rmd
│   ├── 03_integrated_analysis_and_figures.Rmd
│   └── 03_integrated_analysis_and_figures_clean.Rmd
├── R/                    # Reusable R functions and local-data resolver
├── figures/              # Figure-related material
├── legacy/               # Original historical analysis scripts
├── reports/              # Analysis reports and documentation
├── results/              # Generated figures and statistical results
├── scripts/              # Additional analysis scripts
├── .Renviron.example     # Template for the external data location
├── marine-p-deficiency-metat.rproj
└── README.md
```

Research inputs and generated datasets are kept outside the Git repository and
resolved through `R/paths.R`.

## Data availability

The raw sequencing datasets associated with the published study are available from the **European Nucleotide Archive (ENA)**:

- **ENA study:** `PRJEB94162`
- **Sample accessions:** `ERS25306826`–`ERS25306860`
- **Sequence accessions:** `ERR15316847`–`ERR15316881`

ENA study page:

<https://www.ebi.ac.uk/ena/browser/view/PRJEB94162>

Input datasets, metadata, intermediate files, and processed tables used by the R
analysis are not stored in Git. They live under an external data root configured
for each workstation.

Copy `.Renviron.example` to `.Renviron` and set `MARINE_P_DATA_DIR` to the
absolute path of the external `data/` directory:

```bash
cp .Renviron.example .Renviron
```

```text
MARINE_P_DATA_DIR=/absolute/path/to/marine-p-deficiency-metat/data
```

The project-local `.Renviron` file is ignored by Git. If the variable is not
set, the analysis falls back to a local `data/` directory inside the repository;
that directory is also ignored by Git.

Large sequencing files should be retrieved directly from ENA rather than committed to Git.

**PENDIENTE:** document the exact mapping between ENA accessions, experimental treatment, time point, biological replicate, and prokaryotic/eukaryotic library type in a repository sample sheet.

**PENDIENTE:** document which processed outputs from `nf-core/metatdenovo` are required as direct inputs by each R notebook and whether those derived tables can be redistributed through the repository.

## Requirements and validation status

The downstream analysis is implemented primarily in **R**.

Core packages currently used across the reconstructed workflows include:

- `edgeR`
- `vegan`
- `ggplot2`
- `dplyr`
- `tidyr`
- `stringr`
- `readxl`
- `arrow`
- `ggh4x`
- `cowplot`
- `patchwork`
- `magick`
- `plotrix`
- `here`

The repository has already undergone substantial refactoring:

- Repository structure established
- Original R environment reconstructed
- Historical package dependencies recovered
- Original scripts archived
- Prokaryotic preprocessing refactored
- Eukaryotic preprocessing refactored
- Differential-expression analysis refactored
- Multivariate analyses refactored
- Figure-generation workflow separated

Full reproducibility has **not yet been declared complete**.

Remaining validation work includes:

- Freezing exact R/package dependency versions
- Confirming all required input datasets and derived tables
- Running the complete analysis from a clean environment
- Comparing regenerated statistics and figures with the published outputs
- Documenting the exact relationship between ENA raw reads, upstream metatranscriptomic processing, and downstream R inputs

**PENDIENTE:** record a complete software environment (`renv.lock`, container, or equivalent) once dependency validation is finished.

## Reproducibility scope

This repository is a curated reconstruction and refactoring of the downstream computational analyses associated with the published study.

Its purpose is to:

- Preserve the original research analysis
- Separate prokaryotic, eukaryotic, and integrated workflows
- Remove historical hard-coded paths and unnecessary duplicated operations
- Document the relationship between experimental samples, metatranscriptomic data, statistical analyses, and publication figures
- Improve transparency and long-term reproducibility
- Retain the original analytical results before further modernization

The repository does not currently claim complete raw-read-to-publication reproducibility.

Raw reads are publicly available through ENA, while upstream metatranscriptomic processing was performed with `nf-core/metatdenovo`. The current GitHub project focuses primarily on reconstructing and validating the downstream R analysis.

Full end-to-end reproducibility will require linking these layers explicitly:

```text
ENA raw reads
     |
     v
nf-core/metatdenovo
     |
     v
Annotated count / abundance tables
     |
     v
R metatranscriptomic analyses
     |
     v
Statistical results
     |
     v
Publication figures
```

## Citation

If you use material from this repository, please cite the associated publication:

**Delgadillo-Nuño E., Teira E., Fernández E., Justel-Díez M., Di Leo D., Lundin D., Pinhassi J. & Martínez-García S. (2026).**  
*Dealing with phosphorus deficiency: contrasting strategies in marine phytoplankton and bacteria.*  
**ISME Communications**, 6(1), ycag035.  
<https://doi.org/10.1093/ismeco/ycag035>

## Author

**Erick Delgadillo-Nuño**

Marine microbial ecology · Metatranscriptomics · Bioinformatics · R · Reproducible scientific workflows
