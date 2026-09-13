# Derived data manifest

`data/derived/` is the complete input snapshot for `analysis/paper_figures.R`.
No file in this directory requires access to the former raw count, taxonomy,
annotation, or saved-workspace files.

## Files retained

| File | Rows | Used for | Provenance and verification |
| --- | ---: | --- | --- |
| `figure_1_environmental.csv` | 132 | Figure 1 | Environmental means and standard errors exported from the pre-refactor plotting object. |
| `figure_1_map.png` | — | Figure 1 | Map panel extracted losslessly from the final publication Figure 1. |
| `figure_2_taxonomic_composition.csv` | 300 | Figure 2A | Prokaryotes aggregated from the paper-associated final mean-TPM table; eukaryotes reconstructed from the current annotated TPM table and visually checked against the final figure. |
| `figure_2_nmds_scores.csv` | 35 | Figure 2B | Bray–Curtis nMDS scores recalculated from the current full count tables with `vegan::metaMDS` and seed 123 (18 prokaryotic and 17 eukaryotic samples). |
| `figure_2_permanova.csv` | 6 | Figure 2 analysis | PERMANOVA table accompanying the two nMDS analyses. |
| `figure_3_totals.csv` | 72 | Figure 3 | Domain-level phosphorus-gene means, errors, and significance letters from the paper plotting object. |
| `figure_3_taxonomic_contributions.csv` | 1,109 | Figure 3 | Taxonomic contributions to phosphorus-gene categories from the paper plotting object. |
| `figure_4_prokaryote_dge.parquet` | 33,438 | Figure 4 | Prokaryotic DGE results with plotting annotations; visually checked against the final paper panel. |
| `figure_5_eukaryote_dge.parquet` | 120,798 | Figure 5 | Eukaryotic DGE results with plotting annotations; visually checked against the final paper panel. |
| `figure_6_correlations.csv` | 72 | Figure 6 | Six compact correlation tables from the operational workbook; shared sheets were checked cell-for-cell against the December 2024 final workbook. |
| `supplementary_figure_1_proportions.csv` | 72 | Supplementary Figure 1 | Phosphorus-gene category proportions from the paper plotting object. |
| `supplementary_figure_2_dge_overview.parquet` | 1,588,620 | Supplementary Figure 2 | Minimal prokaryotic/eukaryotic DGE overview columns needed for the figure. |

## Source identities

The source filenames and hashes below make the local selection auditable without
embedding workstation-specific paths in the repository.

| Source file | SHA-256 | Note |
| --- | --- | --- |
| `INTERES_Prok_WP2_TPMs_mean_Annotated_1KEGG_ko.parquet` | `ae3b41bcb0656583fe6e15c235054f0126f1baa5c628896d9cf964efbefceb4d` | Paper-associated prokaryotic mean table: 1,220,497 rows. |
| `prok_mean_tpms.parquet` | `10abd0de99380540b6ee534a3c45d0e1b88755d75be3f4f94742e7f857824b84` | Later project copy: 1,220,452 rows; not used for Figure 2A. |
| `euk_tpms_annotated.parquet` | `3af9cb097d9c16d45a772cdfc7bab051e4c4c705742dd4aea3442a3f2fb8cf1f` | Best verified local eukaryotic source; no distinct paper-final mean file was found. |
| `Dataset P y Correlaciones Fnales Dic2024_INTERES.xlsx` | `301518097dff7314e00096596e50096ec1ae8c447fbb97e4d123a9e228b3478c` | Historical publication workbook. |
| `Correlations.xlsx` | `82c05a26fb8b4dbb5997e16926b31e2f6999a88debf4e337799a3a75af31b4e4` | Operational workbook containing the six Figure 6 source sheets. |

The small plotting objects not tied to a separate named source file were
exported from the former `analysis/.RData`. That workspace and its `.Rhistory`
have been verified byte-for-byte against the preserved pre-refactor backup.

## Resolved and unresolved differences

1. **Resolved — prokaryotic Figure 2A:** the paper-associated mean table has 45
   rows absent from the later project table and slightly different TPM values.
   The committed derived table uses the paper-associated source.
2. **Resolved — former Figure 2A overwrite:** the old integrated notebook
   overwrote the full composition object with a Rhodobacterales exploration.
   The committed table restores the full taxonomic composition.
3. **Resolved — Figure 1 map:** the preliminary project map contained an
   informal working title. It was replaced with the final publication panel.
4. **Verified — correlations:** all sheets common to the historical and
   operational workbooks are cell-identical; historical `Hoja1` equals the
   operational `High_affinity` sheet. Five additional operational sheets are
   required for the remaining Figure 6 panels.
5. **Unresolved with certainty — eukaryotic mean source:** no separate final
   eukaryotic mean dataset was located. The best available current source was
   used and its plotted composition matches the paper visually.
6. **Preserved explicitly — DGE comparison names:** Figures 4 and 5 use the
   publication labels, while the stored `comparison` column retains the actual
   model contrast orientation. No sign inversion was applied without evidence.

## Reproducibility notes

- Prokaryotic nMDS stress: `0.03789844`.
- Eukaryotic nMDS stress: `0.05781713`.
- Checksums for every committed derived file are stored in `SHA256SUMS`.
- Derived files should be changed only together with this manifest and the
  regenerated figures.
