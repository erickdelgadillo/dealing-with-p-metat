# Species-labelled exploratory phosphorus plots from INTERES/WP2.
# These are original GTDB assignments, NOT sequence-validated species calls.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

outdir <- file.path("results", "phosphorus_species")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

targets <- c("Planktomarina_temperata", "Amylibacter_sediminis",
             "HTCC2207_sp018607525", "HTCC2207_sp002334915")
input <- file.path("data", "phosphorus_indicators")
d <- read.csv(file.path(input, "transcripts.csv"), stringsAsFactors = FALSE)
tx <- read.csv(file.path(input, "original_taxonomy.csv"), stringsAsFactors = FALSE)
libs <- read.csv(file.path(input, "library_qc.csv"), stringsAsFactors = FALSE)
env <- read.csv(file.path("data", "derived", "figure_6_correlations.csv"),
                stringsAsFactors = FALSE)

# Match taxonomic identity by ORF identifier; do not infer species from genus.
tx <- unique(tx[c("compartment", "feature_id", "Species")])
if (anyDuplicated(tx[c("compartment", "feature_id")])) stop("Conflicting source taxonomy")
names(tx)[3] <- "source_species"
d <- merge(d, tx, by = c("compartment", "feature_id"), all.x = TRUE)
d$source_species <- trimws(d$source_species)
d <- subset(d, compartment == "Prokaryotes" & !ambiguous_ko &
            source_species %in% targets)
if (!nrow(d)) stop("No selected annotated transcripts")
if (any(!is.finite(d$tpm) | d$tpm < 0)) stop("Invalid TPM")
if (anyDuplicated(d[c("feature_id", "treatment", "time", "replic")]))
  stop("Duplicate transcript-library rows")

# Expand the exact KO-derived gene labels already used in the genus analysis.
genes <- strsplit(d$matched_genes, ";", fixed = TRUE)
d <- d[rep(seq_len(nrow(d)), lengths(genes)), , drop = FALSE]
d$gene <- unlist(genes, use.names = FALSE)
d <- subset(d, gene %in% c("pstS", "pit"))
if (!nrow(d)) stop("No target gene annotation")

# A missing feature row is an operational non-detection, not biological absence.
by <- c("source_species", "gene", "treatment", "time", "replic")
sums <- aggregate(d$tpm, by = d[by], FUN = sum)
names(sums)[ncol(sums)] <- "tpm"
lib <- subset(libs, compartment == "Prokaryotes",
              select = c("treatment", "time", "replic"))
lib <- unique(lib)
observed_pairs <- unique(d[c("source_species", "gene")])
grid <- merge(observed_pairs, lib)
z <- merge(grid, sums, by = by, all.x = TRUE)
z$tpm[is.na(z$tpm)] <- 0

# Technical replicate means; one bag per treatment, so no inferential CI.
state <- aggregate(z$tpm, by = z[c("source_species", "gene", "treatment", "time")],
                   FUN = mean)
names(state)[ncol(state)] <- "mean_tpm"
state$sample <- paste(state$treatment, state$time, sep = "_")
e <- subset(env, panel == "High_affinity" & Final_Taxonomy == "Rhodobacterales",
            select = c("sample", "predictor", "predictor_value"))
e <- unique(subset(e, predictor == "PO4", select = c("sample", "predictor_value")))
if (anyDuplicated(e$sample) || nrow(e) != 6L) stop("Expected six historical PO4 states")
names(e)[2] <- "PO4"
state <- merge(state, e, by = "sample", all.x = TRUE)
if (anyNA(state$PO4)) stop("Missing historical phosphate values")

# Counts are annotated feature IDs, not genomes, cells or independent gene copies.
counts <- aggregate(d$feature_id, by = d[c("source_species", "gene")],
                    FUN = function(x) length(unique(x)))
names(counts)[ncol(counts)] <- "n_feature_ids"
write.csv(counts, file.path(outdir, "species_feature_counts.csv"), row.names = FALSE)
write.csv(state, file.path(outdir, "species_gene_states.csv"), row.names = FALSE)
write.csv(z, file.path(outdir, "species_gene_technical_libraries.csv"), row.names = FALSE)

theme_set(theme_bw(base_size = 11))
for (sp in targets) {
  dat <- subset(state, source_species == sp)
  if (!nrow(dat)) next
  dat$gene <- factor(dat$gene, levels = c("pstS", "pit"))
  dat$state <- factor(dat$sample,
                      levels = c("C_0h", "R_0h", "R+P_0h",
                                 "C_72h", "R_72h", "R+P_72h"))
  p1 <- ggplot(dat, aes(PO4, mean_tpm)) +
    geom_point(size = 2.2, color = "#333333") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.45, color = "#777777") +
    scale_x_log10() +
    facet_wrap(~gene, scales = "free_y") +
    labs(x = expression(PO[4]~(mu*M)), y = "Mean TPM",
         subtitle = "Six experimental states; descriptive linear fit") +
    theme(panel.grid.minor = element_blank())
  p2 <- ggplot(dat, aes(state, mean_tpm, group = gene)) +
    geom_point(size = 2.1, color = "#333333") +
    facet_wrap(~gene, scales = "free_y") +
    labs(x = NULL, y = "Mean TPM",
         subtitle = "C, R and R+P at 0 and 72 h") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())
  plot <- (p1 / p2) + plot_annotation(
    title = gsub("_", " ", sp, fixed = TRUE),
    subtitle = "Original GTDB species label; sequence identity unverified",
    caption = paste("INTERES/WP2 | KO-strict | technical-library means |",
                    "one mesocosm per treatment; trends are exploratory"))
  name <- gsub("[^A-Za-z0-9_]", "_", sp)
  ggsave(file.path(outdir, paste0(name, ".png")), plot,
         width = 8.5, height = 8.1, dpi = 300, bg = "white")
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot,
         width = 8.5, height = 8.1, device = cairo_pdf, bg = "white")
  ggsave(file.path(outdir, paste0(name, ".svg")), plot,
         width = 8.5, height = 8.1, device = svg, bg = "white")
}
message("Generated species-labelled summaries in ", outdir)
