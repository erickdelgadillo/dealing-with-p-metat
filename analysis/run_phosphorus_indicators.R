# Run from the repository root. No upstream data or arrow needed for the snapshot.
if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = "UTC")
source("analysis/phosphorus_indicators.R")
source("analysis/audit_historical_regressions.R")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

save_scientific <- function(plot, stem, width = 12, height = 7) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, device = grDevices::svg)
}

pst_pit_states <- function(states) {
  d <- states[states$target_group %in% c("Rhodobacterales", "Pseudomonadales"), ]
  keys <- c("compartment", "taxon", "rank", "target_group", "treatment", "time", "sample", "annotation_mode")
  pst <- d %>% filter(gene %in% c("pstA", "pstB", "pstC", "pstS")) %>%
    group_by(across(all_of(keys))) %>% summarise(pst_sum_tpm = sum(mean_tpm), n_pst_subunits = n_distinct(gene), .groups = "drop")
  pit <- d %>% filter(gene == "pit") %>% select(all_of(keys), pit_tpm = mean_tpm)
  z <- full_join(pst, pit, by = keys)
  z$pst_pit <- safe_ratio(z$pst_sum_tpm, z$pit_tpm)
  z$ratio_status <- ifelse(is.finite(z$pst_pit), "estimable", "missing or zero denominator/numerator")
  z
}

run_phosphorus_indicators <- function(project_dir = ".") {
  out <- file.path(project_dir, "results", "phosphorus_indicators")
  figs <- file.path(out, "figures")
  dir.create(figs, recursive = TRUE, showWarnings = FALSE)
  write <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
  x <- read_indicator_inputs(project_dir)
  env <- environmental_states(project_dir)
  write(env, "environmental_states")
  write(audit_historical_regressions(project_dir), "historical_regression_audit")
  analyses <- lapply(c("strict", "inclusive"), function(mode) {
    a <- aggregate_indicators(x, mode)
    a$associations <- indicator_associations(a$states, env)
    a$candidates <- rank_candidates(a$states, a$associations, a$expanded)
    a
  })
  names(analyses) <- c("strict", "inclusive")
  for (mode in names(analyses)) {
    a <- analyses[[mode]]
    write(a$states, paste0("taxon_gene_states_", mode))
    write(a$libraries, paste0("taxon_gene_libraries_", mode))
    write(a$associations, paste0("associations_", mode))
    write(a$candidates, paste0("candidate_ranking_", mode))
    write(summarise_candidate_taxa(a$candidates), paste0("candidate_taxa_", mode))
    write(indicator_contrasts(a$states), paste0("contrasts_", mode))
    write(pst_pit_states(a$states), paste0("pst_pit_by_state_", mode))
  }
  # Candidate rows include IDs and original KO membership. This records the
  # additional expression admitted by multi-KO annotations without hiding it.
  sensitivity <- full_join(
    analyses$strict$states %>% select(-annotation_mode),
    analyses$inclusive$states %>% select(-annotation_mode),
    by = c("compartment", "taxon", "rank", "target_group", "gene", "treatment", "time", "sample"),
    suffix = c("_strict", "_inclusive"))
  write(sensitivity, "annotation_sensitivity")
  a <- analyses$strict
  colors <- c("C" = "#0072B2", "R" = "#D55E00", "R+P" = "#009E73")
  theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F3F5"), plot.caption = element_text(hjust = 0),
    legend.position = "bottom"))
  for (g in c("Rhodobacterales", "Pseudomonadales", "Dinoflagellates", "Diatoms")) {
    candidates <- summarise_candidate_taxa(a$candidates)
    candidates <- candidates[candidates$target_group == g, ]
    # Select for readability using the documented ranking, not a p-value cut-off.
    taxa <- head(unique(candidates$taxon), 4)
    selected <- a$states[a$states$target_group == g & a$states$taxon %in% taxa & is_target(a$states), ]
    selected$state <- factor(selected$sample, levels = c("C_0h", "C_72h", "R_0h", "R_72h", "R+P_0h", "R+P_72h"))
    slug <- tolower(g)
    cap <- "Means of technical libraries; one bag per treatment. Strict single-KO annotations; exploratory only."
    # Shares always use all annotated taxa in the same compartment and gene.
    shares <- selected %>% select(state, gene, taxon, mean_contribution)
    residual <- shares %>% group_by(state, gene) %>%
      summarise(mean_contribution = pmax(0, 1 - sum(mean_contribution)), .groups = "drop") %>%
      mutate(taxon = "All remaining taxa")
    contribution <- bind_rows(shares, residual)
    p <- ggplot(contribution, aes(state, mean_contribution, fill = taxon)) +
      geom_col(width = 0.8) + facet_wrap(~gene, ncol = 3) +
      scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste(g, "— contribution to each gene's compartment total"),
        x = NULL, y = "Share of annotated gene TPM", fill = NULL, caption = cap) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save_scientific(p, file.path(figs, paste0("contribution_", slug)), 12, if (g %in% c("Diatoms", "Dinoflagellates")) 5.5 else 8)
    selected$label <- paste(selected$taxon, selected$gene, sep = " | ")
    p <- ggplot(selected, aes(state, label, fill = log10(1 + mean_tpm))) +
      geom_tile(colour = "white", linewidth = 0.3) + scale_fill_viridis_c(name = "log10(1+TPM)") +
      labs(title = paste(g, "— taxon–gene expression"), x = NULL, y = NULL, caption = cap) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save_scientific(p, file.path(figs, paste0("heatmap_", slug)), 12, max(4, 2 + length(unique(selected$label)) * 0.27))
    # Show a single prespecified gene per selected taxon: pstS for bacteria,
    # phoA for eukaryotes. Every other tested association remains in the tables.
    plot_gene <- if (g %in% c("Rhodobacterales", "Pseudomonadales")) "pstS" else "phoA"
    z <- selected[selected$gene == plot_gene, ]
    if (nrow(z)) {
      np <- env[env$predictor == "NP", c("sample", "predictor_value")]
      z <- left_join(z, np, by = "sample")
      make_reg <- function(metric, ylabel) {
        ggplot(z, aes(predictor_value, .data[[metric]])) +
          geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#777777", linetype = 2, linewidth = 0.5) +
          geom_point(aes(colour = treatment, shape = time), size = 2.7) +
          facet_wrap(~taxon, scales = "free_y", ncol = 2) + scale_colour_manual(values = colors) +
          labs(x = "N:P (source environmental mean)", y = ylabel, colour = "Treatment", shape = "Time")
      }
      p <- (make_reg("mean_tpm", paste(plot_gene, "(compartment TPM)")) /
        make_reg("mean_within_taxon_fraction", paste(plot_gene, "/ total taxon RNA (fraction)"))) +
        plot_annotation(title = paste(g, "— descriptive association with N:P"), caption = cap)
      save_scientific(p, file.path(figs, paste0("np_regressions_", slug)), 12, 10)
      p <- ggplot(z, aes(time, mean_tpm, colour = treatment, group = treatment)) +
        geom_line() + geom_point(size = 2.6) +
        geom_linerange(aes(ymin = min_tpm, ymax = max_tpm), linewidth = 0.5) +
        facet_wrap(~taxon, scales = "free_y", ncol = 2) + scale_colour_manual(values = colors) +
        labs(title = paste(g, "—", plot_gene, "across experimental states"), x = NULL, y = "Compartment TPM",
          colour = "Treatment", caption = paste(cap, "Ranges show technical min–max, not biological confidence intervals."))
      save_scientific(p, file.path(figs, paste0("conditions_", slug)), 12, 7)
    }
  }
  ratios <- pst_pit_states(a$states)
  ratios <- left_join(ratios, env[env$predictor == "PO4", c("sample", "predictor_value")], by = "sample")
  keep <- ratios %>% group_by(target_group, taxon) %>% summarise(n = sum(is.finite(pst_pit)), .groups = "drop") %>% filter(n == 6)
  ratios <- semi_join(ratios, keep, by = c("target_group", "taxon"))
  if (nrow(ratios)) {
    p <- ggplot(ratios, aes(predictor_value, pst_pit)) + geom_point(aes(colour = treatment, shape = time), size = 2.5) +
      facet_wrap(~target_group + taxon, scales = "free_y", ncol = 3) + scale_colour_manual(values = colors) +
      labs(title = "Taxon-specific pst:pit ratios", x = "PO4 (µmol/L)", y = "Sum(pstA,pstB,pstC,pstS) / pit",
        caption = "Ratios of state mean TPM; zero/absent pit gives NA. Subunits summed explicitly; no pseudocount. No causal inference.")
    save_scientific(p, file.path(figs, "pst_pit_vs_phosphate"), 13, max(6, 2.5 * ceiling(nrow(keep) / 3)))
  }
  writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))
  message("Generated compact tables and PNG/PDF/SVG figures under ", out)
  invisible(analyses)
}

if (sys.nframe() == 0L) run_phosphorus_indicators()
