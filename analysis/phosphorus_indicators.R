suppressPackageStartupMessages(library(dplyr))

safe_ratio <- function(numerator, denominator) {
  ifelse(is.finite(denominator) & denominator > 0, numerator / denominator, NA_real_)
}

read_indicator_inputs <- function(project_dir = ".") {
  p <- file.path(project_dir, "data", "phosphorus_indicators")
  read <- function(f) read.csv(file.path(p, f), stringsAsFactors = FALSE, na.strings = "")
  list(transcripts = read("transcripts.csv"), totals = read("taxon_totals.csv"),
       libraries = read("library_qc.csv"))
}

validate_indicator_inputs <- function(x) {
  d <- x$transcripts; l <- x$libraries; t <- x$totals
  keys <- c("compartment", "feature_id", "treatment", "time", "replic")
  if (anyDuplicated(d[keys])) stop("Duplicate transcript-library key")
  lk <- c("compartment", "treatment", "time", "replic")
  if (anyDuplicated(l[lk])) stop("Duplicate library key")
  if (anyDuplicated(t[c(lk, "taxon", "rank", "target_group")])) stop("Duplicate taxon total")
  if (any(!is.finite(d$tpm) | d$tpm < 0) || any(!is.finite(t$total_taxon_tpm) | t$total_taxon_tpm < 0)) stop("Invalid TPM")
  if (any(abs(l$total_tpm - 1e6) > 1)) stop("Invalid library TPM sum")
  if (nrow(anti_join(d, l, by = lk))) stop("Transcript without a library")
  if (any(!d$treatment %in% c("C", "R", "R+P")) || any(!d$time %in% c("0h", "72h"))) stop("Unexpected experimental state")
  invisible(TRUE)
}

expand_genes <- function(d, mode = c("strict", "inclusive")) {
  mode <- match.arg(mode)
  if (mode == "strict") d <- d[!d$ambiguous_ko, ]
  genes <- strsplit(d$matched_genes, ";", fixed = TRUE)
  out <- d[rep(seq_len(nrow(d)), lengths(genes)), ]
  out$gene <- unlist(genes, use.names = FALSE)
  out
}

aggregate_indicators <- function(x, mode = "strict") {
  validate_indicator_inputs(x)
  d <- expand_genes(x$transcripts, mode)
  tk <- c("compartment", "taxon", "rank", "target_group")
  lk <- c("compartment", "treatment", "time", "replic")
  # Complete only observed libraries. R+P 0h euk R3 is never invented.
  panels <- distinct(d, across(all_of(c(tk, "gene"))))
  grid <- merge(panels, x$libraries[lk], by = "compartment")
  sums <- d %>% group_by(across(all_of(c(tk, "gene", "treatment", "time", "replic")))) %>%
    summarise(tpm = sum(tpm), n_features_detected = n_distinct(feature_id), .groups = "drop")
  lib <- left_join(grid, sums, by = c(tk, "gene", "treatment", "time", "replic"))
  lib$not_detected <- is.na(lib$tpm)
  lib$tpm[lib$not_detected] <- 0
  lib$n_features_detected[lib$not_detected] <- 0L
  # Sparse positive-TPM source: absent rows are operational non-detections.
  lib <- left_join(lib, x$totals, by = c(tk, "treatment", "time", "replic"))
  lib$total_taxon_tpm[is.na(lib$total_taxon_tpm)] <- 0
  lib$within_taxon_fraction <- safe_ratio(lib$tpm, lib$total_taxon_tpm)
  totals <- lib %>% group_by(compartment, gene, treatment, time, replic) %>%
    summarise(all_taxa_gene_tpm = sum(tpm), .groups = "drop")
  lib <- left_join(lib, totals, by = c("compartment", "gene", "treatment", "time", "replic"))
  lib$contribution <- safe_ratio(lib$tpm, lib$all_taxa_gene_tpm)
  states <- lib %>% group_by(across(all_of(c(tk, "gene", "treatment", "time")))) %>%
    summarise(mean_tpm = mean(tpm), min_tpm = min(tpm), max_tpm = max(tpm),
      n_technical = n(), n_detected = sum(tpm > 0),
      mean_taxon_tpm = mean(total_taxon_tpm),
      mean_within_taxon_fraction = if (all(is.na(within_taxon_fraction))) NA_real_ else mean(within_taxon_fraction, na.rm = TRUE),
      mean_contribution = if (all(is.na(contribution))) NA_real_ else mean(contribution, na.rm = TRUE), .groups = "drop")
  states$sample <- paste(states$treatment, states$time, sep = "_")
  states$annotation_mode <- mode
  list(libraries = lib, states = states, expanded = d)
}

environmental_states <- function(project_dir = ".") {
  d <- read.csv(file.path(project_dir, "data", "derived", "figure_6_correlations.csv"))
  env <- unique(d[c("sample", "predictor", "predictor_value")])
  if (anyDuplicated(env[c("sample", "predictor")])) stop("Inconsistent environmental predictor")
  # Figure 1 includes other available environmental means at matching states.
  e <- read.csv(file.path(project_dir, "data", "derived", "figure_1_environmental.csv"))
  e$sample <- paste(e$treatment, e$time, sep = "_")
  e <- e[e$time %in% c("0h", "72h"), ]
  names(e)[names(e) == "variable"] <- "predictor"
  names(e)[names(e) == "mean"] <- "predictor_value"
  e <- e[c("sample", "predictor", "predictor_value")]
  # Preserve Figure 6 values for its four predictors, even if Figure 1 differs.
  e <- anti_join(e, env[c("sample", "predictor")], by = c("sample", "predictor"))
  env <- bind_rows(env, e)
  if (anyDuplicated(env[c("sample", "predictor")])) stop("Environmental merge is not unique")
  env
}

fit_descriptive <- function(y, x) {
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]
  empty <- data.frame(n_states = length(y), slope = NA_real_, intercept = NA_real_, r2 = NA_real_,
    rho_spearman = NA_real_, p_independence_assumed = NA_real_,
    loo_slope_min = NA_real_, loo_slope_max = NA_real_, loo_direction_stable = FALSE)
  if (length(y) < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) return(empty)
  m <- lm(y ~ x)
  loo <- vapply(seq_along(y), function(i) {
    if (length(unique(x[-i])) < 2) return(NA_real_)
    unname(coef(lm(y[-i] ~ x[-i]))[2])
  }, numeric(1))
  empty$slope <- unname(coef(m)[2]); empty$intercept <- unname(coef(m)[1])
  empty$r2 <- summary(m)$r.squared
  empty$rho_spearman <- suppressWarnings(cor(y, x, method = "spearman"))
  empty$p_independence_assumed <- coef(summary(m))[2, 4]
  if (any(is.finite(loo))) {
    empty$loo_slope_min <- min(loo, na.rm = TRUE); empty$loo_slope_max <- max(loo, na.rm = TRUE)
    empty$loo_direction_stable <- all(is.finite(loo)) && all(sign(loo) == sign(empty$slope))
  }
  empty
}

is_target <- function(d) {
  (d$target_group %in% c("Rhodobacterales", "Pseudomonadales") & d$gene %in% c("pstA", "pstB", "pstC", "pstS", "pit")) |
    (d$target_group %in% c("Dinoflagellates", "Diatoms") & d$gene == "phoA")
}

indicator_associations <- function(states, env) {
  d <- states[is_target(states), ]
  d <- merge(d, env, by = "sample")
  keys <- c("compartment", "target_group", "taxon", "rank", "gene", "annotation_mode", "predictor")
  groups <- split(d, interaction(d[keys], drop = TRUE, lex.order = TRUE))
  results <- bind_rows(lapply(groups, function(z) {
    bind_rows(lapply(c("mean_tpm", "mean_within_taxon_fraction"), function(metric) {
      cbind(z[1, keys, drop = FALSE], metric = metric,
            fit_descriptive(z[[metric]], z$predictor_value))
    }))
  }))
  # One family: every taxon x gene x available predictor x normalization in this mode.
  results$q_BH_exploratory <- p.adjust(results$p_independence_assumed, "BH")
  results$formula <- "untransformed_response ~ untransformed_predictor"
  results$inference <- "descriptive only; 1 bag/treatment; dependent times"
  results
}

indicator_contrasts <- function(states) {
  d <- states[is_target(states), ]
  keys <- c("compartment", "target_group", "taxon", "rank", "gene", "annotation_mode", "time")
  r <- d[d$treatment == "R", ]; p <- d[d$treatment == "R+P", ]
  z <- merge(r, p, by = keys, suffixes = c("_R", "_RP"))
  z$tpm_R_over_RP <- safe_ratio(z$mean_tpm_R, z$mean_tpm_RP)
  z$within_R_over_RP <- safe_ratio(z$mean_within_taxon_fraction_R, z$mean_within_taxon_fraction_RP)
  z$interpretation <- "descriptive between-bag contrast; no biological replication"
  z
}

rank_candidates <- function(states, associations, expanded) {
  d <- states[is_target(states), ]
  keys <- c("compartment", "target_group", "taxon", "rank", "gene", "annotation_mode")
  groups <- split(d, interaction(d[keys], drop = TRUE, lex.order = TRUE))
  out <- bind_rows(lapply(groups, function(z) {
    a <- associations[associations$taxon == z$taxon[1] & associations$gene == z$gene[1] &
      associations$compartment == z$compartment[1], ]
    get <- function(pred, metric, field) {
      v <- a[[field]][a$predictor == pred & a$metric == metric]; if (length(v)) v[1] else NA
    }
    e <- expanded[expanded$compartment == z$compartment[1] & expanded$taxon == z$taxon[1] & expanded$gene == z$gene[1], ]
    contrast <- function(metric) {
      vapply(c("0h", "72h"), function(tm) {
        rr <- z[[metric]][z$treatment == "R" & z$time == tm]
        pp <- z[[metric]][z$treatment == "R+P" & z$time == tm]
        length(rr) == 1 && length(pp) == 1 && is.finite(rr) && is.finite(pp) && rr > pp
      }, logical(1))
    }
    poison <- z$rank[1] %in% c("family", "class", "order", "lineage_group", "unresolved")
    coherent <- isTRUE(get("PO4", "mean_within_taxon_fraction", "slope") < 0) &&
      isTRUE(get("PO4", "mean_within_taxon_fraction", "loo_direction_stable")) && all(contrast("mean_within_taxon_fraction"))
    detected <- sum(z$n_detected) / sum(z$n_technical)
    status <- if (poison) "group only; organism unresolved" else if (coherent && detected >= 0.8) "exploratory priority for validation" else "detected; insufficient coherent P-stress pattern"
    cbind(z[1, keys], data.frame(n_features = n_distinct(e$feature_id),
      feature_ids = paste(sort(unique(e$feature_id)), collapse = ";"),
      kos = paste(sort(unique(e$matched_kos)), collapse = ";"),
      detection_fraction = detected, mean_tpm = mean(z$mean_tpm),
      mean_gene_contribution = mean(z$mean_contribution, na.rm = TRUE),
      NP_slope_within = get("NP", "mean_within_taxon_fraction", "slope"),
      PO4_slope_within = get("PO4", "mean_within_taxon_fraction", "slope"),
      PO4_loo_stable = get("PO4", "mean_within_taxon_fraction", "loo_direction_stable"),
      PO4_q_exploratory = get("PO4", "mean_within_taxon_fraction", "q_BH_exploratory"),
      R_gt_RP_within_both_times = all(contrast("mean_within_taxon_fraction")),
      R_gt_RP_tpm_both_times = all(contrast("mean_tpm")),
      provisional_name = any(e$taxon_provisional), status = status,
      limitations = "1 bag/treatment; no paired DNA; orthology unvalidated; no species claim"))
  }))
  # Transparent priority order; high transcript contribution is only a last tie-breaker.
  out <- out %>% group_by(target_group) %>% arrange(
    match(status, c("exploratory priority for validation", "detected; insufficient coherent P-stress pattern", "group only; organism unresolved")),
    desc(detection_fraction), desc(R_gt_RP_within_both_times), desc(PO4_loo_stable),
    desc(mean_gene_contribution), taxon, gene, .by_group = TRUE) %>%
    mutate(rank_within_group = row_number()) %>% ungroup()
  out
}

summarise_candidate_taxa <- function(candidates) {
  keys <- c("compartment", "target_group", "taxon", "rank", "annotation_mode")
  groups <- split(candidates, interaction(candidates[keys], drop = TRUE, lex.order = TRUE))
  out <- bind_rows(lapply(groups, function(z) {
    bacterial <- z$compartment[1] == "Prokaryotes"
    primary <- if (bacterial) c("pstS", "pstB", "pstC", "pstA", "pit") else "phoA"
    z <- z[order(match(z$gene, primary)), ]
    marker <- z[1, ]
    resolved <- z$rank[1] == "genus"
    supported <- z$status == "exploratory priority for validation"
    n_pst <- sum(supported & z$gene %in% c("pstA", "pstB", "pstC", "pstS"))
    status <- if (!resolved) "organism unresolved" else if (isTRUE(supported[1]) && (!bacterial || marker$gene != "pit"))
      "exploratory candidate; validate experimentally" else "detected taxon; insufficient marker evidence"
    cbind(z[1, keys], data.frame(primary_gene = marker$gene,
      genes_detected = paste(z$gene, collapse = ";"), n_coherent_pst_subunits = n_pst,
      primary_detection_fraction = marker$detection_fraction,
      primary_mean_tpm = marker$mean_tpm, primary_gene_contribution = marker$mean_gene_contribution,
      primary_NP_slope_within = marker$NP_slope_within,
      primary_PO4_slope_within = marker$PO4_slope_within,
      primary_PO4_loo_stable = marker$PO4_loo_stable,
      primary_PO4_q_exploratory = marker$PO4_q_exploratory,
      primary_R_gt_RP_both_times = marker$R_gt_RP_within_both_times,
      feature_ids = paste(sort(unique(unlist(strsplit(z$feature_ids, ";", fixed = TRUE)))), collapse = ";"),
      provisional_name = any(z$provisional_name), status = status,
      limitation = "no species identification, no independent validation, no paired DNA"))
  }))
  out %>% group_by(target_group) %>% arrange(
    match(status, c("exploratory candidate; validate experimentally", "detected taxon; insufficient marker evidence", "organism unresolved")),
    desc(n_coherent_pst_subunits), desc(primary_detection_fraction), desc(primary_gene_contribution),
    taxon, .by_group = TRUE) %>% mutate(rank_within_group = row_number()) %>% ungroup()
}
