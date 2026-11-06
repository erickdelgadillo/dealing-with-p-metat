# Run from the repository root with Rscript --vanilla.
# Synthetic observations below are TEST FIXTURES ONLY, never analysis inputs.
source("analysis/phosphorus_indicators.R")
source("analysis/prepare_phosphorus_inputs.R")
source("analysis/audit_historical_regressions.R")

expect_error <- function(expr) {
  caught <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  stopifnot(caught)
}
near <- function(x, y, tolerance = 1e-9) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
}

# Exact KO tokens: reject substrings, keep independent KO assignments distinct.
stopifnot(identical(ko_tokens(c("ko:K01077; K03306,K01077", NA, "K010770 fooK02040")),
                    c("K01077", "K03306")))
# TEST FIXTURE: ambiguous annotation is excluded only in strict mode.
synthetic_annotation <- data.frame(feature_id = c("test_single", "test_multiple"),
  matched_genes = c("phoA", "phoA;pit"), ambiguous_ko = c(FALSE, TRUE))
stopifnot(nrow(expand_genes(synthetic_annotation, "strict")) == 1L,
          nrow(expand_genes(synthetic_annotation, "inclusive")) == 3L,
          !any(expand_genes(synthetic_annotation, "strict")$feature_id == "test_multiple"))

# TEST FIXTURE: shifted eukaryote headings must never produce a species claim
# or promote arbitrary values from the Genus/Species columns to a genus.
synthetic_taxonomy <- data.frame(
  Domain = rep("Eukaryota", 4), Phylum = rep("test_lineage", 4),
  Class = c("Dinoflagellata", "Dinoflagellata", "test_class", "test_class"),
  Order = c("test_order", "test_order", "Bacillariophyta", "Bacillariophyta"),
  Family = c("test_family", "Suessiales", "test_family", "test_family"),
  Genus = rep("Do_not_accept_this_column", 4),
  Species = c("unclassified", "unclassified", "Leptocylindrus", "Unreviewed genus"),
  Final_Taxonomy = c("Dinoflagellata", "Suessiales", "Bacillariophyta", "Bacillariophyta"))
assigned <- assign_taxa(synthetic_taxonomy, "Eukaryotes")
stopifnot(identical(assigned$rank, c("class", "order", "genus", "lineage_group")),
          assigned$taxon[3] == "Leptocylindrus",
          !any(assigned$rank == "species"),
          !any(grepl("Do_not_accept|Unreviewed", assigned$taxon)),
          all(assigned$target_group == c("Dinoflagellates", "Dinoflagellates", "Diatoms", "Diatoms")))

stopifnot(near(safe_ratio(c(2, 0), c(4, 4)), c(0.5, 0)),
          all(is.na(safe_ratio(c(2, 0, 1), c(0, 0, NA)))))

# Real committed snapshot: read-only checks of the experimental design.
x <- read_indicator_inputs()
validate_indicator_inputs(x)
stopifnot(nrow(x$libraries) == 35L,
          sum(x$libraries$compartment == "Prokaryotes") == 18L,
          sum(x$libraries$compartment == "Eukaryotes") == 17L)
bad <- x
bad$transcripts <- rbind(bad$transcripts, bad$transcripts[1, ])
expect_error(validate_indicator_inputs(bad))
bad <- x
bad$transcripts$tpm[1] <- -1
expect_error(validate_indicator_inputs(bad))

a <- aggregate_indicators(x, "strict")
stopifnot(!any(a$expanded$ambiguous_ko),
          !any(a$libraries$compartment == "Eukaryotes" & a$libraries$treatment == "R+P" &
                     a$libraries$time == "0h" & a$libraries$replic == "R3"),
          all(a$states$n_technical[a$states$compartment == "Eukaryotes" &
                    a$states$treatment == "R+P" & a$states$time == "0h"] == 2L),
          all(a$states$n_technical[a$states$compartment == "Prokaryotes"] == 3L))

# All-taxa denominators must include taxa outside the four requested groups.
library_gene_keys <- c("compartment", "gene", "treatment", "time", "replic")
expected_gene_totals <- a$expanded %>% group_by(across(all_of(library_gene_keys))) %>%
  summarise(expected_total = sum(tpm), .groups = "drop")
denominators <- distinct(a$libraries, across(all_of(library_gene_keys)), all_taxa_gene_tpm)
denominators <- inner_join(denominators, expected_gene_totals, by = library_gene_keys)
stopifnot(nrow(denominators) == nrow(expected_gene_totals),
          near(denominators$all_taxa_gene_tpm, denominators$expected_total),
          any(a$expanded$target_group == "Other"))
shares <- a$libraries %>% group_by(across(all_of(library_gene_keys))) %>%
  summarise(total = first(all_taxa_gene_tpm), share_sum = sum(contribution), .groups = "drop")
stopifnot(all(abs(shares$share_sum[shares$total > 0] - 1) < 1e-10))

# Taxon denominator is all RNA for the assigned taxon, not only P-target genes.
taxon_keys <- c("compartment", "taxon", "rank", "target_group", "treatment", "time", "replic")
within <- left_join(a$libraries, x$totals, by = taxon_keys, suffix = c("", "_source"))
stopifnot(all(abs(within$total_taxon_tpm - within$total_taxon_tpm_source) < 1e-9 |
                (within$total_taxon_tpm == 0 & is.na(within$total_taxon_tpm_source))),
          all(is.na(within$within_taxon_fraction[within$total_taxon_tpm == 0])),
          all(within$within_taxon_fraction[within$total_taxon_tpm > 0] <= 1 + 1e-10))

# Observed sparse omissions contribute operational zeros before averaging.
stopifnot(any(a$libraries$not_detected),
          all(a$libraries$tpm[a$libraries$not_detected] == 0))
state_keys <- c("compartment", "taxon", "rank", "target_group", "gene", "treatment", "time")
means <- a$libraries %>% group_by(across(all_of(state_keys))) %>%
  summarise(expected_mean = sum(tpm) / n(), expected_n = n(), .groups = "drop")
means <- inner_join(a$states, means, by = state_keys)
stopifnot(near(means$mean_tpm, means$expected_mean),
          identical(means$n_technical, means$expected_n))

env <- environmental_states()
stopifnot(!anyDuplicated(env[c("sample", "predictor")]))
associations <- indicator_associations(a$states, env)
candidates <- rank_candidates(a$states, associations, a$expanded)
stopifnot(nrow(candidates) > 0L,
          all(c("PO4_slope_within", "PO4_loo_stable", "status", "feature_ids") %in% names(candidates)),
          !any(grepl("species", candidates$rank)),
          all(nzchar(candidates$feature_ids)),
          !any(candidates$status[candidates$rank != "genus"] == "exploratory priority for validation"))
# A real genus with a consistent P contrast must remain eligible after tibble
# extraction; regression/get() bugs previously suppressed every coherent status.
lepto <- candidates[candidates$taxon == "Leptocylindrus" & candidates$gene == "phoA", ]
stopifnot(nrow(lepto) == 1L, isTRUE(lepto$PO4_loo_stable[[1]]),
          isTRUE(lepto$R_gt_RP_within_both_times[[1]]),
          lepto$status[[1]] == "exploratory priority for validation")

# Independently known historical values anchor transformation/model choices.
historical <- audit_historical_regressions()
stopifnot(nrow(historical) == 12L, all(historical$n == 6L),
          near(historical$q_BH_exploratory, p.adjust(historical$p_independence_assumed, "BH")),
          near(historical$p_historical_factor,
               historical$p_independence_assumed * sqrt(6 / 100)))
h <- historical[historical$panel == "Organic_P" & historical$taxon == "Diatom", ]
stopifnot(near(h$r2, 0.5598076, 1e-6), near(h$p_independence_assumed, 0.087120465, 1e-6))
h <- historical[historical$panel == "High_affinity" & historical$taxon == "Rhodobacterales", ]
stopifnot(near(h$r2, 0.8929074, 1e-6), near(h$p_independence_assumed, 0.004464306, 1e-6))
h <- historical[historical$panel == "Organic_p_and_Phosphonate" & historical$response == "phoA", ]
stopifnot(near(h$r2, 0.4852289, 1e-6), near(h$p_independence_assumed, 0.124125653, 1e-6))
message("All phosphorus-indicator checks passed; synthetic fixtures were used only in tests.")
