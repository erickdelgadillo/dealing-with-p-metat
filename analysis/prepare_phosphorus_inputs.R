# Extract a compact, auditable snapshot. Original files are read only.
suppressPackageStartupMessages(library(dplyr))

ko_map <- c(K02038 = "pstA", K02036 = "pstB", K02037 = "pstC",
            K02040 = "pstS", K03306 = "pit", K14640 = "pit", K01077 = "phoA")

ko_tokens <- function(x) {
  z <- unlist(strsplit(paste(x[!is.na(x)], collapse = ";"), "[,;[:space:]]+"))
  unique(sub("^ko:", "", z[grepl("^(ko:)?K[0-9]{5}$", z)]))
}

clean_taxon <- function(x) {
  !is.na(x) & nzchar(trimws(x)) & !grepl("uncl|unclassified|^No |^unknown$", x, ignore.case = TRUE)
}

assign_taxa <- function(d, compartment) {
  for (nm in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Final_Taxonomy"))
    d[[nm]] <- trimws(d[[nm]])
  d$compartment <- compartment
  d$target_group <- "Other"
  if (compartment == "Prokaryotes") {
    d$taxon <- "Unresolved prokaryote"; d$rank <- "unresolved"
    for (rank in c("Order", "Family", "Genus")) {
      ok <- clean_taxon(d[[rank]])
      d$taxon[ok] <- d[[rank]][ok]; d$rank[ok] <- tolower(rank)
    }
    for (g in c("Rhodobacterales", "Pseudomonadales"))
      d$target_group[!is.na(d$Order) & d$Order == g] <- g
    d$taxon_provisional <- grepl("[0-9]|^[A-Z]{3,}", d$taxon)
  } else {
    # The source's eukaryote rank headers are shifted. Never take Genus/Species
    # literally. Only these reviewed lineage labels are accepted in this study.
    d$taxon <- ifelse(clean_taxon(d$Final_Taxonomy), d$Final_Taxonomy, "Unresolved eukaryote")
    d$rank <- "unresolved"
    dino <- !is.na(d$Class) & d$Class == "Dinoflagellata"
    diatom <- !is.na(d$Order) & d$Order == "Bacillariophyta"
    d$target_group[dino] <- "Dinoflagellates"
    d$target_group[diatom] <- "Diatoms"
    d$taxon[dino] <- "Dinophyceae (genus unresolved)"; d$rank[dino] <- "class"
    suess <- dino & !is.na(d$Family) & d$Family == "Suessiales"
    d$taxon[suess] <- "Suessiales (genus unresolved)"; d$rank[suess] <- "order"
    d$taxon[diatom] <- "Bacillariophyta (genus unresolved)"; d$rank[diatom] <- "lineage_group"
    lepto <- diatom & !is.na(d$Species) & d$Species == "Leptocylindrus"
    d$taxon[lepto] <- "Leptocylindrus"; d$rank[lepto] <- "genus"
    d$taxon_provisional <- FALSE
  }
  d
}

sha256 <- function(path) {
  out <- system2("sha256sum", shQuote(normalizePath(path)), stdout = TRUE)
  if (!length(out) || !is.null(attr(out, "status"))) stop("Cannot hash ", path)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

read_selected_tsv <- function(path, ids = NULL, eggnog = FALSE) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  selected <- character(); header <- NULL
  repeat {
    lines <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    if (is.null(header)) {
      h <- if (eggnog) which(startsWith(lines, "#query\t")) else 1L
      if (!length(h)) stop("Missing annotation header: ", path)
      header <- sub("^#", "", lines[h[1]])
    }
    keys <- sub("\t.*$", "", lines)
    if (is.null(ids)) {
      # Fast candidate scan followed by exact validation of the KO field below.
      hit <- grepl(paste(paste0("ko:", names(ko_map)), collapse = "|"), lines)
    } else hit <- keys %in% ids
    selected <- c(selected, lines[hit])
  }
  if (!length(selected)) stop("No original annotations matched: ", path)
  z <- read.delim(text = paste(c(header, selected), collapse = "\n"), quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(ids)) z <- z[vapply(z$KEGG_ko, function(k) any(ko_tokens(k) %in% names(ko_map)), logical(1)), ]
  z
}

extract_original_annotations <- function(source_root, project_dir = ".") {
  out <- file.path(project_dir, "data", "phosphorus_indicators")
  d <- read.csv(file.path(out, "transcripts.csv"), stringsAsFactors = FALSE)
  sources <- read.csv(file.path(out, "sources.csv"), stringsAsFactors = FALSE)
  hits <- tax <- list()
  for (comp in c("Prokaryotes", "Eukaryotes")) {
    ids <- unique(d$feature_id[d$compartment == comp])
    rel <- file.path("raw", tolower(comp), "all_samples_emapper_annotations.tsv")
    p <- file.path(source_root, rel)
    z <- read_selected_tsv(p, ids, eggnog = TRUE)
    names(z)[names(z) == "query"] <- "feature_id"
    z$compartment <- comp
    if (anyDuplicated(z$feature_id)) stop("Duplicated original functional annotation")
    if (!all(ids %in% z$feature_id)) stop("Missing original functional annotations")
    hits[[comp]] <- z
    sources <- bind_rows(sources, data.frame(source_relative_path = rel, sha256 = sha256(p), rows_wp2 = nrow(z), role = "selected_original_eggnog_hits"))
    rel <- file.path("raw", tolower(comp), if (comp == "Prokaryotes") "taxonomy_classification_gtdb.tsv" else "taxonomy_classification.tsv.gz")
    p <- file.path(source_root, rel)
    z <- read_selected_tsv(p, ids)
    names(z)[names(z) == "orf"] <- "feature_id"
    z$compartment <- comp
    if (anyDuplicated(z$feature_id)) stop("Duplicated original taxonomy")
    # Unclassified transcripts can be absent from the original classifier table.
    tax[[comp]] <- z
    sources <- bind_rows(sources, data.frame(source_relative_path = rel, sha256 = sha256(p), rows_wp2 = nrow(z), role = "selected_original_taxonomy"))
  }
  write.csv(bind_rows(hits), file.path(out, "original_annotation_hits.csv"), row.names = FALSE, na = "")
  write.csv(bind_rows(tax), file.path(out, "original_taxonomy.csv"), row.names = FALSE, na = "")
  sources <- distinct(sources, source_relative_path, .keep_all = TRUE)
  write.csv(sources, file.path(out, "sources.csv"), row.names = FALSE, na = "")
  files <- sort(list.files(out, pattern = "\\.csv$", full.names = TRUE))
  writeLines(paste(vapply(files, sha256, ""), basename(files), sep = "  "), file.path(out, "SHA256SUMS"))
}

prepare_phosphorus_inputs <- function(source_root, project_dir = ".") {
  out <- file.path(project_dir, "data", "phosphorus_indicators")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  all_targets <- all_totals <- all_qc <- sources <- metadata <- list()
  for (comp in c("Prokaryotes", "Eukaryotes")) {
    prefix <- if (comp == "Prokaryotes") "prok" else "euk"
    subdir <- tolower(comp)
    rel <- file.path("processed", subdir, paste0(prefix, "_tpms_annotated.parquet"))
    path <- file.path(source_root, rel)
    d <- arrow::read_parquet(path)
    d <- d[d$Workpackage == "WP2", ]
    id <- if (comp == "Prokaryotes") "orf" else "geneid"
    names(d)[names(d) == id] <- "feature_id"
    key <- c("feature_id", "treatment", "time", "replic")
    if (anyDuplicated(d[key])) stop("Duplicated expression key in ", rel)
    if (any(!is.finite(d$tpm) | d$tpm < 0)) stop("Invalid TPM in ", rel)
    d <- assign_taxa(d, comp)
    all_qc[[comp]] <- d %>% group_by(compartment, treatment, time, replic) %>%
      summarise(total_tpm = sum(tpm), n_features = n(), n_zero = sum(tpm == 0), .groups = "drop")
    if (any(abs(all_qc[[comp]]$total_tpm - 1e6) > 1)) stop("TPM library sum differs from 1e6")
    totals <- d %>% group_by(compartment, taxon, rank, target_group, treatment, time, replic) %>%
      summarise(total_taxon_tpm = sum(tpm), .groups = "drop")
    # The processed tables retain only three KO slots. Consult the full raw
    # EggNOG KO list so a fourth or later match is neither missed nor hidden.
    ann <- read_selected_tsv(file.path(source_root, "raw", subdir, "all_samples_emapper_annotations.tsv"), eggnog = TRUE)
    ann <- ann[c("query", "KEGG_ko")]
    names(ann) <- c("feature_id", "original_kos")
    if (anyDuplicated(ann$feature_id)) stop("Duplicate functional annotation key")
    tokens <- lapply(ann$original_kos, ko_tokens)
    ann$matched_kos <- vapply(tokens, function(z) paste(intersect(z, names(ko_map)), collapse = ";"), "")
    ann$matched_genes <- vapply(tokens, function(z) paste(unique(unname(ko_map[intersect(z, names(ko_map))])), collapse = ";"), "")
    ann$ambiguous_ko <- lengths(tokens) > 1L
    ann$pit_legacy_K14640 <- vapply(tokens, function(z) "K14640" %in% z, logical(1))
    ann <- ann[nzchar(ann$matched_genes), ]
    target <- inner_join(d, ann, by = "feature_id")
    all_totals[[comp]] <- semi_join(totals, distinct(target, compartment, taxon, rank, target_group),
      by = c("compartment", "taxon", "rank", "target_group"))
    all_targets[[comp]] <- target %>% select(compartment, feature_id, treatment, time, replic,
      tpm, nReplics, Domain, Phylum, Class, Order, Family, Genus, Species, Final_Taxonomy,
      KEGG_ko, KEGG_ko2, KEGG_ko3, original_kos, PFAMs, matched_kos, matched_genes, ambiguous_ko,
      pit_legacy_K14640, taxon, rank, target_group, taxon_provisional)
    sources[[rel]] <- data.frame(source_relative_path = rel, sha256 = sha256(path),
      rows_wp2 = nrow(d), role = "expression_annotation_taxonomy")
    mrel <- file.path("metadata", subdir, paste0("INTERES_", if (prefix == "prok") "Prok" else "Euk", "_samples_tags.csv"))
    mp <- file.path(source_root, mrel)
    m <- read.csv(mp, check.names = FALSE, fileEncoding = "UTF-8-BOM")
    m <- m[grepl("^WP2_", m$group), ]
    m$compartment <- comp
    metadata[[comp]] <- m
    sources[[mrel]] <- data.frame(source_relative_path = mrel, sha256 = sha256(mp), rows_wp2 = nrow(m), role = "sample_labels")
    rm(d); gc(verbose = FALSE)
  }
  write.csv(bind_rows(all_targets), file.path(out, "transcripts.csv"), row.names = FALSE, na = "")
  write.csv(bind_rows(all_totals), file.path(out, "taxon_totals.csv"), row.names = FALSE, na = "")
  write.csv(bind_rows(all_qc), file.path(out, "library_qc.csv"), row.names = FALSE, na = "")
  write.csv(bind_rows(sources), file.path(out, "sources.csv"), row.names = FALSE, na = "")
  write.csv(bind_rows(metadata), file.path(out, "sample_metadata.csv"), row.names = FALSE, na = "")
  files <- sort(list.files(out, pattern = "\\.csv$", full.names = TRUE))
  writeLines(paste(vapply(files, sha256, ""), basename(files), sep = "  "), file.path(out, "SHA256SUMS"))
  extract_original_annotations(source_root, project_dir)
  invisible(out)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript analysis/prepare_phosphorus_inputs.R SOURCE_DATA_ROOT")
  prepare_phosphorus_inputs(args[1])
}
