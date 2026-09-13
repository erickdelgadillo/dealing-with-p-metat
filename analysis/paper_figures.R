suppressPackageStartupMessages({
  library(arrow)
  library(cowplot)
  library(dplyr)
  library(ggh4x)
  library(ggplot2)
  library(grid)
  library(magick)
  library(patchwork)
  library(readr)
  library(scales)
})

project_dir <- normalizePath(
  if (file.exists(file.path("analysis", "paper_figures.R"))) "." else "..",
  winslash = "/",
  mustWork = TRUE
)
data_dir <- file.path(project_dir, "data", "derived")
figures_dir <- file.path(project_dir, "results", "figures")
supplementary_dir <- file.path(figures_dir, "supplementary")
dir.create(supplementary_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- function(name) file.path(data_dir, name)

taxon_colors <- c(
  "Alphaproteobacteria" = "blue4", "Rhodobacterales" = "deepskyblue4",
  "Puniceispirillales" = "royalblue2", "Rhizobiales" = "turquoise2",
  "Gammaproteobacteria" = "purple4", "Beggiatoales" = "mediumpurple2",
  "Burkholderiales" = "deeppink2", "Enterobacterales" = "magenta2",
  "Pseudomonadales" = "plum3", "SAR86" = "purple2",
  "Bacteroidota" = "peachpuff4", "Flavobacteriales" = "peachpuff1",
  "Firmicutes" = "lightcyan2", "Cyanobacteria" = "chartreuse2",
  "Actinobacteriota" = "darkred", "Planctomycetota" = "aquamarine2",
  "Archaea" = "darkorange2", "Other prokaryotes" = "grey33",
  "Spirotrichea" = "slateblue4", "Choreotrichida" = "slateblue",
  "Tintinnida" = "slateblue1", "Dinoflagellata" = "navy",
  "Gonyaulacales" = "turquoise3", "Alexandrium" = "turquoise1",
  "Gymnodiniales" = "lightskyblue", "Peridiniales" = "royalblue",
  "Suessiales" = "skyblue3", "Symbiodinium" = "skyblue1",
  "Syndiniales" = "lightsteelblue1", "Stramenopiles" = "darkgreen",
  "Ochrophyta" = "springgreen4", "Bacillariophyta" = "lightgreen",
  "Thalassiosira" = "green3", "Chaetoceros" = "darkolivegreen1",
  "Skeletonema" = "springgreen4", "Pseudo-nitzschia" = "green1",
  "Glaucophyta" = "orangered3", "Chlorophyta" = "aquamarine2",
  "Amoebozoa" = "lightpink4", "Hacrobia" = "khaki3",
  "Rhizaria" = "thistle3", "Other eukaryotes" = "grey33"
)

prok_taxa <- names(taxon_colors)[seq_len(match("Other prokaryotes", names(taxon_colors)))]
euk_taxa <- names(taxon_colors)[(match("Other prokaryotes", names(taxon_colors)) + 1):length(taxon_colors)]
subcategory_levels <- c(
  "P starvation response", "High affinity transporters",
  "Low affinity transporters", "Organic P hydrolysis",
  "Phosphonate metabolism", "Pi mobilization"
)
comparison_levels <- c("R vs C", "R+P vs C", "R+P vs R")
comparison_labels <- c("C vs R", "C vs R+P", "R vs R+P")

base_theme <- theme_bw(base_size = 8) +
  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines")
  )

# Figure 1 -----------------------------------------------------------------
environment <- read_csv(data_file("figure_1_environmental.csv"), show_col_types = FALSE) %>%
  mutate(
    treatment = factor(treatment, levels = c("C", "R", "R+P")),
    time = factor(time, levels = c("0h", "24h", "48h", "72h"))
  )

environment_1 <- environment %>%
  filter(!variable %in% c("PO4", "SiO4H4", "ab", "BP", "chla", "PP")) %>%
  mutate(variable = factor(
    variable,
    levels = c("DOC", "DON", "DIN", "NH4", "NO3", "NP"),
    labels = c("DOC\n(µmol L⁻¹)", "DON\n(µmol L⁻¹)", "DIN\n(µmol L⁻¹)",
               "NH₄⁺\n(µmol L⁻¹)", "NO₃⁻\n(µmol L⁻¹)", "N:P\nratio")
  ))

environment_2 <- environment %>%
  filter(!variable %in% c("DOC", "DON", "DIN", "NH4", "NO3", "NP")) %>%
  mutate(variable = factor(
    variable,
    levels = c("PO4", "SiO4H4", "ab", "BP", "chla", "PP"),
    labels = c("PO₄³⁻\n(µmol L⁻¹)", "Si(OH)₄\n(µmol L⁻¹)",
               "Bacterial\nabundance\n(cells mL⁻¹)",
               "Bacterial\nproduction\n(µg C L⁻¹ d⁻¹)",
               "Chlorophyll a\n(mg m⁻³)",
               "Primary\nproduction\n(µg C L⁻¹ h⁻¹)")
  ))

environment_plot <- function(data, colors, y_scales) {
  ggplot(data, aes(time, mean, fill = variable)) +
    geom_errorbar(aes(ymin = mean - errorstandar, ymax = mean + errorstandar), width = 0.4) +
    geom_col() +
    facet_nested(variable ~ treatment, scales = "free_y") +
    facetted_pos_scales(y = y_scales) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

environment_plot_1 <- environment_plot(
  environment_1,
  c("#006BC9", "#0079C3", "#0086BC", "#0094B5", "#00A1AE", "#00AEA8"),
  list(
    scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)),
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 5)),
    scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)),
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 1)),
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1.5)),
    scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20))
  )
)
environment_plot_2 <- environment_plot(
  environment_2,
  c("#00C99A", "#00D794", "#00E48D", "#00F286", "#00FF80", "#00FF72"),
  list(
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 1)),
    scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)),
    scale_y_continuous(limits = c(0, 3500000), breaks = seq(0, 3500000, 1000000)),
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 25)),
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1.5)),
    scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))
  )
)
environment_panels <- environment_plot_1 | environment_plot_2
map_grob <- rasterGrob(image_read(data_file("figure_1_map.png")), interpolate = TRUE)
figure_1 <- ggdraw() +
  draw_plot(map_grob, x = 0, y = 0.57, width = 1, height = 0.43) +
  draw_plot(environment_panels, x = 0, y = 0, width = 1, height = 0.58) +
  draw_label("b)", x = 0.01, y = 0.58, hjust = 0, fontface = "bold", size = 16)
ggsave(file.path(figures_dir, "Figure_1.png"), figure_1, width = 18, height = 24, units = "cm", dpi = 600, bg = "white")

# Figure 2 -----------------------------------------------------------------
taxonomy <- read_csv(data_file("figure_2_taxonomic_composition.csv"), show_col_types = FALSE) %>%
  mutate(
    domain = factor(domain, levels = c("Prokaryotes", "Eukaryotes")),
    treatment = factor(treatment, levels = c("C", "R", "R+P")),
    time = factor(time, levels = c("0h", "72h")),
    Final_Taxonomy = factor(Final_Taxonomy, levels = c(prok_taxa, euk_taxa))
  )
taxonomy_plot <- ggplot(taxonomy, aes(treatment, tpm, fill = Final_Taxonomy)) +
  geom_col(position = "fill") +
  facet_nested(~domain + time, scales = "free_x") +
  scale_fill_manual(values = taxon_colors, drop = TRUE) +
  scale_y_continuous(labels = percent) +
  labs(title = "a) Metatranscriptomic Taxonomy Distribution", x = NULL, y = NULL, fill = NULL) +
  base_theme +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")

nmds <- read_csv(data_file("figure_2_nmds_scores.csv"), show_col_types = FALSE) %>%
  mutate(treatment_time = factor(
    treatment_time,
    levels = c("Control 0h", "Control 72h", "River 0h", "River 72h", "River+P 0h", "River+P 72h")
  ))

ellipse_points <- function(data, npoints = 100) {
  if (nrow(data) < 3) return(NULL)
  covariance <- cov.wt(data[c("NMDS1", "NMDS2")], wt = rep(1 / nrow(data), nrow(data)))$cov
  center <- colMeans(data[c("NMDS1", "NMDS2")])
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  points <- t(center + t(circle %*% chol(covariance)))
  data.frame(NMDS1 = points[, 1], NMDS2 = points[, 2])
}

nmds_plot <- function(domain_name, colors, shapes) {
  points <- nmds %>% filter(domain == domain_name)
  ellipses <- bind_rows(lapply(split(points, points$treatment_time), function(group) {
    outline <- ellipse_points(group)
    if (is.null(outline)) return(NULL)
    mutate(outline, treatment_time = group$treatment_time[[1]])
  }))
  ggplot() +
    geom_polygon(data = ellipses, aes(NMDS1, NMDS2, fill = treatment_time), alpha = 0.2, show.legend = FALSE) +
    geom_point(data = points, aes(NMDS1, NMDS2, colour = treatment_time, shape = treatment_time), size = 3) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(title = domain_name, colour = NULL, shape = NULL, x = "NMDS1", y = "NMDS2") +
    base_theme +
    theme(plot.title = element_text(face = "bold", hjust = 0.04), legend.position = "top")
}

prok_nmds <- nmds_plot(
  "Prokaryotes",
  c("Control 0h" = "deepskyblue", "Control 72h" = "deepskyblue4", "River 0h" = "deeppink",
    "River 72h" = "deeppink4", "River+P 0h" = "aquamarine", "River+P 72h" = "aquamarine4"),
  c("Control 0h" = 19, "Control 72h" = 19, "River 0h" = 17, "River 72h" = 17, "River+P 0h" = 15, "River+P 72h" = 15)
)
euk_nmds <- nmds_plot(
  "Eukaryotes",
  c("Control 0h" = "slateblue", "Control 72h" = "slateblue4", "River 0h" = "green1",
    "River 72h" = "green4", "River+P 0h" = "orangered1", "River+P 72h" = "orangered4"),
  c("Control 0h" = 19, "Control 72h" = 19, "River 0h" = 17, "River 72h" = 17, "River+P 0h" = 15, "River+P 72h" = 15)
)
nmds_title <- ggdraw() + draw_label("b) nMDS Analysis of Prokaryotic and Eukaryotic Metatranscriptome", fontface = "bold", size = 10)
nmds_panel <- plot_grid(nmds_title, plot_grid(prok_nmds, euk_nmds, ncol = 2), ncol = 1, rel_heights = c(0.05, 1))
figure_2 <- plot_grid(taxonomy_plot, nmds_panel, ncol = 1, rel_heights = c(1.5, 1))
ggsave(file.path(figures_dir, "Figure_2.pdf"), figure_2, width = 18, height = 24, units = "cm", device = cairo_pdf)

# Figure 3 and Supplementary Figure 1 --------------------------------------
figure_3_totals <- read_csv(data_file("figure_3_totals.csv"), show_col_types = FALSE) %>%
  mutate(
    Domain = factor(Domain, levels = c("Prokaryotes", "Eukaryotes")),
    Subcategory = factor(Subcategory, levels = subcategory_levels),
    treatment = factor(treatment, levels = c("C", "R", "R+P")),
    time = factor(time, levels = c("0h", "72h"))
  )
figure_3_taxa <- read_csv(data_file("figure_3_taxonomic_contributions.csv"), show_col_types = FALSE) %>%
  mutate(
    Domain = factor(Domain, levels = c("Prokaryotes", "Eukaryotes")),
    Subcategory = factor(Subcategory, levels = subcategory_levels),
    treatment = factor(treatment, levels = c("C", "R", "R+P")),
    time = factor(time, levels = c("0h", "72h")),
    Final_Taxonomy = factor(Final_Taxonomy, levels = c(prok_taxa, euk_taxa))
  )
figure_3 <- ggplot(figure_3_totals, aes(treatment, mean)) +
  geom_errorbar(aes(ymin = mean - errorstandar, ymax = mean + errorstandar), width = 0.4) +
  geom_col(data = figure_3_taxa, aes(y = tpm, fill = Final_Taxonomy)) +
  geom_text(aes(y = mean + errorstandar * 0.85, label = label), vjust = -0.5) +
  facet_nested(Subcategory ~ Domain + time, scales = "free_y") +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 2000), breaks = seq(0, 2000, 400)),
    scale_y_continuous(limits = c(0, 7000), breaks = seq(0, 7000, 1400)),
    scale_y_continuous(limits = c(0, 8500), breaks = seq(0, 8000, 1600)),
    scale_y_continuous(limits = c(0, 1250), breaks = seq(0, 1200, 240)),
    scale_y_continuous(limits = c(0, 210), breaks = seq(0, 200, 40)),
    scale_y_continuous(limits = c(0, 4000), breaks = seq(0, 4000, 800))
  )) +
  scale_fill_manual(values = taxon_colors, drop = TRUE) +
  labs(x = NULL, y = "TPM", fill = NULL) +
  base_theme +
  theme(legend.position = "right")
ggsave(file.path(figures_dir, "Figure_3.pdf"), figure_3, width = 18, height = 24, units = "cm", device = cairo_pdf)

subcategory_colors <- c(
  "P starvation\nresponse" = "aquamarine1", "High affinity\ntransporters" = "plum3",
  "Low affinity\ntransporters" = "slateblue1", "Organic P\nhydrolysis" = "turquoise3",
  "Phosphonate\nmetabolism" = "aquamarine3", "Pi\nmobilization" = "lightsteelblue1"
)
supplementary_1 <- read_csv(data_file("supplementary_figure_1_proportions.csv"), show_col_types = FALSE) %>%
  mutate(
    Domain = factor(Domain, levels = c("Prokaryotes", "Eukaryotes")),
    treatment = factor(treatment, levels = c("C", "R", "R+P")),
    time = factor(time, levels = c("0h", "72h")),
    Subcategory = factor(Subcategory, levels = names(subcategory_colors))
  )
supplementary_figure_1 <- ggplot(supplementary_1, aes(treatment, Proportion, fill = Subcategory)) +
  geom_col(position = "fill") +
  geom_text(aes(label = ifelse(Proportion > 1, sprintf("%.2f", Proportion), "")), position = position_fill(vjust = 0.5), size = 2) +
  facet_nested(~Domain + time) +
  scale_fill_manual(values = subcategory_colors) +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = NULL, fill = NULL) +
  base_theme
ggsave(file.path(supplementary_dir, "Supplementary_Figure_1.png"), supplementary_figure_1, width = 18, height = 14, units = "cm", dpi = 600)

# DGE figures ---------------------------------------------------------------
overview <- read_parquet(data_file("supplementary_figure_2_dge_overview.parquet")) %>%
  mutate(
    significant = FDR <= 0.05,
    display_comparison = factor(comparison, levels = comparison_levels, labels = comparison_labels),
    time = factor(time, levels = c("0h", "72h")),
    transparency = ifelse(PValue <= 0.05 & between(logFC, -1, 1), 0.1,
                          ifelse(PValue >= 0.05 & (logFC <= 1 | logFC >= -1), 0.1, 1))
  )

overview_plot <- function(domain_name, show_fdr = TRUE) {
  data <- overview %>% filter(domain == domain_name)
  plot <- ggplot(data, aes(logFC, logCPM, colour = FDR, shape = significant)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "dimgrey", alpha = 0.3) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "snow3", alpha = 0.3) +
    geom_point(alpha = data$transparency) +
    facet_nested(time ~ domain + display_comparison) +
    scale_colour_viridis_c("FDR") +
    guides(shape = "none") +
    labs(x = "logFC", y = "logCPM") +
    base_theme
  if (!show_fdr) plot <- plot + guides(colour = "none")
  plot
}
supplementary_figure_2 <- overview_plot("Prokaryotes") / overview_plot("Eukaryotes", FALSE)
ggsave(file.path(supplementary_dir, "Supplementary_Figure_2.png"), supplementary_figure_2, width = 18, height = 24, units = "cm", dpi = 600)

plot_selected_dge <- function(path, domain_name, id_col, genes, threshold) {
  data <- read_parquet(data_file(path)) %>%
    filter(Category == "Phosphorus metabolism", gene_protein %in% genes) %>%
    mutate(
      Time = factor(Time, levels = c("0h", "72h")),
      display_comparison = factor(comparation, levels = comparison_levels, labels = comparison_labels),
      Subcategory = factor(Subcategory, levels = subcategory_levels),
      Final_Taxonomy = factor(Final_Taxonomy, levels = if (domain_name == "Prokaryotes") prok_taxa else euk_taxa),
      transparency = ifelse(logCPM <= 1, 0.1,
        ifelse(PValue < 0.05 & between(logFC, -threshold, threshold), 0.1,
          ifelse(PValue > 0.05 & (logFC < threshold | logFC > -threshold), 0.05, 1)))
    )
  data[[id_col]] <- forcats::fct_inorder(data[[id_col]])
  ggplot(data, aes(x = logFC, y = .data[[id_col]], size = logCPM, fill = Final_Taxonomy, alpha = transparency)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "dimgrey", alpha = 0.3) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "snow3", alpha = 0.3) +
    geom_point(shape = 21) +
    geom_vline(xintercept = 0) +
    facet_nested(Subcategory + gene_protein ~ Time + display_comparison, space = "free", strip = strip_nested(size = "variable")) +
    scale_fill_manual(values = taxon_colors, drop = TRUE) +
    scale_size(range = c(1, 4), breaks = c(1, 3, 5)) +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_alpha_identity() +
    labs(x = "logFC", y = NULL, fill = NULL, size = "logCPM") +
    base_theme +
    theme(axis.text.y = element_blank(), legend.position = "top")
}

figure_4 <- plot_selected_dge(
  "figure_4_prokaryote_dge.parquet", "Prokaryotes", "orf",
  c("pstA", "pstB", "pstC", "pstS", "nptA", "pit", "acp5", "phoA", "phoD",
    "phoB", "phoR", "phoU", "phnC", "phnD", "phnE", "ppA", "ppK1"),
  2
)
ggsave(file.path(figures_dir, "Figure_4.pdf"), figure_4, width = 18, height = 24, units = "cm", device = cairo_pdf)

figure_5 <- plot_selected_dge(
  "figure_5_eukaryote_dge.parquet", "Eukaryotes", "geneid",
  c("pho84", "ugpQ", "nptA", "pit", "acp5", "phoA", "pho-3", "phoD",
    "pho80", "pho81", "phoR", "pepM", "phnX", "phnY", "ppd", "ppA", "EPT1", "PCYT2"),
  1.9
)
ggsave(file.path(figures_dir, "Figure_5.pdf"), figure_5, width = 18, height = 24, units = "cm", device = cairo_pdf)

# Figure 6 -----------------------------------------------------------------
correlations <- read_csv(data_file("figure_6_correlations.csv"), show_col_types = FALSE)
panel_data <- function(panel, response = NULL, taxon = NULL) {
  data <- correlations %>% filter(.data$panel == .env$panel)
  if (!is.null(response)) data <- data %>% filter(.data$response == .env$response)
  if (!is.null(taxon)) data <- data %>% filter(.data$Final_Taxonomy == .env$taxon)
  data
}
good_p <- function(model, n = 6) coef(summary(model))[2, 4] * sqrt(n / 100)

high <- panel_data("High_affinity")
fit_power <- function(data) {
  model <- lm(log(response_value) ~ log(predictor_value), data = data)
  list(a = exp(coef(model)[1]), b = coef(model)[2], r2 = summary(model)$r.squared, p = good_p(model))
}
high_fits <- lapply(split(high, high$Final_Taxonomy), fit_power)
high_plot <- ggplot(high, aes(predictor_value, response_value, colour = Final_Taxonomy)) +
  geom_point(size = 3) +
  stat_smooth(method = "nls", formula = y ~ a * x^b, method.args = list(start = list(a = 1, b = -1)), se = FALSE, linetype = "dashed") +
  scale_colour_manual(values = c("Pseudomonadales" = "plum3", "Rhodobacterales" = "deepskyblue4")) +
  annotate("text", x = 0.15, y = 38,
           label = "Rhodobacterales: R² = 0.893   p = 0.001",
           hjust = 0, colour = "deepskyblue4", fontface = "bold", size = 3) +
  annotate("text", x = 0.15, y = 34,
           label = "Pseudomonadales: R² = 0.750   p = 0.006",
           hjust = 0, colour = "plum3", fontface = "bold", size = 3) +
  labs(title = "a) High affinity transport vs PO₄³⁻", x = "PO₄³⁻ (µmol L⁻¹)", y = "pst:pit ratio", colour = NULL) +
  base_theme + theme(legend.position = "none")

low <- panel_data("Low_affinity", "nptA", "Diatom")
low_fit <- fit_power(low)
low_plot <- ggplot(low, aes(predictor_value, response_value)) +
  geom_point(colour = "green4", size = 3) +
  stat_function(fun = function(x) low_fit$a * x^low_fit$b, colour = "green4", linetype = "dotted") +
  annotate("text", x = 0.18, y = 15000,
           label = "Diatoms: R² = 0.702   p = 0.009",
           hjust = 0, colour = "green4", fontface = "bold", size = 3) +
  labs(title = "b) Low affinity transport vs PO₄³⁻", x = "PO₄³⁻ (µmol L⁻¹)", y = "nptA (TPM)") + base_theme

pi_data <- panel_data("P_mobilization", "ppA", "Diatom")
pi_model <- lm(log(response_value) ~ predictor_value, data = pi_data)
pi_plot <- ggplot(pi_data, aes(predictor_value, response_value)) +
  geom_point(colour = "green4", size = 3) +
  stat_smooth(method = "nls", formula = y ~ a * exp(b * x), method.args = list(start = list(a = 10, b = 0.5)), se = FALSE, colour = "green4", linetype = "dashed") +
  annotate("text", x = 0.6, y = 3000,
           label = "Diatoms: R² = 0.72   p = 0.008",
           hjust = 0, colour = "green4", fontface = "bold", size = 3) +
  scale_y_log10(labels = comma) +
  labs(title = "c) Internal P mobilization vs DIN", x = "DIN (µmol L⁻¹)", y = "ppA (TPM)") + base_theme

organic <- panel_data("Organic_P")
organic_plot <- ggplot(organic, aes(predictor_value, response_value, colour = Final_Taxonomy)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_colour_manual(values = c("Pseudomonadales" = "plum3", "Diatom" = "green4")) +
  annotate("text", x = 5, y = 2000,
           label = "Diatoms: R² = 0.560   p = 0.021",
           hjust = 0, colour = "green4", fontface = "bold", size = 3) +
  annotate("text", x = 5, y = 1700,
           label = "Pseudomonadales: R² = 0.483   p = 0.031",
           hjust = 0, colour = "plum3", fontface = "bold", size = 3) +
  scale_y_continuous(labels = comma) +
  labs(title = "d) Organic P utilization vs N:P", x = "N:P ratio", y = "phoA (TPM)", colour = NULL) +
  base_theme + theme(legend.position = "none")

organic_phosphonate <- panel_data("Organic_p_and_Phosphonate")
organic_phosphonate_plot <- ggplot(organic_phosphonate, aes(predictor_value, response_value, colour = response)) +
  geom_point(size = 3) +
  geom_smooth(method = "nls", formula = y ~ a * exp(b * x), method.args = list(start = list(a = 100, b = -0.1)), se = FALSE, linetype = "dashed") +
  scale_colour_manual(values = c("acp5" = "royalblue4", "pho_3" = "royalblue2", "phoA" = "lightskyblue4", "pepM" = "lightskyblue3")) +
  annotate("text", x = 6.6, y = 925, label = "Dinoflagellates:",
           hjust = 0, colour = "black", fontface = "bold", size = 3) +
  annotate("text", x = 6.6, y = 845, label = "acp5:   R² = 0.942   p = 0.001",
           hjust = 0, colour = "royalblue4", fontface = "bold", size = 3) +
  annotate("text", x = 6.6, y = 765, label = "pho-3: R² = 0.920   p = 0.001",
           hjust = 0, colour = "royalblue2", fontface = "bold", size = 3) +
  annotate("text", x = 6.6, y = 685, label = "phoA:  R² = 0.775   p = 0.005",
           hjust = 0, colour = "lightskyblue4", fontface = "bold", size = 3) +
  annotate("text", x = 6.6, y = 605, label = "pepM: R² = 0.806   p = 0.004",
           hjust = 0, colour = "lightskyblue3", fontface = "bold", size = 3) +
  scale_y_continuous(labels = comma) +
  labs(title = "e) Organic P and phosphonate utilization vs PP", x = "Primary production (µg C L⁻¹ h⁻¹)", y = "TPM", colour = NULL) + base_theme

phosphonate <- panel_data("Phosphonate", "phnC")
phosphonate_plot <- ggplot(phosphonate, aes(predictor_value, response_value)) +
  geom_point(colour = "deepskyblue4", size = 3) + geom_smooth(method = "lm", se = FALSE, colour = "deepskyblue4", linetype = "dotted") +
  annotate("text", x = 4.5, y = 16.5,
           label = "Rhodobacterales: R² = 0.620   p = 0.015",
           hjust = 0, colour = "deepskyblue4", fontface = "bold", size = 3) +
  labs(title = "f) Phosphonate utilization vs PP", x = "Primary production (µg C L⁻¹ h⁻¹)", y = "phnC (TPM)") + base_theme

figure_6 <- (high_plot | low_plot) / (pi_plot | organic_plot) / (organic_phosphonate_plot | phosphonate_plot)
ggsave(file.path(figures_dir, "Figure_6.png"), figure_6, width = 18, height = 24, units = "cm", dpi = 600)

message("Generated Figures 1-6 and Supplementary Figures 1-2 from data/derived.")
