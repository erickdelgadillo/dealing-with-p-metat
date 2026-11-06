# Recompute metrics; retain the published artwork and script unchanged.
audit_historical_regressions <- function(project_dir = ".") {
  d <- read.csv(file.path(project_dir, "data", "derived", "figure_6_correlations.csv"))
  groups <- split(d, interaction(d$panel, d$Final_Taxonomy, d$response, drop = TRUE))
  ans <- do.call(rbind, lapply(groups, function(z) {
    form <- if (z$panel[1] %in% c("High_affinity", "Low_affinity")) "log(y) ~ log(x)" else
      if (z$panel[1] %in% c("P_mobilization", "Organic_p_and_Phosphonate")) "log(y) ~ x" else "y ~ x"
    x <- z$predictor_value; y <- z$response_value
    m <- lm(as.formula(form))
    data.frame(panel = z$panel[1], taxon = z$Final_Taxonomy[1], response = z$response[1],
      predictor = z$predictor[1], formula = form, n = nrow(z), intercept = unname(coef(m)[1]),
      slope = unname(coef(m)[2]), r2 = summary(m)$r.squared,
      r2_scale = if (grepl("log(y)", form, fixed = TRUE)) "log_response" else "original_response",
      p_independence_assumed = coef(summary(m))[2, 4],
      p_historical_factor = coef(summary(m))[2, 4] * sqrt(6 / 100),
      limitation = "conventional p assumes independent states; not valid treatment inference")
  }))
  rownames(ans) <- NULL
  ans$q_BH_exploratory <- p.adjust(ans$p_independence_assumed, "BH")
  ans
}
if (sys.nframe() == 0L) {
  out <- file.path("results", "phosphorus_indicators")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(audit_historical_regressions(), file.path(out, "historical_regression_audit.csv"), row.names = FALSE)
}
