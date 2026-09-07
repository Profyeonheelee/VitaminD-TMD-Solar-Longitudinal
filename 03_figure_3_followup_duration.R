# =============================================================================
# Figure 3. Associations of follow-up duration with changes in serum 25(OH)D
# and pain intensity among vitamin D-prescribed patients
# -----------------------------------------------------------------------------
# Figure itself contains:
#   - panel labels "A" and "B" only
#   - axis labels
#   - individual observations
#   - adjusted marginal trajectories
#   - HC3 robust 95% confidence bands
#   - statistical annotations
#
# The full panel descriptions are provided in the manuscript legend.
#
# Population:
#   Patients with recorded vitamin D prescription (expected n = 149)
#
# Panel A:
#   x = observed clinical follow-up interval (months)
#   y = change in serum 25(OH)D (follow-up - baseline)
#
# Panel B:
#   x = observed clinical follow-up interval (months)
#   y = pain reduction, VAS (baseline - follow-up)
#
# Color scheme
# -----------------------------------------------------------------------------
# Predominantly grayscale.
# Muted teal is used only when the adjusted association has p < 0.05.
#
# Interpretation
# -----------------------------------------------------------------------------
# Clinical_followup_months is the observed baseline-to-follow-up interval.
# It is not prescription duration, adherence, or cumulative exposure.
#
# Outputs
# -----------------------------------------------------------------------------
# Figure_3_followup_duration.png
# Figure_3_followup_duration.tiff
# Figure_3_followup_duration.pdf
# Figure_3_followup_duration_statistics.xlsx
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "patchwork",
  "sandwich",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sandwich)
library(openxlsx)


# =============================================================================
# 1. File path
# =============================================================================

xlsx_path <- file.path(
  "data",
  "Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx"
)

sheet_name <- "Analysis_Data"
output_dir <- file.path("outputs", "Figure_3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(xlsx_path)) {
  stop("Data file not found: ", xlsx_path)
}


# =============================================================================
# 2. Import and validation
# =============================================================================

dat0 <- read_excel(
  xlsx_path,
  sheet = sheet_name
)

required_vars <- c(
  "Study_ID",
  "Sex",
  "Age_years",
  "Symptom_duration_months",
  "Baseline_25OHD_ng_mL",
  "Followup_25OHD_ng_mL",
  "Delta_25OHD_ng_mL",
  "Vitamin_D_prescription",
  "Clinical_followup_months",
  "Baseline_pain_intensity_VAS",
  "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS",
  "Delta_SolarRad_60d_MJm2"
)

missing_vars <- setdiff(required_vars, names(dat0))

if (length(missing_vars) > 0) {
  stop(
    "Missing required variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

dat <- dat0 %>%
  mutate(
    Sex = factor(Sex),
    Rx = as.integer(Vitamin_D_prescription),
    log_symptom_duration = log1p(Symptom_duration_months)
  )

if (nrow(dat) != 241) {
  stop("Expected full cohort N=241; found N=", nrow(dat))
}

if (anyDuplicated(dat$Study_ID) > 0) {
  stop("Duplicate Study_ID detected.")
}

if (!all(na.omit(unique(dat$Rx)) %in% c(0, 1))) {
  stop("Vitamin_D_prescription must be coded 0/1.")
}

rx <- dat %>%
  filter(Rx == 1)

if (nrow(rx) != 149) {
  stop(
    "Expected prescribed cohort n=149; found n=",
    nrow(rx)
  )
}


# =============================================================================
# 3. Complete-case datasets
# =============================================================================

bio_vars <- c(
  "Study_ID",
  "Clinical_followup_months",
  "Delta_25OHD_ng_mL",
  "Baseline_25OHD_ng_mL",
  "Age_years",
  "Sex",
  "Delta_SolarRad_60d_MJm2"
)

pain_vars <- c(
  "Study_ID",
  "Clinical_followup_months",
  "Pain_intensity_reduction_VAS",
  "Baseline_pain_intensity_VAS",
  "Baseline_25OHD_ng_mL",
  "Age_years",
  "Sex",
  "log_symptom_duration",
  "Delta_SolarRad_60d_MJm2"
)

bio_dat <- rx %>%
  select(all_of(bio_vars)) %>%
  filter(complete.cases(.))

pain_dat <- rx %>%
  select(all_of(pain_vars)) %>%
  filter(complete.cases(.))


# =============================================================================
# 4. Formatting helpers
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_num <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits)
}

fmt_ci <- function(est, lo, hi, digits = 2) {
  paste0(
    fmt_num(est, digits),
    " (",
    fmt_num(lo, digits),
    " to ",
    fmt_num(hi, digits),
    ")"
  )
}


# =============================================================================
# 5. HC3 robust coefficient extractor
# =============================================================================

hc3_term <- function(model, term_name) {

  V <- sandwich::vcovHC(model, type = "HC3")

  beta <- coef(model)[term_name]
  se <- sqrt(diag(V))[term_name]
  df_res <- df.residual(model)

  crit <- qt(0.975, df = df_res)

  stat <- beta / se

  p <- 2 * pt(
    abs(stat),
    df = df_res,
    lower.tail = FALSE
  )

  tibble(
    beta = unname(beta),
    se = unname(se),
    ci_low = unname(beta - crit * se),
    ci_high = unname(beta + crit * se),
    p = unname(p)
  )
}


# =============================================================================
# 6. Marginally standardized HC3 prediction curve
# =============================================================================

marginal_robust_curve <- function(
  model,
  data,
  xvar,
  xgrid
) {

  V <- sandwich::vcovHC(model, type = "HC3")

  beta <- coef(model)
  df_res <- df.residual(model)
  crit <- qt(0.975, df = df_res)

  tt <- delete.response(terms(model))

  out <- lapply(
    xgrid,
    function(xval) {

      nd <- data
      nd[[xvar]] <- xval

      X <- model.matrix(
        tt,
        data = nd,
        contrasts.arg = model$contrasts
      )

      X <- X[, names(beta), drop = FALSE]

      xbar <- colMeans(X)

      fit <- as.numeric(xbar %*% beta)

      se <- sqrt(
        as.numeric(
          t(xbar) %*% V %*% xbar
        )
      )

      data.frame(
        x = xval,
        fit = fit,
        lower = fit - crit * se,
        upper = fit + crit * se
      )
    }
  )

  bind_rows(out)
}


# =============================================================================
# 7. Panel A analysis
# =============================================================================

bio_spearman <- suppressWarnings(
  cor.test(
    bio_dat$Clinical_followup_months,
    bio_dat$Delta_25OHD_ng_mL,
    method = "spearman",
    exact = FALSE
  )
)

bio_model <- lm(
  Delta_25OHD_ng_mL ~
    Clinical_followup_months +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex +
    Delta_SolarRad_60d_MJm2,
  data = bio_dat
)

bio_adj <- hc3_term(
  bio_model,
  "Clinical_followup_months"
)

bio_xgrid <- seq(
  min(bio_dat$Clinical_followup_months),
  max(bio_dat$Clinical_followup_months),
  length.out = 160
)

bio_curve <- marginal_robust_curve(
  model = bio_model,
  data = bio_dat,
  xvar = "Clinical_followup_months",
  xgrid = bio_xgrid
)


# =============================================================================
# 8. Panel B analysis
# =============================================================================

pain_spearman <- suppressWarnings(
  cor.test(
    pain_dat$Clinical_followup_months,
    pain_dat$Pain_intensity_reduction_VAS,
    method = "spearman",
    exact = FALSE
  )
)

pain_model <- lm(
  Pain_intensity_reduction_VAS ~
    Clinical_followup_months +
    Baseline_pain_intensity_VAS +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex +
    log_symptom_duration +
    Delta_SolarRad_60d_MJm2,
  data = pain_dat
)

pain_adj <- hc3_term(
  pain_model,
  "Clinical_followup_months"
)

pain_xgrid <- seq(
  min(pain_dat$Clinical_followup_months),
  max(pain_dat$Clinical_followup_months),
  length.out = 160
)

pain_curve <- marginal_robust_curve(
  model = pain_model,
  data = pain_dat,
  xvar = "Clinical_followup_months",
  xgrid = pain_xgrid
)


# =============================================================================
# 9. Visual system
# =============================================================================

COL_POINT       <- "#A9A9A9"
COL_POINT_EDGE  <- "#6F6F6F"
COL_NONSIG      <- "#4A4A4A"
COL_NONSIG_BAND <- "#BDBDBD"

COL_SIG         <- "#137C8B"
COL_SIG_BAND    <- "#8FC8CE"

COL_ZERO        <- "#8A8A8A"
COL_TEXT        <- "#252525"
COL_PANEL_BG    <- "#FCFCFC"

bio_sig <- bio_adj$p < 0.05
pain_sig <- pain_adj$p < 0.05

bio_line_col <- if (bio_sig) COL_SIG else COL_NONSIG
bio_band_col <- if (bio_sig) COL_SIG_BAND else COL_NONSIG_BAND

pain_line_col <- if (pain_sig) COL_SIG else COL_NONSIG
pain_band_col <- if (pain_sig) COL_SIG_BAND else COL_NONSIG_BAND


# =============================================================================
# 10. Statistical annotation text
# =============================================================================

bio_annotation <- paste0(
  "n = ", nrow(bio_dat),
  "\nSpearman \u03c1 = ",
  sprintf("%.3f", unname(bio_spearman$estimate)),
  "; p ",
  ifelse(
    bio_spearman$p.value < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(bio_spearman$p.value))
  ),
  "\nAdjusted \u03b2/month = ",
  fmt_ci(
    bio_adj$beta,
    bio_adj$ci_low,
    bio_adj$ci_high,
    digits = 2
  ),
  "; p ",
  ifelse(
    bio_adj$p < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(bio_adj$p))
  )
)

pain_annotation <- paste0(
  "n = ", nrow(pain_dat),
  "\nSpearman \u03c1 = ",
  sprintf("%.3f", unname(pain_spearman$estimate)),
  "; p ",
  ifelse(
    pain_spearman$p.value < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(pain_spearman$p.value))
  ),
  "\nAdjusted \u03b2/month = ",
  fmt_ci(
    pain_adj$beta,
    pain_adj$ci_low,
    pain_adj$ci_high,
    digits = 2
  ),
  "; p ",
  ifelse(
    pain_adj$p < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(pain_adj$p))
  )
)


# =============================================================================
# 11. Shared theme
# =============================================================================

theme_fig2 <- theme_classic(
  base_family = "Arial",
  base_size = 11
) +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = COL_PANEL_BG,
      colour = NA
    ),
    axis.line = element_line(
      colour = "#505050",
      linewidth = 0.55
    ),
    axis.ticks = element_line(
      colour = "#505050",
      linewidth = 0.45
    ),
    axis.text = element_text(
      colour = COL_TEXT,
      size = 9.5
    ),
    axis.title = element_text(
      colour = COL_TEXT,
      size = 10.5,
      face = "bold"
    ),
    plot.margin = margin(
      10, 12, 10, 10
    )
  )


# =============================================================================
# 12. Panel A
# =============================================================================

panel_a <- ggplot(
  bio_dat,
  aes(
    x = Clinical_followup_months,
    y = Delta_25OHD_ng_mL
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_ZERO
  ) +

  geom_point(
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    fill = COL_POINT,
    colour = COL_POINT_EDGE,
    alpha = 0.48
  ) +

  geom_ribbon(
    data = bio_curve,
    aes(
      x = x,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = bio_band_col,
    alpha = if (bio_sig) 0.30 else 0.22
  ) +

  geom_line(
    data = bio_curve,
    aes(
      x = x,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = bio_line_col,
    linewidth = 1.25
  ) +

  # Panel label only
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "A",
    hjust = -0.35,
    vjust = 1.10,
    size = 5.6,
    fontface = "bold",
    family = "Arial",
    colour = COL_TEXT
  ) +

  # Larger statistical annotation box
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = bio_annotation,
    hjust = 1.03,
    vjust = 1.12,
    size = 3.0,
    family = "Arial",
    lineheight = 1.10,
    colour = if (bio_sig) COL_SIG else COL_TEXT,
    fill = "white",
    label.size = 0.28,
    label.padding = grid::unit(0.46, "lines"),
    label.r = grid::unit(0.05, "lines")
  ) +

  labs(
    x = "Observed clinical follow-up, months",
    y = expression(Delta*"25(OH)D, ng/mL")
  ) +

  theme_fig2


# =============================================================================
# 13. Panel B
# =============================================================================

panel_b <- ggplot(
  pain_dat,
  aes(
    x = Clinical_followup_months,
    y = Pain_intensity_reduction_VAS
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_ZERO
  ) +

  geom_point(
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    fill = COL_POINT,
    colour = COL_POINT_EDGE,
    alpha = 0.48
  ) +

  geom_ribbon(
    data = pain_curve,
    aes(
      x = x,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = pain_band_col,
    alpha = if (pain_sig) 0.30 else 0.22
  ) +

  geom_line(
    data = pain_curve,
    aes(
      x = x,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = pain_line_col,
    linewidth = 1.25
  ) +

  # Panel label only
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "B",
    hjust = -0.35,
    vjust = 1.10,
    size = 5.6,
    fontface = "bold",
    family = "Arial",
    colour = COL_TEXT
  ) +

  # Larger statistical annotation box
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = pain_annotation,
    hjust = 1.03,
    vjust = 1.12,
    size = 3.0,
    family = "Arial",
    lineheight = 1.10,
    colour = if (pain_sig) COL_SIG else COL_TEXT,
    fill = "white",
    label.size = 0.28,
    label.padding = grid::unit(0.46, "lines"),
    label.r = grid::unit(0.05, "lines")
  ) +

  labs(
    x = "Observed clinical follow-up, months",
    y = "Pain reduction, VAS"
  ) +

  theme_fig2


# =============================================================================
# 14. Combine panels
# =============================================================================
# No title and no subtitle inside the figure.

figure2 <- panel_a + panel_b +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

print(figure2)


# =============================================================================
# 15. Statistical audit tables
# =============================================================================

main_stats <- bind_rows(
  tibble(
    Panel = "A",
    Outcome = "Delta serum 25(OH)D, ng/mL",
    N = nrow(bio_dat),
    Exposure = "Clinical follow-up interval, months",
    Spearman_rho = unname(bio_spearman$estimate),
    Spearman_p = bio_spearman$p.value,
    Adjusted_beta_per_month = bio_adj$beta,
    HC3_CI_low = bio_adj$ci_low,
    HC3_CI_high = bio_adj$ci_high,
    Adjusted_p = bio_adj$p,
    Significant_adjusted_p_lt_0_05 = bio_sig
  ),
  tibble(
    Panel = "B",
    Outcome = "Pain reduction, VAS",
    N = nrow(pain_dat),
    Exposure = "Clinical follow-up interval, months",
    Spearman_rho = unname(pain_spearman$estimate),
    Spearman_p = pain_spearman$p.value,
    Adjusted_beta_per_month = pain_adj$beta,
    HC3_CI_low = pain_adj$ci_low,
    HC3_CI_high = pain_adj$ci_high,
    Adjusted_p = pain_adj$p,
    Significant_adjusted_p_lt_0_05 = pain_sig
  )
)

bio_model_audit <- {
  V <- sandwich::vcovHC(bio_model, type = "HC3")
  se <- sqrt(diag(V))
  df_res <- df.residual(bio_model)
  crit <- qt(0.975, df = df_res)
  b <- coef(bio_model)
  p <- 2 * pt(abs(b / se), df = df_res, lower.tail = FALSE)

  tibble(
    Term = names(b),
    Estimate = unname(b),
    HC3_SE = unname(se),
    CI_low = unname(b - crit * se),
    CI_high = unname(b + crit * se),
    HC3_p = unname(p)
  )
}

pain_model_audit <- {
  V <- sandwich::vcovHC(pain_model, type = "HC3")
  se <- sqrt(diag(V))
  df_res <- df.residual(pain_model)
  crit <- qt(0.975, df = df_res)
  b <- coef(pain_model)
  p <- 2 * pt(abs(b / se), df = df_res, lower.tail = FALSE)

  tibble(
    Term = names(b),
    Estimate = unname(b),
    HC3_SE = unname(se),
    CI_low = unname(b - crit * se),
    CI_high = unname(b + crit * se),
    HC3_p = unname(p)
  )
}

analysis_info <- tibble(
  Item = c(
    "Figure population",
    "Expected prescribed cohort N",
    "Panel labels",
    "Figure title/subtitle",
    "Time axis",
    "Time-axis interpretation",
    "Panel A outcome",
    "Panel A adjusted covariates",
    "Panel B outcome",
    "Panel B adjusted covariates",
    "Adjusted uncertainty",
    "Adjusted trajectory",
    "Color rule",
    "Prescription duration",
    "Missing data",
    "R version requested for manuscript"
  ),
  Value = c(
    "Vitamin D-prescribed patients",
    "149",
    "A and B only",
    "None inside figure",
    "Clinical_followup_months",
    paste(
      "Observed baseline-to-follow-up interval only;",
      "not prescription duration, adherence, or cumulative dose"
    ),
    "Delta_25OHD_ng_mL = follow-up minus baseline",
    paste(
      "baseline 25(OH)D + age + sex +",
      "Delta 60-day solar radiation"
    ),
    "Pain_intensity_reduction_VAS = baseline minus follow-up VAS",
    paste(
      "baseline VAS + baseline 25(OH)D + age + sex +",
      "log1p symptom duration + Delta 60-day solar radiation"
    ),
    "HC3 heteroskedasticity-robust standard errors and 95% CI",
    "Marginal standardization over observed prescribed-patient covariate distribution",
    "Muted teal only when adjusted p<0.05; otherwise grayscale",
    "Excluded from Figure 3 and final Paper 1 inference",
    "Complete-case analysis; no imputation",
    "R 4.5.1"
  )
)


# =============================================================================
# 16. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Figure_3_followup_duration.png"
)

tiff_file <- file.path(
  output_dir,
  "Figure_3_followup_duration.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Figure_3_followup_duration.pdf"
)

ggsave(
  filename = png_file,
  plot = figure2,
  width = 10.8,
  height = 5.1,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure2,
  width = 10.8,
  height = 5.1,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = figure2,
  width = 10.8,
  height = 5.1,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 17. Export statistics workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Figure_3_followup_duration_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "Figure3_Main_Stats")
addWorksheet(wb, "Biochemical_Model")
addWorksheet(wb, "Pain_Model")
addWorksheet(wb, "Biochemical_Curve")
addWorksheet(wb, "Pain_Curve")
addWorksheet(wb, "Analysis_Info")

writeData(wb, "Figure3_Main_Stats", main_stats)
writeData(wb, "Biochemical_Model", bio_model_audit)
writeData(wb, "Pain_Model", pain_model_audit)
writeData(wb, "Biochemical_Curve", bio_curve)
writeData(wb, "Pain_Curve", pain_curve)
writeData(wb, "Analysis_Info", analysis_info)

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sh in names(wb)) {
  temp <- readWorkbook(wb, sheet = sh)

  if (ncol(temp) > 0) {
    addStyle(
      wb,
      sh,
      header_style,
      rows = 1,
      cols = seq_len(ncol(temp)),
      gridExpand = TRUE
    )

    freezePane(
      wb,
      sh,
      firstRow = TRUE
    )

    setColWidths(
      wb,
      sh,
      cols = seq_len(ncol(temp)),
      widths = "auto"
    )
  }
}

setColWidths(
  wb,
  "Analysis_Info",
  cols = 1:2,
  widths = c(38, 95)
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 18. Console output
# =============================================================================

cat("\n============================================================\n")
cat("FIGURE 2 COMPLETED\n")
cat("============================================================\n\n")

cat("Panel A: Delta 25(OH)D over observed follow-up\n")
cat("N =", nrow(bio_dat), "\n")
cat(
  "Spearman rho =",
  sprintf("%.3f", unname(bio_spearman$estimate)),
  "; p =", fmt_p(bio_spearman$p.value), "\n"
)
cat(
  "Adjusted beta/month =",
  fmt_ci(
    bio_adj$beta,
    bio_adj$ci_low,
    bio_adj$ci_high,
    2
  ),
  "; p =", fmt_p(bio_adj$p), "\n\n"
)

cat("Panel B: Pain reduction over observed follow-up\n")
cat("N =", nrow(pain_dat), "\n")
cat(
  "Spearman rho =",
  sprintf("%.3f", unname(pain_spearman$estimate)),
  "; p =", fmt_p(pain_spearman$p.value), "\n"
)
cat(
  "Adjusted beta/month =",
  fmt_ci(
    pain_adj$beta,
    pain_adj$ci_low,
    pain_adj$ci_high,
    2
  ),
  "; p =", fmt_p(pain_adj$p), "\n\n"
)

cat("Saved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("============================================================\n")
