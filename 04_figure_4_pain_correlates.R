# =============================================================================
# Figure 4. Clinical, biochemical, and environmental correlates of pain intensity
# R 4.5.1
#
# Analysis A: Baseline concurrent correlates of baseline VAS
# Analysis B: Follow-up concurrent correlates of follow-up VAS
# Analysis C: Baseline predictors of follow-up VAS adjusted for baseline VAS
#
# Analysis notes
# - Pain-location variables are intentionally excluded from Paper 1.
# - Exploratory analyses only.
# - Spearman correlations use BH-FDR within each prespecified family.
# - Adjusted predictor models use HC3 robust standard errors and BH-FDR
#   across the predictor family.
# - Associations are interpreted as exploratory and noncausal.
# =============================================================================

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2",
  "patchwork", "sandwich", "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with install.packages()."
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sandwich)
library(openxlsx)

message("Running Figure 4: pain-intensity correlates")

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------

data_path <- file.path(
  "data",
  "Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx"
)
output_dir <- file.path("outputs", "Figure_4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

dat0 <- read_xlsx(data_path, sheet = "Analysis_Data")

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------

to_binary01 <- function(x, var_name) {
  if (is.logical(x)) return(as.integer(x))

  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(na.omit(x)))
    if (all(ux %in% c(0, 1))) return(as.integer(x))
    if (all(ux %in% c(1, 2))) return(as.integer(x == 2))
  }

  z <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(z))

  no_values <- c("0","no","n","absent","negative","none","false")
  yes_values <- c("1","yes","y","present","positive","true")

  out[z %in% no_values] <- 0L
  out[z %in% yes_values] <- 1L

  unresolved <- unique(z[!is.na(z) & is.na(out)])
  if (length(unresolved) > 0) {
    stop(
      "Could not convert ", var_name,
      " to 0/1. Values: ",
      paste(unresolved, collapse = ", ")
    )
  }
  out
}

to_female01 <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(na.omit(x)))
    if (all(ux %in% c(0, 1))) return(as.integer(x))
    if (all(ux %in% c(1, 2))) return(as.integer(x == 2))
  }

  z <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("male","m","man","남","남성")] <- 0L
  out[z %in% c("female","f","woman","여","여성")] <- 1L

  unresolved <- unique(z[!is.na(z) & is.na(out)])
  if (length(unresolved) > 0) {
    stop(
      "Could not convert Sex to female=1/male=0. Values: ",
      paste(unresolved, collapse = ", ")
    )
  }
  out
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x), "NA",
    formatC(x, format = "f", digits = digits)
  )
}

fmt_signed <- function(x, digits = 3) {
  ifelse(
    is.na(x), "NA",
    paste0(
      ifelse(x > 0, "+", ""),
      formatC(x, format = "f", digits = digits)
    )
  )
}

# -----------------------------------------------------------------------------
# 3. Analysis dataset
# -----------------------------------------------------------------------------

dat <- dat0 %>%
  mutate(
    Female = to_female01(Sex),

    Bruxism_bin = to_binary01(Bruxism, "Bruxism"),
    Clenching_bin = to_binary01(Clenching, "Clenching"),
    Parafunction_bin = to_binary01(Parafunction, "Parafunction"),
    TMJ_noise_bin = to_binary01(TMJ_noise, "TMJ_noise"),
    Self_reported_TMJ_noise_bin =
      to_binary01(Self_reported_TMJ_noise, "Self_reported_TMJ_noise"),
    Objective_locking_bin =
      to_binary01(Objective_locking, "Objective_locking"),
    Self_reported_locking_bin =
      to_binary01(Self_reported_locking, "Self_reported_locking"),
    Stiffness_bin = to_binary01(Stiffness, "Stiffness"),
    Trauma_history_bin =
      to_binary01(Trauma_history, "Trauma_history"),

    log_symptom_duration = log1p(Symptom_duration_months)
  )

# -----------------------------------------------------------------------------
# 4. Baseline concurrent family
# -----------------------------------------------------------------------------

baseline_spec <- tibble(
  Variable = c(
    "Baseline_MMO_mm",
    "Baseline_CMO_mm",
    "Self_reported_TMJ_noise_bin",
    "Bruxism_bin",
    "Age_years",
    "Female",
    "Trauma_history_bin",
    "Parafunction_bin",
    "TMJ_noise_bin",
    "Symptom_duration_months",
    "Stiffness_bin",
    "Clenching_bin",
    "SolarRad_preBL_90d_mean_MJm2",
    "Objective_locking_bin",
    "Baseline_25OHD_ng_mL",
    "Self_reported_locking_bin",
    "Sunshine_preBL_90d_mean_hr",
    "Sunshine_preBL_60d_mean_hr",
    "SolarRad_preBL_60d_mean_MJm2"
  ),
  Label = c(
    "MUO",
    "PFO",
    "Self-reported TMJ noise",
    "Self-reported bruxism",
    "Age",
    "Female sex",
    "Trauma history",
    "Parafunction",
    "Clinically detected TMJ noise",
    "Symptom duration",
    "Stiffness",
    "Clenching",
    "90-day solar radiation",
    "Clinically detected jaw locking",
    "Serum 25(OH)D",
    "Self-reported jaw locking",
    "90-day sunshine duration",
    "60-day sunshine duration",
    "60-day solar radiation"
  )
)

run_spearman_family <- function(data, outcome, spec) {
  bind_rows(
    lapply(seq_len(nrow(spec)), function(i) {
      v <- spec$Variable[i]

      d <- data %>%
        select(all_of(c(outcome, v))) %>%
        filter(complete.cases(.))

      ct <- suppressWarnings(
        cor.test(
          d[[v]],
          d[[outcome]],
          method = "spearman",
          exact = FALSE
        )
      )

      tibble(
        Variable = v,
        Label = spec$Label[i],
        N = nrow(d),
        Spearman_rho = unname(ct$estimate),
        p_value = ct$p.value
      )
    })
  ) %>%
    mutate(
      q_BH = p.adjust(p_value, method = "BH"),
      FDR_significant = q_BH < 0.05
    ) %>%
    arrange(p_value)
}

baseline_corr <- run_spearman_family(
  dat,
  "Baseline_pain_intensity_VAS",
  baseline_spec
)

# -----------------------------------------------------------------------------
# 5. Follow-up concurrent family
# -----------------------------------------------------------------------------
# These are variables that are meaningfully interpretable at the follow-up visit.
# Age and sex are retained for comparison with baseline; biochemical, function,
# and environmental measures use follow-up values.

followup_spec <- tibble(
  Variable = c(
    "Age_years",
    "Female",
    "Followup_25OHD_ng_mL",
    "Followup_CMO_mm",
    "Followup_MMO_mm",
    "SolarRad_preFU_60d_mean_MJm2",
    "SolarRad_preFU_90d_mean_MJm2",
    "Sunshine_preFU_60d_mean_hr",
    "Sunshine_preFU_90d_mean_hr"
  ),
  Label = c(
    "Age",
    "Female sex",
    "Serum 25(OH)D",
    "PFO",
    "MUO",
    "60-day solar radiation",
    "90-day solar radiation",
    "60-day sunshine duration",
    "90-day sunshine duration"
  )
)

missing_fu <- setdiff(
  followup_spec$Variable,
  names(dat)
)

if (length(missing_fu) > 0) {
  stop(
    "Missing follow-up variable(s): ",
    paste(missing_fu, collapse = ", ")
  )
}

followup_corr <- run_spearman_family(
  dat,
  "Followup_pain_intensity_VAS",
  followup_spec
)

# -----------------------------------------------------------------------------
# 6. Baseline predictors of follow-up VAS, adjusted for baseline VAS
# -----------------------------------------------------------------------------
# Candidate family mirrors the baseline exploratory set.
# Each predictor is assessed in a separate model.
#
# Core adjustment:
#   baseline VAS + clinical follow-up interval
# Plus age and sex unless the predictor itself is age or sex.
#
# Continuous predictor effects are expressed per 1 SD increase.
# Binary predictor effects are present vs absent.
# Outcome remains follow-up VAS in original 0-10 units.

predictor_spec <- baseline_spec %>%
  mutate(
    Binary = Variable %in% c(
      "Self_reported_TMJ_noise_bin",
      "Bruxism_bin",
      "Female",
      "Trauma_history_bin",
      "Parafunction_bin",
      "TMJ_noise_bin",
      "Stiffness_bin",
      "Clenching_bin",
      "Objective_locking_bin",
      "Self_reported_locking_bin"
    )
  )

fit_predictor <- function(v, label, binary_flag) {

  covars <- c(
    "Baseline_pain_intensity_VAS",
    "Clinical_followup_months"
  )

  if (v != "Age_years") covars <- c(covars, "Age_years")
  if (v != "Female") covars <- c(covars, "Female")

  needed <- unique(c(
    "Followup_pain_intensity_VAS",
    v,
    covars
  ))

  d <- dat %>%
    select(all_of(needed)) %>%
    filter(complete.cases(.))

  if (!binary_flag) {
    sd_v <- sd(d[[v]], na.rm = TRUE)
    if (!is.finite(sd_v) || sd_v == 0) {
      stop("Zero/invalid SD for predictor: ", v)
    }
    newv <- paste0(v, "_perSD")
    d[[newv]] <- d[[v]] / sd_v
    predictor_term <- newv
    scale_note <- paste0("Per 1 SD (", signif(sd_v, 4), " original units)")
  } else {
    predictor_term <- v
    sd_v <- NA_real_
    scale_note <- "Present vs absent"
  }

  f <- reformulate(
    termlabels = c(predictor_term, covars),
    response = "Followup_pain_intensity_VAS"
  )

  fit <- lm(f, data = d)
  V <- sandwich::vcovHC(fit, type = "HC3")

  b <- coef(fit)[predictor_term]
  se <- sqrt(diag(V))[predictor_term]

  df_res <- df.residual(fit)
  crit <- qt(0.975, df = df_res)

  p <- 2 * pt(
    abs(b / se),
    df = df_res,
    lower.tail = FALSE
  )

  tibble(
    Variable = v,
    Label = label,
    Binary = binary_flag,
    N = nrow(d),
    Predictor_SD = sd_v,
    Scale = scale_note,
    Beta = unname(b),
    CI_low = unname(b - crit * se),
    CI_high = unname(b + crit * se),
    p_value = unname(p)
  )
}

baseline_predictors <- bind_rows(
  lapply(seq_len(nrow(predictor_spec)), function(i) {
    fit_predictor(
      predictor_spec$Variable[i],
      predictor_spec$Label[i],
      predictor_spec$Binary[i]
    )
  })
) %>%
  mutate(
    q_BH = p.adjust(p_value, method = "BH"),
    FDR_significant = q_BH < 0.05
  ) %>%
  arrange(p_value)

# -----------------------------------------------------------------------------
# 7. Positive-association summaries
# -----------------------------------------------------------------------------

baseline_positive_raw <- baseline_corr %>%
  filter(
    Spearman_rho > 0,
    p_value < 0.05
  )

baseline_positive_fdr <- baseline_corr %>%
  filter(
    Spearman_rho > 0,
    q_BH < 0.05
  )

followup_positive_raw <- followup_corr %>%
  filter(
    Spearman_rho > 0,
    p_value < 0.05
  )

followup_positive_fdr <- followup_corr %>%
  filter(
    Spearman_rho > 0,
    q_BH < 0.05
  )

predictor_positive_raw <- baseline_predictors %>%
  filter(
    Beta > 0,
    p_value < 0.05
  )

predictor_positive_fdr <- baseline_predictors %>%
  filter(
    Beta > 0,
    q_BH < 0.05
  )

# -----------------------------------------------------------------------------
# 8. Figure data
# -----------------------------------------------------------------------------

plot_corr <- function(
  df,
  xlab,
  highlight_label = NULL,
  highlight_color = "#7A3E2E",
  significant_color = "#C23B6B"
) {

  pdat <- df %>%
    arrange(Spearman_rho) %>%
    mutate(
      Label = factor(Label, levels = Label),

      txt = paste0(
        "\u03c1 = ", fmt_signed(Spearman_rho, 3),
        "; p ", ifelse(
          p_value < 0.001,
          "<0.001",
          paste0("= ", fmt_p(p_value))
        ),
        "; q ", ifelse(
          q_BH < 0.001,
          "<0.001",
          paste0("= ", fmt_p(q_BH))
        )
      )
    )

  # Base-R color classification outside mutate().
  # This is deliberately simple and cannot create a zero-length mutate column.
  pdat$Color_group <- rep("p ≥ 0.05", nrow(pdat))
  pdat$Color_group[pdat$p_value < 0.05] <- "p < 0.05"

  if (!is.null(highlight_label) && length(highlight_label) == 1) {
    idx_highlight <- (
      as.character(pdat$Label) == highlight_label &
      pdat$p_value < 0.05
    )
    pdat$Color_group[idx_highlight] <- "Highlighted"
  }

  xmin <- min(c(pdat$Spearman_rho, 0), na.rm = TRUE)
  xmax <- max(c(pdat$Spearman_rho, 0), na.rm = TRUE)
  span <- xmax - xmin
  if (!is.finite(span) || span <= 0) span <- 1

  # Compact positive x-axis, per final figure preference.
  text_x <- max(0.22, xmax + 0.05)
  upper <- 0.50
  lower <- xmin - 0.08 * span

  ggplot(pdat, aes(x = Spearman_rho, y = Label)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "#8A8A8A",
      linewidth = 0.55
    ) +
    geom_segment(
      aes(
        x = 0,
        xend = Spearman_rho,
        y = Label,
        yend = Label,
        colour = Color_group
      ),
      linewidth = 0.75
    ) +
    geom_point(
      aes(colour = Color_group),
      size = 3.0
    ) +
    geom_text(
      aes(
        x = text_x,
        label = txt,
        colour = Color_group
      ),
      hjust = 0,
      family = "Arial",
      size = 2.55,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Highlighted" = highlight_color,
        "p < 0.05" = significant_color,
        "p ≥ 0.05" = "#6D6D6D"
      )
    ) +
    scale_x_continuous(
      limits = c(lower, upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = xlab,
      y = NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#FCFCFC", colour = NA),
      axis.text = element_text(colour = "#252525", size = 9.0),
      axis.text.y = element_text(size = 9.0),
      axis.title = element_text(face = "bold", size = 10.2),
      legend.position = "none",
      plot.margin = margin(12, 18, 12, 12)
    )
}

# Panel A:
# - Age: dark reddish-brown
# - Other nominally significant associations (raw p < 0.05): rose-pink
# - p >= 0.05: gray
panel_a <- plot_corr(
  baseline_corr,
  "Association with baseline pain intensity (Spearman \u03c1)",
  highlight_label = "Age",
  highlight_color = "#7A3E2E",
  significant_color = "#C23B6B"
)

# Panel B:
# - All nominally significant associations (raw p < 0.05): refined green
# - p >= 0.05: gray
panel_b <- plot_corr(
  followup_corr,
  "Association with follow-up pain intensity (Spearman \u03c1)",
  highlight_label = NULL,
  highlight_color = "#2F7D5A",
  significant_color = "#2F7D5A"
)

efig5 <- panel_a / panel_b +
  plot_layout(heights = c(1.45, 1.0))

print(efig5)

# -----------------------------------------------------------------------------
# 9. Export
# -----------------------------------------------------------------------------

xlsx_out <- file.path(
  output_dir,
  "Figure_4_pain_correlates_statistics.xlsx"
)

wb <- createWorkbook()

sheets <- list(
  "Baseline_correlates" = baseline_corr,
  "Followup_correlates" = followup_corr,
  "Baseline_predictors_FU_VAS" = baseline_predictors,
  "BL_positive_rawP" = baseline_positive_raw,
  "BL_positive_FDR" = baseline_positive_fdr,
  "FU_positive_rawP" = followup_positive_raw,
  "FU_positive_FDR" = followup_positive_fdr,
  "Predictor_positive_rawP" = predictor_positive_raw,
  "Predictor_positive_FDR" = predictor_positive_fdr
)

for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]])
  freezePane(wb, nm, firstRow = TRUE)
  setColWidths(wb, nm, cols = 1:ncol(sheets[[nm]]), widths = "auto")
}

saveWorkbook(wb, xlsx_out, overwrite = TRUE)

png_out <- file.path(
  output_dir,
  "Figure_4_pain_correlates.png"
)

tiff_out <- file.path(
  output_dir,
  "Figure_4_pain_correlates.tiff"
)

pdf_out <- file.path(
  output_dir,
  "Figure_4_pain_correlates.pdf"
)

ggsave(
  png_out, efig5,
  width = 12, height = 11,
  units = "in", dpi = 400, bg = "white"
)

ggsave(
  tiff_out, efig5,
  width = 12, height = 11,
  units = "in", dpi = 600,
  compression = "lzw", bg = "white"
)

ggsave(
  pdf_out, efig5,
  width = 12, height = 11,
  units = "in", device = cairo_pdf, bg = "white"
)

# -----------------------------------------------------------------------------
# 10. Console report
# -----------------------------------------------------------------------------

cat("\n====================================================\n")
cat("BASELINE VAS CORRELATES\n")
cat("====================================================\n")
print(baseline_corr)

cat("\nPositive associations with raw p < 0.05:\n")
print(baseline_positive_raw)

cat("\nPositive associations surviving BH-FDR:\n")
print(baseline_positive_fdr)

cat("\n====================================================\n")
cat("FOLLOW-UP VAS CONCURRENT CORRELATES\n")
cat("====================================================\n")
print(followup_corr)

cat("\nPositive associations with raw p < 0.05:\n")
print(followup_positive_raw)

cat("\nPositive associations surviving BH-FDR:\n")
print(followup_positive_fdr)

cat("\n====================================================\n")
cat("BASELINE PREDICTORS OF FOLLOW-UP VAS\n")
cat("Adjusted for baseline VAS + follow-up interval + age/sex as applicable\n")
cat("====================================================\n")
print(baseline_predictors)

cat("\nPositive predictors with raw p < 0.05:\n")
print(predictor_positive_raw)

cat("\nPositive predictors surviving BH-FDR:\n")
print(predictor_positive_fdr)

cat("\nSaved:\n")
cat(xlsx_out, "\n")
cat(png_out, "\n")
cat(tiff_out, "\n")
cat(pdf_out, "\n")
