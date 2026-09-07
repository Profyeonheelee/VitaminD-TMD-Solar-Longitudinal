# =============================================================================
# Supplementary Figure S5. Environmental-exposure sensitivity analyses
# Distribution + paired-window coefficient display
#
# Panels
#   (A) Environmental change -> change in serum 25(OH)D
#   (B) Environmental change -> pain reduction
#
# Visual concept
#   - Each exposure family occupies ONE horizontal row:
#       Solar radiation
#       Sunshine duration
#   - 60-day and 90-day adjusted effect estimates are displayed as two
#     colored centers on the same row and CONNECTED to emphasize the
#     60-day -> 90-day change.
#   - Small gray points show the bootstrap distribution of adjusted effects.
#   - Two gray horizontal dashed guide lines mark the two exposure-family rows.
#   - 60-day = muted orange
#   - 90-day = muted green
#   - Analytic point estimates and 95% CIs remain HC3 robust estimates.
#   - Bootstrap dots are used for VISUALIZATION of the sampling distribution,
#     not to replace the HC3 confidence intervals or p-values.
#
# Effect scale
#   Adjusted beta per +1 SD increase in each environmental-change exposure.
#   This standardization permits direct visual comparison between solar
#   radiation and sunshine duration, which have different physical units.
#
# Panel A models
#   Delta 25(OH)D ~ environmental change
#                   + baseline 25(OH)D
#                   + age + sex
#                   + clinical follow-up interval
#
# Panel B models
#   Pain reduction ~ environmental change
#                    + baseline VAS
#                    + age + sex
#                    + clinical follow-up interval
#
# Multiplicity
#   BH q-values are calculated separately across the four biochemical
#   sensitivity comparisons and the four pain sensitivity comparisons.
#
# Bootstrap
#   Default B = 500 per model.
#   A fixed original-sample exposure SD is used to express bootstrap
#   coefficients on the same per-1-SD scale as the analytic estimate.
#
# Figure rules
#   - No overall title inside figure
#   - No panel letters inside figure
#   - Orange circle = 60-day
#   - Green triangle = 90-day
#   - Small gray points = bootstrap coefficient distribution
#   - Colored horizontal bars = HC3 95% CI
#   - Gray connector = 60-day -> 90-day center shift
#   - Vertical dashed line = beta 0
#
# Outputs
#   Supplementary_Figure_S5_environmental_sensitivity.png
#   Supplementary_Figure_S5_environmental_sensitivity.tiff
#   Supplementary_Figure_S5_environmental_sensitivity.pdf
#   Supplementary_Figure_S5_environmental_sensitivity_statistics.xlsx
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
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
library(tidyr)
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
output_dir <- file.path("outputs", "Supplementary_Figure_S5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(xlsx_path)) {
  stop("Data file not found: ", xlsx_path)
}


# =============================================================================
# 2. Bootstrap settings
# =============================================================================

BOOT_B <- 500
BOOT_SEED <- 20260825

set.seed(BOOT_SEED)


# =============================================================================
# 3. Import and validation
# =============================================================================

dat0 <- read_excel(
  xlsx_path,
  sheet = sheet_name
)

required_vars <- c(
  "Study_ID",
  "Sex",
  "Age_years",
  "Baseline_25OHD_ng_mL",
  "Delta_25OHD_ng_mL",
  "Baseline_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS",
  "Clinical_followup_months",
  "Delta_SolarRad_60d_MJm2",
  "Delta_SolarRad_90d_MJm2",
  "Delta_Sunshine_60d_hr",
  "Delta_Sunshine_90d_hr"
)

missing_vars <- setdiff(
  required_vars,
  names(dat0)
)

if (length(missing_vars) > 0) {
  stop(
    "Missing required variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

if (anyDuplicated(dat0$Study_ID) > 0) {
  stop("Duplicate Study_ID detected.")
}

dat <- dat0 %>%
  mutate(
    Sex = factor(Sex)
  )


# =============================================================================
# 4. Formatting helpers
# =============================================================================

fmt_p <- function(p) {

  out <- rep(
    NA_character_,
    length(p)
  )

  out[is.na(p)] <- "NA"

  idx_lt <- !is.na(p) & p < 0.001
  out[idx_lt] <- "<0.001"

  idx_other <- !is.na(p) & p >= 0.001
  out[idx_other] <- sprintf(
    "%.3f",
    p[idx_other]
  )

  out
}


fmt_num <- function(
  x,
  digits = 2
) {

  out <- rep(
    NA_character_,
    length(x)
  )

  out[is.na(x)] <- "NA"

  idx <- !is.na(x)

  out[idx] <- formatC(
    x[idx],
    format = "f",
    digits = digits
  )

  out
}


fmt_signed <- function(
  x,
  digits = 2
) {

  out <- rep(
    NA_character_,
    length(x)
  )

  out[is.na(x)] <- "NA"

  idx <- !is.na(x)

  out[idx] <- paste0(
    ifelse(
      x[idx] > 0,
      "+",
      ""
    ),
    formatC(
      x[idx],
      format = "f",
      digits = digits
    )
  )

  out
}


# =============================================================================
# 5. Exposure specification
# =============================================================================

exposure_spec <- tibble(
  Exposure_variable = c(
    "Delta_SolarRad_60d_MJm2",
    "Delta_SolarRad_90d_MJm2",
    "Delta_Sunshine_60d_hr",
    "Delta_Sunshine_90d_hr"
  ),

  Exposure_family = c(
    "Solar radiation",
    "Solar radiation",
    "Sunshine duration",
    "Sunshine duration"
  ),

  Window = c(
    "60-day",
    "90-day",
    "60-day",
    "90-day"
  ),

  Raw_unit = c(
    "MJ/m²/day",
    "MJ/m²/day",
    "hr/day",
    "hr/day"
  )
)


# =============================================================================
# 6. Build model dataset + formula
# =============================================================================

make_model_data <- function(
  data,
  exposure,
  outcome_type
) {

  if (outcome_type == "biochemical") {

    vars <- c(
      "Delta_25OHD_ng_mL",
      exposure,
      "Baseline_25OHD_ng_mL",
      "Age_years",
      "Sex",
      "Clinical_followup_months"
    )

    d <- data %>%
      select(
        all_of(vars)
      ) %>%
      filter(
        complete.cases(.)
      )

    form <- reformulate(
      termlabels = c(
        exposure,
        "Baseline_25OHD_ng_mL",
        "Age_years",
        "Sex",
        "Clinical_followup_months"
      ),
      response = "Delta_25OHD_ng_mL"
    )

  } else if (outcome_type == "pain") {

    vars <- c(
      "Pain_intensity_reduction_VAS",
      exposure,
      "Baseline_pain_intensity_VAS",
      "Age_years",
      "Sex",
      "Clinical_followup_months"
    )

    d <- data %>%
      select(
        all_of(vars)
      ) %>%
      filter(
        complete.cases(.)
      )

    form <- reformulate(
      termlabels = c(
        exposure,
        "Baseline_pain_intensity_VAS",
        "Age_years",
        "Sex",
        "Clinical_followup_months"
      ),
      response = "Pain_intensity_reduction_VAS"
    )

  } else {

    stop(
      "outcome_type must be 'biochemical' or 'pain'."
    )
  }

  list(
    data = d,
    formula = form
  )
}


# =============================================================================
# 7. Analytic HC3 estimate
# =============================================================================

fit_hc3_effect <- function(
  data,
  exposure,
  outcome_type
) {

  obj <- make_model_data(
    data,
    exposure,
    outcome_type
  )

  d <- obj$data
  form <- obj$formula

  model <- lm(
    form,
    data = d
  )

  V <- sandwich::vcovHC(
    model,
    type = "HC3"
  )

  b <- coef(model)[exposure]
  se <- sqrt(diag(V))[exposure]

  df_res <- df.residual(model)
  crit <- qt(
    0.975,
    df = df_res
  )

  p <- 2 * pt(
    abs(b / se),
    df = df_res,
    lower.tail = FALSE
  )

  ci_low <- b - crit * se
  ci_high <- b + crit * se

  exposure_sd <- sd(
    d[[exposure]],
    na.rm = TRUE
  )

  tibble(
    Exposure_variable = exposure,
    Outcome_type = outcome_type,
    N = nrow(d),
    Exposure_SD = exposure_sd,

    Beta_raw = unname(b),
    CI_low_raw = unname(ci_low),
    CI_high_raw = unname(ci_high),

    Beta_per_SD = unname(
      b * exposure_sd
    ),

    CI_low_per_SD = unname(
      ci_low * exposure_sd
    ),

    CI_high_per_SD = unname(
      ci_high * exposure_sd
    ),

    p_value = unname(p)
  )
}


# =============================================================================
# 8. Bootstrap distribution
# =============================================================================

bootstrap_effect <- function(
  data,
  exposure,
  outcome_type,
  fixed_exposure_sd,
  B = BOOT_B
) {

  obj <- make_model_data(
    data,
    exposure,
    outcome_type
  )

  d <- obj$data
  form <- obj$formula

  n <- nrow(d)

  out <- rep(
    NA_real_,
    B
  )

  for (b in seq_len(B)) {

    idx <- sample.int(
      n,
      size = n,
      replace = TRUE
    )

    db <- d[
      idx,
      ,
      drop = FALSE
    ]

    fit_b <- try(
      lm(
        form,
        data = db
      ),
      silent = TRUE
    )

    if (inherits(fit_b, "try-error")) {
      next
    }

    co <- coef(fit_b)

    if (
      !exposure %in% names(co) ||
      !is.finite(co[exposure])
    ) {
      next
    }

    out[b] <- unname(
      co[exposure] *
        fixed_exposure_sd
    )
  }

  tibble(
    Bootstrap_iteration = seq_len(B),
    Beta_per_SD_boot = out
  ) %>%
    filter(
      is.finite(Beta_per_SD_boot)
    )
}


# =============================================================================
# 9. Run analytic models
# =============================================================================

biochem_results <- bind_rows(
  lapply(
    exposure_spec$Exposure_variable,
    function(v) {
      fit_hc3_effect(
        dat,
        v,
        "biochemical"
      )
    }
  )
) %>%
  left_join(
    exposure_spec,
    by = "Exposure_variable"
  ) %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )


pain_results <- bind_rows(
  lapply(
    exposure_spec$Exposure_variable,
    function(v) {
      fit_hc3_effect(
        dat,
        v,
        "pain"
      )
    }
  )
) %>%
  left_join(
    exposure_spec,
    by = "Exposure_variable"
  ) %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )


# =============================================================================
# 10. Bootstrap all 8 models
# =============================================================================

bootstrap_all <- function(
  analytic_results,
  outcome_type
) {

  rows <- lapply(
    seq_len(nrow(analytic_results)),
    function(i) {

      r <- analytic_results[i, ]

      cat(
        "Bootstrap ",
        outcome_type,
        ": ",
        r$Exposure_family,
        " / ",
        r$Window,
        "\n",
        sep = ""
      )

      bootstrap_effect(
        data = dat,
        exposure = r$Exposure_variable,
        outcome_type = outcome_type,
        fixed_exposure_sd = r$Exposure_SD,
        B = BOOT_B
      ) %>%
        mutate(
          Exposure_variable = r$Exposure_variable,
          Exposure_family = r$Exposure_family,
          Window = r$Window,
          Outcome_type = outcome_type
        )
    }
  )

  bind_rows(rows)
}


boot_biochem <- bootstrap_all(
  biochem_results,
  "biochemical"
)

boot_pain <- bootstrap_all(
  pain_results,
  "pain"
)


# =============================================================================
# 11. Positioning
# =============================================================================

row_position <- function(family) {

  ifelse(
    family == "Solar radiation",
    2,
    1
  )
}


window_offset <- function(window) {

  ifelse(
    window == "60-day",
    0.12,
    -0.12
  )
}


prepare_analytic_plot <- function(x) {

  x %>%
    mutate(
      Row_center = row_position(
        Exposure_family
      ),

      Y = Row_center +
        window_offset(
          Window
        ),

      Window = factor(
        Window,
        levels = c(
          "60-day",
          "90-day"
        )
      ),

      Stat_label = paste0(
        Window,
        ": β = ",
        fmt_signed(
          Beta_per_SD,
          2
        ),
        " (",
        fmt_num(
          CI_low_per_SD,
          2
        ),
        " to ",
        fmt_num(
          CI_high_per_SD,
          2
        ),
        "); p ",
        ifelse(
          p_value < 0.001,
          "<0.001",
          paste0(
            "= ",
            fmt_p(p_value)
          )
        ),
        "; q ",
        ifelse(
          q_value < 0.001,
          "<0.001",
          paste0(
            "= ",
            fmt_p(q_value)
          )
        )
      )
    )
}


prepare_boot_plot <- function(x) {

  set.seed(
    BOOT_SEED + 1
  )

  x %>%
    mutate(
      Row_center = row_position(
        Exposure_family
      ),

      Base_Y = Row_center +
        window_offset(
          Window
        ),

      Y_jitter = Base_Y +
        runif(
          n(),
          min = -0.055,
          max = 0.055
        )
    )
}


plot_biochem <- prepare_analytic_plot(
  biochem_results
)

plot_pain <- prepare_analytic_plot(
  pain_results
)

boot_biochem_plot <- prepare_boot_plot(
  boot_biochem
)

boot_pain_plot <- prepare_boot_plot(
  boot_pain
)


# =============================================================================
# 12. 60-day -> 90-day connector data
# =============================================================================

make_connector_data <- function(plot_dat) {

  plot_dat %>%
    select(
      Exposure_family,
      Window,
      Beta_per_SD,
      Y
    ) %>%
    pivot_wider(
      names_from = Window,
      values_from = c(
        Beta_per_SD,
        Y
      )
    ) %>%
    transmute(
      Exposure_family = Exposure_family,

      x60 = `Beta_per_SD_60-day`,
      y60 = `Y_60-day`,

      x90 = `Beta_per_SD_90-day`,
      y90 = `Y_90-day`
    )
}


conn_biochem <- make_connector_data(
  plot_biochem
)

conn_pain <- make_connector_data(
  plot_pain
)


# =============================================================================
# 13. Colors and theme
# =============================================================================

COL_TEXT <- "#252525"
COL_AXIS <- "#505050"
COL_ZERO <- "#8B8B8B"
COL_GUIDE <- "#C2C2C2"

# Refined 60-day amber/orange
COL_60 <- "#D97706"

# Refined 90-day green
COL_90 <- "#2F7D5A"

COL_BOOT <- "#A9A9A9"
COL_CONNECT <- "#737373"


theme_efig3 <- theme_classic(
  base_family = "Arial",
  base_size = 11
) +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),

    panel.background = element_rect(
      fill = "#FCFCFC",
      colour = NA
    ),

    axis.line = element_line(
      colour = COL_AXIS,
      linewidth = 0.55
    ),

    axis.ticks = element_line(
      colour = COL_AXIS,
      linewidth = 0.45
    ),

    axis.text = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 9.3
    ),

    axis.text.y = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 9.6
    ),

    axis.title = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.3,
      face = "bold"
    ),

    legend.title = element_blank(),

    legend.text = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 9.1
    ),

    plot.margin = margin(
      14,
      18,
      14,
      14
    )
  )


# =============================================================================
# 14. Refined panel builder
# =============================================================================

make_refined_panel <- function(
  analytic_dat,
  boot_dat,
  connector_dat,
  x_label
) {

  # Keep the numerical annotation outside the coefficient cloud.
  all_x <- c(
    boot_dat$Beta_per_SD_boot,
    analytic_dat$CI_low_per_SD,
    analytic_dat$CI_high_per_SD,
    0
  )

  x_min_data <- min(
    all_x,
    na.rm = TRUE
  )

  x_max_data <- max(
    all_x,
    na.rm = TRUE
  )

  span <- x_max_data - x_min_data

  if (
    !is.finite(span) ||
    span <= 0
  ) {
    span <- 1
  }

  text_x <- x_max_data +
    0.16 * span

  x_lower <- x_min_data -
    0.08 * span

  x_upper <- x_max_data +
    1.10 * span

  p <- ggplot() +

    # Two horizontal gray dotted guide lines, one per exposure family.
    geom_hline(
      yintercept = c(
        1,
        2
      ),
      linetype = "dashed",
      linewidth = 0.55,
      colour = COL_GUIDE
    ) +

    # Null effect.
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.65,
      colour = COL_ZERO
    ) +

    # Bootstrap distribution: all small gray points.
    geom_point(
      data = boot_dat,
      aes(
        x = Beta_per_SD_boot,
        y = Y_jitter
      ),
      shape = 16,
      size = 0.90,
      alpha = 0.16,
      colour = COL_BOOT
    ) +

    # Connect the two centers within each exposure family:
    # 60-day -> 90-day.
    geom_segment(
      data = connector_dat,
      aes(
        x = x60,
        y = y60,
        xend = x90,
        yend = y90
      ),
      linewidth = 1.05,
      colour = COL_CONNECT,
      alpha = 0.80
    ) +

    # HC3 95% CIs.
    geom_segment(
      data = analytic_dat,
      aes(
        x = CI_low_per_SD,
        xend = CI_high_per_SD,
        y = Y,
        yend = Y,
        colour = Window
      ),
      linewidth = 1.50,
      lineend = "round"
    ) +

    # CI end caps.
    geom_segment(
      data = analytic_dat,
      aes(
        x = CI_low_per_SD,
        xend = CI_low_per_SD,
        y = Y - 0.035,
        yend = Y + 0.035,
        colour = Window
      ),
      linewidth = 1.0
    ) +

    geom_segment(
      data = analytic_dat,
      aes(
        x = CI_high_per_SD,
        xend = CI_high_per_SD,
        y = Y - 0.035,
        yend = Y + 0.035,
        colour = Window
      ),
      linewidth = 1.0
    ) +

    # Analytic centers.
    geom_point(
      data = analytic_dat,
      aes(
        x = Beta_per_SD,
        y = Y,
        colour = Window,
        shape = Window
      ),
      size = 4.0,
      stroke = 1.0
    ) +

    # Numeric effect text.
    geom_text(
      data = analytic_dat,
      aes(
        x = text_x,
        y = Y,
        label = Stat_label,
        colour = Window
      ),
      hjust = 0,
      family = "Arial",
      size = 2.82,
      show.legend = FALSE
    ) +

    scale_colour_manual(
      values = c(
        "60-day" = COL_60,
        "90-day" = COL_90
      )
    ) +

    scale_shape_manual(
      values = c(
        "60-day" = 16,
        "90-day" = 17
      )
    ) +

    scale_y_continuous(
      breaks = c(
        1,
        2
      ),
      labels = c(
        "Sunshine duration",
        "Solar radiation"
      ),
      limits = c(
        0.55,
        2.45
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +

    scale_x_continuous(
      limits = c(
        x_lower,
        x_upper
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +

    labs(
      x = x_label,
      y = NULL
    ) +

    theme_efig3 +

    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )

  p
}


# =============================================================================
# 15. Panels
# =============================================================================

panel_a <- make_refined_panel(
  analytic_dat = plot_biochem,
  boot_dat = boot_biochem_plot,
  connector_dat = conn_biochem,
  x_label = "Adjusted change in Δ25(OH)D per 1-SD increase in exposure change, ng/mL"
)


panel_b <- make_refined_panel(
  analytic_dat = plot_pain,
  boot_dat = boot_pain_plot,
  connector_dat = conn_pain,
  x_label = "Adjusted change in pain reduction per 1-SD increase in exposure change, VAS"
)


# =============================================================================
# 16. Combine
# =============================================================================

efigure3 <- panel_a / panel_b +
  plot_layout(
    heights = c(
      1,
      1
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

print(efigure3)


# =============================================================================
# 17. Audit tables
# =============================================================================

all_effects <- bind_rows(
  biochem_results %>%
    mutate(
      Outcome = "Δ25(OH)D"
    ),

  pain_results %>%
    mutate(
      Outcome = "Pain reduction"
    )
) %>%
  select(
    Outcome,
    Exposure_family,
    Window,
    Exposure_variable,
    Raw_unit,
    N,
    Exposure_SD,
    Beta_raw,
    CI_low_raw,
    CI_high_raw,
    Beta_per_SD,
    CI_low_per_SD,
    CI_high_per_SD,
    p_value,
    q_value
  )


bootstrap_summary <- bind_rows(
  boot_biochem %>%
    mutate(
      Outcome = "Δ25(OH)D"
    ),

  boot_pain %>%
    mutate(
      Outcome = "Pain reduction"
    )
) %>%
  group_by(
    Outcome,
    Exposure_family,
    Window,
    Exposure_variable
  ) %>%
  summarise(
    Bootstrap_successful_N = n(),
    Bootstrap_mean = mean(
      Beta_per_SD_boot
    ),
    Bootstrap_SD = sd(
      Beta_per_SD_boot
    ),
    Bootstrap_median = median(
      Beta_per_SD_boot
    ),
    Bootstrap_2.5pct = quantile(
      Beta_per_SD_boot,
      0.025
    ),
    Bootstrap_97.5pct = quantile(
      Beta_per_SD_boot,
      0.975
    ),
    .groups = "drop"
  )


analysis_info <- tibble(
  Item = c(
    "Figure purpose",
    "Panel A",
    "Panel B",
    "60-day color",
    "90-day color",
    "60-day symbol",
    "90-day symbol",
    "Gray dots",
    "Colored horizontal bars",
    "Center connector",
    "Horizontal guide lines",
    "Vertical reference line",
    "Displayed effect scale",
    "Panel A model",
    "Panel B model",
    "Multiplicity",
    "Bootstrap repetitions",
    "Bootstrap role",
    "Panel letters",
    "Solar-radiation interpretation",
    "R version requested for manuscript"
  ),

  Value = c(
    paste(
      "Compare 60-day versus 90-day environmental windows while showing",
      "the sampling distribution of each adjusted effect"
    ),

    "Environmental change versus change in serum 25(OH)D",

    "Environmental change versus pain reduction",

    COL_60,

    COL_90,

    "Circle",

    "Triangle",

    paste(
      "Bootstrap distributions of adjusted beta estimates;",
      "all successful bootstrap estimates displayed"
    ),

    "Analytic HC3 robust 95% confidence intervals",

    "Gray line connects analytic 60-day and 90-day centers within each exposure family",

    "Two gray dashed horizontal lines mark solar-radiation and sunshine-duration rows",

    "Dashed vertical line at beta = 0",

    "Adjusted beta per 1-SD increase in environmental-change exposure",

    paste(
      "Delta 25(OH)D ~ environmental change + baseline 25(OH)D",
      "+ age + sex + clinical follow-up interval"
    ),

    paste(
      "Pain reduction ~ environmental change + baseline VAS",
      "+ age + sex + clinical follow-up interval"
    ),

    paste(
      "Benjamini-Hochberg q-values calculated separately across",
      "the four biochemical and four pain comparisons"
    ),

    as.character(
      BOOT_B
    ),

    paste(
      "Bootstrap dots are visual distribution aids only;",
      "HC3 confidence intervals and p-values remain the inferential results"
    ),

    "Not drawn inside figure; add (A) and (B) during manuscript layout",

    "Ambient solar-energy proxy; not direct individual UVB exposure",

    "R 4.5.1"
  )
)


# =============================================================================
# 18. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Supplementary_Figure_S5_environmental_sensitivity.png"
)

tiff_file <- file.path(
  output_dir,
  "Supplementary_Figure_S5_environmental_sensitivity.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Figure_S5_environmental_sensitivity.pdf"
)

ggsave(
  filename = png_file,
  plot = efigure3,
  width = 12.0,
  height = 8.2,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = efigure3,
  width = 12.0,
  height = 8.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = efigure3,
  width = 12.0,
  height = 8.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 19. Export audit workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Supplementary_Figure_S5_environmental_sensitivity_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "All_Analytic_Effects"
)

addWorksheet(
  wb,
  "Biochemical"
)

addWorksheet(
  wb,
  "Pain"
)

addWorksheet(
  wb,
  "Bootstrap_Summary"
)

addWorksheet(
  wb,
  "Bootstrap_Biochemical"
)

addWorksheet(
  wb,
  "Bootstrap_Pain"
)

addWorksheet(
  wb,
  "Analysis_Info"
)

writeData(
  wb,
  "All_Analytic_Effects",
  all_effects
)

writeData(
  wb,
  "Biochemical",
  biochem_results
)

writeData(
  wb,
  "Pain",
  pain_results
)

writeData(
  wb,
  "Bootstrap_Summary",
  bootstrap_summary
)

writeData(
  wb,
  "Bootstrap_Biochemical",
  boot_biochem
)

writeData(
  wb,
  "Bootstrap_Pain",
  boot_pain
)

writeData(
  wb,
  "Analysis_Info",
  analysis_info
)

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sh in names(wb)) {

  temp <- readWorkbook(
    wb,
    sheet = sh
  )

  if (ncol(temp) > 0) {

    addStyle(
      wb,
      sh,
      header_style,
      rows = 1,
      cols = seq_len(
        ncol(temp)
      ),
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
      cols = seq_len(
        ncol(temp)
      ),
      widths = "auto"
    )
  }
}

setColWidths(
  wb,
  "Analysis_Info",
  cols = 1:2,
  widths = c(
    40,
    110
  )
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 20. Console output
# =============================================================================

cat("\n========================================================\n")
cat("Supplementary Figure S5 sensitivity analysis completed\n")
cat("========================================================\n\n")

cat("Analytic effects:\n")
print(
  all_effects
)

cat("\nBootstrap summary:\n")
print(
  bootstrap_summary
)

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("========================================================\n")
