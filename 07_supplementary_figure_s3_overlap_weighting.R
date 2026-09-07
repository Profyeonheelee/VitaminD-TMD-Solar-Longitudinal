# =============================================================================
# Supplementary Figure S3. Propensity-score overlap-weighting diagnostics
#
# Panels
#   A. Propensity-score distributions in the low-vitamin-D cohort
#   B. Covariate balance before and after overlap weighting (Love plot)
#
# Final propensity-score model
#   Vitamin D prescription ~
#     age
#     + sex
#     + baseline serum 25(OH)D
#     + log1p(symptom duration)
#     + baseline pain intensity VAS
#     + clinically detected TMJ noise
#     + trauma history
#     + baseline 60-day mean solar radiation
#
# Overlap weights
#   Rx = 1: weight = 1 - PS
#   Rx = 0: weight = PS
#
# Expected low-vitamin-D cohort:
#   Total N = 200
#   No prescription = 51
#   Prescription = 149
#
# Figure style
#   - No overall title inside the image
#   - No panel letters inside the image
#   - Grayscale + muted orange; post-weighting balance markers are vivid blue
#   - Absolute standardized mean difference threshold = 0.10
#
# Outputs
#   Supplementary_Figure_S3_overlap_weighting.png
#   Supplementary_Figure_S3_overlap_weighting.tiff
#   Supplementary_Figure_S3_overlap_weighting.pdf
#   Supplementary_Figure_S3_overlap_weighting_statistics.xlsx
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
library(openxlsx)


# =============================================================================
# 1. File path
# =============================================================================

xlsx_path <- file.path(
  "data",
  "Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx"
)

sheet_name <- "Analysis_Data"
output_dir <- file.path("outputs", "Supplementary_Figure_S3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(xlsx_path)) {
  stop("Data file not found: ", xlsx_path)
}


# =============================================================================
# 2. Import
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
  "Vitamin_D_prescription",
  "Low_Vitamin_D_cohort",
  "Baseline_pain_intensity_VAS",
  "TMJ_noise",
  "Trauma_history",
  "SolarRad_preBL_60d_mean_MJm2"
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


# =============================================================================
# 3. Binary-variable conversion helper
# =============================================================================

to_binary01 <- function(x, var_name) {

  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x) || is.integer(x)) {

    ux <- sort(unique(na.omit(x)))

    if (all(ux %in% c(0, 1))) {
      return(as.integer(x))
    }

    if (all(ux %in% c(1, 2))) {
      return(as.integer(x == 2))
    }

    stop(
      var_name,
      " is numeric but not coded as 0/1 or 1/2. Unique values: ",
      paste(ux, collapse = ", ")
    )
  }

  z <- trimws(
    tolower(
      as.character(x)
    )
  )

  out <- rep(
    NA_integer_,
    length(z)
  )

  no_values <- c(
    "0",
    "no",
    "n",
    "absent",
    "negative",
    "none",
    "false"
  )

  yes_values <- c(
    "1",
    "yes",
    "y",
    "present",
    "positive",
    "true"
  )

  out[z %in% no_values] <- 0L
  out[z %in% yes_values] <- 1L

  unresolved <- unique(
    z[
      !is.na(z) &
        is.na(out)
    ]
  )

  if (length(unresolved) > 0) {
    stop(
      "Could not convert ",
      var_name,
      " to 0/1. Unresolved value(s): ",
      paste(unresolved, collapse = ", ")
    )
  }

  out
}


# =============================================================================
# 4. Sex conversion
# =============================================================================

to_female01 <- function(x) {

  if (is.numeric(x) || is.integer(x)) {

    ux <- sort(unique(na.omit(x)))

    if (all(ux %in% c(0, 1))) {
      # Preserve existing 0/1 coding.
      return(as.integer(x))
    }

    if (all(ux %in% c(1, 2))) {
      # Common coding: 1 male, 2 female.
      return(as.integer(x == 2))
    }
  }

  z <- trimws(
    tolower(
      as.character(x)
    )
  )

  out <- rep(
    NA_integer_,
    length(z)
  )

  male_values <- c(
    "male",
    "m",
    "man",
    "남",
    "남성"
  )

  female_values <- c(
    "female",
    "f",
    "woman",
    "여",
    "여성"
  )

  out[z %in% male_values] <- 0L
  out[z %in% female_values] <- 1L

  unresolved <- unique(
    z[
      !is.na(z) &
        is.na(out)
    ]
  )

  if (length(unresolved) > 0) {
    stop(
      "Could not convert Sex to female=1/male=0. Unresolved value(s): ",
      paste(unresolved, collapse = ", ")
    )
  }

  out
}


# =============================================================================
# 5. Low-vitamin-D cohort and final PS variables
# =============================================================================

low_flag <- to_binary01(
  dat0$Low_Vitamin_D_cohort,
  "Low_Vitamin_D_cohort"
)

rx_all <- to_binary01(
  dat0$Vitamin_D_prescription,
  "Vitamin_D_prescription"
)

tmj_noise_all <- to_binary01(
  dat0$TMJ_noise,
  "TMJ_noise"
)

trauma_all <- to_binary01(
  dat0$Trauma_history,
  "Trauma_history"
)

female_all <- to_female01(
  dat0$Sex
)

dat <- dat0 %>%
  mutate(
    Low_flag = low_flag,
    Rx = rx_all,
    Female = female_all,
    TMJ_noise_bin = tmj_noise_all,
    Trauma_bin = trauma_all,
    log_symptom_duration = log1p(
      Symptom_duration_months
    )
  ) %>%
  filter(
    Low_flag == 1
  ) %>%
  transmute(
    Study_ID = Study_ID,
    Rx = Rx,
    Age_years = Age_years,
    Female = Female,
    Baseline_25OHD_ng_mL = Baseline_25OHD_ng_mL,
    log_symptom_duration = log_symptom_duration,
    Baseline_pain_intensity_VAS = Baseline_pain_intensity_VAS,
    TMJ_noise = TMJ_noise_bin,
    Trauma_history = Trauma_bin,
    SolarRad_preBL_60d_mean_MJm2 = SolarRad_preBL_60d_mean_MJm2
  )


# =============================================================================
# 6. Complete-case PS dataset
# =============================================================================

ps_vars <- c(
  "Rx",
  "Age_years",
  "Female",
  "Baseline_25OHD_ng_mL",
  "log_symptom_duration",
  "Baseline_pain_intensity_VAS",
  "TMJ_noise",
  "Trauma_history",
  "SolarRad_preBL_60d_mean_MJm2"
)

dat_ps <- dat %>%
  filter(
    complete.cases(
      across(
        all_of(ps_vars)
      )
    )
  )

cat("\n===========================================\n")
cat("Supplementary Figure S3: PS overlap-weighting diagnostics\n")
cat("===========================================\n")
cat("Low-vitamin-D cohort before CC restriction: ", nrow(dat), "\n", sep = "")
cat("PS complete-case N: ", nrow(dat_ps), "\n", sep = "")
cat(
  "No Rx: ",
  sum(dat_ps$Rx == 0),
  "; Rx: ",
  sum(dat_ps$Rx == 1),
  "\n",
  sep = ""
)


# =============================================================================
# 7. Expected-N audit
# =============================================================================

if (nrow(dat) != 200) {
  warning(
    "Expected low-vitamin-D cohort N=200; found N=",
    nrow(dat)
  )
}

if (
  sum(dat$Rx == 0, na.rm = TRUE) != 51 ||
  sum(dat$Rx == 1, na.rm = TRUE) != 149
) {
  warning(
    "Expected low-vitamin-D group counts 51 no-Rx / 149 Rx. Found ",
    sum(dat$Rx == 0, na.rm = TRUE),
    " / ",
    sum(dat$Rx == 1, na.rm = TRUE),
    "."
  )
}


# =============================================================================
# 8. Propensity-score model
# =============================================================================

ps_model <- glm(
  Rx ~
    Age_years +
    Female +
    Baseline_25OHD_ng_mL +
    log_symptom_duration +
    Baseline_pain_intensity_VAS +
    TMJ_noise +
    Trauma_history +
    SolarRad_preBL_60d_mean_MJm2,
  data = dat_ps,
  family = binomial(
    link = "logit"
  )
)

dat_ps <- dat_ps %>%
  mutate(
    PS = predict(
      ps_model,
      type = "response"
    ),

    OW = ifelse(
      Rx == 1,
      1 - PS,
      PS
    ),

    Exposure = factor(
      Rx,
      levels = c(0, 1),
      labels = c(
        "Low VitD / no Rx",
        "Low VitD / Rx"
      )
    )
  )


# =============================================================================
# 9. Positivity audit
# =============================================================================

if (
  any(dat_ps$PS <= 0) ||
  any(dat_ps$PS >= 1)
) {
  warning(
    "At least one estimated propensity score is exactly 0 or 1."
  )
}


# =============================================================================
# 10. Weighted summary helpers
# =============================================================================

weighted_mean <- function(x, w) {

  ok <- is.finite(x) &
    is.finite(w) &
    !is.na(x) &
    !is.na(w)

  sum(
    w[ok] * x[ok]
  ) /
    sum(w[ok])
}


weighted_var <- function(x, w) {

  ok <- is.finite(x) &
    is.finite(w) &
    !is.na(x) &
    !is.na(w)

  x <- x[ok]
  w <- w[ok]

  m <- sum(w * x) /
    sum(w)

  sum(
    w * (x - m)^2
  ) /
    sum(w)
}


smd_continuous <- function(
  x,
  treatment,
  weights = NULL
) {

  if (is.null(weights)) {
    weights <- rep(
      1,
      length(x)
    )
  }

  x0 <- x[
    treatment == 0
  ]

  x1 <- x[
    treatment == 1
  ]

  w0 <- weights[
    treatment == 0
  ]

  w1 <- weights[
    treatment == 1
  ]

  m0 <- weighted_mean(
    x0,
    w0
  )

  m1 <- weighted_mean(
    x1,
    w1
  )

  v0 <- weighted_var(
    x0,
    w0
  )

  v1 <- weighted_var(
    x1,
    w1
  )

  pooled_sd <- sqrt(
    (v0 + v1) / 2
  )

  if (
    !is.finite(pooled_sd) ||
    pooled_sd == 0
  ) {
    return(0)
  }

  (m1 - m0) / pooled_sd
}


smd_binary <- function(
  x,
  treatment,
  weights = NULL
) {

  if (is.null(weights)) {
    weights <- rep(
      1,
      length(x)
    )
  }

  x0 <- x[
    treatment == 0
  ]

  x1 <- x[
    treatment == 1
  ]

  w0 <- weights[
    treatment == 0
  ]

  w1 <- weights[
    treatment == 1
  ]

  p0 <- weighted_mean(
    x0,
    w0
  )

  p1 <- weighted_mean(
    x1,
    w1
  )

  pooled_sd <- sqrt(
    (
      p0 * (1 - p0) +
        p1 * (1 - p1)
    ) / 2
  )

  if (
    !is.finite(pooled_sd) ||
    pooled_sd == 0
  ) {
    return(0)
  }

  (p1 - p0) / pooled_sd
}


# =============================================================================
# 11. Covariate balance
# =============================================================================

balance_spec <- tibble(
  Variable = c(
    "Age_years",
    "Female",
    "Baseline_25OHD_ng_mL",
    "log_symptom_duration",
    "Baseline_pain_intensity_VAS",
    "TMJ_noise",
    "Trauma_history",
    "SolarRad_preBL_60d_mean_MJm2"
  ),

  Label = c(
    "Age",
    "Female sex",
    "Baseline 25(OH)D",
    "Log symptom duration",
    "Baseline pain intensity",
    "Clinically detected TMJ noise",
    "Trauma history",
    "Baseline 60-day solar radiation"
  ),

  Type = c(
    "continuous",
    "binary",
    "continuous",
    "continuous",
    "continuous",
    "binary",
    "binary",
    "continuous"
  )
)


balance_rows <- lapply(
  seq_len(nrow(balance_spec)),
  function(i) {

    v <- balance_spec$Variable[i]
    type <- balance_spec$Type[i]

    x <- dat_ps[[v]]

    if (type == "binary") {

      smd_before <- smd_binary(
        x,
        dat_ps$Rx,
        weights = rep(
          1,
          nrow(dat_ps)
        )
      )

      smd_after <- smd_binary(
        x,
        dat_ps$Rx,
        weights = dat_ps$OW
      )

    } else {

      smd_before <- smd_continuous(
        x,
        dat_ps$Rx,
        weights = rep(
          1,
          nrow(dat_ps)
        )
      )

      smd_after <- smd_continuous(
        x,
        dat_ps$Rx,
        weights = dat_ps$OW
      )
    }

    tibble(
      Variable = v,
      Label = balance_spec$Label[i],
      Type = type,
      SMD_before = smd_before,
      Abs_SMD_before = abs(
        smd_before
      ),
      SMD_after = smd_after,
      Abs_SMD_after = abs(
        smd_after
      )
    )
  }
)

balance <- bind_rows(
  balance_rows
)


# =============================================================================
# 12. Effective sample size
# =============================================================================

ess <- function(w) {

  w <- w[
    is.finite(w) &
      !is.na(w)
  ]

  (sum(w)^2) /
    sum(w^2)
}

ess_table <- tibble(
  Group = c(
    "Low VitD / no Rx",
    "Low VitD / Rx",
    "Overall"
  ),

  Unweighted_N = c(
    sum(dat_ps$Rx == 0),
    sum(dat_ps$Rx == 1),
    nrow(dat_ps)
  ),

  Overlap_weight_ESS = c(
    ess(
      dat_ps$OW[
        dat_ps$Rx == 0
      ]
    ),
    ess(
      dat_ps$OW[
        dat_ps$Rx == 1
      ]
    ),
    ess(dat_ps$OW)
  )
)


# =============================================================================
# 13. PS distribution summary
# =============================================================================

ps_summary <- dat_ps %>%
  group_by(
    Exposure
  ) %>%
  summarise(
    N = n(),
    PS_mean = mean(PS),
    PS_SD = sd(PS),
    PS_median = median(PS),
    PS_Q1 = quantile(
      PS,
      0.25
    ),
    PS_Q3 = quantile(
      PS,
      0.75
    ),
    PS_min = min(PS),
    PS_max = max(PS),
    .groups = "drop"
  )



# =============================================================================
# 13B. Concise numeric annotations for the figure
# =============================================================================

ps_annot_no_rx <- ps_summary %>%
  filter(
    Exposure == "Low VitD / no Rx"
  )

ps_annot_rx <- ps_summary %>%
  filter(
    Exposure == "Low VitD / Rx"
  )

panel_a_annotation <- paste0(
  "No Rx: n = ",
  ps_annot_no_rx$N,
  ", mean PS = ",
  sprintf(
    "%.3f",
    ps_annot_no_rx$PS_mean
  ),
  "\nRx: n = ",
  ps_annot_rx$N,
  ", mean PS = ",
  sprintf(
    "%.3f",
    ps_annot_rx$PS_mean
  )
)

max_smd_before <- max(
  balance$Abs_SMD_before,
  na.rm = TRUE
)

max_smd_after <- max(
  balance$Abs_SMD_after,
  na.rm = TRUE
)

n_imbalanced_before <- sum(
  balance$Abs_SMD_before > 0.10,
  na.rm = TRUE
)

n_imbalanced_after <- sum(
  balance$Abs_SMD_after > 0.10,
  na.rm = TRUE
)

panel_b_annotation <- paste0(
  "Max |SMD|: ",
  sprintf(
    "%.3f",
    max_smd_before
  ),
  " \u2192 ",
  sprintf(
    "%.3f",
    max_smd_after
  ),
  "\n|SMD| > 0.10: ",
  n_imbalanced_before,
  " \u2192 ",
  n_imbalanced_after,
  " covariates"
)


# =============================================================================
# 14. Visual system
# =============================================================================

COL_TEXT <- "#252525"
COL_AXIS <- "#505050"

COL_GRAY <- "#5B5B5B"
COL_GRAY_LIGHT <- "#BEBEBE"

COL_ORANGE <- "#C96A1B"
COL_ORANGE_LIGHT <- "#E8B07A"

COL_BLUE <- "#0066FF"

theme_efig <- theme_classic(
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
      size = 9.4
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
      size = 9.2,
      colour = COL_TEXT
    ),

    plot.margin = margin(
      10, 12, 10, 12
    )
  )


# =============================================================================
# 15. Panel A: propensity-score distributions
# =============================================================================

panel_a <- ggplot(
  dat_ps,
  aes(
    x = PS,
    colour = Exposure,
    fill = Exposure
  )
) +

  geom_density(
    aes(
      y = after_stat(density)
    ),
    linewidth = 1.15,
    alpha = 0.16,
    adjust = 1.0
  ) +

  geom_rug(
    aes(
      colour = Exposure
    ),
    sides = "b",
    alpha = 0.22,
    linewidth = 0.35
  ) +

  annotate(
    "label",
    x = 0.02,
    y = Inf,
    label = panel_a_annotation,
    hjust = 0,
    vjust = 1.10,
    size = 3.0,
    family = "Arial",
    lineheight = 1.08,
    colour = COL_TEXT,
    fill = "white",
    label.size = 0.28,
    label.padding = grid::unit(
      0.35,
      "lines"
    ),
    label.r = grid::unit(
      0.05,
      "lines"
    )
  ) +

  scale_colour_manual(
    values = c(
      "Low VitD / no Rx" = COL_GRAY,
      "Low VitD / Rx" = COL_ORANGE
    )
  ) +

  scale_fill_manual(
    values = c(
      "Low VitD / no Rx" = COL_GRAY_LIGHT,
      "Low VitD / Rx" = COL_ORANGE_LIGHT
    )
  ) +

  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +

  labs(
    x = "Propensity score",
    y = "Density"
  ) +

  theme_efig +
  theme(
    legend.position = c(
      0.02,
      0.82
    ),
    legend.justification = c(
      0,
      1
    ),
    legend.background = element_rect(
      fill = "white",
      colour = "#D0D0D0",
      linewidth = 0.30
    ),
    legend.key = element_blank()
  )


# =============================================================================
# 16. Panel B: Love plot
# =============================================================================

love_dat <- balance %>%
  select(
    Label,
    Abs_SMD_before,
    Abs_SMD_after
  ) %>%
  pivot_longer(
    cols = c(
      Abs_SMD_before,
      Abs_SMD_after
    ),
    names_to = "Stage",
    values_to = "Abs_SMD"
  ) %>%
  mutate(
    Stage = factor(
      Stage,
      levels = c(
        "Abs_SMD_before",
        "Abs_SMD_after"
      ),
      labels = c(
        "Before weighting",
        "After overlap weighting"
      )
    ),

    Label = factor(
      Label,
      levels = rev(
        balance$Label[
          order(
            balance$Abs_SMD_before,
            decreasing = FALSE
          )
        ]
      )
    )
  )

max_x <- max(
  0.35,
  ceiling(
    max(
      love_dat$Abs_SMD,
      na.rm = TRUE
    ) * 20
  ) / 20 + 0.05
)

panel_b <- ggplot(
  love_dat,
  aes(
    x = Abs_SMD,
    y = Label,
    colour = Stage,
    shape = Stage
  )
) +

  geom_vline(
    xintercept = 0.10,
    linetype = "dashed",
    linewidth = 0.65,
    colour = "#888888"
  ) +

  geom_segment(
    data = balance,
    aes(
      x = Abs_SMD_after,
      xend = Abs_SMD_before,
      y = Label,
      yend = Label
    ),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "#C8C8C8"
  ) +

  geom_point(
    size = 3.1,
    stroke = 0.85
  ) +

  scale_colour_manual(
    values = c(
      "Before weighting" = COL_GRAY,
      "After overlap weighting" = COL_BLUE
    )
  ) +

  scale_shape_manual(
    values = c(
      "Before weighting" = 16,
      "After overlap weighting" = 17
    )
  ) +

  scale_x_continuous(
    limits = c(
      0,
      max_x
    ),
    breaks = pretty(
      c(
        0,
        max_x
      ),
      n = 6
    ),
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +

  labs(
    x = "Absolute standardized mean difference",
    y = NULL,
    subtitle = panel_b_annotation
  ) +

  theme_efig +
  theme(
    axis.text.y = element_text(
      family = "Arial",
      size = 9.1,
      colour = COL_TEXT
    ),

    plot.subtitle = element_text(
      family = "Arial",
      size = 9.3,
      face = "bold",
      colour = COL_TEXT,
      hjust = 1,
      lineheight = 1.08,
      margin = margin(
        b = 7
      )
    ),

    legend.position = "bottom",

    legend.direction = "horizontal",

    legend.box = "horizontal"
  )


# =============================================================================
# 17. Combine
# =============================================================================

efigure1 <- panel_a + panel_b +
  plot_layout(
    widths = c(
      1.00,
      1.15
    ),
    guides = "keep"
  )

print(efigure1)


# =============================================================================
# 18. Diagnostics / console audit
# =============================================================================

cat("\nPropensity-score model coefficients:\n")
print(
  summary(ps_model)$coefficients
)

cat("\nCovariate balance:\n")
print(
  balance %>%
    select(
      Label,
      SMD_before,
      Abs_SMD_before,
      SMD_after,
      Abs_SMD_after
    )
)

cat(
  "\nMaximum absolute SMD before weighting = ",
  sprintf(
    "%.4f",
    max(
      balance$Abs_SMD_before
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum absolute SMD after overlap weighting = ",
  sprintf(
    "%.4f",
    max(
      balance$Abs_SMD_after
    )
  ),
  "\n",
  sep = ""
)

cat("\nEffective sample size:\n")
print(ess_table)


# =============================================================================
# 19. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Supplementary_Figure_S3_overlap_weighting.png"
)

tiff_file <- file.path(
  output_dir,
  "Supplementary_Figure_S3_overlap_weighting.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Figure_S3_overlap_weighting.pdf"
)

ggsave(
  filename = png_file,
  plot = efigure1,
  width = 11.0,
  height = 5.6,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = efigure1,
  width = 11.0,
  height = 5.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = efigure1,
  width = 11.0,
  height = 5.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 20. Export audit workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Supplementary_Figure_S3_overlap_weighting_statistics.xlsx"
)

ps_model_table <- as.data.frame(
  summary(ps_model)$coefficients
) %>%
  tibble::rownames_to_column(
    "Term"
  )

analysis_info <- tibble(
  Item = c(
    "Analysis cohort",
    "Expected cohort N",
    "Exposure",
    "PS model",
    "Overlap weight for Rx",
    "Overlap weight for no Rx",
    "Balance metric",
    "Balance threshold",
    "Panel A",
    "Panel A annotation",
    "Panel B",
    "Panel B annotation",
    "Missing data",
    "R version requested for manuscript"
  ),

  Value = c(
    "Low-vitamin-D cohort",
    "200 (51 no Rx; 149 Rx)",
    "Recorded vitamin D prescription",
    paste(
      "Age + sex + baseline 25(OH)D + log1p symptom duration",
      "+ baseline VAS + clinically detected TMJ noise",
      "+ trauma history + baseline 60-day solar radiation"
    ),
    "1 - propensity score",
    "Propensity score",
    "Absolute standardized mean difference",
    "0.10",
    "Unweighted propensity-score density distributions",
    "Group n and mean propensity score",
    "Absolute SMD before and after overlap weighting",
    "Maximum absolute SMD and number of covariates with |SMD| > 0.10 before and after weighting; post-weighting markers shown in vivid blue",
    "Complete-case PS estimation; no imputation",
    "R 4.5.1"
  )
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "PS_Model"
)

addWorksheet(
  wb,
  "PS_Distribution"
)

addWorksheet(
  wb,
  "Covariate_Balance"
)

addWorksheet(
  wb,
  "Effective_Sample_Size"
)

addWorksheet(
  wb,
  "Patient_Level_Diagnostics"
)

addWorksheet(
  wb,
  "Analysis_Info"
)

writeData(
  wb,
  "PS_Model",
  ps_model_table
)

writeData(
  wb,
  "PS_Distribution",
  ps_summary
)

writeData(
  wb,
  "Covariate_Balance",
  balance
)

writeData(
  wb,
  "Effective_Sample_Size",
  ess_table
)

writeData(
  wb,
  "Patient_Level_Diagnostics",
  dat_ps %>%
    select(
      Study_ID,
      Exposure,
      PS,
      OW
    )
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
    36,
    105
  )
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 21. Final console output
# =============================================================================

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("===========================================\n")
