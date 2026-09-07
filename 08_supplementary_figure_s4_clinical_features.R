# =============================================================================
# Supplementary Figure S4. Baseline serum 25(OH)D and TMD-related clinical features
#
# Panels
#   A. Unadjusted associations across 12 prespecified TMD-related features
#      -> Spearman rho with BH-FDR q-values
#
#   B. Adjusted associations with continuous clinical outcomes
#      -> beta per +10-ng/mL higher baseline serum 25(OH)D
#
#   C. Adjusted associations with binary clinical features
#      -> odds ratio per +10-ng/mL higher baseline serum 25(OH)D
#
# Final feature set
#   Continuous:
#     Pain intensity VAS
#     PFO
#     MUO
#
#   Binary:
#     Self-reported bruxism
#     Clenching
#     Parafunction
#     Clinically detected TMJ noise
#     Self-reported TMJ noise
#     Clinically detected jaw locking
#     Self-reported jaw locking
#     Stiffness
#     Trauma history
#
# Adjusted models
#   Continuous outcome ~ baseline 25(OH)D per 10 ng/mL
#                        + age + sex
#                        + log1p(symptom duration)
#                        + baseline 60-day mean solar radiation
#
#   Binary outcome ~ baseline 25(OH)D per 10 ng/mL
#                    + age + sex
#                    + log1p(symptom duration)
#                    + baseline 60-day mean solar radiation
#
# HC3 heteroskedasticity-robust standard errors are used for all adjusted
# linear and logistic models.
#
# Multiplicity
#   - Unadjusted BH q-values across all 12 prespecified features.
#   - Adjusted BH q-values across all 12 prespecified adjusted models.
#
# Figure style
#   - No overall title inside figure
#   - No panel letters inside figure
#   - Panel A: muted violet highlights q < 0.05; otherwise gray
#   - Panel B/C: muted teal highlights adjusted p < 0.05; otherwise gray
#   - Numeric annotation includes effect estimate, 95% CI, p, and q
#
# Outputs
#   Supplementary_Figure_S4_clinical_features.png
#   Supplementary_Figure_S4_clinical_features.tiff
#   Supplementary_Figure_S4_clinical_features.pdf
#   Supplementary_Figure_S4_clinical_features_statistics.xlsx
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
output_dir <- file.path("outputs", "Supplementary_Figure_S4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(xlsx_path)) {
  stop("Data file not found: ", xlsx_path)
}


# =============================================================================
# 2. Import and validation
# =============================================================================

dat0 <- read_xlsx(
  xlsx_path,
  sheet = sheet_name
)

required_vars <- c(
  "Study_ID",
  "Sex",
  "Age_years",
  "Symptom_duration_months",
  "Baseline_25OHD_ng_mL",
  "Baseline_pain_intensity_VAS",
  "Baseline_CMO_mm",
  "Baseline_MMO_mm",
  "Bruxism",
  "Clenching",
  "Parafunction",
  "TMJ_noise",
  "Self_reported_TMJ_noise",
  "Objective_locking",
  "Self_reported_locking",
  "Stiffness",
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

if (nrow(dat0) != 241) {
  warning(
    "Expected N=241; found N=",
    nrow(dat0)
  )
}

if (anyDuplicated(dat0$Study_ID) > 0) {
  stop("Duplicate Study_ID detected.")
}


# =============================================================================
# 3. Variable conversion helpers
# =============================================================================

to_binary01 <- function(x, var_name) {

  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x) || is.integer(x)) {

    ux <- sort(
      unique(
        na.omit(x)
      )
    )

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


to_female01 <- function(x) {

  if (is.numeric(x) || is.integer(x)) {

    ux <- sort(
      unique(
        na.omit(x)
      )
    )

    if (all(ux %in% c(0, 1))) {
      return(as.integer(x))
    }

    if (all(ux %in% c(1, 2))) {
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
# 4. Analysis dataset
# =============================================================================

dat <- dat0 %>%
  mutate(
    Female = to_female01(
      Sex
    ),

    Baseline_25OHD_per10 =
      Baseline_25OHD_ng_mL / 10,

    log_symptom_duration =
      log1p(
        Symptom_duration_months
      ),

    Bruxism_bin = to_binary01(
      Bruxism,
      "Bruxism"
    ),

    Clenching_bin = to_binary01(
      Clenching,
      "Clenching"
    ),

    Parafunction_bin = to_binary01(
      Parafunction,
      "Parafunction"
    ),

    TMJ_noise_bin = to_binary01(
      TMJ_noise,
      "TMJ_noise"
    ),

    Self_reported_TMJ_noise_bin = to_binary01(
      Self_reported_TMJ_noise,
      "Self_reported_TMJ_noise"
    ),

    Objective_locking_bin = to_binary01(
      Objective_locking,
      "Objective_locking"
    ),

    Self_reported_locking_bin = to_binary01(
      Self_reported_locking,
      "Self_reported_locking"
    ),

    Stiffness_bin = to_binary01(
      Stiffness,
      "Stiffness"
    ),

    Trauma_history_bin = to_binary01(
      Trauma_history,
      "Trauma_history"
    )
  )


# =============================================================================
# 5. Feature specification
# =============================================================================

feature_spec <- tibble(
  Variable = c(
    "Baseline_pain_intensity_VAS",
    "Baseline_CMO_mm",
    "Baseline_MMO_mm",
    "Bruxism_bin",
    "Clenching_bin",
    "Parafunction_bin",
    "TMJ_noise_bin",
    "Self_reported_TMJ_noise_bin",
    "Objective_locking_bin",
    "Self_reported_locking_bin",
    "Stiffness_bin",
    "Trauma_history_bin"
  ),

  Label = c(
    "Pain intensity VAS",
    "PFO",
    "MUO",
    "Self-reported bruxism",
    "Clenching",
    "Parafunction",
    "Clinically detected TMJ noise",
    "Self-reported TMJ noise",
    "Clinically detected jaw locking",
    "Self-reported jaw locking",
    "Stiffness",
    "Trauma history"
  ),

  Type = c(
    "continuous",
    "continuous",
    "continuous",
    rep(
      "binary",
      9
    )
  )
)


# =============================================================================
# 6. Formatting helpers
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
# 7. Unadjusted Spearman associations
# =============================================================================

unadjusted_rows <- lapply(
  seq_len(
    nrow(feature_spec)
  ),
  function(i) {

    v <- feature_spec$Variable[i]

    d <- dat %>%
      select(
        Baseline_25OHD_ng_mL,
        all_of(v)
      ) %>%
      filter(
        complete.cases(.)
      )

    ct <- suppressWarnings(
      cor.test(
        d$Baseline_25OHD_ng_mL,
        d[[v]],
        method = "spearman",
        exact = FALSE
      )
    )

    tibble(
      Variable = v,
      Label = feature_spec$Label[i],
      Type = feature_spec$Type[i],
      N = nrow(d),
      Spearman_rho = unname(
        ct$estimate
      ),
      p_unadjusted = ct$p.value
    )
  }
)

unadjusted <- bind_rows(
  unadjusted_rows
) %>%
  mutate(
    q_unadjusted = p.adjust(
      p_unadjusted,
      method = "BH"
    )
  )


# =============================================================================
# 8. HC3 adjusted continuous models
# =============================================================================

fit_continuous <- function(
  outcome
) {

  vars <- c(
    outcome,
    "Baseline_25OHD_per10",
    "Age_years",
    "Female",
    "log_symptom_duration",
    "SolarRad_preBL_60d_mean_MJm2"
  )

  d <- dat %>%
    select(
      all_of(vars)
    ) %>%
    filter(
      complete.cases(.)
    )

  form <- reformulate(
    termlabels = c(
      "Baseline_25OHD_per10",
      "Age_years",
      "Female",
      "log_symptom_duration",
      "SolarRad_preBL_60d_mean_MJm2"
    ),
    response = outcome
  )

  fit <- lm(
    form,
    data = d
  )

  V <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )

  term <- "Baseline_25OHD_per10"

  b <- coef(fit)[term]
  se <- sqrt(
    diag(V)
  )[term]

  df_res <- df.residual(fit)

  crit <- qt(
    0.975,
    df = df_res
  )

  p <- 2 * pt(
    abs(b / se),
    df = df_res,
    lower.tail = FALSE
  )

  tibble(
    Variable = outcome,
    N_adjusted = nrow(d),
    Estimate = unname(b),
    CI_low = unname(
      b - crit * se
    ),
    CI_high = unname(
      b + crit * se
    ),
    p_adjusted = unname(p),
    Effect_measure = "Beta"
  )
}


continuous_vars <- feature_spec %>%
  filter(
    Type == "continuous"
  ) %>%
  pull(
    Variable
  )

adjusted_continuous <- bind_rows(
  lapply(
    continuous_vars,
    fit_continuous
  )
)


# =============================================================================
# 9. HC3 adjusted logistic models
# =============================================================================

fit_binary <- function(
  outcome
) {

  vars <- c(
    outcome,
    "Baseline_25OHD_per10",
    "Age_years",
    "Female",
    "log_symptom_duration",
    "SolarRad_preBL_60d_mean_MJm2"
  )

  d <- dat %>%
    select(
      all_of(vars)
    ) %>%
    filter(
      complete.cases(.)
    )

  form <- reformulate(
    termlabels = c(
      "Baseline_25OHD_per10",
      "Age_years",
      "Female",
      "log_symptom_duration",
      "SolarRad_preBL_60d_mean_MJm2"
    ),
    response = outcome
  )

  fit <- glm(
    form,
    data = d,
    family = binomial(
      link = "logit"
    )
  )

  V <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )

  term <- "Baseline_25OHD_per10"

  b <- coef(fit)[term]
  se <- sqrt(
    diag(V)
  )[term]

  z <- b / se

  p <- 2 * pnorm(
    abs(z),
    lower.tail = FALSE
  )

  ci_low_log <- b - 1.96 * se
  ci_high_log <- b + 1.96 * se

  tibble(
    Variable = outcome,
    N_adjusted = nrow(d),
    Estimate = exp(
      unname(b)
    ),
    CI_low = exp(
      unname(ci_low_log)
    ),
    CI_high = exp(
      unname(ci_high_log)
    ),
    p_adjusted = unname(p),
    Effect_measure = "OR"
  )
}


binary_vars <- feature_spec %>%
  filter(
    Type == "binary"
  ) %>%
  pull(
    Variable
  )

adjusted_binary <- bind_rows(
  lapply(
    binary_vars,
    fit_binary
  )
)


# =============================================================================
# 10. Combine adjusted results + BH q-values across all 12
# =============================================================================

adjusted_all <- bind_rows(
  adjusted_continuous,
  adjusted_binary
) %>%
  left_join(
    feature_spec,
    by = "Variable"
  ) %>%
  mutate(
    q_adjusted = p.adjust(
      p_adjusted,
      method = "BH"
    )
  )


# =============================================================================
# 11. Cross-check table
# =============================================================================

audit_table <- unadjusted %>%
  select(
    Variable,
    Label,
    Type,
    N_unadjusted = N,
    Spearman_rho,
    p_unadjusted,
    q_unadjusted
  ) %>%
  left_join(
    adjusted_all %>%
      select(
        Variable,
        N_adjusted,
        Estimate,
        CI_low,
        CI_high,
        p_adjusted,
        q_adjusted,
        Effect_measure
      ),
    by = "Variable"
  )


# =============================================================================
# 12. Plot-data preparation
# =============================================================================

label_order_all <- rev(
  feature_spec$Label
)

plot_unadjusted <- unadjusted %>%
  mutate(
    Label = factor(
      Label,
      levels = label_order_all
    ),

    Significant = factor(
      q_unadjusted < 0.05,
      levels = c(
        FALSE,
        TRUE
      ),
      labels = c(
        "q ≥ 0.05",
        "q < 0.05"
      )
    ),

    Stat_label = paste0(
      "\u03c1 = ",
      fmt_signed(
        Spearman_rho,
        3
      ),
      "; p ",
      ifelse(
        p_unadjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            p_unadjusted
          )
        )
      ),
      "; q ",
      ifelse(
        q_unadjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            q_unadjusted
          )
        )
      )
    )
  )


plot_cont <- adjusted_all %>%
  filter(
    Type == "continuous"
  ) %>%
  mutate(
    Label = factor(
      Label,
      levels = rev(
        feature_spec$Label[
          feature_spec$Type == "continuous"
        ]
      )
    ),

    Significant = factor(
      p_adjusted < 0.05,
      levels = c(
        FALSE,
        TRUE
      ),
      labels = c(
        "p ≥ 0.05",
        "p < 0.05"
      )
    ),

    Stat_label = paste0(
      "\u03b2 = ",
      fmt_signed(
        Estimate,
        2
      ),
      " (",
      fmt_num(
        CI_low,
        2
      ),
      " to ",
      fmt_num(
        CI_high,
        2
      ),
      "); p ",
      ifelse(
        p_adjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            p_adjusted
          )
        )
      ),
      "; q ",
      ifelse(
        q_adjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            q_adjusted
          )
        )
      )
    )
  )


plot_bin <- adjusted_all %>%
  filter(
    Type == "binary"
  ) %>%
  mutate(
    Label = factor(
      Label,
      levels = rev(
        feature_spec$Label[
          feature_spec$Type == "binary"
        ]
      )
    ),

    Significant = factor(
      p_adjusted < 0.05,
      levels = c(
        FALSE,
        TRUE
      ),
      labels = c(
        "p ≥ 0.05",
        "p < 0.05"
      )
    ),

    Stat_label = paste0(
      "OR = ",
      fmt_num(
        Estimate,
        2
      ),
      " (",
      fmt_num(
        CI_low,
        2
      ),
      " to ",
      fmt_num(
        CI_high,
        2
      ),
      "); p ",
      ifelse(
        p_adjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            p_adjusted
          )
        )
      ),
      "; q ",
      ifelse(
        q_adjusted < 0.001,
        "<0.001",
        paste0(
          "= ",
          fmt_p(
            q_adjusted
          )
        )
      )
    )
  )


# =============================================================================
# 13. Visual system
# =============================================================================

COL_TEXT <- "#252525"
COL_AXIS <- "#505050"
COL_NULL <- "#8A8A8A"

COL_GRAY <- "#6D6D6D"
COL_GRAY_LIGHT <- "#BEBEBE"

# Refined rose-pink for unadjusted q < 0.05
COL_VIOLET <- "#C23B6B"

# Refined teal for adjusted p < 0.05
COL_TEAL <- "#237A78"

theme_efig4 <- theme_classic(
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
      size = 9.2
    ),

    axis.text.y = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 9.2
    ),

    axis.title = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.3,
      face = "bold"
    ),

    legend.position = "none",

    plot.margin = margin(
      12,
      18,
      12,
      12
    )
  )


# =============================================================================
# 14. Panel A: unadjusted Spearman forest/dot plot
# =============================================================================

rho_min <- min(
  plot_unadjusted$Spearman_rho,
  0,
  na.rm = TRUE
)

rho_max <- max(
  plot_unadjusted$Spearman_rho,
  0,
  na.rm = TRUE
)

rho_span <- rho_max - rho_min

if (
  !is.finite(rho_span) ||
  rho_span <= 0
) {
  rho_span <- 1
}

rho_text_x <- rho_max +
  0.16 * rho_span

rho_upper <- rho_max +
  1.40 * rho_span

rho_lower <- rho_min -
  0.08 * rho_span

panel_a <- ggplot(
  plot_unadjusted,
  aes(
    x = Spearman_rho,
    y = Label
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.60,
    colour = COL_NULL
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = Spearman_rho,
      y = Label,
      yend = Label,
      colour = Significant
    ),
    linewidth = 0.80,
    alpha = 0.72
  ) +

  geom_point(
    aes(
      colour = Significant
    ),
    size = 3.2
  ) +

  geom_text(
    aes(
      x = rho_text_x,
      label = Stat_label,
      colour = Significant
    ),
    hjust = 0,
    family = "Arial",
    size = 2.70
  ) +

  scale_colour_manual(
    values = c(
      "q ≥ 0.05" = COL_GRAY,
      "q < 0.05" = COL_VIOLET
    )
  ) +

  scale_x_continuous(
    limits = c(
      rho_lower,
      rho_upper
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = "Unadjusted Spearman correlation coefficient (\u03c1)",
    y = NULL
  ) +

  theme_efig4


# =============================================================================
# 15. Panel B: adjusted continuous outcomes
# =============================================================================

cont_min <- min(
  plot_cont$CI_low,
  0,
  na.rm = TRUE
)

cont_max <- max(
  plot_cont$CI_high,
  0,
  na.rm = TRUE
)

cont_span <- cont_max - cont_min

if (
  !is.finite(cont_span) ||
  cont_span <= 0
) {
  cont_span <- 1
}

cont_text_x <- cont_max +
  0.16 * cont_span

cont_upper <- cont_max +
  1.35 * cont_span

cont_lower <- cont_min -
  0.08 * cont_span

panel_b <- ggplot(
  plot_cont,
  aes(
    x = Estimate,
    y = Label
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.60,
    colour = COL_NULL
  ) +

  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      colour = Significant
    ),
    height = 0.12,
    linewidth = 1.0
  ) +

  geom_point(
    aes(
      colour = Significant
    ),
    size = 3.5
  ) +

  geom_text(
    aes(
      x = cont_text_x,
      label = Stat_label,
      colour = Significant
    ),
    hjust = 0,
    family = "Arial",
    size = 2.78
  ) +

  scale_colour_manual(
    values = c(
      "p ≥ 0.05" = COL_GRAY,
      "p < 0.05" = COL_TEAL
    )
  ) +

  scale_x_continuous(
    limits = c(
      cont_lower,
      cont_upper
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = "Adjusted difference per +10-ng/mL higher baseline 25(OH)D",
    y = NULL
  ) +

  theme_efig4


# =============================================================================
# 16. Panel C: adjusted binary outcomes
# =============================================================================

# Use log scale because OR null = 1 and CIs are multiplicative.
bin_min <- min(
  plot_bin$CI_low,
  1,
  na.rm = TRUE
)

bin_max <- max(
  plot_bin$CI_high,
  1,
  na.rm = TRUE
)

# Right annotation region on multiplicative scale.
bin_text_x <- bin_max * 1.20
bin_upper <- bin_max * 3.10
bin_lower <- max(
  0.20,
  bin_min / 1.15
)

panel_c <- ggplot(
  plot_bin,
  aes(
    x = Estimate,
    y = Label
  )
) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.60,
    colour = COL_NULL
  ) +

  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      colour = Significant
    ),
    height = 0.12,
    linewidth = 1.0
  ) +

  geom_point(
    aes(
      colour = Significant
    ),
    size = 3.5
  ) +

  geom_text(
    aes(
      x = bin_text_x,
      label = Stat_label,
      colour = Significant
    ),
    hjust = 0,
    family = "Arial",
    size = 2.70
  ) +

  scale_colour_manual(
    values = c(
      "p ≥ 0.05" = COL_GRAY,
      "p < 0.05" = COL_TEAL
    )
  ) +

  scale_x_log10(
    limits = c(
      bin_lower,
      bin_upper
    ),
    breaks = c(
      0.5,
      0.75,
      1,
      1.5,
      2,
      3
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = "Adjusted odds ratio per +10-ng/mL higher baseline 25(OH)D",
    y = NULL
  ) +

  theme_efig4


# =============================================================================
# 17. Combine
# =============================================================================

efigure4 <- panel_a /
  panel_b /
  panel_c +
  plot_layout(
    heights = c(
      1.65,
      0.75,
      1.30
    )
  )

print(
  efigure4
)


# =============================================================================
# 18. Console audit
# =============================================================================

cat("\n====================================================\n")
cat("Supplementary Figure S4 baseline 25(OH)D clinical-feature analysis\n")
cat("====================================================\n\n")

cat("Unadjusted associations:\n")
print(
  unadjusted %>%
    select(
      Label,
      N,
      Spearman_rho,
      p_unadjusted,
      q_unadjusted
    )
)

cat("\nAdjusted associations:\n")
print(
  adjusted_all %>%
    select(
      Label,
      Type,
      N_adjusted,
      Effect_measure,
      Estimate,
      CI_low,
      CI_high,
      p_adjusted,
      q_adjusted
    )
)


# =============================================================================
# 19. Analysis information
# =============================================================================

analysis_info <- tibble(
  Item = c(
    "Panel A",
    "Panel A multiplicity",
    "Panel B",
    "Panel C",
    "Adjusted exposure scale",
    "Adjusted covariates",
    "Adjusted continuous models",
    "Adjusted binary models",
    "Adjusted multiplicity",
    "Unadjusted color rule",
    "Adjusted color rule",
    "Panel letters",
    "Missing data",
    "R version requested for manuscript"
  ),

  Value = c(
    paste(
      "Spearman correlations of baseline serum 25(OH)D with",
      "12 prespecified TMD-related clinical features"
    ),

    "Benjamini-Hochberg q-values across all 12 unadjusted feature comparisons",

    paste(
      "Adjusted HC3 linear-regression beta coefficients for pain VAS, PFO, and MUO"
    ),

    paste(
      "Adjusted HC3 logistic-regression odds ratios for nine binary clinical features"
    ),

    "Per +10-ng/mL higher baseline serum 25(OH)D",

    paste(
      "Age, sex, log1p symptom duration, and baseline 60-day mean solar radiation"
    ),

    "Multivariable linear regression with HC3 robust standard errors",

    "Multivariable logistic regression with HC3 robust standard errors",

    "Benjamini-Hochberg q-values across all 12 adjusted models",

    paste(
      "Refined rose-pink for q<0.05; gray otherwise"
    ),

    paste(
      "Muted teal for adjusted p<0.05; gray otherwise"
    ),

    "Not drawn inside figure; add (A), (B), and (C) during manuscript layout",

    "Outcome-specific complete-case analysis; no imputation",

    "R 4.5.1"
  )
)


# =============================================================================
# 20. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Supplementary_Figure_S4_clinical_features.png"
)

tiff_file <- file.path(
  output_dir,
  "Supplementary_Figure_S4_clinical_features.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Figure_S4_clinical_features.pdf"
)

ggsave(
  filename = png_file,
  plot = efigure4,
  width = 12.0,
  height = 13.0,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = efigure4,
  width = 12.0,
  height = 13.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = efigure4,
  width = 12.0,
  height = 13.0,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 21. Export audit workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Supplementary_Figure_S4_clinical_features_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Combined_Audit"
)

addWorksheet(
  wb,
  "Unadjusted"
)

addWorksheet(
  wb,
  "Adjusted_Continuous"
)

addWorksheet(
  wb,
  "Adjusted_Binary"
)

addWorksheet(
  wb,
  "Adjusted_All"
)

addWorksheet(
  wb,
  "Analysis_Info"
)

writeData(
  wb,
  "Combined_Audit",
  audit_table
)

writeData(
  wb,
  "Unadjusted",
  unadjusted
)

writeData(
  wb,
  "Adjusted_Continuous",
  adjusted_continuous
)

writeData(
  wb,
  "Adjusted_Binary",
  adjusted_binary
)

writeData(
  wb,
  "Adjusted_All",
  adjusted_all
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
    38,
    110
  )
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 22. Console output
# =============================================================================

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("====================================================\n")
