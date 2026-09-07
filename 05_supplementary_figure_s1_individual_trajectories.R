# =============================================================================
# Supplementary Figure S1. Individual longitudinal changes in serum 25(OH)D
# and pain intensity among vitamin D-prescribed patients
#
# Panels
#   A. Individual baseline-to-follow-up serum 25(OH)D trajectories
#   B. Individual baseline-to-follow-up pain-intensity trajectories
#
# Figure rules
#   - No overall title inside figure
#   - No panel letters inside figure
#   - Individual patients shown as thin light-gray paired trajectories
#   - Group mean trajectory overlaid prominently
#   - Panel A mean trajectory: red
#   - Panel B mean trajectory: dark blue
#   - 95% CI shown for baseline and follow-up means
#   - Compact annotation reports:
#       n
#       baseline mean
#       follow-up mean
#       mean paired change with 95% CI
#       paired two-sided p-value
#
# Interpretation
#   These are within-patient descriptive longitudinal changes among patients
#   with a recorded vitamin D prescription. They do not estimate a causal
#   prescription effect because there is no untreated counterfactual within
#   this panel
#
# Change definitions
#   Delta 25(OH)D = follow-up - baseline
#   Pain reduction = baseline - follow-up VAS
#
# Outputs
#   Supplementary_Figure_S1_individual_trajectories.png
#   Supplementary_Figure_S1_individual_trajectories.tiff
#   Supplementary_Figure_S1_individual_trajectories.pdf
#   Supplementary_Figure_S1_individual_trajectories_statistics.xlsx
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
output_dir <- file.path("outputs", "Supplementary_Figure_S1")
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
  "Vitamin_D_prescription",
  "Baseline_25OHD_ng_mL",
  "Followup_25OHD_ng_mL",
  "Delta_25OHD_ng_mL",
  "Baseline_pain_intensity_VAS",
  "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS"
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
# 4. Prescribed cohort
# =============================================================================

rx01 <- to_binary01(
  dat0$Vitamin_D_prescription,
  "Vitamin_D_prescription"
)

dat_rx <- dat0 %>%
  mutate(
    Rx = rx01
  ) %>%
  filter(
    Rx == 1
  )

cat("\n=============================================\n")
cat("Supplementary Figure S1: prescribed-patient trajectories\n")
cat("=============================================\n")
cat("Recorded vitamin D prescription N = ", nrow(dat_rx), "\n", sep = "")

if (nrow(dat_rx) != 149) {
  warning(
    "Expected prescribed-patient N=149; found N=",
    nrow(dat_rx)
  )
}


# =============================================================================
# 5. Paired complete-case datasets
# =============================================================================

dat_d <- dat_rx %>%
  select(
    Study_ID,
    Baseline_25OHD_ng_mL,
    Followup_25OHD_ng_mL
  ) %>%
  filter(
    complete.cases(.)
  ) %>%
  mutate(
    Delta_25OHD_calc =
      Followup_25OHD_ng_mL -
      Baseline_25OHD_ng_mL
  )

dat_pain <- dat_rx %>%
  select(
    Study_ID,
    Baseline_pain_intensity_VAS,
    Followup_pain_intensity_VAS
  ) %>%
  filter(
    complete.cases(.)
  ) %>%
  mutate(
    Pain_reduction_calc =
      Baseline_pain_intensity_VAS -
      Followup_pain_intensity_VAS
  )


# =============================================================================
# 6. Paired-change summary helper
# =============================================================================

paired_summary <- function(
  baseline,
  followup,
  change,
  change_direction = "FU_minus_BL"
) {

  n <- length(change)

  baseline_mean <- mean(baseline)
  baseline_sd <- sd(baseline)

  followup_mean <- mean(followup)
  followup_sd <- sd(followup)

  change_mean <- mean(change)
  change_sd <- sd(change)
  change_se <- change_sd / sqrt(n)

  crit <- qt(
    0.975,
    df = n - 1
  )

  change_ci_low <- change_mean - crit * change_se
  change_ci_high <- change_mean + crit * change_se

  # Equivalent to paired t-test for the supplied paired-change definition.
  tt <- t.test(
    change,
    mu = 0,
    alternative = "two.sided"
  )

  tibble(
    N = n,
    Baseline_mean = baseline_mean,
    Baseline_SD = baseline_sd,
    Followup_mean = followup_mean,
    Followup_SD = followup_sd,
    Mean_change = change_mean,
    Change_SD = change_sd,
    Change_CI_low = change_ci_low,
    Change_CI_high = change_ci_high,
    Paired_t = unname(tt$statistic),
    df = unname(tt$parameter),
    p_value = tt$p.value,
    Change_definition = change_direction
  )
}


sum_d <- paired_summary(
  baseline = dat_d$Baseline_25OHD_ng_mL,
  followup = dat_d$Followup_25OHD_ng_mL,
  change = dat_d$Delta_25OHD_calc,
  change_direction = "Follow-up minus baseline"
)

sum_pain <- paired_summary(
  baseline = dat_pain$Baseline_pain_intensity_VAS,
  followup = dat_pain$Followup_pain_intensity_VAS,
  change = dat_pain$Pain_reduction_calc,
  change_direction = "Baseline minus follow-up; positive = pain improvement"
)


# =============================================================================
# 7. Mean and 95% CI at each timepoint
# =============================================================================

mean_ci <- function(x) {

  n <- sum(!is.na(x))
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  se <- s / sqrt(n)
  crit <- qt(0.975, df = n - 1)

  tibble(
    N = n,
    Mean = m,
    SD = s,
    CI_low = m - crit * se,
    CI_high = m + crit * se
  )
}


summary_d_time <- bind_rows(
  mean_ci(
    dat_d$Baseline_25OHD_ng_mL
  ) %>%
    mutate(
      Time = "Baseline"
    ),

  mean_ci(
    dat_d$Followup_25OHD_ng_mL
  ) %>%
    mutate(
      Time = "Follow-up"
    )
) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Baseline",
        "Follow-up"
      )
    )
  )


summary_pain_time <- bind_rows(
  mean_ci(
    dat_pain$Baseline_pain_intensity_VAS
  ) %>%
    mutate(
      Time = "Baseline"
    ),

  mean_ci(
    dat_pain$Followup_pain_intensity_VAS
  ) %>%
    mutate(
      Time = "Follow-up"
    )
) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Baseline",
        "Follow-up"
      )
    )
  )


# =============================================================================
# 8. Long-format individual trajectories
# =============================================================================

long_d <- dat_d %>%
  select(
    Study_ID,
    Baseline_25OHD_ng_mL,
    Followup_25OHD_ng_mL
  ) %>%
  pivot_longer(
    cols = c(
      Baseline_25OHD_ng_mL,
      Followup_25OHD_ng_mL
    ),
    names_to = "Time",
    values_to = "Value"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Baseline_25OHD_ng_mL",
        "Followup_25OHD_ng_mL"
      ),
      labels = c(
        "Baseline",
        "Follow-up"
      )
    )
  )


long_pain <- dat_pain %>%
  select(
    Study_ID,
    Baseline_pain_intensity_VAS,
    Followup_pain_intensity_VAS
  ) %>%
  pivot_longer(
    cols = c(
      Baseline_pain_intensity_VAS,
      Followup_pain_intensity_VAS
    ),
    names_to = "Time",
    values_to = "Value"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Baseline_pain_intensity_VAS",
        "Followup_pain_intensity_VAS"
      ),
      labels = c(
        "Baseline",
        "Follow-up"
      )
    )
  )


# =============================================================================
# 9. Formatting helpers
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_num <- function(
  x,
  digits = 2
) {
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


# =============================================================================
# 10. Annotation text
# =============================================================================

ann_d <- paste0(
  "n = ",
  sum_d$N,
  "\nMean: ",
  fmt_num(sum_d$Baseline_mean, 2),
  " \u2192 ",
  fmt_num(sum_d$Followup_mean, 2),
  " ng/mL",
  "\nMean \u0394 = ",
  ifelse(
    sum_d$Mean_change >= 0,
    "+",
    ""
  ),
  fmt_num(sum_d$Mean_change, 2),
  " (95% CI ",
  fmt_num(sum_d$Change_CI_low, 2),
  " to ",
  fmt_num(sum_d$Change_CI_high, 2),
  ")",
  "\nPaired p ",
  ifelse(
    sum_d$p_value < 0.001,
    "< 0.001",
    paste0(
      "= ",
      fmt_p(sum_d$p_value)
    )
  )
)


ann_pain <- paste0(
  "n = ",
  sum_pain$N,
  "\nMean VAS: ",
  fmt_num(sum_pain$Baseline_mean, 2),
  " \u2192 ",
  fmt_num(sum_pain$Followup_mean, 2),
  "\nMean pain reduction = ",
  ifelse(
    sum_pain$Mean_change >= 0,
    "+",
    ""
  ),
  fmt_num(sum_pain$Mean_change, 2),
  " (95% CI ",
  fmt_num(sum_pain$Change_CI_low, 2),
  " to ",
  fmt_num(sum_pain$Change_CI_high, 2),
  ")",
  "\nPaired p ",
  ifelse(
    sum_pain$p_value < 0.001,
    "< 0.001",
    paste0(
      "= ",
      fmt_p(sum_pain$p_value)
    )
  )
)


# =============================================================================
# 11. Visual system
# =============================================================================

COL_TEXT <- "#252525"
COL_AXIS <- "#505050"

COL_INDIVIDUAL <- "#B8B8B8"
COL_INDIVIDUAL_POINT <- "#9E9E9E"

COL_RED <- "#C00000"
COL_RED_LIGHT <- "#E6A0A0"

COL_BLUE <- "#003A8C"
COL_BLUE_LIGHT <- "#9CB7DC"

COL_REF <- "#888888"

theme_efig2 <- theme_classic(
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
      size = 9.5
    ),

    axis.title = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.5,
      face = "bold"
    ),

    axis.title.y = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.5,
      face = "bold"
    ),

    plot.margin = margin(
      12, 14, 12, 14
    )
  )


# =============================================================================
# 12. Panel A: individual serum 25(OH)D trajectories
# =============================================================================

panel_a <- ggplot(
  long_d,
  aes(
    x = Time,
    y = Value,
    group = Study_ID
  )
) +

  geom_hline(
    yintercept = 20,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_REF
  ) +

  geom_hline(
    yintercept = 30,
    linetype = "dotted",
    linewidth = 0.55,
    colour = COL_REF
  ) +

  geom_line(
    colour = COL_INDIVIDUAL,
    linewidth = 0.48,
    alpha = 0.34
  ) +

  geom_point(
    colour = COL_INDIVIDUAL_POINT,
    size = 1.45,
    alpha = 0.38
  ) +

  geom_errorbar(
    data = summary_d_time,
    aes(
      x = Time,
      ymin = CI_low,
      ymax = CI_high,
      group = 1
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.95,
    colour = COL_RED
  ) +

  geom_line(
    data = summary_d_time,
    aes(
      x = Time,
      y = Mean,
      group = 1
    ),
    inherit.aes = FALSE,
    linewidth = 1.65,
    colour = COL_RED
  ) +

  geom_point(
    data = summary_d_time,
    aes(
      x = Time,
      y = Mean
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 4.0,
    stroke = 1.0,
    fill = "white",
    colour = COL_RED
  ) +

  annotate(
    "text",
    x = 1.50,
    y = Inf,
    label = ann_d,
    hjust = 0.5,
    vjust = 1.15,
    size = 2.85,
    family = "Arial",
    lineheight = 1.08,
    colour = COL_RED
  ) +

  annotate(
    "text",
    x = 2.02,
    y = 20.5,
    label = "20 ng/mL",
    hjust = 0,
    vjust = 0,
    size = 2.8,
    family = "Arial",
    colour = "#666666"
  ) +

  annotate(
    "text",
    x = 2.02,
    y = 30.5,
    label = "30 ng/mL",
    hjust = 0,
    vjust = 0,
    size = 2.8,
    family = "Arial",
    colour = "#666666"
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = "Serum 25(OH)D, ng/mL"
  ) +

  theme_efig2


# =============================================================================
# 13. Panel B: individual pain-intensity trajectories
# =============================================================================

panel_b <- ggplot(
  long_pain,
  aes(
    x = Time,
    y = Value,
    group = Study_ID
  )
) +

  geom_line(
    colour = COL_INDIVIDUAL,
    linewidth = 0.48,
    alpha = 0.34
  ) +

  geom_point(
    colour = COL_INDIVIDUAL_POINT,
    size = 1.45,
    alpha = 0.38
  ) +

  geom_errorbar(
    data = summary_pain_time,
    aes(
      x = Time,
      ymin = CI_low,
      ymax = CI_high,
      group = 1
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.95,
    colour = COL_BLUE
  ) +

  geom_line(
    data = summary_pain_time,
    aes(
      x = Time,
      y = Mean,
      group = 1
    ),
    inherit.aes = FALSE,
    linewidth = 1.65,
    colour = COL_BLUE
  ) +

  geom_point(
    data = summary_pain_time,
    aes(
      x = Time,
      y = Mean
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 4.0,
    stroke = 1.0,
    fill = "white",
    colour = COL_BLUE
  ) +

  annotate(
    "text",
    x = 1.50,
    y = Inf,
    label = ann_pain,
    hjust = 0.5,
    vjust = 1.15,
    size = 2.85,
    family = "Arial",
    lineheight = 1.08,
    colour = COL_BLUE
  ) +

  coord_cartesian(
    ylim = c(
      0,
      max(
        10,
        max(
          long_pain$Value,
          na.rm = TRUE
        )
      )
    ),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = "Pain intensity, VAS"
  ) +

  theme_efig2


# =============================================================================
# 14. Combine
# =============================================================================

efigure2 <- panel_a + panel_b +
  plot_layout(
    widths = c(
      1,
      1
    )
  )

print(efigure2)


# =============================================================================
# 15. Audit tables
# =============================================================================

paired_stats <- bind_rows(
  sum_d %>%
    mutate(
      Outcome = "Serum 25(OH)D"
    ),

  sum_pain %>%
    mutate(
      Outcome = "Pain reduction"
    )
) %>%
  select(
    Outcome,
    everything()
  )


analysis_info <- tibble(
  Item = c(
    "Analysis cohort",
    "Panel A",
    "Panel B",
    "Panel A change definition",
    "Panel B change definition",
    "Individual trajectories",
    "Mean trajectory",
    "Mean uncertainty",
    "Panel A reference lines",
    "Paired p-value",
    "Interpretation",
    "Panel letters",
    "R version requested for manuscript"
  ),

  Value = c(
    "Patients with a recorded vitamin D prescription",
    "Individual baseline-to-follow-up serum 25(OH)D trajectories",
    "Individual baseline-to-follow-up pain-intensity trajectories",
    "Follow-up minus baseline serum 25(OH)D",
    "Baseline minus follow-up VAS; positive value indicates pain improvement",
    "Thin light-gray paired lines",
    "Red for serum 25(OH)D; dark blue for pain intensity",
    "95% CI around observed mean at each timepoint",
    "20 and 30 ng/mL",
    "Two-sided one-sample t-test of the paired change against zero",
    paste(
      "Descriptive within-patient longitudinal changes among prescribed patients;",
      "not a causal estimate of vitamin D prescription effect"
    ),
    "Not drawn inside figure; add (A) and (B) during manuscript layout",
    "R 4.5.1"
  )
)


# =============================================================================
# 16. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Supplementary_Figure_S1_individual_trajectories.png"
)

tiff_file <- file.path(
  output_dir,
  "Supplementary_Figure_S1_individual_trajectories.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Figure_S1_individual_trajectories.pdf"
)

ggsave(
  filename = png_file,
  plot = efigure2,
  width = 11.0,
  height = 5.8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = efigure2,
  width = 11.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = efigure2,
  width = 11.0,
  height = 5.8,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 17. Export statistics workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Supplementary_Figure_S1_individual_trajectories_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Paired_Stats"
)

addWorksheet(
  wb,
  "VitD_Timepoint_Summary"
)

addWorksheet(
  wb,
  "Pain_Timepoint_Summary"
)

addWorksheet(
  wb,
  "VitD_Paired_Data"
)

addWorksheet(
  wb,
  "Pain_Paired_Data"
)

addWorksheet(
  wb,
  "Analysis_Info"
)

writeData(
  wb,
  "Paired_Stats",
  paired_stats
)

writeData(
  wb,
  "VitD_Timepoint_Summary",
  summary_d_time
)

writeData(
  wb,
  "Pain_Timepoint_Summary",
  summary_pain_time
)

writeData(
  wb,
  "VitD_Paired_Data",
  dat_d
)

writeData(
  wb,
  "Pain_Paired_Data",
  dat_pain
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
    105
  )
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 18. Console output
# =============================================================================

cat("\nSerum 25(OH)D paired summary:\n")
print(sum_d)

cat("\nPain paired summary:\n")
print(sum_pain)

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("=============================================\n")
