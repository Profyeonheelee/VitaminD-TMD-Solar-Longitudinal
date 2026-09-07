# =============================================================================
# Supplementary Figure S2. Changes in vitamin D status according to baseline
# vitamin D status and prescription group
# R version: 4.5.1
#
# Definition:
#   Low vitamin D status = serum 25(OH)D <30 ng/mL
#   Higher/normalized status = serum 25(OH)D >=30 ng/mL
#
# Analysis notes:
# - The original three exposure groups are defined partly by BASELINE
#   25(OH)D status:
#     Group 1 = baseline >=30 ng/mL, no prescription
#     Group 2 = baseline <30 ng/mL, no prescription
#     Group 3 = baseline <30 ng/mL, prescription
# - Therefore, NO inferential p-value is calculated for baseline differences
#   in <30-ng/mL prevalence across the three groups. That baseline contrast
#   is deterministic by study-group definition.
# - The inferential comparison is made at FOLLOW-UP.
# - Among Groups 2 and 3, normalization to >=30 ng/mL is also compared
#   descriptively/inferentially at follow-up.
#
# Statistical approach:
# - Exact Clopper-Pearson 95% CIs for proportions in the figure.
# - Follow-up three-group comparison: Pearson chi-square unless any expected
#   cell count is <5, in which case Fisher exact test is used.
# - All three pairwise follow-up comparisons are tested similarly and
#   Holm-adjusted, matching the Table 2 binary-outcome approach.
#
# Outputs:
#   Supplementary_Figure_S2_vitamin_D_status_statistics.xlsx
#   Supplementary_Figure_S2_vitamin_D_status.png
#   Supplementary_Figure_S2_vitamin_D_status.tiff
#   Supplementary_Figure_S2_vitamin_D_status.pdf
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
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
library(openxlsx)


# =============================================================================
# 1. File path
# =============================================================================

xlsx_path <- file.path(
  "data",
  "Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx"
)

sheet_name <- "Analysis_Data"
output_dir <- file.path("outputs", "Supplementary_Figure_S2")
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
  "Three_group_code",
  "Baseline_25OHD_ng_mL",
  "Followup_25OHD_ng_mL"
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
  filter(
    Three_group_code %in% c(1, 2, 3)
  ) %>%
  mutate(
    Group = factor(
      Three_group_code,
      levels = c(1, 2, 3),
      labels = c(
        "Group 1",
        "Group 2",
        "Group 3"
      )
    ),

    BL_low30 = if_else(
      is.na(Baseline_25OHD_ng_mL),
      NA_integer_,
      as.integer(
        Baseline_25OHD_ng_mL < 30
      )
    ),

    FU_low30 = if_else(
      is.na(Followup_25OHD_ng_mL),
      NA_integer_,
      as.integer(
        Followup_25OHD_ng_mL < 30
      )
    ),

    BL_ge30 = if_else(
      is.na(Baseline_25OHD_ng_mL),
      NA_integer_,
      as.integer(
        Baseline_25OHD_ng_mL >= 30
      )
    ),

    FU_ge30 = if_else(
      is.na(Followup_25OHD_ng_mL),
      NA_integer_,
      as.integer(
        Followup_25OHD_ng_mL >= 30
      )
    )
  )


# =============================================================================
# 3. Sanity checks
# =============================================================================

cat("\n============================================================\n")
cat("30-ng/mL vitamin D status transition analysis\n")
cat("============================================================\n\n")

cat("Three-group sample sizes:\n")
print(
  table(dat$Group)
)

cat("\nBaseline low-status counts (<30 ng/mL):\n")
print(
  table(
    dat$Group,
    dat$BL_low30,
    useNA = "ifany"
  )
)

cat("\nFollow-up low-status counts (<30 ng/mL):\n")
print(
  table(
    dat$Group,
    dat$FU_low30,
    useNA = "ifany"
  )
)

# Because group assignment uses baseline 25(OH)D status, these should hold.
expected_bl <- dat %>%
  group_by(Group) %>%
  summarise(
    N = sum(
      !is.na(BL_low30)
    ),
    Low_n = sum(
      BL_low30 == 1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (
  expected_bl$Low_n[
    expected_bl$Group == "Group 1"
  ] != 0
) {
  warning(
    "Group 1 contains baseline 25(OH)D <30 ng/mL; check group coding."
  )
}

if (
  expected_bl$Low_n[
    expected_bl$Group == "Group 2"
  ] !=
  expected_bl$N[
    expected_bl$Group == "Group 2"
  ]
) {
  warning(
    "Not all Group 2 participants are baseline <30 ng/mL; check group coding."
  )
}

if (
  expected_bl$Low_n[
    expected_bl$Group == "Group 3"
  ] !=
  expected_bl$N[
    expected_bl$Group == "Group 3"
  ]
) {
  warning(
    "Not all Group 3 participants are baseline <30 ng/mL; check group coding."
  )
}


# =============================================================================
# 4. Formatting helpers
# =============================================================================

fmt_p <- function(p) {

  if (length(p) == 0 || is.na(p)) {
    return("NA")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf(
    "%.3f",
    p
  )
}


fmt_pct <- function(x) {

  sprintf(
    "%.1f%%",
    x
  )
}


# =============================================================================
# 5. Exact binomial 95% CI
# =============================================================================

exact_prop_ci <- function(
  events,
  n
) {

  if (
    is.na(events) ||
    is.na(n) ||
    n == 0
  ) {
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }

  bt <- binom.test(
    events,
    n,
    conf.level = 0.95
  )

  as.numeric(
    bt$conf.int
  ) * 100
}


# =============================================================================
# 6. Long-format status data
# =============================================================================

long_status <- bind_rows(
  dat %>%
    transmute(
      Study_ID,
      Group,
      Visit = "Baseline",
      Low30 = BL_low30,
      GE30 = BL_ge30
    ),

  dat %>%
    transmute(
      Study_ID,
      Group,
      Visit = "Follow-up",
      Low30 = FU_low30,
      GE30 = FU_ge30
    )
) %>%
  mutate(
    Visit = factor(
      Visit,
      levels = c(
        "Baseline",
        "Follow-up"
      )
    )
  )


# =============================================================================
# 7. Group-by-visit status summary
# =============================================================================

status_summary <- long_status %>%
  filter(
    !is.na(Low30)
  ) %>%
  group_by(
    Group,
    Visit
  ) %>%
  summarise(
    N = n(),
    Low_n = sum(
      Low30 == 1
    ),
    GE30_n = sum(
      Low30 == 0
    ),
    Low_pct = 100 * Low_n / N,
    GE30_pct = 100 * GE30_n / N,
    .groups = "drop"
  )

ci_mat <- t(
  mapply(
    exact_prop_ci,
    status_summary$Low_n,
    status_summary$N
  )
)

status_summary$Low_CI_low <- ci_mat[, 1]
status_summary$Low_CI_high <- ci_mat[, 2]

status_summary <- status_summary %>%
  mutate(
    Low_display = paste0(
      Low_n,
      "/",
      N,
      " (",
      sprintf(
        "%.1f",
        Low_pct
      ),
      "%)"
    ),

    GE30_display = paste0(
      GE30_n,
      "/",
      N,
      " (",
      sprintf(
        "%.1f",
        GE30_pct
      ),
      "%)"
    )
  )


# =============================================================================
# 8. Follow-up three-group test + pairwise Holm tests
#    Same logic as Table 2 binary descriptive analyses
# =============================================================================

d_fu <- dat %>%
  select(
    Group,
    FU_low30
  ) %>%
  filter(
    !is.na(Group),
    !is.na(FU_low30)
  )

tab_fu <- table(
  d_fu$Group,
  d_fu$FU_low30
)

chi_global <- suppressWarnings(
  chisq.test(
    tab_fu,
    correct = FALSE
  )
)

if (
  any(
    chi_global$expected < 5
  )
) {

  global_test <- fisher.test(
    tab_fu
  )

  global_method <- "Fisher exact"

} else {

  global_test <- chi_global
  global_method <- "Pearson chi-square"
}

global_p <- global_test$p.value


pairs <- combn(
  levels(d_fu$Group),
  2,
  simplify = FALSE
)

pairwise_fu <- bind_rows(
  lapply(
    pairs,
    function(pr) {

      d2 <- d_fu %>%
        filter(
          Group %in% pr
        ) %>%
        mutate(
          Group = droplevels(
            Group
          )
        )

      tab2 <- table(
        d2$Group,
        d2$FU_low30
      )

      chi2 <- suppressWarnings(
        chisq.test(
          tab2,
          correct = FALSE
        )
      )

      if (
        any(
          chi2$expected < 5
        )
      ) {

        tst <- fisher.test(
          tab2
        )

        method_used <- "Fisher exact"

      } else {

        tst <- chi2
        method_used <- "Pearson chi-square"
      }

      tibble(
        Comparison = paste(
          pr[1],
          "vs",
          pr[2]
        ),
        Method = method_used,
        p_raw = tst$p.value
      )
    }
  )
) %>%
  mutate(
    p_Holm = p.adjust(
      p_raw,
      method = "holm"
    )
  )


# =============================================================================
# 9. Low-baseline cohort normalization comparison: Group 2 vs Group 3
# =============================================================================

low_bl <- dat %>%
  filter(
    Group %in% c(
      "Group 2",
      "Group 3"
    ),
    BL_low30 == 1,
    !is.na(FU_ge30)
  ) %>%
  mutate(
    Group = droplevels(
      Group
    )
  )

normalization_summary <- low_bl %>%
  group_by(
    Group
  ) %>%
  summarise(
    N = n(),
    Normalized_ge30_n = sum(
      FU_ge30 == 1
    ),
    Remained_low_n = sum(
      FU_ge30 == 0
    ),
    Normalized_ge30_pct =
      100 * Normalized_ge30_n / N,
    Remained_low_pct =
      100 * Remained_low_n / N,
    .groups = "drop"
  ) %>%
  mutate(
    Normalized_display = paste0(
      Normalized_ge30_n,
      "/",
      N,
      " (",
      sprintf(
        "%.1f",
        Normalized_ge30_pct
      ),
      "%)"
    )
  )

tab_norm <- table(
  low_bl$Group,
  low_bl$FU_ge30
)

chi_norm <- suppressWarnings(
  chisq.test(
    tab_norm,
    correct = FALSE
  )
)

if (
  any(
    chi_norm$expected < 5
  )
) {

  norm_test <- fisher.test(
    tab_norm
  )

  norm_method <- "Fisher exact"

} else {

  norm_test <- chi_norm
  norm_method <- "Pearson chi-square"
}

normalization_p <- norm_test$p.value


# =============================================================================
# 10. Individual transition categories
# =============================================================================

transition_detail <- dat %>%
  filter(
    !is.na(BL_low30),
    !is.na(FU_low30)
  ) %>%
  mutate(
    Transition = case_when(
      BL_low30 == 0 &
        FU_low30 == 0 ~
        "Remained ≥30 ng/mL",

      BL_low30 == 0 &
        FU_low30 == 1 ~
        "≥30 to <30 ng/mL",

      BL_low30 == 1 &
        FU_low30 == 0 ~
        "<30 to ≥30 ng/mL",

      BL_low30 == 1 &
        FU_low30 == 1 ~
        "Remained <30 ng/mL",

      TRUE ~ NA_character_
    )
  ) %>%
  count(
    Group,
    Transition,
    name = "n"
  ) %>%
  group_by(
    Group
  ) %>%
  mutate(
    N_group = sum(n),
    Percent = 100 * n / N_group
  ) %>%
  ungroup()


# =============================================================================
# 11. Summary table
# =============================================================================

get_display <- function(
  group_name,
  visit_name,
  field = "Low_display"
) {

  tmp <- status_summary %>%
    filter(
      Group == group_name,
      Visit == visit_name
    )

  if (
    nrow(tmp) == 0
  ) {
    return("")
  }

  tmp[[field]][1]
}


g2_norm <- normalization_summary %>%
  filter(
    Group == "Group 2"
  ) %>%
  pull(
    Normalized_display
  )

g3_norm <- normalization_summary %>%
  filter(
    Group == "Group 3"
  ) %>%
  pull(
    Normalized_display
  )


etable5 <- tibble(
  Outcome = c(
    "Low vitamin D status (<30 ng/mL), baseline",
    "Low vitamin D status (<30 ng/mL), follow-up",
    "Normalization to ≥30 ng/mL among patients with low baseline vitamin D"
  ),

  `Group 1` = c(
    get_display(
      "Group 1",
      "Baseline",
      "Low_display"
    ),

    get_display(
      "Group 1",
      "Follow-up",
      "Low_display"
    ),

    "Not applicable"
  ),

  `Group 2` = c(
    get_display(
      "Group 2",
      "Baseline",
      "Low_display"
    ),

    get_display(
      "Group 2",
      "Follow-up",
      "Low_display"
    ),

    g2_norm
  ),

  `Group 3` = c(
    get_display(
      "Group 3",
      "Baseline",
      "Low_display"
    ),

    get_display(
      "Group 3",
      "Follow-up",
      "Low_display"
    ),

    g3_norm
  ),

  `Comparison p` = c(
    "Not tested: baseline status defines the exposure groups",

    paste0(
      fmt_p(
        global_p
      ),
      " (global)"
    ),

    paste0(
      fmt_p(
        normalization_p
      ),
      " (Group 2 vs Group 3)"
    )
  )
)


# =============================================================================
# 12. Supplementary Figure S2
# =============================================================================

group_colors <- c(
  "Group 1" = "#555555",
  "Group 2" = "#496D89",
  "Group 3" = "#C66A2B"
)

group_strip_labels <- c(
  "Group 1" =
    "Group 1\nSufficient / no prescription",

  "Group 2" =
    "Group 2\nLow / no prescription",

  "Group 3" =
    "Group 3\nLow / prescription"
)


plot_dat <- status_summary %>%
  mutate(
    Visit_num = ifelse(
      Visit == "Baseline",
      1,
      2
    ),

    Point_label = paste0(
      Low_n,
      "/",
      N,
      "\n",
      sprintf(
        "%.1f%%",
        Low_pct
      )
    ),

    Label_y = case_when(
      Low_pct >= 94 ~
        Low_pct - 11,

      Low_pct <= 6 ~
        Low_pct + 10,

      TRUE ~
        Low_pct + 9
    )
  )

change_summary <- status_summary %>%
  select(
    Group,
    Visit,
    Low_pct
  ) %>%
  tidyr::pivot_wider(
    names_from = Visit,
    values_from = Low_pct
  ) %>%
  mutate(
    Delta_pp = `Follow-up` - Baseline,
    Delta_label = paste0(
      "Change: ",
      ifelse(
        Delta_pp > 0,
        "+",
        ""
      ),
      sprintf(
        "%.1f",
        Delta_pp
      ),
      " percentage points"
    ),
    x = 1.50,
    y = 6
  )


normalization_annot <- normalization_summary %>%
  filter(
    Group %in% c(
      "Group 2",
      "Group 3"
    )
  ) %>%
  mutate(
    x = 1.50,
    y = 18,
    Label = paste0(
      "Normalized to ≥30 ng/mL: ",
      Normalized_ge30_n,
      "/",
      N,
      " (",
      sprintf(
        "%.1f",
        Normalized_ge30_pct
      ),
      "%)"
    )
  )

normalization_p_annot <- tibble(
  Group = factor(
    "Group 3",
    levels = levels(dat$Group)
  ),
  x = 1.50,
  y = 11.5,
  Label = paste0(
    "G2 vs G3 normalization p = ",
    fmt_p(
      normalization_p
    )
  )
)


g2g3_row <- pairwise_fu %>%
  filter(
    Comparison == "Group 2 vs Group 3"
  )

g1g2_row <- pairwise_fu %>%
  filter(
    Comparison == "Group 1 vs Group 2"
  )

g1g3_row <- pairwise_fu %>%
  filter(
    Comparison == "Group 1 vs Group 3"
  )


test_caption <- paste0(
  "Follow-up prevalence comparison: global ",
  global_method,
  " p ",
  ifelse(
    global_p < 0.001,
    "<0.001",
    paste0(
      "= ",
      fmt_p(
        global_p
      )
    )
  ),
  ". Baseline-to-follow-up within-group p-values were not tested because baseline ",
  "25(OH)D status is part of the exposure-group definition."
)


efigure5 <- ggplot(
  plot_dat,
  aes(
    x = Visit_num,
    y = Low_pct,
    group = Group,
    colour = Group
  )
) +

  geom_hline(
    yintercept = c(
      25,
      50,
      75
    ),
    colour = "#D8D8D8",
    linewidth = 0.45,
    linetype = "dashed"
  ) +

  geom_line(
    linewidth = 1.25,
    lineend = "round"
  ) +

  geom_errorbar(
    aes(
      ymin = Low_CI_low,
      ymax = Low_CI_high
    ),
    width = 0.07,
    linewidth = 0.8
  ) +

  geom_point(
    shape = 21,
    fill = "white",
    size = 4.2,
    stroke = 1.4
  ) +

  geom_text(
    aes(
      y = Label_y,
      label = Point_label
    ),
    family = "Arial",
    fontface = "bold",
    size = 3.35,
    lineheight = 0.95,
    colour = "#3A3A3A",
    show.legend = FALSE
  ) +

  geom_text(
    data = change_summary,
    aes(
      x = x,
      y = y,
      label = Delta_label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.15,
    colour = "#3A3A3A"
  ) +

  geom_text(
    data = normalization_annot,
    aes(
      x = x,
      y = y,
      label = Label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.0,
    colour = "#3A3A3A"
  ) +

  geom_text(
    data = normalization_p_annot,
    aes(
      x = x,
      y = y,
      label = Label
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.0,
    colour = "#3A3A3A"
  ) +

  facet_wrap(
    ~Group,
    nrow = 1,
    labeller = as_labeller(
      group_strip_labels
    )
  ) +

  scale_colour_manual(
    values = group_colors
  ) +

  scale_x_continuous(
    breaks = c(
      1,
      2
    ),
    labels = c(
      "Baseline",
      "Follow-up"
    ),
    limits = c(
      0.78,
      2.22
    )
  ) +

  scale_y_continuous(
    breaks = seq(
      0,
      100,
      by = 25
    ),
    limits = c(
      0,
      110
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    x = NULL,
    y = "Low vitamin D status (<30 ng/mL), %",
    caption = test_caption
  ) +

  theme_classic(
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

    strip.background = element_blank(),

    strip.text = element_text(
      family = "Arial",
      face = "bold",
      size = 10.4,
      colour = "#252525",
      margin = margin(
        b = 8
      )
    ),

    axis.text = element_text(
      family = "Arial",
      colour = "#252525",
      size = 9.5
    ),

    axis.title.y = element_text(
      family = "Arial",
      face = "bold",
      colour = "#252525",
      size = 10.5,
      margin = margin(
        r = 10
      )
    ),

    axis.line = element_line(
      colour = "#555555",
      linewidth = 0.55
    ),

    axis.ticks = element_line(
      colour = "#555555",
      linewidth = 0.45
    ),

    panel.spacing = grid::unit(
      1.25,
      "lines"
    ),

    legend.position = "none",

    plot.caption = element_text(
      family = "Arial",
      colour = "#454545",
      size = 9.0,
      hjust = 0,
      margin = margin(
        t = 12
      )
    ),

    plot.margin = margin(
      15,
      18,
      12,
      15
    )
  )


print(
  efigure5
)


# =============================================================================
# 13. Export Supplementary Figure S2
# =============================================================================

png_file <- file.path(
  output_dir,
  "Supplementary_Figure_S2_vitamin_D_status.png"
)

tiff_file <- file.path(
  output_dir,
  "Supplementary_Figure_S2_vitamin_D_status.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Figure_S2_vitamin_D_status.pdf"
)

ggsave(
  filename = png_file,
  plot = efigure5,
  width = 11.0,
  height = 5.8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = efigure5,
  width = 11.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = efigure5,
  width = 11.0,
  height = 5.8,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 14. Export statistics workbook
# =============================================================================

followup_global_test <- tibble(
  Outcome = "Low vitamin D status (<30 ng/mL) at follow-up",
  Method = global_method,
  p_value = global_p
)


normalization_test_table <- tibble(
  Outcome = "Normalization to ≥30 ng/mL among baseline-low participants",
  Comparison = "Group 2 vs Group 3",
  Method = norm_method,
  p_value = normalization_p
)


analysis_info <- tibble(
  Item = c(
    "Low vitamin D definition",
    "Higher/normalized status definition",
    "Group 1 definition",
    "Group 2 definition",
    "Group 3 definition",
    "Baseline group test",
    "Follow-up global test",
    "Follow-up pairwise tests",
    "Multiplicity",
    "Proportion confidence intervals",
    "Within-group paired p-values",
    "R version"
  ),

  Value = c(
    "Serum 25(OH)D <30 ng/mL",
    "Serum 25(OH)D ≥30 ng/mL",
    "Baseline ≥30 ng/mL; no vitamin D prescription",
    "Baseline <30 ng/mL; no vitamin D prescription",
    "Baseline <30 ng/mL; vitamin D prescription",
    paste(
      "Not performed because baseline vitamin D status",
      "is part of the three-group definition"
    ),
    paste(
      "Pearson chi-square unless expected cell count <5;",
      "otherwise Fisher exact"
    ),
    paste(
      "Pearson chi-square unless expected cell count <5;",
      "otherwise Fisher exact"
    ),
    "Holm adjustment across the three follow-up pairwise comparisons",
    "Exact Clopper-Pearson 95% confidence intervals",
    paste(
      "Not used because baseline status is structurally determined",
      "by exposure-group definition; transitions are summarized descriptively"
    ),
    "R 4.5.1"
  )
)


xlsx_out <- file.path(
  output_dir,
  "Supplementary_Figure_S2_vitamin_D_status_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Status_Summary"
)

addWorksheet(
  wb,
  "Status_by_Group_Visit"
)

addWorksheet(
  wb,
  "Transition_Detail"
)

addWorksheet(
  wb,
  "Followup_Global"
)

addWorksheet(
  wb,
  "Followup_Pairwise"
)

addWorksheet(
  wb,
  "Normalization_G2_G3"
)

addWorksheet(
  wb,
  "Normalization_Test"
)

addWorksheet(
  wb,
  "Analysis_Info"
)


writeData(
  wb,
  "Status_Summary",
  etable5
)

writeData(
  wb,
  "Status_by_Group_Visit",
  status_summary
)

writeData(
  wb,
  "Transition_Detail",
  transition_detail
)

writeData(
  wb,
  "Followup_Global",
  followup_global_test
)

writeData(
  wb,
  "Followup_Pairwise",
  pairwise_fu
)

writeData(
  wb,
  "Normalization_G2_G3",
  normalization_summary
)

writeData(
  wb,
  "Normalization_Test",
  normalization_test_table
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

  tmp <- readWorkbook(
    wb,
    sheet = sh
  )

  if (ncol(tmp) > 0) {

    addStyle(
      wb,
      sh,
      header_style,
      rows = 1,
      cols = seq_len(
        ncol(tmp)
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
        ncol(tmp)
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
  xlsx_out,
  overwrite = TRUE
)


# =============================================================================
# 15. Console report
# =============================================================================

cat("\n============================================================\n")
cat("Status-transition summary\n")
cat("============================================================\n")
print(
  etable5
)

cat("\nFollow-up global comparison:\n")
print(
  followup_global_test
)

cat("\nFollow-up pairwise comparisons:\n")
print(
  pairwise_fu
)

cat("\nNormalization among baseline-low participants:\n")
print(
  normalization_summary
)

cat("\nGroup 2 vs Group 3 normalization test:\n")
print(
  normalization_test_table
)

cat("\nTransition detail:\n")
print(
  transition_detail
)

cat("\nSaved files:\n")
cat(xlsx_out, "\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat("============================================================\n")
