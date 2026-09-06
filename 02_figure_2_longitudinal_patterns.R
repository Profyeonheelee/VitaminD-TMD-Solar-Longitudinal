# =============================================================================
# Figure 2. Longitudinal serum 25(OH)D patterns, ambient solar radiation, and pain improvement
# Panels A-F
#
# Panels
#   A. Serum 25(OH)D trajectories across three exposure groups
#   B. Follow-up serum 25(OH)D distributions across exposure groups
#   C. Pre-baseline 60-day solar radiation -> baseline serum 25(OH)D
#   D. Pre-follow-up 60-day solar radiation -> follow-up serum 25(OH)D
#   E. Change in 60-day solar radiation -> change in serum 25(OH)D
#   F. Change in serum 25(OH)D -> pain reduction
#
# Figure rules
#   - No overall title inside figure
#   - No panel letters inside the figure
#   - Y-axis title typography unified
#   - Grayscale default
#   - Muted orange highlights Group 3 and statistically significant
#     adjusted associations
#   - All p-values shown in the figure are calculated in this script
#
# Statistical principles
#   - Panel B uses baseline-adjusted follow-up 25(OH)D, not raw change-score
#     group inference, to reduce regression-to-the-mean concerns.
#   - HC3 heteroskedasticity-robust covariance is used for adjusted models.
#   - Panel F beta is expressed per +10-ng/mL increase in Delta 25(OH)D.
#
# Outputs
#   Figure_2_longitudinal_patterns.png
#   Figure_2_longitudinal_patterns.tiff
#   Figure_2_longitudinal_patterns.pdf
#   Figure_2_longitudinal_patterns_statistics.xlsx
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
output_dir <- file.path("outputs", "Figure_2")
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
  "Three_group_code",
  "Vitamin_D_prescription",
  "Clinical_followup_months",
  "Baseline_pain_intensity_VAS",
  "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS",
  "SolarRad_preBL_60d_mean_MJm2",
  "SolarRad_preFU_60d_mean_MJm2",
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

    Group = factor(
      Three_group_code,
      levels = c(1, 2, 3),
      labels = c("Group 1", "Group 2", "Group 3")
    ),

    Rx = as.integer(Vitamin_D_prescription),

    log_symptom_duration = log1p(Symptom_duration_months),

    Delta_25OHD_per10 = Delta_25OHD_ng_mL / 10
  )

if (nrow(dat) != 241) {
  stop("Expected N=241; found N=", nrow(dat))
}

if (anyDuplicated(dat$Study_ID) > 0) {
  stop("Duplicate Study_ID detected.")
}

group_counts <- table(dat$Group)

if (!all(group_counts == c(41, 51, 149))) {
  stop(
    "Unexpected group counts. Found: ",
    paste(
      names(group_counts),
      as.integer(group_counts),
      collapse = "; "
    )
  )
}

if (!all(na.omit(unique(dat$Rx)) %in% c(0, 1))) {
  stop("Vitamin_D_prescription must be coded 0/1.")
}


# =============================================================================
# 3. Formatting helpers
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
# 4. HC3 helper functions
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


model_audit <- function(model) {

  V <- sandwich::vcovHC(model, type = "HC3")

  b <- coef(model)
  se <- sqrt(diag(V))
  df_res <- df.residual(model)
  crit <- qt(0.975, df = df_res)

  p <- 2 * pt(
    abs(b / se),
    df = df_res,
    lower.tail = FALSE
  )

  tibble(
    Term = names(b),
    Estimate = unname(b),
    HC3_SE = unname(se),
    CI_low = unname(b - crit * se),
    CI_high = unname(b + crit * se),
    HC3_p = unname(p)
  )
}


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
# 5. Panel A: serum 25(OH)D trajectories by group
# =============================================================================

traj_long <- dat %>%
  select(
    Study_ID,
    Group,
    Baseline_25OHD_ng_mL,
    Followup_25OHD_ng_mL
  ) %>%
  pivot_longer(
    cols = c(
      Baseline_25OHD_ng_mL,
      Followup_25OHD_ng_mL
    ),
    names_to = "Time",
    values_to = "OHD"
  ) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "Baseline_25OHD_ng_mL",
        "Followup_25OHD_ng_mL"
      ),
      labels = c("Baseline", "Follow-up")
    )
  ) %>%
  filter(
    !is.na(OHD),
    !is.na(Group)
  )

traj_summary <- traj_long %>%
  group_by(
    Group,
    Time
  ) %>%
  summarise(
    N = n(),
    Mean = mean(OHD),
    SD = sd(OHD),
    SE = SD / sqrt(N),
    tcrit = qt(0.975, df = N - 1),
    CI_low = Mean - tcrit * SE,
    CI_high = Mean + tcrit * SE,
    .groups = "drop"
  )


# =============================================================================
# 6. Panel B: follow-up serum 25(OH)D distributions by group
#    Robust baseline-adjusted inference
# =============================================================================

dat_b <- dat %>%
  select(
    Group,
    Followup_25OHD_ng_mL,
    Baseline_25OHD_ng_mL,
    Age_years,
    Sex
  ) %>%
  filter(complete.cases(.))

mod_b <- lm(
  Followup_25OHD_ng_mL ~
    Group +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex,
  data = dat_b
)

V_b <- sandwich::vcovHC(
  mod_b,
  type = "HC3"
)

b_coef <- coef(mod_b)
df_b <- df.residual(mod_b)

g2_term <- "GroupGroup 2"
g3_term <- "GroupGroup 3"

if (!all(c(g2_term, g3_term) %in% names(b_coef))) {
  stop("Could not locate expected group coefficients in Panel B model.")
}

# Global robust Wald F test
R_global <- matrix(
  0,
  nrow = 2,
  ncol = length(b_coef)
)

colnames(R_global) <- names(b_coef)

R_global[1, g2_term] <- 1
R_global[2, g3_term] <- 1

Rb <- as.numeric(R_global %*% b_coef)

RVRT <- R_global %*%
  V_b %*%
  t(R_global)

W <- as.numeric(
  t(Rb) %*%
    solve(RVRT) %*%
    Rb
)

q_global <- 2
F_global <- W / q_global

p_global_b <- pf(
  F_global,
  df1 = q_global,
  df2 = df_b,
  lower.tail = FALSE
)


contrast_test <- function(
  contrast_vector,
  label
) {

  contrast_vector <- as.numeric(contrast_vector)
  names(contrast_vector) <- names(b_coef)

  est <- sum(
    contrast_vector * b_coef
  )

  se <- sqrt(
    as.numeric(
      t(contrast_vector) %*%
        V_b %*%
        contrast_vector
    )
  )

  crit <- qt(
    0.975,
    df = df_b
  )

  tval <- est / se

  p <- 2 * pt(
    abs(tval),
    df = df_b,
    lower.tail = FALSE
  )

  tibble(
    Contrast = label,
    Estimate = est,
    SE = se,
    CI_low = est - crit * se,
    CI_high = est + crit * se,
    p_raw = p
  )
}


c_g1_g2 <- rep(0, length(b_coef))
names(c_g1_g2) <- names(b_coef)
c_g1_g2[g2_term] <- -1

c_g1_g3 <- rep(0, length(b_coef))
names(c_g1_g3) <- names(b_coef)
c_g1_g3[g3_term] <- -1

c_g2_g3 <- rep(0, length(b_coef))
names(c_g2_g3) <- names(b_coef)
c_g2_g3[g2_term] <- 1
c_g2_g3[g3_term] <- -1

pair_b <- bind_rows(
  contrast_test(
    c_g1_g2,
    "G1 vs G2"
  ),
  contrast_test(
    c_g1_g3,
    "G1 vs G3"
  ),
  contrast_test(
    c_g2_g3,
    "G2 vs G3"
  )
) %>%
  mutate(
    p_Holm = p.adjust(
      p_raw,
      method = "holm"
    )
  )


get_pair_p <- function(label) {
  pair_b$p_Holm[
    pair_b$Contrast == label
  ]
}


pair_b_annotation <- paste0(
  "Adjusted global p ",
  ifelse(
    p_global_b < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(p_global_b))
  ),
  "\nG1 vs G2: Holm p ",
  ifelse(
    get_pair_p("G1 vs G2") < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(get_pair_p("G1 vs G2")))
  ),
  "\nG1 vs G3: Holm p ",
  ifelse(
    get_pair_p("G1 vs G3") < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(get_pair_p("G1 vs G3")))
  ),
  "\nG2 vs G3: Holm p ",
  ifelse(
    get_pair_p("G2 vs G3") < 0.001,
    "< 0.001",
    paste0("= ", fmt_p(get_pair_p("G2 vs G3")))
  )
)


# =============================================================================
# 7. Panel C: pre-BL solar radiation -> baseline 25(OH)D
# =============================================================================

dat_c <- dat %>%
  select(
    SolarRad_preBL_60d_mean_MJm2,
    Baseline_25OHD_ng_mL,
    Age_years,
    Sex
  ) %>%
  filter(complete.cases(.))

sp_c <- suppressWarnings(
  cor.test(
    dat_c$SolarRad_preBL_60d_mean_MJm2,
    dat_c$Baseline_25OHD_ng_mL,
    method = "spearman",
    exact = FALSE
  )
)

mod_c <- lm(
  Baseline_25OHD_ng_mL ~
    SolarRad_preBL_60d_mean_MJm2 +
    Age_years +
    Sex,
  data = dat_c
)

adj_c <- hc3_term(
  mod_c,
  "SolarRad_preBL_60d_mean_MJm2"
)

grid_c <- seq(
  min(dat_c$SolarRad_preBL_60d_mean_MJm2),
  max(dat_c$SolarRad_preBL_60d_mean_MJm2),
  length.out = 180
)

curve_c <- marginal_robust_curve(
  mod_c,
  dat_c,
  "SolarRad_preBL_60d_mean_MJm2",
  grid_c
)


# =============================================================================
# 8. Panel D: pre-FU solar radiation -> follow-up 25(OH)D
# =============================================================================

dat_d <- dat %>%
  select(
    SolarRad_preFU_60d_mean_MJm2,
    Followup_25OHD_ng_mL,
    Baseline_25OHD_ng_mL,
    Age_years,
    Sex,
    Clinical_followup_months
  ) %>%
  filter(complete.cases(.))

sp_d <- suppressWarnings(
  cor.test(
    dat_d$SolarRad_preFU_60d_mean_MJm2,
    dat_d$Followup_25OHD_ng_mL,
    method = "spearman",
    exact = FALSE
  )
)

mod_d <- lm(
  Followup_25OHD_ng_mL ~
    SolarRad_preFU_60d_mean_MJm2 +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex +
    Clinical_followup_months,
  data = dat_d
)

adj_d <- hc3_term(
  mod_d,
  "SolarRad_preFU_60d_mean_MJm2"
)

grid_d <- seq(
  min(dat_d$SolarRad_preFU_60d_mean_MJm2),
  max(dat_d$SolarRad_preFU_60d_mean_MJm2),
  length.out = 180
)

curve_d <- marginal_robust_curve(
  mod_d,
  dat_d,
  "SolarRad_preFU_60d_mean_MJm2",
  grid_d
)


# =============================================================================
# 9. Panel E: Delta solar radiation -> Delta 25(OH)D
# =============================================================================

dat_e <- dat %>%
  select(
    Delta_SolarRad_60d_MJm2,
    Delta_25OHD_ng_mL,
    Baseline_25OHD_ng_mL,
    Age_years,
    Sex,
    Clinical_followup_months
  ) %>%
  filter(complete.cases(.))

sp_e <- suppressWarnings(
  cor.test(
    dat_e$Delta_SolarRad_60d_MJm2,
    dat_e$Delta_25OHD_ng_mL,
    method = "spearman",
    exact = FALSE
  )
)

mod_e <- lm(
  Delta_25OHD_ng_mL ~
    Delta_SolarRad_60d_MJm2 +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex +
    Clinical_followup_months,
  data = dat_e
)

adj_e <- hc3_term(
  mod_e,
  "Delta_SolarRad_60d_MJm2"
)

grid_e <- seq(
  min(dat_e$Delta_SolarRad_60d_MJm2),
  max(dat_e$Delta_SolarRad_60d_MJm2),
  length.out = 180
)

curve_e <- marginal_robust_curve(
  mod_e,
  dat_e,
  "Delta_SolarRad_60d_MJm2",
  grid_e
)


# =============================================================================
# 10. Panel F: Delta 25(OH)D -> pain reduction
# =============================================================================

dat_f <- dat %>%
  select(
    Delta_25OHD_ng_mL,
    Delta_25OHD_per10,
    Pain_intensity_reduction_VAS,
    Baseline_pain_intensity_VAS,
    Baseline_25OHD_ng_mL,
    Age_years,
    Sex,
    log_symptom_duration,
    Clinical_followup_months,
    Delta_SolarRad_60d_MJm2,
    Rx
  ) %>%
  filter(complete.cases(.))

sp_f <- suppressWarnings(
  cor.test(
    dat_f$Delta_25OHD_ng_mL,
    dat_f$Pain_intensity_reduction_VAS,
    method = "spearman",
    exact = FALSE
  )
)

mod_f <- lm(
  Pain_intensity_reduction_VAS ~
    Delta_25OHD_per10 +
    Baseline_pain_intensity_VAS +
    Baseline_25OHD_ng_mL +
    Age_years +
    Sex +
    log_symptom_duration +
    Clinical_followup_months +
    Delta_SolarRad_60d_MJm2 +
    Rx,
  data = dat_f
)

adj_f <- hc3_term(
  mod_f,
  "Delta_25OHD_per10"
)

grid_f_raw <- seq(
  min(dat_f$Delta_25OHD_ng_mL),
  max(dat_f$Delta_25OHD_ng_mL),
  length.out = 180
)

curve_f <- marginal_robust_curve(
  mod_f,
  dat_f,
  "Delta_25OHD_per10",
  grid_f_raw / 10
) %>%
  mutate(
    x_raw = x * 10
  )


# =============================================================================
# 11. Visual system
# =============================================================================

COL_TEXT <- "#252525"
COL_AXIS <- "#505050"
COL_POINT <- "#A6A6A6"
COL_POINT_EDGE <- "#666666"
COL_ZERO <- "#888888"
COL_PANEL_BG <- "#FCFCFC"

COL_ORANGE <- "#C96A1B"
COL_ORANGE_LIGHT <- "#E8B07A"

COL_G1 <- "#303030"
COL_G2 <- "#777777"
COL_G3 <- COL_ORANGE

COL_NONSIG <- "#4A4A4A"
COL_NONSIG_BAND <- "#BDBDBD"

sig_c <- adj_c$p < 0.05
sig_d <- adj_d$p < 0.05
sig_e <- adj_e$p < 0.05
sig_f <- adj_f$p < 0.05

line_c <- if (sig_c) COL_ORANGE else COL_NONSIG
band_c <- if (sig_c) COL_ORANGE_LIGHT else COL_NONSIG_BAND

line_d <- if (sig_d) COL_ORANGE else COL_NONSIG
band_d <- if (sig_d) COL_ORANGE_LIGHT else COL_NONSIG_BAND

line_e <- if (sig_e) COL_ORANGE else COL_NONSIG
band_e <- if (sig_e) COL_ORANGE_LIGHT else COL_NONSIG_BAND

line_f <- if (sig_f) COL_ORANGE else COL_NONSIG
band_f <- if (sig_f) COL_ORANGE_LIGHT else COL_NONSIG_BAND


# =============================================================================
# 12. Annotation helpers
# =============================================================================

make_env_annotation <- function(
  n,
  sp,
  adj,
  beta_label = "Adjusted β"
) {

  paste0(
    "n = ", n,
    "\nSpearman ρ = ",
    sprintf("%.3f", unname(sp$estimate)),
    "; p ",
    ifelse(
      sp$p.value < 0.001,
      "< 0.001",
      paste0("= ", fmt_p(sp$p.value))
    ),
    "\n",
    beta_label,
    " = ",
    fmt_ci(
      adj$beta,
      adj$ci_low,
      adj$ci_high,
      2
    ),
    "; p ",
    ifelse(
      adj$p < 0.001,
      "< 0.001",
      paste0("= ", fmt_p(adj$p))
    )
  )
}

ann_c <- make_env_annotation(
  nrow(dat_c),
  sp_c,
  adj_c
)

ann_d <- make_env_annotation(
  nrow(dat_d),
  sp_d,
  adj_d
)

ann_e <- make_env_annotation(
  nrow(dat_e),
  sp_e,
  adj_e
)

ann_f <- make_env_annotation(
  nrow(dat_f),
  sp_f,
  adj_f,
  beta_label = "Adjusted β / +10 ng/mL"
)


# =============================================================================
# 13. Shared theme
# =============================================================================

theme_fig3 <- theme_classic(
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

    axis.title = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.2,
      face = "bold"
    ),

    axis.title.y = element_text(
      family = "Arial",
      colour = COL_TEXT,
      size = 10.2,
      face = "bold"
    ),

    legend.text = element_text(
      family = "Arial",
      size = 8.8,
      colour = COL_TEXT
    ),

    legend.title = element_blank(),

    plot.margin = margin(
      10, 12, 10, 14
    )
  )


# =============================================================================
# 14. Panel A
# =============================================================================

panel_a <- ggplot(
  traj_summary,
  aes(
    x = Time,
    y = Mean,
    group = Group,
    colour = Group,
    shape = Group,
    linetype = Group
  )
) +

  geom_hline(
    yintercept = 20,
    linewidth = 0.45,
    linetype = "dashed",
    colour = "#8A8A8A"
  ) +

  geom_hline(
    yintercept = 30,
    linewidth = 0.45,
    linetype = "dotted",
    colour = "#8A8A8A"
  ) +

  geom_errorbar(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    width = 0.07,
    linewidth = 0.70
  ) +

  geom_line(
    linewidth = 1.05
  ) +

  geom_point(
    size = 3.1,
    stroke = 0.75,
    fill = "white"
  ) +

  annotate(
    "text",
    x = 2.03,
    y = 20.6,
    label = "20 ng/mL",
    hjust = 0,
    size = 2.9,
    family = "Arial",
    colour = "#666666"
  ) +

  annotate(
    "text",
    x = 2.03,
    y = 30.6,
    label = "30 ng/mL",
    hjust = 0,
    size = 2.9,
    family = "Arial",
    colour = "#666666"
  ) +

  scale_colour_manual(
    values = c(
      "Group 1" = COL_G1,
      "Group 2" = COL_G2,
      "Group 3" = COL_G3
    ),
    labels = c(
      "Group 1  VitD sufficient / no Rx",
      "Group 2  Low VitD / no Rx",
      "Group 3  Low VitD / Rx"
    )
  ) +

  scale_shape_manual(
    values = c(
      "Group 1" = 21,
      "Group 2" = 24,
      "Group 3" = 22
    ),
    labels = c(
      "Group 1  VitD sufficient / no Rx",
      "Group 2  Low VitD / no Rx",
      "Group 3  Low VitD / Rx"
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Group 1" = "solid",
      "Group 2" = "longdash",
      "Group 3" = "dotdash"
    ),
    labels = c(
      "Group 1  VitD sufficient / no Rx",
      "Group 2  Low VitD / no Rx",
      "Group 3  Low VitD / Rx"
    )
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = "Serum 25(OH)D, ng/mL"
  ) +

  theme_fig3 +
  theme(
    legend.position = "bottom"
  )


# =============================================================================
# 15. Panel B
# =============================================================================

panel_b <- ggplot(
  dat_b,
  aes(
    x = Group,
    y = Followup_25OHD_ng_mL
  )
) +

  geom_violin(
    aes(
      fill = Group,
      colour = Group
    ),
    width = 0.88,
    alpha = 0.18,
    linewidth = 0.65,
    trim = FALSE
  ) +

  geom_boxplot(
    aes(
      colour = Group
    ),
    width = 0.22,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.65
  ) +

  geom_point(
    position = position_jitter(
      width = 0.10,
      height = 0
    ),
    shape = 21,
    size = 1.55,
    stroke = 0.25,
    fill = "#B2B2B2",
    colour = "#606060",
    alpha = 0.58
  ) +

  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = pair_b_annotation,
    hjust = 1.03,
    vjust = 1.10,
    size = 2.55,
    family = "Arial",
    lineheight = 1.08,
    colour = if (p_global_b < 0.05) COL_ORANGE else COL_TEXT,
    fill = "white",
    label.size = 0.28,
    label.padding = grid::unit(0.38, "lines"),
    label.r = grid::unit(0.05, "lines")
  ) +

  scale_fill_manual(
    values = c(
      "Group 1" = "#B8B8B8",
      "Group 2" = "#D0D0D0",
      "Group 3" = COL_ORANGE_LIGHT
    )
  ) +

  scale_colour_manual(
    values = c(
      "Group 1" = COL_G1,
      "Group 2" = COL_G2,
      "Group 3" = COL_G3
    )
  ) +

  scale_x_discrete(
    labels = c(
      "Group 1\nVitD sufficient / no Rx",
      "Group 2\nLow VitD / no Rx",
      "Group 3\nLow VitD / Rx"
    )
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = "Follow-up serum 25(OH)D, ng/mL"
  ) +

  guides(
    fill = "none",
    colour = "none"
  ) +

  theme_fig3 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      family = "Arial",
      size = 8.5,
      lineheight = 0.95
    )
  )


# =============================================================================
# 16. Scatter-panel builder
# =============================================================================

make_scatter_panel <- function(
  data,
  xvar,
  yvar,
  curve,
  x_curve_var = "x",
  panel_letter,
  x_label,
  y_label,
  annotation_text,
  significant,
  line_colour,
  band_colour,
  add_zero_x = FALSE,
  add_zero_y = FALSE
) {

  p <- ggplot(
    data,
    aes(
      x = .data[[xvar]],
      y = .data[[yvar]]
    )
  )

  if (add_zero_y) {
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.50,
        colour = COL_ZERO
      )
  }

  if (add_zero_x) {
    p <- p +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.50,
        colour = COL_ZERO
      )
  }

  p +
    geom_point(
      shape = 21,
      size = 1.85,
      stroke = 0.30,
      fill = COL_POINT,
      colour = COL_POINT_EDGE,
      alpha = 0.52
    ) +

    geom_ribbon(
      data = curve,
      aes(
        x = .data[[x_curve_var]],
        ymin = lower,
        ymax = upper
      ),
      inherit.aes = FALSE,
      fill = band_colour,
      alpha = if (significant) 0.30 else 0.22
    ) +

    geom_line(
      data = curve,
      aes(
        x = .data[[x_curve_var]],
        y = fit
      ),
      inherit.aes = FALSE,
      colour = line_colour,
      linewidth = 1.18
    ) +

    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.03,
      vjust = 1.10,
      size = 2.60,
      family = "Arial",
      lineheight = 1.10,
      colour = if (significant) COL_ORANGE else COL_TEXT,
      fill = "white",
      label.size = 0.28,
      label.padding = grid::unit(0.38, "lines"),
      label.r = grid::unit(0.05, "lines")
    ) +

    coord_cartesian(
      clip = "off"
    ) +

    labs(
      x = x_label,
      y = y_label
    ) +

    theme_fig3 +
    theme(
      legend.position = "none"
    )
}


# =============================================================================
# 17. Panels C-F
# =============================================================================

panel_c <- make_scatter_panel(
  data = dat_c,
  xvar = "SolarRad_preBL_60d_mean_MJm2",
  yvar = "Baseline_25OHD_ng_mL",
  curve = curve_c,
  panel_letter = "C",
  x_label = "Pre-baseline 60-day solar radiation, MJ/m²/day",
  y_label = "Baseline 25(OH)D, ng/mL",
  annotation_text = ann_c,
  significant = sig_c,
  line_colour = line_c,
  band_colour = band_c
)

panel_d <- make_scatter_panel(
  data = dat_d,
  xvar = "SolarRad_preFU_60d_mean_MJm2",
  yvar = "Followup_25OHD_ng_mL",
  curve = curve_d,
  panel_letter = "D",
  x_label = "Pre-follow-up 60-day solar radiation, MJ/m²/day",
  y_label = "Follow-up 25(OH)D, ng/mL",
  annotation_text = ann_d,
  significant = sig_d,
  line_colour = line_d,
  band_colour = band_d
)

panel_e <- make_scatter_panel(
  data = dat_e,
  xvar = "Delta_SolarRad_60d_MJm2",
  yvar = "Delta_25OHD_ng_mL",
  curve = curve_e,
  panel_letter = "E",
  x_label = "Δ60-day mean solar radiation, MJ/m²/day",
  y_label = "Δ25(OH)D, ng/mL",
  annotation_text = ann_e,
  significant = sig_e,
  line_colour = line_e,
  band_colour = band_e,
  add_zero_x = TRUE,
  add_zero_y = TRUE
)

panel_f <- make_scatter_panel(
  data = dat_f,
  xvar = "Delta_25OHD_ng_mL",
  yvar = "Pain_intensity_reduction_VAS",
  curve = curve_f,
  x_curve_var = "x_raw",
  panel_letter = "F",
  x_label = "Δ25(OH)D, ng/mL",
  y_label = "Pain reduction, VAS",
  annotation_text = ann_f,
  significant = sig_f,
  line_colour = line_f,
  band_colour = band_f,
  add_zero_x = TRUE,
  add_zero_y = TRUE
)


# =============================================================================
# 18. Combine 6 panels
# =============================================================================

figure3 <- (
  panel_a + panel_b
) / (
  panel_c + panel_d
) / (
  panel_e + panel_f
) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

print(figure3)


# =============================================================================
# 19. Statistical audit tables
# =============================================================================

panel_b_global <- tibble(
  Analysis = "HC3 robust baseline-adjusted global group test",
  N = nrow(dat_b),
  F_statistic = F_global,
  df1 = q_global,
  df2 = df_b,
  p_value = p_global_b
)

main_stats <- bind_rows(
  tibble(
    Panel = "C",
    Exposure = "Pre-baseline 60-day solar radiation, MJ/m^2/day",
    Outcome = "Baseline 25(OH)D, ng/mL",
    N = nrow(dat_c),
    Spearman_rho = unname(sp_c$estimate),
    Spearman_p = sp_c$p.value,
    Adjusted_beta = adj_c$beta,
    HC3_CI_low = adj_c$ci_low,
    HC3_CI_high = adj_c$ci_high,
    Adjusted_p = adj_c$p,
    Effect_unit = "Per +1 MJ/m^2/day"
  ),
  tibble(
    Panel = "D",
    Exposure = "Pre-follow-up 60-day solar radiation, MJ/m^2/day",
    Outcome = "Follow-up 25(OH)D, ng/mL",
    N = nrow(dat_d),
    Spearman_rho = unname(sp_d$estimate),
    Spearman_p = sp_d$p.value,
    Adjusted_beta = adj_d$beta,
    HC3_CI_low = adj_d$ci_low,
    HC3_CI_high = adj_d$ci_high,
    Adjusted_p = adj_d$p,
    Effect_unit = "Per +1 MJ/m^2/day"
  ),
  tibble(
    Panel = "E",
    Exposure = "Delta 60-day solar radiation, MJ/m^2/day",
    Outcome = "Delta 25(OH)D, ng/mL",
    N = nrow(dat_e),
    Spearman_rho = unname(sp_e$estimate),
    Spearman_p = sp_e$p.value,
    Adjusted_beta = adj_e$beta,
    HC3_CI_low = adj_e$ci_low,
    HC3_CI_high = adj_e$ci_high,
    Adjusted_p = adj_e$p,
    Effect_unit = "Per +1 MJ/m^2/day"
  ),
  tibble(
    Panel = "F",
    Exposure = "Delta 25(OH)D",
    Outcome = "Pain reduction, VAS",
    N = nrow(dat_f),
    Spearman_rho = unname(sp_f$estimate),
    Spearman_p = sp_f$p.value,
    Adjusted_beta = adj_f$beta,
    HC3_CI_low = adj_f$ci_low,
    HC3_CI_high = adj_f$ci_high,
    Adjusted_p = adj_f$p,
    Effect_unit = "Per +10 ng/mL Delta 25(OH)D"
  )
)

analysis_info <- tibble(
  Item = c(
    "Panel A",
    "Panel B",
    "Panel B inference",
    "Panel C",
    "Panel D",
    "Panel E",
    "Panel F",
    "Primary environmental window",
    "Solar-radiation interpretation",
    "Adjusted uncertainty",
    "Color rule",
    "Panel labels",
    "Figure title/subtitle",
    "Missing data",
    "R version requested for manuscript"
  ),
  Value = c(
    "Observed serum 25(OH)D group means at baseline and follow-up with 95% CI",
    "Follow-up serum 25(OH)D distributions: violin + boxplot + jitter",
    paste(
      "HC3 robust model adjusted for baseline 25(OH)D, age, and sex;",
      "global Wald F test plus 3 Holm-adjusted pairwise contrasts"
    ),
    "Pre-baseline 60-day solar radiation vs baseline 25(OH)D; adjusted for age and sex",
    paste(
      "Pre-follow-up 60-day solar radiation vs follow-up 25(OH)D;",
      "adjusted for baseline 25(OH)D, age, sex, and clinical follow-up interval"
    ),
    paste(
      "Delta 60-day solar radiation vs Delta 25(OH)D;",
      "adjusted for baseline 25(OH)D, age, sex, and clinical follow-up interval"
    ),
    paste(
      "Delta 25(OH)D vs pain reduction;",
      "adjusted for baseline VAS, baseline 25(OH)D, age, sex,",
      "log1p symptom duration, clinical follow-up interval,",
      "Delta 60-day solar radiation, and vitamin D prescription status"
    ),
    "60-day",
    "Ambient solar-energy proxy; not direct individual UVB exposure",
    "HC3 heteroskedasticity-robust standard errors and 95% CI",
    "Muted orange highlights Group 3 and statistically significant adjusted associations",
    "No panel letters inside the figure",
    "None inside figure",
    "Outcome-specific complete-case analysis; no imputation",
    "R 4.5.1"
  )
)


# =============================================================================
# 20. Export figure
# =============================================================================

png_file <- file.path(
  output_dir,
  "Figure_2_longitudinal_patterns.png"
)

tiff_file <- file.path(
  output_dir,
  "Figure_2_longitudinal_patterns.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Figure_2_longitudinal_patterns.pdf"
)

ggsave(
  filename = png_file,
  plot = figure3,
  width = 11.4,
  height = 13.0,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure3,
  width = 11.4,
  height = 13.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = figure3,
  width = 11.4,
  height = 13.0,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


# =============================================================================
# 21. Export statistics workbook
# =============================================================================

stats_file <- file.path(
  output_dir,
  "Figure_2_longitudinal_patterns_statistics.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "PanelA_Trajectory")
addWorksheet(wb, "PanelB_Global")
addWorksheet(wb, "PanelB_Pairwise")
addWorksheet(wb, "PanelB_Model")
addWorksheet(wb, "PanelsCtoF_Main")
addWorksheet(wb, "PanelC_Model")
addWorksheet(wb, "PanelD_Model")
addWorksheet(wb, "PanelE_Model")
addWorksheet(wb, "PanelF_Model")
addWorksheet(wb, "PanelC_Curve")
addWorksheet(wb, "PanelD_Curve")
addWorksheet(wb, "PanelE_Curve")
addWorksheet(wb, "PanelF_Curve")
addWorksheet(wb, "Analysis_Info")

writeData(
  wb,
  "PanelA_Trajectory",
  traj_summary
)

writeData(
  wb,
  "PanelB_Global",
  panel_b_global
)

writeData(
  wb,
  "PanelB_Pairwise",
  pair_b
)

writeData(
  wb,
  "PanelB_Model",
  model_audit(mod_b)
)

writeData(
  wb,
  "PanelsCtoF_Main",
  main_stats
)

writeData(
  wb,
  "PanelC_Model",
  model_audit(mod_c)
)

writeData(
  wb,
  "PanelD_Model",
  model_audit(mod_d)
)

writeData(
  wb,
  "PanelE_Model",
  model_audit(mod_e)
)

writeData(
  wb,
  "PanelF_Model",
  model_audit(mod_f)
)

writeData(
  wb,
  "PanelC_Curve",
  curve_c
)

writeData(
  wb,
  "PanelD_Curve",
  curve_d
)

writeData(
  wb,
  "PanelE_Curve",
  curve_e
)

writeData(
  wb,
  "PanelF_Curve",
  curve_f
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
  widths = c(42, 110)
)

saveWorkbook(
  wb,
  stats_file,
  overwrite = TRUE
)


# =============================================================================
# 22. Console output
# =============================================================================

cat("\n============================================================\n")
cat("FIGURE 3 SIX-PANEL ANALYSIS COMPLETED\n")
cat("============================================================\n\n")

cat("Panel B adjusted group comparison\n")
cat(
  "Global p = ",
  fmt_p(p_global_b),
  "\n",
  sep = ""
)

print(
  pair_b %>%
    select(
      Contrast,
      Estimate,
      CI_low,
      CI_high,
      p_raw,
      p_Holm
    )
)

cat("\nPanels C-F\n")

for (i in seq_len(nrow(main_stats))) {

  r <- main_stats[i, ]

  cat(
    "\nPanel ", r$Panel, "\n",
    sep = ""
  )

  cat(
    "Spearman rho = ",
    sprintf("%.3f", r$Spearman_rho),
    "; p = ",
    fmt_p(r$Spearman_p),
    "\n",
    sep = ""
  )

  cat(
    "Adjusted beta = ",
    fmt_ci(
      r$Adjusted_beta,
      r$HC3_CI_low,
      r$HC3_CI_high,
      2
    ),
    "; p = ",
    fmt_p(r$Adjusted_p),
    "\n",
    sep = ""
  )
}

cat("\nSaved files:\n")
cat(png_file, "\n")
cat(tiff_file, "\n")
cat(pdf_file, "\n")
cat(stats_file, "\n")
cat("============================================================\n")
