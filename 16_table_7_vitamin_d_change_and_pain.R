# Table 7. Association between longitudinal change in serum 25(OH)D and pain improvement

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Sex", "Age_years", "Symptom_duration_months", "Clinical_followup_months",
  "Baseline_25OHD_ng_mL", "Delta_25OHD_ng_mL", "Vitamin_D_prescription",
  "Baseline_pain_intensity_VAS", "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS", "Delta_SolarRad_60d_MJm2"
)

dat <- read_analysis_data(required_vars) |>
  mutate(
    Sex = factor(Sex),
    Rx = to_binary01(Vitamin_D_prescription, "Vitamin_D_prescription"),
    log_symptom_duration = log1p(Symptom_duration_months),
    Delta_25OHD_per10 = Delta_25OHD_ng_mL / 10
  )

unadjusted <- spearman_row(dat$Delta_25OHD_ng_mL, dat$Pain_intensity_reduction_VAS)

model <- lm(
  Pain_intensity_reduction_VAS ~ Delta_25OHD_per10 + Baseline_pain_intensity_VAS +
    Baseline_25OHD_ng_mL + Age_years + Sex + log_symptom_duration +
    Clinical_followup_months + Delta_SolarRad_60d_MJm2 + Rx,
  data = dat
)
adjusted <- coef_row(model, "Delta_25OHD_per10")

sensitivity_model <- lm(
  Followup_pain_intensity_VAS ~ Delta_25OHD_per10 + Baseline_pain_intensity_VAS +
    Baseline_25OHD_ng_mL + Age_years + Sex + log_symptom_duration +
    Clinical_followup_months + Delta_SolarRad_60d_MJm2 + Rx,
  data = dat
)
sensitivity <- coef_row(sensitivity_model, "Delta_25OHD_per10")

table_7 <- data.frame(
  Analysis = c("Unadjusted association", "Multivariable-adjusted association"),
  N = c(as.integer(unadjusted["n"]), stats::nobs(model)),
  Exposure = c("Δ25(OH)D", "Δ25(OH)D, per 10-ng/mL increase"),
  Outcome = c("Pain reduction, VAS", "Pain reduction, VAS"),
  `Effect measure` = c("Spearman ρ", "β"),
  `Estimate (95% CI)` = c(
    fmt_num(unadjusted["rho"], 3),
    paste0(fmt_num(adjusted["estimate"], 2), " (", fmt_num(adjusted["lower"], 2), " to ", fmt_num(adjusted["upper"], 2), ")")
  ),
  `p-value` = c(fmt_p(unadjusted["p"]), fmt_p(adjusted["p"])),
  check.names = FALSE
)

audit <- data.frame(
  model = c("Pain reduction", "Follow-up VAS sensitivity"),
  estimate = c(adjusted["estimate"], sensitivity["estimate"]),
  lower = c(adjusted["lower"], sensitivity["lower"]),
  upper = c(adjusted["upper"], sensitivity["upper"]),
  p = c(adjusted["p"], sensitivity["p"])
)

write_table_outputs(table_7, list(Sensitivity_model = audit), "Table_7_vitamin_d_change_and_pain")
