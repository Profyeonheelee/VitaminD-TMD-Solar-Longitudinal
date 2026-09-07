# Table 6. Associations of ambient solar radiation with serum 25(OH)D and pain outcomes

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Sex", "Age_years", "Symptom_duration_months", "Clinical_followup_months",
  "Baseline_25OHD_ng_mL", "Followup_25OHD_ng_mL", "Delta_25OHD_ng_mL",
  "Baseline_pain_intensity_VAS", "Followup_pain_intensity_VAS", "Pain_intensity_reduction_VAS",
  "SolarRad_preBL_60d_mean_MJm2", "SolarRad_preFU_60d_mean_MJm2", "Delta_SolarRad_60d_MJm2"
)

dat <- read_analysis_data(required_vars) |>
  mutate(Sex = factor(Sex), log_symptom_duration = log1p(Symptom_duration_months))

analysis_spec <- list(
  BL_vitd = list(
    section = "Biochemical outcomes", label = "Pre-BL 60-day solar radiation → baseline 25(OH)D",
    exposure = "SolarRad_preBL_60d_mean_MJm2", outcome = "Baseline_25OHD_ng_mL",
    covariates = c("Age_years", "Sex")
  ),
  FU_vitd = list(
    section = "Biochemical outcomes", label = "Pre-FU 60-day solar radiation → follow-up 25(OH)D",
    exposure = "SolarRad_preFU_60d_mean_MJm2", outcome = "Followup_25OHD_ng_mL",
    covariates = c("Baseline_25OHD_ng_mL", "Age_years", "Sex", "Clinical_followup_months")
  ),
  Delta_vitd = list(
    section = "Biochemical outcomes", label = "Δ60-day solar radiation → Δ25(OH)D",
    exposure = "Delta_SolarRad_60d_MJm2", outcome = "Delta_25OHD_ng_mL",
    covariates = c("Baseline_25OHD_ng_mL", "Age_years", "Sex", "Clinical_followup_months")
  ),
  BL_pain = list(
    section = "Pain outcomes", label = "Pre-BL 60-day solar radiation → baseline VAS",
    exposure = "SolarRad_preBL_60d_mean_MJm2", outcome = "Baseline_pain_intensity_VAS",
    covariates = c("Age_years", "Sex", "log_symptom_duration")
  ),
  FU_pain = list(
    section = "Pain outcomes", label = "Pre-FU 60-day solar radiation → follow-up VAS",
    exposure = "SolarRad_preFU_60d_mean_MJm2", outcome = "Followup_pain_intensity_VAS",
    covariates = c("Baseline_pain_intensity_VAS", "Age_years", "Sex", "Clinical_followup_months")
  ),
  Delta_pain = list(
    section = "Pain outcomes", label = "Δ60-day solar radiation → pain reduction",
    exposure = "Delta_SolarRad_60d_MJm2", outcome = "Pain_intensity_reduction_VAS",
    covariates = c("Baseline_pain_intensity_VAS", "Age_years", "Sex", "Clinical_followup_months")
  )
)

results <- lapply(analysis_spec, function(spec) {
  unadjusted <- spearman_row(dat[[spec$exposure]], dat[[spec$outcome]])
  model <- lm(reformulate(c(spec$exposure, spec$covariates), response = spec$outcome), data = dat)
  adjusted <- coef_row(model, spec$exposure)
  data.frame(
    section = spec$section, association = spec$label, n = unadjusted["n"],
    rho = unadjusted["rho"], p = unadjusted["p"],
    beta = adjusted["estimate"], lower = adjusted["lower"], upper = adjusted["upper"],
    adjusted_p = adjusted["p"]
  )
}) |> bind_rows()

results$fdr_q <- p.adjust(results$p, method = "BH")
results$adjusted_fdr_q <- p.adjust(results$adjusted_p, method = "BH")

section_row <- function(label) data.frame(
  `Exposure–outcome association` = label, N = "", `Spearman ρ` = "", `p-value` = "",
  `FDR q-value` = "", `Adjusted β per +1 MJ/m²/day (95% CI)` = "", `Adjusted p` = "",
  check.names = FALSE
)

format_rows <- function(x) x |>
  transmute(
    `Exposure–outcome association` = association, N = as.integer(n),
    `Spearman ρ` = fmt_num(rho, 3), `p-value` = fmt_p(p), `FDR q-value` = fmt_p(fdr_q),
    `Adjusted β per +1 MJ/m²/day (95% CI)` = paste0(fmt_num(beta, 2), " (", fmt_num(lower, 2), " to ", fmt_num(upper, 2), ")"),
    `Adjusted p` = fmt_p(adjusted_p)
  )

table_6 <- bind_rows(
  section_row("Biochemical outcomes"),
  format_rows(filter(results, section == "Biochemical outcomes")),
  section_row("Pain outcomes"),
  format_rows(filter(results, section == "Pain outcomes"))
)

write_table_outputs(table_6, list(Adjusted_FDR_audit = results), "Table_6_solar_radiation_associations")
