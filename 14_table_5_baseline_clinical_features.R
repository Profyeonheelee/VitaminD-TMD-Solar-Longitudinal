# Table 5. Associations of baseline serum 25(OH)D with TMD-related clinical features

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Sex", "Age_years", "Symptom_duration_months", "Baseline_25OHD_ng_mL",
  "Baseline_pain_intensity_VAS", "Baseline_CMO_mm", "Baseline_MMO_mm", "Bruxism",
  "Clenching", "Parafunction", "TMJ_noise", "Self_reported_TMJ_noise",
  "Objective_locking", "Self_reported_locking", "Stiffness", "Trauma_history",
  "SolarRad_preBL_60d_mean_MJm2"
)

dat <- read_analysis_data(required_vars) |>
  mutate(
    Sex = factor(Sex),
    log_symptom_duration = log1p(Symptom_duration_months),
    Bruxism_bin = to_binary01(Bruxism, "Bruxism"),
    Clenching_bin = to_binary01(Clenching, "Clenching"),
    Parafunction_bin = to_binary01(Parafunction, "Parafunction"),
    TMJ_noise_bin = to_binary01(TMJ_noise, "TMJ_noise"),
    Self_reported_TMJ_noise_bin = to_binary01(Self_reported_TMJ_noise, "Self_reported_TMJ_noise"),
    Objective_locking_bin = to_binary01(Objective_locking, "Objective_locking"),
    Self_reported_locking_bin = to_binary01(Self_reported_locking, "Self_reported_locking"),
    Stiffness_bin = to_binary01(Stiffness, "Stiffness"),
    Trauma_history_bin = to_binary01(Trauma_history, "Trauma_history"),
    Baseline_25OHD_per10 = Baseline_25OHD_ng_mL / 10
  )

feature_spec <- data.frame(
  label = c(
    "Present TMD-related pain intensity, VAS", "PFO, mm", "MUO, mm",
    "Self-reported bruxism", "Clenching", "Parafunction",
    "Clinically detected TMJ noise", "Self-reported TMJ noise",
    "Clinically detected jaw locking", "Self-reported jaw locking",
    "Stiffness", "Trauma history"
  ),
  variable = c(
    "Baseline_pain_intensity_VAS", "Baseline_CMO_mm", "Baseline_MMO_mm",
    "Bruxism_bin", "Clenching_bin", "Parafunction_bin", "TMJ_noise_bin",
    "Self_reported_TMJ_noise_bin", "Objective_locking_bin",
    "Self_reported_locking_bin", "Stiffness_bin", "Trauma_history_bin"
  ),
  type = c(rep("continuous", 3), rep("binary", 9)),
  stringsAsFactors = FALSE
)

results <- lapply(seq_len(nrow(feature_spec)), function(i) {
  variable <- feature_spec$variable[i]
  unadjusted <- spearman_row(dat$Baseline_25OHD_ng_mL, dat[[variable]])
  formula <- reformulate(
    c("Baseline_25OHD_per10", "Age_years", "Sex", "log_symptom_duration", "SolarRad_preBL_60d_mean_MJm2"),
    response = variable
  )
  model <- if (feature_spec$type[i] == "binary") glm(formula, family = binomial(), data = dat) else lm(formula, data = dat)
  adjusted <- coef_row(model, "Baseline_25OHD_per10", exponentiate = feature_spec$type[i] == "binary")
  data.frame(
    feature = feature_spec$label[i], type = feature_spec$type[i],
    n = unadjusted["n"], rho = unadjusted["rho"], p = unadjusted["p"],
    effect = ifelse(feature_spec$type[i] == "binary", "OR", "β"),
    estimate = adjusted["estimate"], lower = adjusted["lower"], upper = adjusted["upper"],
    adjusted_p = adjusted["p"]
  )
}) |> bind_rows()

results$fdr_q <- p.adjust(results$p, method = "BH")
results$adjusted_fdr_q <- p.adjust(results$adjusted_p, method = "BH")

table_5 <- results |>
  transmute(
    `Clinical feature` = feature,
    N = as.integer(n),
    `Spearman ρ` = fmt_num(rho, 3),
    `p-value` = fmt_p(p),
    `FDR q-value` = fmt_p(fdr_q),
    Effect = effect,
    `Adjusted estimate per 10-ng/mL higher 25(OH)D (95% CI)` = paste0(
      fmt_num(estimate, 2), " (", fmt_num(lower, 2), " to ", fmt_num(upper, 2), ")"
    ),
    `Adjusted p` = fmt_p(adjusted_p)
  )

write_table_outputs(table_5, list(Adjusted_FDR_audit = results), "Table_5_baseline_clinical_features")
