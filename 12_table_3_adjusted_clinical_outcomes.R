# Table 3. Adjusted longitudinal clinical outcomes across vitamin D exposure groups

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Three_group_code", "Sex", "Age_years",
  "Baseline_pain_intensity_VAS", "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS", "Baseline_CMO_mm", "Followup_CMO_mm",
  "Baseline_MMO_mm", "Followup_MMO_mm"
)

dat <- read_analysis_data(required_vars) |>
  mutate(
    Group = make_group(Three_group_code), Sex = factor(Sex),
    Response_ge2 = ifelse(Baseline_pain_intensity_VAS >= 2, as.integer(Pain_intensity_reduction_VAS >= 2), NA_integer_),
    Response_ge30 = ifelse(Baseline_pain_intensity_VAS > 0, as.integer(Pain_intensity_reduction_VAS / Baseline_pain_intensity_VAS >= 0.30), NA_integer_)
  )

fit_outcome <- function(data, formula, exponentiate = FALSE) {
  model <- if (exponentiate) glm(formula, family = poisson(link = "log"), data = data) else lm(formula, data = data)
  vc <- hc3_vcov(model)
  terms <- grep("^Group", names(coef(model)), value = TRUE)
  list(
    model = model, vcov = vc, n = stats::nobs(model),
    global_p = robust_wald(model, terms, vc),
    contrasts = group_contrasts(model, exponentiate)
  )
}

models <- list(
  Pain = fit_outcome(dat, Followup_pain_intensity_VAS ~ Group + Baseline_pain_intensity_VAS + Age_years + Sex),
  PFO = fit_outcome(dat, Followup_CMO_mm ~ Group + Baseline_CMO_mm + Age_years + Sex),
  MUO = fit_outcome(dat, Followup_MMO_mm ~ Group + Baseline_MMO_mm + Age_years + Sex),
  Response2 = fit_outcome(dat |> filter(!is.na(Response_ge2)), Response_ge2 ~ Group + Baseline_pain_intensity_VAS + Age_years + Sex, TRUE),
  Response30 = fit_outcome(dat |> filter(!is.na(Response_ge30)), Response_ge30 ~ Group + Baseline_pain_intensity_VAS + Age_years + Sex, TRUE)
)

make_row <- function(label, fit, effect_label) {
  cts <- fit$contrasts
  data.frame(
    Outcome = label,
    `Analysis N` = fit$n,
    `Adjusted overall p-value` = fmt_p(fit$global_p),
    `Group 1 vs 2` = format_effect(cts[1, ], effect_label, cts[1, "p_holm"]),
    `Group 1 vs 3` = format_effect(cts[2, ], effect_label, cts[2, "p_holm"]),
    `Group 2 vs 3` = format_effect(cts[3, ], effect_label, cts[3, "p_holm"]),
    check.names = FALSE
  )
}

table_3 <- bind_rows(
  make_row("Follow-up present TMD-related pain intensity, VAS", models$Pain, "MD"),
  make_row("Follow-up PFO, mm", models$PFO, "MD"),
  make_row("Follow-up MUO, mm", models$MUO, "MD"),
  make_row("VAS reduction ≥2 points", models$Response2, "RR"),
  make_row("VAS reduction ≥30%", models$Response30, "RR")
)

audit <- bind_rows(lapply(names(models), function(nm) transform(as.data.frame(models[[nm]]$contrasts), outcome = nm)))
write_table_outputs(table_3, list(Contrasts = audit), "Table_3_adjusted_clinical_outcomes")
