# Table 4. Overlap-weighted associations of vitamin D prescription with
# longitudinal outcomes among patients with low vitamin D

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Sex", "Age_years", "Symptom_duration_months",
  "Baseline_25OHD_ng_mL", "Followup_25OHD_ng_mL", "Vitamin_D_prescription",
  "Low_Vitamin_D_cohort", "Clinical_followup_months", "Baseline_pain_intensity_VAS",
  "Followup_pain_intensity_VAS", "Pain_intensity_reduction_VAS",
  "Baseline_CMO_mm", "Followup_CMO_mm", "Baseline_MMO_mm", "Followup_MMO_mm",
  "TMJ_noise", "Trauma_history", "SolarRad_preBL_60d_mean_MJm2"
)

dat <- read_analysis_data(required_vars) |>
  mutate(
    Rx = to_binary01(Vitamin_D_prescription, "Vitamin_D_prescription"),
    Female = to_female01(Sex),
    TMJ_noise_bin = to_binary01(TMJ_noise, "TMJ_noise"),
    Trauma_bin = to_binary01(Trauma_history, "Trauma_history"),
    log_symptom_duration = log1p(Symptom_duration_months),
    Response_ge2 = ifelse(Baseline_pain_intensity_VAS >= 2, as.integer(Pain_intensity_reduction_VAS >= 2), NA_integer_),
    Response_ge30 = ifelse(Baseline_pain_intensity_VAS > 0, as.integer(Pain_intensity_reduction_VAS / Baseline_pain_intensity_VAS >= 0.30), NA_integer_)
  ) |>
  filter(Baseline_25OHD_ng_mL < 30)

ps_formula <- Rx ~ Age_years + Female + Baseline_25OHD_ng_mL +
  log_symptom_duration + Baseline_pain_intensity_VAS + TMJ_noise_bin +
  Trauma_bin + SolarRad_preBL_60d_mean_MJm2

ps_model <- glm(ps_formula, family = binomial(), data = dat)
dat$propensity_score <- predict(ps_model, newdata = dat, type = "response")
dat$overlap_weight <- ifelse(dat$Rx == 1L, 1 - dat$propensity_score, dat$propensity_score)

fit_weighted <- function(data, formula, family = gaussian()) {
  vars <- all.vars(formula)
  d <- data[complete.cases(data[, c(vars, "overlap_weight")]), ]
  design <- survey::svydesign(ids = ~1, weights = ~overlap_weight, data = d)
  model <- survey::svyglm(formula, design = design, family = family)
  term <- "Rx"
  b <- unname(coef(model)[term])
  se <- unname(sqrt(diag(vcov(model)))[term])
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  lo <- b - qnorm(0.975) * se
  hi <- b + qnorm(0.975) * se
  exponentiate <- identical(family$link, "log")
  estimate <- if (exponentiate) exp(c(estimate = b, lower = lo, upper = hi)) else c(estimate = b, lower = lo, upper = hi)
  list(data = d, model = model, result = c(estimate, p = p))
}

models <- list(
  VitaminD = fit_weighted(dat, Followup_25OHD_ng_mL ~ Rx + Baseline_25OHD_ng_mL + Age_years + Female + Clinical_followup_months),
  Pain = fit_weighted(dat, Followup_pain_intensity_VAS ~ Rx + Baseline_pain_intensity_VAS + Age_years + Female + Clinical_followup_months),
  PFO = fit_weighted(dat, Followup_CMO_mm ~ Rx + Baseline_CMO_mm + Age_years + Female + Clinical_followup_months),
  MUO = fit_weighted(dat, Followup_MMO_mm ~ Rx + Baseline_MMO_mm + Age_years + Female + Clinical_followup_months),
  Response2 = fit_weighted(dat |> filter(!is.na(Response_ge2)), Response_ge2 ~ Rx + Baseline_pain_intensity_VAS + Age_years + Female + Clinical_followup_months, quasipoisson(link = "log")),
  Response30 = fit_weighted(dat |> filter(!is.na(Response_ge30)), Response_ge30 ~ Rx + Baseline_pain_intensity_VAS + Age_years + Female + Clinical_followup_months, quasipoisson(link = "log"))
)

raw_values <- function(fit, outcome, binary = FALSE) {
  d <- fit$data
  fun <- if (binary) fmt_fraction_pct else fmt_mean_sd
  c(no_rx = fun(d[[outcome]][d$Rx == 0L]), rx = fun(d[[outcome]][d$Rx == 1L]))
}

make_row <- function(label, fit, outcome, measure, binary = FALSE) {
  raw <- raw_values(fit, outcome, binary)
  r <- fit$result
  data.frame(
    Outcome = label,
    `Low vitamin D / no Rx (n=51)` = raw[1],
    `Low vitamin D / Rx (n=149)` = raw[2],
    `Analysis N` = nrow(fit$data),
    `Effect measure` = measure,
    `Adjusted estimate (95% CI)` = format_effect(r, measure),
    `p-value` = fmt_p(r["p"]),
    check.names = FALSE
  )
}

section_row <- function(label) data.frame(
  Outcome = label, `Low vitamin D / no Rx (n=51)` = "", `Low vitamin D / Rx (n=149)` = "", `Analysis N` = "",
  `Effect measure` = "", `Adjusted estimate (95% CI)` = "", `p-value` = "", check.names = FALSE
)

table_4 <- bind_rows(
  section_row("Biochemical outcome"),
  make_row("Follow-up serum 25(OH)D, ng/mL", models$VitaminD, "Followup_25OHD_ng_mL", "MD"),
  section_row("Clinical outcomes"),
  make_row("Follow-up pain intensity, VAS", models$Pain, "Followup_pain_intensity_VAS", "MD"),
  make_row("Follow-up PFO, mm", models$PFO, "Followup_CMO_mm", "MD"),
  make_row("Follow-up MUO, mm", models$MUO, "Followup_MMO_mm", "MD"),
  make_row("VAS reduction ≥2 points", models$Response2, "Response_ge2", "RR", TRUE),
  make_row("VAS reduction ≥30%", models$Response30, "Response_ge30", "RR", TRUE)
)

audit <- bind_rows(lapply(names(models), function(nm) {
  r <- as.data.frame(t(models[[nm]]$result))
  r$outcome <- nm
  r$n <- nrow(models[[nm]]$data)
  r
}))
write_table_outputs(table_4, list(Weighted_models = audit), "Table_4_overlap_weighted_outcomes")
