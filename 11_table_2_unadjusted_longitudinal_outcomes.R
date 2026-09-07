# Table 2. Unadjusted longitudinal biochemical and clinical outcomes according
# to vitamin D status and prescription group

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Three_group_code", "Baseline_25OHD_ng_mL", "Followup_25OHD_ng_mL",
  "Delta_25OHD_ng_mL", "Baseline_pain_intensity_VAS", "Followup_pain_intensity_VAS",
  "Pain_intensity_reduction_VAS", "Baseline_CMO_mm", "Followup_CMO_mm",
  "Change_CMO_mm", "Baseline_MMO_mm", "Followup_MMO_mm", "Change_MMO_mm"
)

dat <- read_analysis_data(required_vars) |> mutate(Group = make_group(Three_group_code))
stopifnot(nrow(dat) == 241L)

continuous_test <- function(data, outcome) {
  model <- lm(reformulate("Group", outcome), data = data)
  terms <- grep("^Group", names(coef(model)), value = TRUE)
  list(model = model, global_p = robust_wald(model, terms), contrasts = group_contrasts(model))
}

binary_group_test <- function(data, outcome) {
  test <- safe_categorical_test(data[[outcome]], data$Group)
  list(global_p = as.numeric(test["p"]), contrasts = pairwise_categorical(data[[outcome]], data$Group))
}

make_continuous_row <- function(label, variable, data = dat, signed = FALSE, test = NULL) {
  summaries <- vapply(levels(data$Group), function(g) fmt_mean_sd(data[[variable]][data$Group == g], signed = signed), character(1))
  data.frame(
    Outcome = label, `Group 1 (n=41)` = summaries[1], `Group 2 (n=51)` = summaries[2],
    `Group 3 (n=149)` = summaries[3], `Unadjusted p-value` = if (is.null(test)) "—" else fmt_p(test$global_p),
    `Post hoc comparison` = "—", check.names = FALSE
  )
}

make_binary_row <- function(label, outcome, eligible = rep(TRUE, nrow(dat)), modified_poisson = FALSE, show_denominator = FALSE) {
  d <- dat[eligible & !is.na(dat[[outcome]]), ]
  if (modified_poisson) {
    model <- glm(reformulate("Group", outcome), family = poisson(link = "log"), data = d)
    vc <- hc3_vcov(model)
    terms <- grep("^Group", names(coef(model)), value = TRUE)
    test <- list(global_p = robust_wald(model, terms, vc), contrasts = group_contrasts(model, TRUE))
  } else {
    test <- binary_group_test(d, outcome)
  }
  summary_fun <- if (show_denominator) fmt_fraction_pct else fmt_n_pct
  summaries <- vapply(levels(dat$Group), function(g) summary_fun(d[[outcome]][d$Group == g]), character(1))
  sig <- test$contrasts[test$contrasts$p_holm < 0.05, , drop = FALSE]
  comparison <- if ("comparison" %in% names(sig)) sig$comparison else rownames(sig)
  post <- if (nrow(sig) == 0L) "—" else paste0(comparison, ", Holm-adjusted p", ifelse(sig$p_holm < 0.001, "<0.001", paste0("=", fmt_p(sig$p_holm))), collapse = "; ")
  data.frame(
    Outcome = label, `Group 1 (n=41)` = summaries[1], `Group 2 (n=51)` = summaries[2],
    `Group 3 (n=149)` = summaries[3], `Unadjusted p-value` = fmt_p(test$global_p),
    `Post hoc comparison` = post, check.names = FALSE
  )
}

fu_vitd <- continuous_test(dat, "Followup_25OHD_ng_mL")
fu_pain <- continuous_test(dat, "Followup_pain_intensity_VAS")
pfo_dat <- dat |> filter(complete.cases(Baseline_CMO_mm, Followup_CMO_mm))
muo_dat <- dat |> filter(complete.cases(Baseline_MMO_mm, Followup_MMO_mm))
fu_pfo <- continuous_test(pfo_dat, "Followup_CMO_mm")
fu_muo <- continuous_test(muo_dat, "Followup_MMO_mm")

add_posthoc <- function(row, result) {
  sig <- result$contrasts[result$contrasts[, "p_holm"] < 0.05, , drop = FALSE]
  row[["Post hoc comparison"]] <- if (nrow(sig) == 0L) "—" else paste(
    vapply(seq_len(nrow(sig)), function(i) format_effect(sig[i, ], "MD", sig[i, "p_holm"]), character(1)),
    collapse = "; "
  )
  row
}

header_row <- function(label) data.frame(
  Outcome = label, `Group 1 (n=41)` = "", `Group 2 (n=51)` = "", `Group 3 (n=149)` = "",
  `Unadjusted p-value` = "", `Post hoc comparison` = "", check.names = FALSE
)

vitd_low20 <- as.integer(dat$Followup_25OHD_ng_mL >= 20)
vitd_low30 <- as.integer(dat$Followup_25OHD_ng_mL >= 30)
vitd_inc5 <- as.integer(dat$Delta_25OHD_ng_mL >= 5)
vitd_inc10 <- as.integer(dat$Delta_25OHD_ng_mL >= 10)
dat <- dat |> mutate(FU_ge20 = vitd_low20, FU_ge30 = vitd_low30, Increase_ge5 = vitd_inc5, Increase_ge10 = vitd_inc10,
                     Response_ge2 = ifelse(Baseline_pain_intensity_VAS >= 2, as.integer(Pain_intensity_reduction_VAS >= 2), NA_integer_),
                     Response_ge30 = ifelse(Baseline_pain_intensity_VAS > 0, as.integer(Pain_intensity_reduction_VAS / Baseline_pain_intensity_VAS >= 0.30), NA_integer_))

table_2 <- bind_rows(
  header_row("Serum 25(OH)D, ng/mL"),
  make_continuous_row("  Baseline", "Baseline_25OHD_ng_mL"),
  add_posthoc(make_continuous_row("  Follow-up", "Followup_25OHD_ng_mL", test = fu_vitd), fu_vitd),
  make_continuous_row("  Δ25(OH)D", "Delta_25OHD_ng_mL", signed = TRUE),
  make_binary_row("Follow-up 25(OH)D ≥20 ng/mL", "FU_ge20"),
  make_binary_row("Follow-up 25(OH)D ≥30 ng/mL", "FU_ge30"),
  make_binary_row("Increase in 25(OH)D ≥5 ng/mL", "Increase_ge5"),
  make_binary_row("Increase in 25(OH)D ≥10 ng/mL", "Increase_ge10"),
  header_row("Present TMD-related pain intensity, VAS"),
  make_continuous_row("  Baseline", "Baseline_pain_intensity_VAS"),
  add_posthoc(make_continuous_row("  Follow-up", "Followup_pain_intensity_VAS", test = fu_pain), fu_pain),
  make_continuous_row("Pain reduction, VAS points", "Pain_intensity_reduction_VAS", signed = TRUE),
  make_binary_row("VAS reduction ≥2 points", "Response_ge2", modified_poisson = TRUE, show_denominator = TRUE),
  make_binary_row("VAS reduction ≥30%", "Response_ge30", modified_poisson = TRUE, show_denominator = TRUE),
  header_row("PFO, mm"),
  make_continuous_row("  Baseline", "Baseline_CMO_mm", pfo_dat),
  add_posthoc(make_continuous_row("  Follow-up", "Followup_CMO_mm", pfo_dat, test = fu_pfo), fu_pfo),
  make_continuous_row("  Change in PFO", "Change_CMO_mm", pfo_dat, signed = TRUE),
  header_row("MUO, mm"),
  make_continuous_row("  Baseline", "Baseline_MMO_mm", muo_dat),
  add_posthoc(make_continuous_row("  Follow-up", "Followup_MMO_mm", muo_dat, test = fu_muo), fu_muo),
  make_continuous_row("  Change in MUO", "Change_MMO_mm", muo_dat, signed = TRUE)
)

audit <- bind_rows(
  transform(as.data.frame(fu_vitd$contrasts), outcome = "Follow-up 25(OH)D"),
  transform(as.data.frame(fu_pain$contrasts), outcome = "Follow-up VAS"),
  transform(as.data.frame(fu_pfo$contrasts), outcome = "Follow-up PFO"),
  transform(as.data.frame(fu_muo$contrasts), outcome = "Follow-up MUO")
)
write_table_outputs(table_2, list(Continuous_contrasts = audit), "Table_2_unadjusted_longitudinal_outcomes")
