# Table 1. Baseline characteristics and follow-up timing according to vitamin D
# status and prescription group

source(file.path("R", "table_utils.R"))

required_vars <- c(
  "Study_ID", "Sex", "Age_years", "Symptom_duration_months",
  "Baseline_25OHD_ng_mL", "Three_group_code", "VitD_draw_season_BL",
  "SolarRad_preBL_60d_mean_MJm2", "Baseline_pain_intensity_VAS",
  "Baseline_CMO_mm", "Baseline_MMO_mm", "Bruxism", "Clenching",
  "TMJ_noise", "Objective_locking", "Stiffness", "Trauma_history",
  "Clinical_followup_months"
)

dat <- read_analysis_data(required_vars) |>
  mutate(
    Group = make_group(Three_group_code),
    Female = to_female01(Sex),
    Deficient = as.integer(Baseline_25OHD_ng_mL < 20),
    Self_reported_bruxism = to_binary01(Bruxism, "Bruxism"),
    Clenching_bin = to_binary01(Clenching, "Clenching"),
    TMJ_noise_bin = to_binary01(TMJ_noise, "TMJ_noise"),
    Locking_bin = to_binary01(Objective_locking, "Objective_locking"),
    Stiffness_bin = to_binary01(Stiffness, "Stiffness"),
    Trauma_bin = to_binary01(Trauma_history, "Trauma_history")
  )

stopifnot(nrow(dat) == 241L, table(dat$Group) == c(41L, 51L, 149L))

continuous_spec <- list(
  Age = c("Age, years", "Age_years", "welch"),
  Symptom = c("Symptom duration, months", "Symptom_duration_months", "kruskal"),
  VitaminD = c("Baseline 25(OH)D, ng/mL", "Baseline_25OHD_ng_mL", "welch"),
  Solar = c("Pre-BL 60-day solar radiation, MJ/m²/day", "SolarRad_preBL_60d_mean_MJm2", "welch"),
  Pain = c("Baseline present TMD-related pain intensity, VAS", "Baseline_pain_intensity_VAS", "welch"),
  PFO = c("Baseline PFO, mm", "Baseline_CMO_mm", "welch"),
  MUO = c("Baseline MUO, mm", "Baseline_MMO_mm", "welch"),
  Followup = c("Follow-up interval, months", "Clinical_followup_months", "kruskal")
)

continuous_rows <- lapply(continuous_spec, function(spec) {
  label <- spec[1]
  variable <- spec[2]
  method <- spec[3]
  p <- switch(
    method,
    welch = oneway.test(dat[[variable]] ~ dat$Group, var.equal = FALSE)$p.value,
    kruskal = kruskal.test(dat[[variable]] ~ dat$Group)$p.value,
    descriptive = NA_real_
  )
  data.frame(
    Characteristic = label,
    Total = fmt_mean_sd(dat[[variable]]),
    `Group 1` = group_summary(dat, variable)[1],
    `Group 2` = group_summary(dat, variable)[2],
    `Group 3` = group_summary(dat, variable)[3],
    `Overall p-value` = fmt_p(p),
    `Post hoc comparison` = "—",
    check.names = FALSE
  )
})
names(continuous_rows) <- names(continuous_spec)

categorical_spec <- list(
  Female = c("Female sex", "Female", "pearson"),
  Deficient = c("Vitamin D deficiency (<20 ng/mL)", "Deficient", "pearson"),
  Bruxism = c("Self-reported bruxism", "Self_reported_bruxism", "pearson"),
  Clenching = c("Clenching", "Clenching_bin", "pearson"),
  Noise = c("Clinically detected TMJ noise", "TMJ_noise_bin", "pearson"),
  Locking = c("Clinically detected jaw locking", "Locking_bin", "pearson"),
  Stiffness = c("Jaw stiffness", "Stiffness_bin", "pearson"),
  Trauma = c("Trauma history", "Trauma_bin", "fisher")
)

categorical_rows <- lapply(categorical_spec, function(spec) {
  label <- spec[1]
  variable <- spec[2]
  method <- spec[3]
  p <- if (method == "descriptive") NA_real_ else as.numeric(safe_categorical_test(dat[[variable]], dat$Group)["p"])
  data.frame(
    Characteristic = label,
    Total = fmt_n_pct(dat[[variable]]),
    `Group 1` = fmt_n_pct(dat[[variable]][dat$Group == "Group 1"]),
    `Group 2` = fmt_n_pct(dat[[variable]][dat$Group == "Group 2"]),
    `Group 3` = fmt_n_pct(dat[[variable]][dat$Group == "Group 3"]),
    `Overall p-value` = fmt_p(p),
    `Post hoc comparison` = "—",
    check.names = FALSE
  )
})
names(categorical_rows) <- names(categorical_spec)

season_levels <- c("Winter", "Spring", "Summer", "Autumn")
season <- factor(dat$VitD_draw_season_BL, levels = season_levels)
season_p <- safe_categorical_test(season, dat$Group)["p"] |> as.numeric()
season_header <- data.frame(
  Characteristic = "Blood-draw season", Total = "", `Group 1` = "", `Group 2` = "", `Group 3` = "",
  `Overall p-value` = fmt_p(season_p), `Post hoc comparison` = "—", check.names = FALSE
)
season_rows <- lapply(season_levels, function(level) {
  indicator <- as.integer(season == level)
  data.frame(
    Characteristic = paste0("  ", level), Total = fmt_n_pct(indicator),
    `Group 1` = fmt_n_pct(indicator[dat$Group == "Group 1"]),
    `Group 2` = fmt_n_pct(indicator[dat$Group == "Group 2"]),
    `Group 3` = fmt_n_pct(indicator[dat$Group == "Group 3"]),
    `Overall p-value` = "", `Post hoc comparison` = "", check.names = FALSE
  )
})

posthoc <- list(
  Age = games_howell(dat$Age_years, dat$Group),
  Symptom = dunn_holm(dat$Symptom_duration_months, dat$Group),
  VitaminD = games_howell(dat$Baseline_25OHD_ng_mL, dat$Group),
  Solar = games_howell(dat$SolarRad_preBL_60d_mean_MJm2, dat$Group),
  Pain = games_howell(dat$Baseline_pain_intensity_VAS, dat$Group),
  PFO = games_howell(dat$Baseline_CMO_mm, dat$Group),
  MUO = games_howell(dat$Baseline_MMO_mm, dat$Group),
  Followup = dunn_holm(dat$Clinical_followup_months, dat$Group),
  Female = pairwise_categorical(dat$Female, dat$Group),
  Deficient = pairwise_categorical(dat$Deficient, dat$Group),
  Bruxism = pairwise_categorical(dat$Self_reported_bruxism, dat$Group),
  Clenching = pairwise_categorical(dat$Clenching_bin, dat$Group),
  Noise = pairwise_categorical(dat$TMJ_noise_bin, dat$Group),
  Locking = pairwise_categorical(dat$Locking_bin, dat$Group),
  Stiffness = pairwise_categorical(dat$Stiffness_bin, dat$Group),
  Trauma = pairwise_categorical(dat$Trauma_bin, dat$Group)
)

format_posthoc <- function(x, p_column) {
  significant <- x[x[[p_column]] < 0.05, , drop = FALSE]
  if (nrow(significant) == 0L) return("—")
  paste0(significant$comparison, ", p_adj", ifelse(significant[[p_column]] < 0.001, "<0.001", paste0("=", fmt_p(significant[[p_column]]))), collapse = "; ")
}

for (nm in names(posthoc)) {
  p_column <- if ("p_adjusted" %in% names(posthoc[[nm]])) "p_adjusted" else "p_holm"
  if (nm %in% names(continuous_rows)) continuous_rows[[nm]][["Post hoc comparison"]] <- format_posthoc(posthoc[[nm]], p_column)
  if (nm %in% names(categorical_rows)) categorical_rows[[nm]][["Post hoc comparison"]] <- format_posthoc(posthoc[[nm]], p_column)
}

table_1 <- bind_rows(
  continuous_rows$Age,
  categorical_rows$Female,
  continuous_rows$Symptom,
  continuous_rows$VitaminD,
  categorical_rows$Deficient,
  season_header,
  bind_rows(season_rows),
  continuous_rows$Solar,
  continuous_rows$Pain,
  continuous_rows$PFO,
  continuous_rows$MUO,
  categorical_rows$Bruxism,
  categorical_rows$Clenching,
  categorical_rows$Noise,
  categorical_rows$Locking,
  categorical_rows$Stiffness,
  categorical_rows$Trauma,
  continuous_rows$Followup
)

names(table_1)[2:5] <- c("Total (N=241)", "Group 1 (n=41)", "Group 2 (n=51)", "Group 3 (n=149)")
audit <- lapply(names(posthoc), function(nm) transform(posthoc[[nm]], variable = nm)) |> bind_rows()
write_table_outputs(table_1, list(Post_hoc_tests = audit), "Table_1_baseline_characteristics")
