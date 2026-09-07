# VitaminD-TMD-Solar-Longitudinal
R code for a longitudinal study of vitamin D and clinical outcomes in temporomandibular disorders.
# Vitamin D, Solar Radiation, and Longitudinal TMD Outcomes

R code accompanying the study:
> **Vitamin D Prescription, Solar Radiation, and Longitudinal Pain Outcomes in Temporomandibular Disorders**

This repository contains the analysis code used to examine longitudinal
relationships among recorded vitamin D prescription, ambient solar radiation,
serum 25-hydroxyvitamin D [25(OH)D], present TMD-related pain intensity, and
mandibular opening in patients with temporomandibular disorders (TMD).

## Software

The analyses were developed for R 4.5.1. Required packages are:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "ggplot2", "patchwork",
  "sandwich", "openxlsx", "survey"
))
```

## Data

Place the analysis workbook at:

```text
data/Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx
```

The workbook must contain a worksheet named `Analysis_Data`. The analysis
dataset is not included in this repository.

## Scripts

### Tables

| Script | Output |
|---|---|
| `R/10_table_1_baseline_characteristics.R` | Table 1: baseline characteristics and follow-up timing |
| `R/11_table_2_unadjusted_longitudinal_outcomes.R` | Table 2: unadjusted longitudinal biochemical and clinical outcomes |
| `R/12_table_3_adjusted_clinical_outcomes.R` | Table 3: adjusted longitudinal clinical outcomes |
| `R/13_table_4_overlap_weighted_outcomes.R` | Table 4: overlap-weighted prescription-associated outcomes |
| `R/14_table_5_baseline_clinical_features.R` | Table 5: baseline serum 25(OH)D and TMD-related clinical features |
| `R/15_table_6_solar_radiation_associations.R` | Table 6: solar-radiation associations with biochemical and pain outcomes |
| `R/16_table_7_vitamin_d_change_and_pain.R` | Table 7: longitudinal change in serum 25(OH)D and pain improvement |

### Figures

| Script | Output |
|---|---|
| `R/01_figure_1_study_design.R` | Figure 1: study design and environmental-exposure framework |
| `R/02_figure_2_longitudinal_patterns.R` | Figure 2: longitudinal 25(OH)D patterns, solar radiation, and pain improvement |
| `R/03_figure_3_followup_duration.R` | Figure 3: follow-up duration and longitudinal changes |
| `R/04_figure_4_pain_correlates.R` | Figure 4: correlates of pain intensity |
| `R/05_supplementary_figure_s1_individual_trajectories.R` | Supplementary Figure S1: individual trajectories |
| `R/06_supplementary_figure_s2_vitamin_d_status.R` | Supplementary Figure S2: vitamin D status transitions |
| `R/07_supplementary_figure_s3_overlap_weighting.R` | Supplementary Figure S3: overlap-weighting diagnostics |
| `R/08_supplementary_figure_s4_clinical_features.R` | Supplementary Figure S4: baseline 25(OH)D and TMD-related clinical features |
| `R/09_supplementary_figure_s5_environmental_sensitivity.R` | Supplementary Figure S5: environmental-exposure sensitivity analyses |

Run scripts from the repository root. Each script checks its required columns
before analysis. Table scripts write manuscript-formatted Excel and CSV files,
together with analysis-audit sheets, to `outputs/tables/`. Figure scripts write
PNG, TIFF, PDF, and analysis-audit files to their own subdirectories under
`outputs/`.

To regenerate all tables, run:

```r
source("run_all_tables.R")
```

To regenerate all figures, run:

```r
source("run_all_figures.R")
```

To regenerate all tables and figures, run:

```r
source("run_all.R")
```

Figure 1 uses the manuscript-verified cohort counts and does not read the
patient-level workbook. Supplementary Figure S5 uses 500 bootstrap samples per
model and therefore takes longer than the other scripts.

## Terminology

**PFO** and **MUO** denote pain-free opening and maximum unassisted opening,
respectively, in accordance with DC/TMD terminology. **Self-reported bruxism**
denotes bruxism-related behavior assessed using a single patient-reported
yes/no item rather than a formal clinical or instrumental diagnosis.

## Statistical notes

The scripts preserve the specifications described in the manuscript,
including HC3 heteroskedasticity-robust standard errors, Holm correction for
prespecified pairwise comparisons, Benjamini–Hochberg false-discovery-rate
correction for exploratory analysis families, and outcome-specific complete
case analysis without imputation.

