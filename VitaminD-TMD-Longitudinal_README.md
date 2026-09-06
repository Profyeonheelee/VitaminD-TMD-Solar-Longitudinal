# Vitamin D, Solar Radiation, and Longitudinal TMD Outcomes

R code accompanying the study:

> **Vitamin D Prescription, Solar Radiation, and Longitudinal Pain Outcomes in Temporomandibular Disorders**

This repository contains the analysis and figure-generation code used to examine longitudinal relationships among recorded vitamin D prescription, ambient solar radiation, serum 25-hydroxyvitamin D [25(OH)D], present TMD-related pain intensity, and mandibular opening in patients with temporomandibular disorders (TMD).

## Study overview

The retrospective longitudinal cohort included 241 patients with TMD and repeated serum 25(OH)D measurements. Participants were classified according to baseline vitamin D status and recorded vitamin D prescription:

- **Group 1:** vitamin D-sufficient, no prescription (`n = 41`)
- **Group 2:** low vitamin D, no prescription (`n = 51`)
- **Group 3:** low vitamin D, prescribed vitamin D (`n = 149`)

The primary prescription-associated analysis compared Groups 2 and 3 using propensity-score overlap weighting. Additional analyses evaluated longitudinal biochemical and clinical outcomes, date-linked ambient solar radiation, vitamin D status transitions, and clinical correlates of pain intensity.

The principal finding was a clear prescription-associated biochemical signal without corresponding clear evidence of improvement in TMD-related pain or mandibular opening. Ambient solar radiation was consistently associated with serum 25(OH)D but not with pain.

## Repository contents

```text
VitaminD-TMD-Longitudinal/
├── R/
│   ├── 01_figure_1_study_design.R
│   ├── 02_figure_2_longitudinal_patterns.R
│   ├── 03_figure_3_followup_duration.R
│   ├── 04_figure_4_pain_correlates.R
│   ├── 05_supplementary_figure_s1_individual_trajectories.R
│   ├── 06_supplementary_figure_s2_vitamin_d_status.R
│   ├── 07_supplementary_figure_s3_overlap_weighting.R
│   ├── 08_supplementary_figure_s4_clinical_features.R
│   └── 09_supplementary_figure_s5_environmental_sensitivity.R
├── data/
│   └── README.md
├── outputs/
│   └── .gitkeep
├── run_all_figures.R
└── README.md
```

| Script | Output |
|---|---|
| `01_figure_1_study_design.R` | Study design and longitudinal environmental-exposure framework |
| `02_figure_2_longitudinal_patterns.R` | Longitudinal serum 25(OH)D patterns, ambient solar radiation, and pain improvement |
| `03_figure_3_followup_duration.R` | Follow-up duration and longitudinal changes among prescribed patients |
| `04_figure_4_pain_correlates.R` | Clinical, biochemical, and environmental correlates of pain intensity |
| `05_supplementary_figure_s1_individual_trajectories.R` | Individual serum 25(OH)D and pain trajectories |
| `06_supplementary_figure_s2_vitamin_d_status.R` | Vitamin D status transitions |
| `07_supplementary_figure_s3_overlap_weighting.R` | Propensity-score distribution and covariate balance |
| `08_supplementary_figure_s4_clinical_features.R` | Baseline serum 25(OH)D and TMD-related clinical features |
| `09_supplementary_figure_s5_environmental_sensitivity.R` | Sensitivity analyses of environmental-exposure windows |

## Software requirements

The analyses were developed using **R 4.5.1**. Required packages are:

```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "sandwich",
  "openxlsx"
))
```

## Data requirements

The analysis dataset is not included in this repository. To run the code, place the analysis workbook at:

```text
data/Paper1_VitaminD_TMD_Longitudinal_KMA_Solar_Merged.xlsx
```

The workbook must contain a worksheet named `Analysis_Data`. Each script checks the variables required for its analysis before proceeding.

Meteorological exposures were derived from publicly available Korea Meteorological Administration data for Seoul Automated Synoptic Observing System station 108. The primary environmental exposure was mean daily total solar radiation during the 60 days preceding each serum sampling date; 90-day windows and sunshine-duration measures were evaluated in sensitivity analyses.

## Running the code

Clone or download the repository, open R in the repository root, and install the required packages. To generate every figure in manuscript order, run:

```r
source("run_all_figures.R")
```

Individual figures can be generated separately. For example:

```r
source("R/02_figure_2_longitudinal_patterns.R")
```

Each script writes publication-resolution PNG, TIFF, and PDF files, together with the corresponding statistical audit workbook where applicable, to a dedicated directory under `outputs/`.

Figure 1 uses manuscript-verified cohort counts and does not require access to patient-level data. Supplementary Figure S5 uses bootstrap resampling and may take longer to complete than the other scripts.

## Terminology and variable mapping

- **VAS** refers to present TMD-related pain intensity recorded using the clinic's 0–10 visual analog scale.
- **PFO** refers to pain-free opening, following DC/TMD terminology.
- **MUO** refers to maximum unassisted opening, following DC/TMD terminology.
- **Self-reported bruxism** refers to bruxism-related behavior assessed using a single patient-reported yes/no item and not to a formal clinical or instrumental diagnosis.

## Statistical framework

The scripts reproduce the specifications described in the manuscript, including:

- propensity-score overlap weighting for the primary comparison within the baseline-low-vitamin-D cohort;
- HC3 heteroskedasticity-robust standard errors for adjusted regression models;
- modified Poisson regression with robust standard errors for binary pain-response outcomes;
- Holm correction for prespecified pairwise group comparisons;
- Benjamini-Hochberg false-discovery-rate correction for exploratory analysis families; and
- outcome-specific complete-case analysis without imputation.

Longitudinal change in serum 25(OH)D and environmental exposure was calculated as follow-up minus baseline. Pain reduction was calculated as baseline minus follow-up VAS, so positive values indicate improvement.

## Reproducibility notes

This repository is organized to preserve a direct correspondence between the manuscript figures and R scripts.

## Data and code availability

Meteorological data are publicly available through the [Korea Meteorological Administration Open MET Data Portal](https://data.kma.go.kr/data/grnd/selectAsosRltmList.do?pgmNo=36).

Clinical and laboratory data may be made available by the corresponding author upon reasonable request, subject to institutional and ethical requirements.

## Citation

If you use this code, please cite the associated article. Full citation details will be added following publication.

## Contact

For questions regarding the code or study, please contact the corresponding author through the contact information provided in the article.
