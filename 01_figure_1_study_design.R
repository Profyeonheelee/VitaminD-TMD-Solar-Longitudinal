# Figure 1. Study design and longitudinal environmental-exposure framework
#
# Output:
#   Figure_1_study_design.png
#   Figure_1_study_design.tiff
#   Figure_1_study_design.pdf

# ============================================================
# 0. Packages
# ============================================================
required_pkgs <- c("ggplot2", "patchwork")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install these packages before running the script."
  )
}

library(ggplot2)
library(patchwork)

# ============================================================
# 1. File paths
# ============================================================
# Figure 1 uses the cohort counts reported in the manuscript and does not
# require the analysis workbook.
output_dir <- file.path("outputs", "Figure_1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Manuscript-verified cohort/group counts
# ============================================================
n_total <- 241L
n_g1 <- 41L
n_g2 <- 51L
n_g3 <- 149L
n_low <- 200L

message("Total N = ", n_total)
message("Group 1 = ", n_g1)
message("Group 2 = ", n_g2)
message("Group 3 = ", n_g3)
message("Low-vitamin-D cohort = ", n_low)

# Internal consistency checks
stopifnot(n_total == 241)
stopifnot(n_g1 == 41)
stopifnot(n_g2 == 51)
stopifnot(n_g3 == 149)
stopifnot(n_low == 200)

# ============================================================
# 3. Common figure theme
# ============================================================
base_family <- "Arial"

theme_framework <- theme_void(base_family = base_family) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0,
      margin = margin(b = 8)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 4. Panel A: Cohort flow and exposure groups
# ============================================================
panel_a <- ggplot() +

  # ----------------------------------------------------------
  # Parent cohort
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 3.25, xmax = 6.75,
    ymin = 7.35, ymax = 8.85,
    fill = "white",
    colour = "black",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 5.0, y = 8.37,
    label = paste0(
      "TMD patients with ≥2 serum 25(OH)D measurements\n",
      "N = ", n_total
    ),
    size = 4.0,
    fontface = "bold",
    family = base_family,
    lineheight = 1.05
  ) +
  annotate(
    "text",
    x = 5.0, y = 7.72,
    label = "Retrospective longitudinal cohort",
    size = 3.35,
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Branching connectors
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 5.0, xend = 5.0,
    y = 7.35, yend = 6.60,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 1.65, xend = 8.35,
    y = 6.60, yend = 6.60,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 1.65, xend = 1.65,
    y = 6.60, yend = 5.87,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 5.0, xend = 5.0,
    y = 6.60, yend = 5.87,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 8.35, xend = 8.35,
    y = 6.60, yend = 5.87,
    linewidth = 0.65,
    colour = "black"
  ) +

  # ----------------------------------------------------------
  # Group 1
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 0.25, xmax = 3.05,
    ymin = 3.15, ymax = 5.87,
    fill = "grey95",
    colour = "black",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = 1.65, y = 5.43,
    label = "Group 1",
    size = 3.8,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 1.65, y = 4.67,
    label = "Vitamin D sufficient\nBaseline 25(OH)D ≥30 ng/mL",
    size = 3.25,
    family = base_family,
    lineheight = 1.05
  ) +
  annotate(
    "text",
    x = 1.65, y = 3.88,
    label = "No vitamin D prescription",
    size = 3.15,
    family = base_family
  ) +
  annotate(
    "text",
    x = 1.65, y = 3.42,
    label = paste0("n = ", n_g1),
    size = 3.6,
    fontface = "bold",
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Group 2
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 3.60, xmax = 6.40,
    ymin = 3.15, ymax = 5.87,
    fill = "grey88",
    colour = "black",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = 5.0, y = 5.43,
    label = "Group 2",
    size = 3.8,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 5.0, y = 4.67,
    label = "Low vitamin D\nBaseline 25(OH)D <30 ng/mL",
    size = 3.25,
    family = base_family,
    lineheight = 1.05
  ) +
  annotate(
    "text",
    x = 5.0, y = 3.88,
    label = "No vitamin D prescription",
    size = 3.15,
    family = base_family
  ) +
  annotate(
    "text",
    x = 5.0, y = 3.42,
    label = paste0("n = ", n_g2),
    size = 3.6,
    fontface = "bold",
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Group 3
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 6.95, xmax = 9.75,
    ymin = 3.15, ymax = 5.87,
    fill = "grey78",
    colour = "black",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = 8.35, y = 5.43,
    label = "Group 3",
    size = 3.8,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 8.35, y = 4.67,
    label = "Low vitamin D\nBaseline 25(OH)D <30 ng/mL",
    size = 3.25,
    family = base_family,
    lineheight = 1.05
  ) +
  annotate(
    "text",
    x = 8.35, y = 3.88,
    label = "Vitamin D prescribed",
    size = 3.15,
    family = base_family
  ) +
  annotate(
    "text",
    x = 8.35, y = 3.42,
    label = paste0("n = ", n_g3),
    size = 3.6,
    fontface = "bold",
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Bracket: primary low-vitamin-D comparison
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 3.85, xend = 9.50,
    y = 2.82, yend = 2.82,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 3.85, xend = 3.85,
    y = 2.82, yend = 3.02,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 9.50, xend = 9.50,
    y = 2.82, yend = 3.02,
    linewidth = 0.65,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 6.675, y = 2.33,
    label = paste0(
      "Primary prescription-associated comparison: low-vitamin-D cohort (n = ",
      n_low, ")\n",
      "No Rx (n = ", n_g2, ") vs Rx (n = ", n_g3,
      "), propensity-score overlap weighting"
    ),
    size = 2.90,
    family = base_family,
    fontface = "bold",
    lineheight = 1.08
  ) +

  # ----------------------------------------------------------
  # Group-definition note
  # ----------------------------------------------------------
  annotate(
    "text",
    x = 5.0, y = 1.55,
    label = paste0(
      "Exposure groups were defined using baseline serum 25(OH)D status ",
      "and vitamin D prescription."
    ),
    size = 2.95,
    family = base_family
  ) +

  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(1.15, 9.15),
    clip = "off"
  ) +
  labs(
    title = "A  Study cohort and exposure groups"
  ) +
  theme_framework


# ============================================================
# 5. Panel B: Longitudinal assessments and environmental
#             exposure windows
# ============================================================
panel_b <- ggplot() +

  # ----------------------------------------------------------
  # Main time axis
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 0.0, xend = 10.0,
    y = 5.0, yend = 5.0,
    linewidth = 0.8,
    colour = "black",
    arrow = arrow(
      length = grid::unit(0.18, "cm"),
      type = "closed"
    )
  ) +

  # ----------------------------------------------------------
  # Baseline and follow-up time points
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 2.5, xend = 2.5,
    y = 4.75, yend = 5.25,
    linewidth = 0.9,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = 7.5, xend = 7.5,
    y = 4.75, yend = 5.25,
    linewidth = 0.9,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 2.5, y = 5.58,
    label = "Baseline",
    size = 3.8,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 7.5, y = 5.58,
    label = "Follow-up",
    size = 3.8,
    fontface = "bold",
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Baseline measurement box
  # Pain-location count intentionally removed.
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 1.22, xmax = 3.78,
    ymin = 2.78, ymax = 4.25,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2.5, y = 3.90,
    label = "Serum 25(OH)D",
    size = 3.2,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 2.5, y = 3.42,
    label = "VAS",
    size = 3.0,
    family = base_family
  ) +
  annotate(
    "text",
    x = 2.5, y = 3.02,
    label = "PFO • MUO • TMD features",
    size = 3.0,
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Follow-up measurement box
  # Pain-location count intentionally removed.
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 6.22, xmax = 8.78,
    ymin = 2.78, ymax = 4.25,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 7.5, y = 3.90,
    label = "Serum 25(OH)D",
    size = 3.2,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 7.5, y = 3.42,
    label = "VAS",
    size = 3.0,
    family = base_family
  ) +
  annotate(
    "text",
    x = 7.5, y = 3.02,
    label = "PFO • MUO",
    size = 3.0,
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Baseline environmental windows
  # 60-day = PRIMARY
  # 90-day = SENSITIVITY
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 0.55, xend = 2.5,
    y = 7.05, yend = 7.05,
    linewidth = 4.0,
    colour = "grey70",
    lineend = "butt"
  ) +
  annotate(
    "segment",
    x = 1.20, xend = 2.5,
    y = 6.42, yend = 6.42,
    linewidth = 4.0,
    colour = "grey35",
    lineend = "butt"
  ) +
  annotate(
    "text",
    x = 1.52, y = 7.48,
    label = "Pre-BL 90-day exposure\n(sensitivity)",
    size = 2.70,
    family = base_family,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 1.84, y = 5.98,
    label = "Pre-BL 60-day exposure (primary)",
    size = 2.75,
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Follow-up environmental windows
  # ----------------------------------------------------------
  annotate(
    "segment",
    x = 5.55, xend = 7.5,
    y = 7.05, yend = 7.05,
    linewidth = 4.0,
    colour = "grey70",
    lineend = "butt"
  ) +
  annotate(
    "segment",
    x = 6.20, xend = 7.5,
    y = 6.42, yend = 6.42,
    linewidth = 4.0,
    colour = "grey35",
    lineend = "butt"
  ) +
  annotate(
    "text",
    x = 6.52, y = 7.48,
    label = "Pre-FU 90-day exposure\n(sensitivity)",
    size = 2.70,
    family = base_family,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = 6.84, y = 5.98,
    label = "Pre-FU 60-day exposure (primary)",
    size = 2.75,
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Prescription STATUS only.
  # No duration bar, because prescription duration is excluded
  # from the final Paper 1 analysis.
  # ----------------------------------------------------------
  annotate(
    "rect",
    xmin = 3.55, xmax = 6.45,
    ymin = 1.35, ymax = 2.05,
    fill = "grey92",
    colour = "grey35",
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = 5.0, y = 1.70,
    label = "Vitamin D prescription status (Group 3)",
    size = 3.0,
    family = base_family,
    fontface = "bold"
  ) +

  # ----------------------------------------------------------
  # Environmental data source
  # ----------------------------------------------------------
  annotate(
    "text",
    x = 5.0, y = 8.13,
    label = paste0(
      "KMA Seoul ASOS station 108: daily solar radiation ",
      "and sunshine duration"
    ),
    size = 3.15,
    fontface = "bold",
    family = base_family
  ) +

  # ----------------------------------------------------------
  # Change-score definitions
  # ----------------------------------------------------------
  annotate(
    "text",
    x = 5.0, y = 0.50,
    label = paste0(
      "Longitudinal change: follow-up minus baseline; ",
      "pain reduction: baseline VAS minus follow-up VAS"
    ),
    size = 2.75,
    family = base_family
  ) +

  coord_cartesian(
    xlim = c(0, 10.35),
    ylim = c(0.15, 8.55),
    clip = "off"
  ) +
  labs(
    title = "B  Longitudinal assessments and environmental exposure windows"
  ) +
  theme_framework


# ============================================================
# 6. Combine panels
# ============================================================
# No overall title inside the figure.
figure1 <- panel_a / panel_b +
  plot_layout(
    heights = c(1.00, 1.05)
  )

# Preview in RStudio
print(figure1)


# ============================================================
# 7. Export
# ============================================================
png_file <- file.path(
  output_dir,
  "Figure_1_study_design.png"
)

tiff_file <- file.path(
  output_dir,
  "Figure_1_study_design.tiff"
)

pdf_file <- file.path(
  output_dir,
  "Figure_1_study_design.pdf"
)

ggsave(
  filename = png_file,
  plot = figure1,
  width = 8.3,
  height = 9.4,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure1,
  width = 8.3,
  height = 9.4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = figure1,
  width = 8.3,
  height = 9.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Saved:")
message(png_file)
message(tiff_file)
message(pdf_file)


# ============================================================
# 8. Updated manuscript-style Figure 1 legend
# ============================================================
figure_legend <- paste0(
  "Figure 1. Study design and longitudinal environmental-exposure framework. ",
  "(A) Participants with temporomandibular disorder and at least two serum ",
  "25-hydroxyvitamin D [25(OH)D] measurements were classified into three ",
  "groups according to baseline vitamin D status and vitamin D prescription: ",
  "vitamin D-sufficient without prescription (Group 1), low vitamin D without ",
  "prescription (Group 2), and low vitamin D with prescription (Group 3). ",
  "Groups 2 and 3 constituted the low-vitamin-D cohort used for the primary ",
  "prescription-associated comparison with propensity-score overlap weighting. ",
  "(B) Serum 25(OH)D, pain intensity, and mandibular opening were assessed at ",
  "baseline and follow-up; baseline TMD-related clinical features were also ",
  "recorded. Daily total solar radiation and sunshine duration from Seoul ASOS ",
  "station 108 of the Korea Meteorological Administration were linked to serum ",
  "sampling dates. The 60-day mean solar-radiation window was the primary ",
  "environmental exposure, whereas 90-day windows and sunshine duration were ",
  "used for sensitivity or comparison analyses. Vitamin D prescription was ",
  "treated as an observed clinical exposure status rather than a measure of ",
  "adherence or cumulative dose. PFO, pain-free opening; ",
  "KMA, Korea Meteorological Administration; MUO, maximum unassisted ",
  "opening; TMD, temporomandibular disorder; VAS, visual analog scale."
)

cat("\n\n", figure_legend, "\n")
