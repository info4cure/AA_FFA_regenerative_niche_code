
# ==============================================================================
# FIGURE 03
# BULGE/ISTHMUS REGENERATIVE-NICHE ARCHITECTURE
#
# A  Patient-level Bulge abundance
# B  Patient-level functional phenotype
# C  K3 functional territories in Bulge UMAP
# D  Patient-level disease occupancy of K3 territories
# E  Molecular/functional phenotype of K3 territories
#
# IMPORTANT
# - Patient is the statistical replicate.
# - Frozen statistics are reused.
# - No new hypothesis testing is performed.
# - Primary display threshold: FDR < 0.05.
# - HC = grey; FFA = green; AA = red.
# - SVG + PDF + PNG are generated.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "scales",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(svglite)
})


# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------

ROOT <- normalizePath(
  ".",
  mustWork = TRUE
)

BASE <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/",
  "22_MULTICOMPARTMENT_DISEASE_AXES/",
  "22F_BULGE_STEM_NICHE_INTEGRITY"
)

F1_METRICS <- file.path(
  BASE,
  "22F1_CELL_ARCHITECTURE/tables/",
  "22F1_Bulge_patient_metrics.csv"
)

F1_STATS <- file.path(
  BASE,
  "22F1_CELL_ARCHITECTURE/tables/",
  "22F1_Bulge_patient_statistics.csv"
)

F2_PATIENT <- file.path(
  BASE,
  "22F2_FUNCTIONAL_STATE_DENSITY/",
  "22F2A_SCORE_FREEZE/tables/",
  "22F2A_Bulge_patient_functional_scores.csv"
)

F2_GLOBAL <- file.path(
  BASE,
  "22F2_FUNCTIONAL_STATE_DENSITY/",
  "22F2A_SCORE_FREEZE/tables/",
  "22F2A_Bulge_global_Kruskal_Wallis.csv"
)

F2_PAIRWISE <- file.path(
  BASE,
  "22F2_FUNCTIONAL_STATE_DENSITY/",
  "22F2A_SCORE_FREEZE/tables/",
  "22F2A_Bulge_pairwise_Wilcoxon.csv"
)

F3_LEDGER <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3A3_KMEANS/tables/",
  "22F3A3_Bulge_spatial_territory_ledger.csv"
)

F3_COMPOSITION <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F2_PATIENT_DA/tables/",
  "22F3F2_patient_K3_complete_composition.csv"
)

F3_GLOBAL <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F2_PATIENT_DA/tables/",
  "22F3F2_K3_global_KW.csv"
)

F3_PAIRWISE <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F2_PATIENT_DA/tables/",
  "22F3F2_K3_pairwise_Wilcoxon.csv"
)

F3_EFFECT <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F2_PATIENT_DA/tables/",
  "22F3F2_K3_effect_summary.csv"
)

F3_PHENOTYPE <- file.path(
  BASE,
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F3_INTEGRATED_PHENOTYPE/tables/",
  "22F3F3_global_territory_phenotype.csv"
)


# ------------------------------------------------------------------------------
# Frozen Bulge trajectory geometry for disease-split occupancy panel
# ------------------------------------------------------------------------------

F4F5_BASE <- file.path(
  BASE,
  "22F4_BULGE_MONOCLE3_TRAJECTORIES/",
  "22F4F_PUBLICATION_FIGURES/",
  "22F4F5_TRAJECTORY_DYNAMIC_GENES/tables"
)

F4F5_CELLS <- file.path(
  F4F5_BASE,
  "22F4F5_cell_UMAP_trajectory_source.csv"
)

F4F5_NODES <- file.path(
  F4F5_BASE,
  "22F4F5_principal_graph_nodes.csv"
)

F4F5_EDGES <- file.path(
  F4F5_BASE,
  "22F4F5_principal_graph_edges.csv"
)

required_files <- c(
  F1_METRICS,
  F1_STATS,
  F2_PATIENT,
  F2_GLOBAL,
  F2_PAIRWISE,
  F3_LEDGER,
  F3_COMPOSITION,
  F3_GLOBAL,
  F3_PAIRWISE,
  F3_EFFECT,
  F3_PHENOTYPE,
  F4F5_CELLS,
  F4F5_NODES,
  F4F5_EDGES
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0) {
  stop(
    "Missing required inputs:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------------------------
# 2. OUTPUT
# ------------------------------------------------------------------------------

OUT_ROOT <- file.path(
  ROOT,
  "PIPELINE_V2/results/publication_figures/",
  "FIG03_BULGE_REGENERATIVE_NICHE"
)

PDF_DIR <- file.path(
  OUT_ROOT,
  "PDF"
)

SVG_DIR <- file.path(
  OUT_ROOT,
  "SVG"
)

PNG_DIR <- file.path(
  OUT_ROOT,
  "PNG"
)

SOURCE_DIR <- file.path(
  OUT_ROOT,
  "SOURCE_DATA"
)

AUDIT_DIR <- file.path(
  OUT_ROOT,
  "AUDIT"
)

for (x in c(
  PDF_DIR,
  SVG_DIR,
  PNG_DIR,
  SOURCE_DIR,
  AUDIT_DIR
)) {
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------------------------
# 3. PALETTES
# ------------------------------------------------------------------------------

COL_HC <- "#8F8D8E"
COL_FFA <- "#4D9E4F"
COL_AA <- "#A8142A"

DISEASE_COLS <- c(
  "HC" = COL_HC,
  "FFA" = COL_FFA,
  "AA" = COL_AA
)

DISEASE_ORDER <- c(
  "HC",
  "FFA",
  "AA"
)

# K3 colours are figure-specific functional-state colours.
# They deliberately do not reuse disease colours.

K3_COLS <- c(
  "T1" = "#6C757D",
  "T2" = "#D89B3C",
  "T3" = "#3973AC"
)

COL_LOW <- "#2166AC"
COL_MID <- "#F7F7F7"
COL_HIGH <- "#B2182B"

COL_DARK <- "#303030"
COL_GRID <- "#EAEAEA"


# ------------------------------------------------------------------------------
# 4. GENERAL THEME
# ------------------------------------------------------------------------------

theme_pub <- function(
  base_size = 8
) {

  theme_classic(
    base_size = base_size
  ) +
    theme(
      plot.background =
        element_rect(
          fill = "white",
          colour = NA
        ),

      panel.background =
        element_blank(),

      plot.title =
        element_text(
          face = "bold",
          size = base_size + 1.4,
          hjust = 0
        ),

      plot.subtitle =
        element_text(
          size = base_size,
          colour = "#444444",
          hjust = 0
        ),

      axis.title =
        element_text(
          size = base_size
        ),

      axis.text =
        element_text(
          size = base_size - 0.4,
          colour = COL_DARK
        ),

      legend.title =
        element_text(
          face = "bold",
          size = base_size - 0.2
        ),

      legend.text =
        element_text(
          size = base_size - 0.4
        ),

      strip.background =
        element_rect(
          fill = "#F4F4F4",
          colour = "#666666",
          linewidth = 0.3
        ),

      strip.text =
        element_text(
          face = "bold",
          size = base_size
        )
    )
}


save_all <- function(
  plot,
  stem,
  width,
  height
) {

  ggsave(
    file.path(
      SVG_DIR,
      paste0(
        stem,
        ".svg"
      )
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    device = svglite,
    bg = "white"
  )

  ggsave(
    file.path(
      PDF_DIR,
      paste0(
        stem,
        ".pdf"
      )
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  ggsave(
    file.path(
      PNG_DIR,
      paste0(
        stem,
        ".png"
      )
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white"
  )
}


# ------------------------------------------------------------------------------
# 5. READ
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "READ INPUTS"
)
message(
  "============================================================"
)

f1_metrics <- read_csv(
  F1_METRICS,
  show_col_types = FALSE
)

f1_stats <- read_csv(
  F1_STATS,
  show_col_types = FALSE
)

f2_patient <- read_csv(
  F2_PATIENT,
  show_col_types = FALSE
)

f2_global <- read_csv(
  F2_GLOBAL,
  show_col_types = FALSE
)

f2_pairwise <- read_csv(
  F2_PAIRWISE,
  show_col_types = FALSE
)

ledger <- read_csv(
  F3_LEDGER,
  show_col_types = FALSE
)

k3_comp <- read_csv(
  F3_COMPOSITION,
  show_col_types = FALSE
)

k3_global <- read_csv(
  F3_GLOBAL,
  show_col_types = FALSE
)

k3_pairwise <- read_csv(
  F3_PAIRWISE,
  show_col_types = FALSE
)

k3_pheno <- read_csv(
  F3_PHENOTYPE,
  show_col_types = FALSE
)


# ------------------------------------------------------------------------------
# 6. FORMAL AUDIT
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "FORMAL INPUT AUDIT"
)
message(
  "============================================================"
)


if (nrow(f1_metrics) != 16) {
  stop(
    "AUDIT FAIL - expected 16 patients in 22F1."
  )
}

if (nrow(f2_patient) != 16) {
  stop(
    "AUDIT FAIL - expected 16 patients in 22F2."
  )
}

if (nrow(ledger) != 5420) {
  stop(
    "AUDIT FAIL - expected 5,420 Bulge/Isthmus cells."
  )
}

if (
  !setequal(
    unique(
      ledger$disease
    ),
    DISEASE_ORDER
  )
) {
  stop(
    "AUDIT FAIL - disease labels differ from HC/FFA/AA."
  )
}

if (
  !setequal(
    unique(
      ledger$territory_K3
    ),
    c(
      "T1",
      "T2",
      "T3"
    )
  )
) {
  stop(
    "AUDIT FAIL - expected K3 territories T1/T2/T3."
  )
}

if (nrow(k3_comp) != 48) {
  stop(
    "AUDIT FAIL - expected 16 patients x 3 territories = 48 rows."
  )
}

if (nrow(k3_pheno) != 3) {
  stop(
    "AUDIT FAIL - expected three K3 territory phenotype rows."
  )
}


# Frozen key results.

stem_AA_FFA <- f2_pairwise |>
  filter(
    metric == "stem_niche_mean",
    group1 == "AA",
    group2 == "FFA"
  )

cd200_AA_FFA <- f2_pairwise |>
  filter(
    metric == "CD200_mean",
    group1 == "AA",
    group2 == "FFA"
  )

cycle_AA_FFA <- f2_pairwise |>
  filter(
    metric == "cycling_mean",
    group1 == "AA",
    group2 == "FFA"
  )

mhcii_FFA_HC <- f2_pairwise |>
  filter(
    metric == "MHCII_mean",
    group1 == "FFA",
    group2 == "HC"
  )

if (
  nrow(stem_AA_FFA) != 1 ||
  abs(stem_AA_FFA$FDR - 0.0294) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected AA-vs-FFA stem-niche result."
  )
}

if (
  nrow(cd200_AA_FFA) != 1 ||
  abs(cd200_AA_FFA$FDR - 0.0294) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected AA-vs-FFA CD200 result."
  )
}

if (
  nrow(cycle_AA_FFA) != 1 ||
  abs(cycle_AA_FFA$FDR - 0.0375) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected AA-vs-FFA cycling result."
  )
}

if (
  nrow(mhcii_FFA_HC) != 1 ||
  abs(mhcii_FFA_HC$FDR - 0.0294) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected FFA-vs-HC MHCII result."
  )
}


T2_global <- k3_global |>
  filter(
    territory == "T2"
  )

T3_global <- k3_global |>
  filter(
    territory == "T3"
  )

if (
  abs(
    T2_global$FDR - 0.0129
  ) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected T2 global FDR."
  )
}

if (
  abs(
    T3_global$FDR - 0.0129
  ) > 0.002
) {
  stop(
    "AUDIT FAIL - unexpected T3 global FDR."
  )
}


message(
  "PASS - 16 patients"
)

message(
  "PASS - 5,420 Bulge/Isthmus cells"
)

message(
  "PASS - K3 = T1/T2/T3"
)

message(
  "PASS - AA vs FFA stem-niche FDR = ",
  signif(
    stem_AA_FFA$FDR,
    3
  )
)

message(
  "PASS - AA vs FFA CD200 FDR = ",
  signif(
    cd200_AA_FFA$FDR,
    3
  )
)

message(
  "PASS - AA vs FFA cycling FDR = ",
  signif(
    cycle_AA_FFA$FDR,
    3
  )
)

message(
  "PASS - FFA vs HC MHCII FDR = ",
  signif(
    mhcii_FFA_HC$FDR,
    3
  )
)

message(
  "PASS - T2 global FDR = ",
  signif(
    T2_global$FDR,
    3
  )
)

message(
  "PASS - T3 global FDR = ",
  signif(
    T3_global$FDR,
    3
  )
)


# ------------------------------------------------------------------------------
# 7. HELPER FOR BRACKETS
# ------------------------------------------------------------------------------

add_bracket <- function(
  p,
  x1,
  x2,
  y,
  label,
  tip = NULL
) {

  if (is.null(tip)) {
    tip <- 0.025 * abs(y)
  }

  if (
    !is.finite(tip) ||
    tip == 0
  ) {
    tip <- 0.02
  }

  p +
    annotate(
      "segment",
      x = x1,
      xend = x2,
      y = y,
      yend = y,
      linewidth = 0.35
    ) +
    annotate(
      "segment",
      x = x1,
      xend = x1,
      y = y,
      yend = y - tip,
      linewidth = 0.35
    ) +
    annotate(
      "segment",
      x = x2,
      xend = x2,
      y = y,
      yend = y - tip,
      linewidth = 0.35
    ) +
    annotate(
      "text",
      x = mean(
        c(
          x1,
          x2
        )
      ),
      y = y + tip,
      label = label,
      size = 2.5
    )
}


# ------------------------------------------------------------------------------
# Divergent functional-score palette
# Shared publication convention:
# blue = relatively lower; white = intermediate; red = relatively higher
# ------------------------------------------------------------------------------

COL_DOWN <- "#2166AC"
COL_ZERO <- "#F7F7F7"
COL_UP   <- "#B2182B"


# ------------------------------------------------------------------------------
# 8. PANEL A - BULGE ABUNDANCE
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "PANEL A - BULGE ABUNDANCE"
)
message(
  "============================================================"
)

A <- f1_metrics |>
  mutate(
    disease =
      factor(
        disease,
        levels = DISEASE_ORDER
      )
  )

write_csv(
  A |>
    select(
      patient_id,
      disease,
      batch,
      n_bulge
    ),
  file.path(
    SOURCE_DIR,
    "Fig03A_Bulge_abundance.csv"
  )
)

pA <- ggplot(
  A,
  aes(
    x = disease,
    y = n_bulge,
    fill = disease
  )
) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    linewidth = 0.45,
    alpha = 0.72
  ) +
  geom_jitter(
    aes(
      colour = disease
    ),
    width = 0.10,
    height = 0,
    size = 2.2,
    alpha = 0.95,
    shape = 21,
    fill = "white",
    stroke = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = DISEASE_COLS,
    guide = "none"
  ) +
  scale_colour_manual(
    values = DISEASE_COLS,
    guide = "none"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.02,
        0.10
      )
    )
  ) +
  labs(
    title = "A  Bulge/Isthmus abundance",
    subtitle = "Patient-level Bulge/Isthmus cell recovery",
    x = NULL,
    y = "Bulge/Isthmus cells per patient"
  ) +
  theme_pub(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        face = "bold"
      ),

    legend.position =
      "none",

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )

save_all(
  pA,
  "Fig03A_Bulge_abundance",
  width = 3.8,
  height = 4.2
)


# ------------------------------------------------------------------------------
# 9. PANEL B - PATIENT-LEVEL BULGE FUNCTIONAL PHENOTYPE
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "PANEL B - BULGE FUNCTIONAL PHENOTYPE"
)
message(
  "============================================================"
)

B_metrics <- c(
  "stem_niche_mean",
  "CD200_mean",
  "cycling_mean",
  "MHCII_mean"
)

B_labels <- c(
  "stem_niche_mean" =
    "Stem-niche score",

  "CD200_mean" =
    "CD200 expression",

  "cycling_mean" =
    "Cycling score",

  "MHCII_mean" =
    "MHC-II score"
)

B <- f2_patient |>
  select(
    patient_id,
    disease,
    batch,
    all_of(
      B_metrics
    )
  ) |>
  pivot_longer(
    cols =
      all_of(
        B_metrics
      ),
    names_to =
      "metric",
    values_to =
      "value"
  ) |>
  mutate(
    disease =
      factor(
        disease,
        levels =
          DISEASE_ORDER
      ),

    metric =
      factor(
        metric,
        levels =
          B_metrics,
        labels =
          unname(
            B_labels[
              B_metrics
            ]
          )
      )
  )

write_csv(
  B,
  file.path(
    SOURCE_DIR,
    "Fig03B_Bulge_functional_scores.csv"
  )
)


# Frozen pairwise results used only for annotations.

B_stats <- f2_pairwise |>
  filter(
    metric %in%
      B_metrics,
    FDR < 0.05
  ) |>
  mutate(
    metric_label =
      unname(
        B_labels[
          metric
        ]
      )
  )

write_csv(
  B_stats,
  file.path(
    SOURCE_DIR,
    "Fig03B_frozen_pairwise_statistics.csv"
  )
)


pB <- ggplot(
  B,
  aes(
    x = disease,
    y = value,
    fill = disease
  )
) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    linewidth = 0.4,
    alpha = 0.68
  ) +
  geom_jitter(
    aes(
      colour = disease
    ),
    width = 0.10,
    height = 0,
    size = 1.75,
    alpha = 0.95,
    shape = 21,
    fill = "white",
    stroke = 0.55,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values =
      DISEASE_COLS,
    guide =
      "none"
  ) +
  scale_colour_manual(
    values =
      DISEASE_COLS,
    guide =
      "none"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0.05,
          0.30
        )
      )
  ) +
  labs(
    title = "A  Functional remodeling of the Bulge/Isthmus niche",
    subtitle = "Patient-level functional scores",
    x = NULL,
    y = NULL
  ) +
  theme_pub(
    base_size = 7.5
  ) +
  theme(
    axis.text.x =
      element_text(
        face = "bold",
        size = 6.8
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 7
      ),

    legend.position =
      "none",

    panel.spacing.x =
      unit(
        0.8,
        "lines"
      ),

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )


# ------------------------------------------------------------------------------
# Add only frozen FDR < 0.05 comparisons.
#
# Stem niche:
# AA > FFA, FDR = 0.0294
#
# CD200:
# AA > FFA, FDR = 0.0294
# FFA < HC, FDR = 0.0294
# AA > HC,  FDR = 0.0375
#
# Cycling:
# AA < FFA, FDR = 0.0375
#
# MHC-II:
# FFA > HC, FDR = 0.0294
# ------------------------------------------------------------------------------

B_ranges <- B |>
  group_by(
    metric
  ) |>
  summarise(
    ymin =
      min(
        value,
        na.rm = TRUE
      ),

    ymax =
      max(
        value,
        na.rm = TRUE
      ),

    yrange =
      ymax - ymin,

    .groups =
      "drop"
  ) |>
  mutate(
    yrange =
      if_else(
        yrange == 0,
        1,
        yrange
      )
  )

get_B_y <- function(
  metric_name,
  multiplier
) {

  rr <- B_ranges |>
    filter(
      as.character(
        metric
      ) ==
        B_labels[
          metric_name
        ]
    )

  rr$ymax +
    multiplier *
    rr$yrange
}


# We use geom segments/text explicitly because facet-specific brackets
# are easier to audit than adding significance packages.

B_brackets <- tibble::tribble(
  ~metric,             ~x1, ~x2, ~level, ~label,
  "Stem-niche score",     2,   3,   0.08, "FDR=0.029",
  "CD200 expression",     1,   2,   0.08, "FDR=0.029",
  "CD200 expression",     2,   3,   0.20, "FDR=0.029",
  "CD200 expression",     1,   3,   0.34, "FDR=0.038",
  "Cycling score",        2,   3,   0.08, "FDR=0.038",
  "MHC-II score",         1,   2,   0.08, "FDR=0.029"
) |>
  left_join(
    B_ranges,
    by =
      "metric"
  ) |>
  mutate(
    y =
      ymax +
      level *
      yrange,

    tip =
      0.025 *
      yrange,

    text_y =
      y +
      0.025 *
      yrange
  )

pB <- pB +
  geom_segment(
    data =
      B_brackets,
    aes(
      x = x1,
      xend = x2,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data =
      B_brackets,
    aes(
      x = x1,
      xend = x1,
      y = y,
      yend = y - tip
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data =
      B_brackets,
    aes(
      x = x2,
      xend = x2,
      y = y,
      yend = y - tip
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_text(
    data =
      B_brackets,
    aes(
      x =
        (
          x1 + x2
        ) / 2,
      y =
        text_y,
      label =
        label
    ),
    inherit.aes = FALSE,
    size = 2.15
  )

save_all(
  pB,
  "Fig03B_Bulge_functional_phenotype",
  width = 8.8,
  height = 4.2
)


# ------------------------------------------------------------------------------
# 10. PANEL C - DISEASE-SPECIFIC OCCUPANCY OF SHARED BULGE TRAJECTORY
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)

message(
  "PANEL C - DISEASE-SPECIFIC OCCUPANCY OF SHARED BULGE TRAJECTORY"
)

message(
  "============================================================"
)

trajectory_cells <- read_csv(
  F4F5_CELLS,
  show_col_types = FALSE
)

trajectory_nodes <- read_csv(
  F4F5_NODES,
  show_col_types = FALSE
)

trajectory_edges <- read_csv(
  F4F5_EDGES,
  show_col_types = FALSE
)


# ------------------------------------------------------------------------------
# Formal geometry audit
# ------------------------------------------------------------------------------

if (nrow(trajectory_cells) != 5420) {
  stop(
    "AUDIT FAIL - expected 5,420 cells in frozen trajectory source."
  )
}

if (n_distinct(trajectory_cells$patient_id) != 16) {
  stop(
    "AUDIT FAIL - expected 16 patients in frozen trajectory source."
  )
}

if (nrow(trajectory_nodes) != 40) {
  stop(
    "AUDIT FAIL - expected 40 frozen principal-graph nodes."
  )
}

if (nrow(trajectory_edges) != 39) {
  stop(
    "AUDIT FAIL - expected 39 frozen principal-graph edges."
  )
}

if (!setequal(
  unique(trajectory_cells$territory_K3),
  c("T1", "T2", "T3")
)) {
  stop(
    "AUDIT FAIL - frozen trajectory source does not contain exactly T1/T2/T3."
  )
}

if (!setequal(
  unique(trajectory_cells$disease),
  c("HC", "FFA", "AA")
)) {
  stop(
    "AUDIT FAIL - frozen trajectory source does not contain HC/FFA/AA."
  )
}


# Exact cell-level correspondence with the Bulge ledger.

ledger_cells <- if ("cell" %in% names(ledger)) {
  ledger$cell
} else if ("cell_id" %in% names(ledger)) {
  ledger$cell_id
} else {
  stop(
    "AUDIT FAIL - unable to identify cell ID column in Bulge ledger."
  )
}

if (!setequal(
  trajectory_cells$cell_id,
  ledger_cells
)) {
  stop(
    "AUDIT FAIL - F4F5 trajectory cells do not match the frozen 5,420-cell Bulge ledger."
  )
}


message(
  "PASS - 5,420 trajectory cells"
)

message(
  "PASS - 16 patients"
)

message(
  "PASS - 40 principal-graph nodes"
)

message(
  "PASS - 39 principal-graph edges"
)

message(
  "PASS - trajectory geometry matches frozen Bulge cell universe"
)


# ------------------------------------------------------------------------------
# Prepare disease-split cell data
# ------------------------------------------------------------------------------

C <- trajectory_cells |>
  mutate(
    disease =
      factor(
        disease,
        levels = DISEASE_ORDER
      ),

    territory_K3 =
      factor(
        territory_K3,
        levels = c(
          "T1",
          "T2",
          "T3"
        )
      )
  )


# Replicate the SAME principal graph across all three disease facets.

C_edges <- tidyr::crossing(
  disease =
    factor(
      DISEASE_ORDER,
      levels = DISEASE_ORDER
    ),
  trajectory_edges
)

C_nodes <- tidyr::crossing(
  disease =
    factor(
      DISEASE_ORDER,
      levels = DISEASE_ORDER
    ),
  trajectory_nodes
)


# ------------------------------------------------------------------------------
# Source data
# ------------------------------------------------------------------------------

write_csv(
  C,
  file.path(
    SOURCE_DIR,
    "Fig03C_disease_specific_trajectory_cells.csv"
  )
)

write_csv(
  C_edges,
  file.path(
    SOURCE_DIR,
    "Fig03C_shared_principal_graph_edges.csv"
  )
)

write_csv(
  C_nodes,
  file.path(
    SOURCE_DIR,
    "Fig03C_shared_principal_graph_nodes.csv"
  )
)


# ------------------------------------------------------------------------------
# Plot
#
# Important:
#   graph first
#   cells second
#
# This preserves the shared topology while allowing disease-specific
# occupancy to remain visually dominant.
# ------------------------------------------------------------------------------

pC <- ggplot() +

  geom_segment(
    data = C_edges,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    colour = "#252525",
    linewidth = 0.48,
    alpha = 0.90,
    lineend = "round"
  ) +

  geom_point(
    data = C,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      fill = territory_K3
    ),
    shape = 21,
    size = 1.10,
    stroke = 0.14,
    colour = "#3A3A3A",
    alpha = 0.92
  ) +

  facet_wrap(
    ~ disease,
    nrow = 1,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = K3_COLS,
    breaks = c(
      "T3",
      "T2",
      "T1"
    ),
    name = "State"
  ) +

  coord_equal() +

  labs(
    title =
      "B  Disease-specific occupancy of the shared Bulge trajectory landscape",

    subtitle =
      "Identical UMAP geometry, K3 states and principal graph across groups",

    x =
      "UMAP 1",

    y =
      "UMAP 2"
  ) +

  theme_pub(
    base_size = 8
  ) +

  theme(
    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face = "bold",
        size = 9.2,
        colour = "#202020"
      ),

    legend.position =
      "right",

    legend.title =
      element_text(
        size = 7.8
      ),

    legend.text =
      element_text(
        size = 7.2
      ),

    axis.title =
      element_text(
        size = 7.6
      ),

    axis.text =
      element_text(
        size = 6.7
      ),

    panel.spacing.x =
      unit(
        0.55,
        "lines"
      ),

    plot.margin =
      margin(
        5,
        7,
        5,
        5
      )
  )


save_all(
  pC,
  "Fig03C_disease_specific_shared_Bulge_trajectory",
  width = 11.2,
  height = 4.7
)


# ------------------------------------------------------------------------------
# 11. PANEL D - DISEASE OCCUPANCY OF K3 TERRITORIES
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "PANEL D - K3 DISEASE OCCUPANCY"
)
message(
  "============================================================"
)

D <- k3_comp |>
  mutate(
    disease =
      factor(
        disease,
        levels =
          DISEASE_ORDER
      ),

    territory =
      factor(
        territory,
        levels = c(
          "T1",
          "T2",
          "T3"
        )
      )
  )

write_csv(
  D,
  file.path(
    SOURCE_DIR,
    "Fig03D_K3_patient_occupancy.csv"
  )
)

D_stats <- k3_pairwise |>
  filter(
    FDR < 0.05
  )

write_csv(
  D_stats,
  file.path(
    SOURCE_DIR,
    "Fig03D_K3_pairwise_FDR005.csv"
  )
)

pD <- ggplot(
  D,
  aes(
    x = disease,
    y = pct_bulge,
    fill = disease
  )
) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    linewidth = 0.42,
    alpha = 0.70
  ) +
  geom_jitter(
    aes(
      colour = disease
    ),
    width = 0.10,
    size = 1.7,
    shape = 21,
    fill = "white",
    stroke = 0.55,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ territory,
    nrow = 1
  ) +
  scale_fill_manual(
    values =
      DISEASE_COLS,
    guide =
      "none"
  ) +
  scale_colour_manual(
    values =
      DISEASE_COLS,
    guide =
      "none"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      112
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "D  Disease-specific occupancy of Bulge territories",
    subtitle = "Percentage of each patient's Bulge/Isthmus cells",
    x = NULL,
    y = "Bulge/Isthmus cells (%)"
  ) +
  theme_pub(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        face = "bold"
      ),

    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "none",

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )


# Significant frozen comparisons:
#
# T2:
# FFA > HC  FDR = 0.0485
# AA  < FFA FDR = 0.0365
#
# T3:
# FFA < HC  FDR = 0.0485
# AA  > FFA FDR = 0.0365

D_brackets <- tibble::tribble(
  ~territory, ~x1, ~x2, ~y, ~label,
  "T2",          1,   2,  91, "FDR=0.049",
  "T2",          2,   3, 102, "FDR=0.037",
  "T3",          1,   2,  91, "FDR=0.049",
  "T3",          2,   3, 102, "FDR=0.037"
) |>
  mutate(
    territory =
      factor(
        territory,
        levels = c(
          "T1",
          "T2",
          "T3"
        )
      )
  )

pD <- pD +
  geom_segment(
    data =
      D_brackets,
    aes(
      x = x1,
      xend = x2,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data =
      D_brackets,
    aes(
      x = x1,
      xend = x1,
      y = y,
      yend = y - 2
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data =
      D_brackets,
    aes(
      x = x2,
      xend = x2,
      y = y,
      yend = y - 2
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_text(
    data =
      D_brackets,
    aes(
      x =
        (
          x1 + x2
        ) / 2,
      y =
        y + 2,
      label =
        label
    ),
    inherit.aes = FALSE,
    size = 2.15
  )

save_all(
  pD,
  "Fig03D_K3_disease_occupancy",
  width = 6.8,
  height = 4.8
)


# ------------------------------------------------------------------------------
# Frozen K3 territory functional phenotype
# ------------------------------------------------------------------------------

TERRITORY_PHENO_FILE <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/",
  "22_MULTICOMPARTMENT_DISEASE_AXES/",
  "22F_BULGE_STEM_NICHE_INTEGRITY/",
  "22F3_BULGE_SPATIAL_TERRITORIES/",
  "22F3F_INTEGRATED_TERRITORY_PHENOTYPE/",
  "22F3F3_INTEGRATED_PHENOTYPE/",
  "tables/",
  "22F3F3_global_territory_phenotype.csv"
)

if (!file.exists(TERRITORY_PHENO_FILE)) {
  stop(
    "AUDIT FAIL - missing frozen K3 territory phenotype file: ",
    TERRITORY_PHENO_FILE
  )
}


# ------------------------------------------------------------------------------
# 12. PANEL E - K3 FUNCTIONAL PHENOTYPE
#
# ------------------------------------------------------------------------------
# 12. PANEL E - K3 FUNCTIONAL PHENOTYPE
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "PANEL E - K3 FUNCTIONAL PHENOTYPE"
)
message(
  "============================================================"
)

# Load frozen K3 global functional phenotype.
# This table contains the precomputed functional characterization of
# territories T1-T3 and is reused without additional statistical testing.

if (!exists("territory_pheno")) {

  if (!exists("TERRITORY_PHENO_FILE")) {
    stop(
      "AUDIT FAIL - TERRITORY_PHENO_FILE is not defined."
    )
  }

  if (!file.exists(TERRITORY_PHENO_FILE)) {
    stop(
      "AUDIT FAIL - missing K3 territory phenotype file: ",
      TERRITORY_PHENO_FILE
    )
  }

  territory_pheno <- read_csv(
    TERRITORY_PHENO_FILE,
    show_col_types = FALSE
  )
}

required_E_cols <- c(
  "territory",
  "n_cells",
  "n_patients",
  "pct_stem_high",
  "stem_niche_mean",
  "cycling_mean",
  "MHCII_mean",
  "CD200_z_mean",
  "CD200_detected_pct"
)

missing_E_cols <- setdiff(
  required_E_cols,
  names(territory_pheno)
)

if (length(missing_E_cols) > 0) {
  stop(
    "AUDIT FAIL - missing columns in K3 territory phenotype table: ",
    paste(
      missing_E_cols,
      collapse = "; "
    )
  )
}

if (!setequal(
  territory_pheno$territory,
  c("T1", "T2", "T3")
)) {
  stop(
    "AUDIT FAIL - expected exactly K3 territories T1/T2/T3."
  )
}

message(
  "PASS - loaded frozen K3 functional phenotype: ",
  nrow(territory_pheno),
  " territories"
)

E <- territory_pheno |>
  select(
    territory,
    n_cells,
    n_patients,
    pct_stem_high,
    stem_niche_mean,
    cycling_mean,
    MHCII_mean,
    CD200_z_mean,
    CD200_detected_pct
  ) |>
  mutate(
    territory =
      factor(
        territory,
        levels = c(
          "T1",
          "T2",
          "T3"
        )
      )
  )

write_csv(
  E,
  file.path(
    SOURCE_DIR,
    "Fig03E_K3_functional_phenotype.csv"
  )
)


# ------------------------------------------------------------------------------
# For visualization, standardize continuous functional variables across the
# three territories. This allows quantities on different original scales to
# be displayed in a common heatmap without changing their relative ordering.
#
# The two percentage variables are not included simultaneously because
# pct_stem_high and CD200_detected_pct partly recapitulate the continuous
# stem-niche/CD200 axes. We retain the four mechanistically distinct axes:
#
#   Stem-niche programme
#   CD200
#   Cycling
#   MHC-II
#
# Raw values remain available in the source-data table above.
# ------------------------------------------------------------------------------

E_long <- E |>
  select(
    territory,
    stem_niche_mean,
    CD200_z_mean,
    cycling_mean,
    MHCII_mean
  ) |>
  pivot_longer(
    cols = -territory,
    names_to = "programme",
    values_to = "value"
  ) |>
  group_by(
    programme
  ) |>
  mutate(
    z =
      as.numeric(
        scale(
          value
        )
      )
  ) |>
  ungroup() |>
  mutate(
    programme =
      recode(
        programme,

        "stem_niche_mean" =
          "Stem-niche",

        "CD200_z_mean" =
          "CD200",

        "cycling_mean" =
          "Cycling",

        "MHCII_mean" =
          "MHC-II"
      ),

    programme =
      factor(
        programme,
        levels = c(
          "Stem-niche",
          "CD200",
          "Cycling",
          "MHC-II"
        )
      )
  )

write_csv(
  E_long,
  file.path(
    SOURCE_DIR,
    "Fig03E_K3_functional_phenotype_standardized.csv"
  )
)


# Raw values to print inside the heatmap.
#
# These retain biological scale while fill represents the relative
# territory-specific enrichment/depletion of each programme.

E_text <- E_long |>
  mutate(
    value_label =
      case_when(
        programme == "Stem-niche" ~
          sprintf(
            "%.2f",
            value
          ),

        programme == "CD200" ~
          sprintf(
            "%.2f",
            value
          ),

        programme == "Cycling" ~
          sprintf(
            "%.2f",
            value
          ),

        programme == "MHC-II" ~
          sprintf(
            "%.2f",
            value
          ),

        TRUE ~
          sprintf(
            "%.2f",
            value
          )
      ),

    text_colour =
      if_else(
        abs(z) >= 0.9,
        "white",
        "#303030"
      )
  )


pE <- ggplot(
  E_long,
  aes(
    x = programme,
    y = territory,
    fill = z
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  geom_text(
    data = E_text,
    aes(
      label = value_label,
      colour = text_colour
    ),
    size = 2.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = COL_DOWN,
    mid = COL_ZERO,
    high = COL_UP,
    midpoint = 0,
    limits = c(
      -1.5,
      1.5
    ),
    oob = scales::squish,
    name = "Relative\nprogramme\nscore"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    title = "E  Functional identity of Bulge territories",
    subtitle = "Fill: standardized across territories; numbers: original mean values",
    x = NULL,
    y = NULL
  ) +
  theme_pub(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1,
        vjust = 1,
        face = "bold"
      ),

    axis.text.y =
      element_text(
        face = "bold",
        size = 8.5
      ),

    axis.line =
      element_blank(),

    axis.ticks =
      element_blank(),

    panel.grid =
      element_blank(),

    legend.position =
      "right",

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )

save_all(
  pE,
  "Fig03E_K3_functional_phenotype",
  width = 5.3,
  height = 4.8
)



# ------------------------------------------------------------------------------
# LOAD FROZEN K3 STATISTICS FOR PANELS F AND FIGURE SUMMARY
# ------------------------------------------------------------------------------
# FROZEN K3 EFFECT TABLE
#
# 22F3F2 already contains the frozen manuscript contrasts:
#   FFA_vs_HC
#   AA_vs_HC
#   AA_vs_FFA
#
# Therefore contrast direction is NOT reconstructed from group1/group2.
# ------------------------------------------------------------------------------

f3_pairwise <- read_csv(
  F3_PAIRWISE,
  show_col_types = FALSE
)

if (nrow(f3_pairwise) != 9) {
  stop(
    "AUDIT FAIL - expected exactly 9 frozen K3 pairwise comparisons; found ",
    nrow(f3_pairwise),
    "."
  )
}

required_k3_cols <- c(
  "contrast",
  "group1",
  "group2",
  "median1",
  "median2",
  "p",
  "FDR"
)

missing_k3_cols <- setdiff(
  required_k3_cols,
  names(f3_pairwise)
)

if (length(missing_k3_cols) > 0) {
  stop(
    "AUDIT FAIL - missing columns in frozen K3 pairwise table: ",
    paste(
      missing_k3_cols,
      collapse = ", "
    )
  )
}

# Standardize territory column name while supporting both frozen
# table conventions used in the pipeline.
if ("territory" %in% names(f3_pairwise)) {

  k3_effect <- f3_pairwise

} else if ("territory_K3" %in% names(f3_pairwise)) {

  k3_effect <- f3_pairwise |>
    rename(
      territory = territory_K3
    )

} else {

  stop(
    "AUDIT FAIL - frozen K3 table contains neither territory nor territory_K3."
  )
}

required_contrasts <- c(
  "FFA_vs_HC",
  "AA_vs_HC",
  "AA_vs_FFA"
)

observed_contrasts <- sort(
  unique(
    k3_effect$contrast
  )
)

if (!setequal(
  observed_contrasts,
  required_contrasts
)) {

  message(
    "Observed contrasts: ",
    paste(
      observed_contrasts,
      collapse = ", "
    )
  )

  stop(
    "AUDIT FAIL - frozen K3 contrast set differs from expected manuscript contrasts."
  )
}

required_territories <- c(
  "T1",
  "T2",
  "T3"
)

observed_territories <- sort(
  unique(
    k3_effect$territory
  )
)

if (!setequal(
  observed_territories,
  required_territories
)) {

  message(
    "Observed territories: ",
    paste(
      observed_territories,
      collapse = ", "
    )
  )

  stop(
    "AUDIT FAIL - frozen K3 territory set differs from T1/T2/T3."
  )
}

# Exact structural audit:
# 3 territories x 3 contrasts = 9 rows.

expected_k3_grid <- tidyr::expand_grid(
  territory = required_territories,
  contrast = required_contrasts
)

observed_k3_grid <- k3_effect |>
  distinct(
    territory,
    contrast
  )

missing_k3_grid <- anti_join(
  expected_k3_grid,
  observed_k3_grid,
  by = c(
    "territory",
    "contrast"
  )
)

duplicate_k3_grid <- k3_effect |>
  count(
    territory,
    contrast,
    name = "n"
  ) |>
  filter(
    n != 1
  )

if (nrow(missing_k3_grid) > 0) {

  print(
    missing_k3_grid,
    n = Inf
  )

  stop(
    "AUDIT FAIL - incomplete K3 territory x contrast grid."
  )
}

if (nrow(duplicate_k3_grid) > 0) {

  print(
    duplicate_k3_grid,
    n = Inf
  )

  stop(
    "AUDIT FAIL - duplicated K3 territory x contrast rows."
  )
}

message(
  "PASS - frozen K3 pairwise statistics loaded: ",
  nrow(k3_effect),
  " comparisons"
)

message(
  "PASS - frozen K3 contrasts reused directly: ",
  paste(
    required_contrasts,
    collapse = ", "
  )
)

message(
  "PASS - complete T1/T2/T3 x contrast grid"
)


# ------------------------------------------------------------------------------
# 13. PANEL F - STEM-HIGH FRACTION BY K3 TERRITORY
#
# This is deliberately shown separately from Panel E because the proportion
# of stem-high cells has an immediately interpretable biological scale.
#
# Frozen values:
#
# T1 = 19.5%
# T2 = 43.9%
# T3 = 76.6%
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "PANEL F - STEM-HIGH FRACTION BY K3 TERRITORY"
)
message(
  "============================================================"
)

F <- E |>
  select(
    territory,
    n_cells,
    n_patients,
    pct_stem_high,
    CD200_detected_pct
  )

write_csv(
  F,
  file.path(
    SOURCE_DIR,
    "Fig03F_K3_stem_high_fraction.csv"
  )
)

pF <- ggplot(
  F,
  aes(
    x = territory,
    y = pct_stem_high,
    fill = territory
  )
) +
  geom_col(
    width = 0.64,
    colour = "#3A3A3A",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label =
        sprintf(
          "%.1f%%",
          pct_stem_high
        )
    ),
    vjust = -0.45,
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values =
      K3_COLS,
    guide =
      "none"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      90
    ),
    breaks = seq(
      0,
      80,
      by = 20
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "C  Stem-high state across Bulge territories",
    subtitle = "Cell-level classification within the shared Bulge/Isthmus manifold",
    x = "Bulge territory",
    y = "Stem-high cells"
  ) +
  theme_pub(
    base_size = 8
  ) +
  theme(
    axis.text.x =
      element_text(
        face = "bold"
      ),

    legend.position =
      "none",

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )

save_all(
  pF,
  "Fig03F_K3_stem_high_fraction",
  width = 4.1,
  height = 4.8
)




# ------------------------------------------------------------------------------
# 14. FINAL FIGURE ASSEMBLY
#
# Main Figure 3 layout:
#
# B  Patient-level functional remodeling of the Bulge/Isthmus niche
#
# C  Disease-specific occupancy of the shared Bulge trajectory landscape
# F  Stem-high enrichment across K3 territories
#
# D  Patient-level differential occupancy of K3 territories
# E  Functional identity of K3 territories
#
# Panel A (absolute Bulge abundance) remains generated as an individual panel
# and in the source-data outputs, but is omitted from the main Figure 3.
# No underlying analysis or statistical result is removed.
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)

message(
  "FINAL FIGURE ASSEMBLY"
)

message(
  "============================================================"
)

message(
  "Main Figure 3: B + C + F + D + E"
)

message(
  "Panel A retained as standalone output but omitted from main assembly"
)


# ------------------------------------------------------------------------------
# ROW 1
# Patient-level functional remodeling
# ------------------------------------------------------------------------------

top_row <-
  pB


# ------------------------------------------------------------------------------
# ROW 2
# Shared trajectory landscape + stem-high interpretation
#
# Panel C receives 80% of the available width.
# Panel F occupies the remaining 20%.
# ------------------------------------------------------------------------------

middle_row <-
  pC +
  pF +
  patchwork::plot_layout(
    widths = c(
      0.80,
      0.20
    )
  )


# ------------------------------------------------------------------------------
# ROW 3
# Quantitative K3 occupancy + functional identity
#
# Panel D requires more width because it contains the three K3 territories
# with patient-level disease distributions and statistical annotations.
# ------------------------------------------------------------------------------

bottom_row <-
  pD +
  pE +
  patchwork::plot_layout(
    widths = c(
      0.64,
      0.36
    )
  )


# ------------------------------------------------------------------------------
# COMPLETE FIGURE
# ------------------------------------------------------------------------------

Figure03_FINAL <-
  top_row /
  middle_row /
  bottom_row +
  patchwork::plot_layout(
    heights = c(
      0.82,
      1.02,
      0.96
    )
  ) +
  patchwork::plot_annotation(
    title =
      "AA and FFA occupy distinct functional states of the Bulge/Isthmus regenerative niche"
  ) &
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11,
        hjust = 0
      )
  )


save_all(
  Figure03_FINAL,
  "Figure03_FINAL",
  width = 15.8,
  height = 12.0
)


message(
  "PASS - Panel A omitted from main Figure 3"
)

message(
  "PASS - Panel B spans complete top row"
)

message(
  "PASS - Panels C/F assembled at 80/20 width"
)

message(
  "PASS - Panels D/E assembled at 64/36 width"
)

# ------------------------------------------------------------------------------
# 15. KEY STATISTICS
# ------------------------------------------------------------------------------

figure_summary <- tibble::tribble(
  ~metric, ~value,

  "Bulge_cells_total",
  nrow(ledger),

  "Bulge_patients_total",
  n_distinct(ledger$patient_id),

  "HC_patients",
  sum(
    f2_patient$disease == "HC"
  ),

  "FFA_patients",
  sum(
    f2_patient$disease == "FFA"
  ),

  "AA_patients",
  sum(
    f2_patient$disease == "AA"
  ),

  "T1_cells",
  sum(
    ledger$territory_K3 == "T1"
  ),

  "T2_cells",
  sum(
    ledger$territory_K3 == "T2"
  ),

  "T3_cells",
  sum(
    ledger$territory_K3 == "T3"
  ),

  "T1_stem_high_pct",
  E$pct_stem_high[
    E$territory == "T1"
  ],

  "T2_stem_high_pct",
  E$pct_stem_high[
    E$territory == "T2"
  ],

  "T3_stem_high_pct",
  E$pct_stem_high[
    E$territory == "T3"
  ],

  "AA_median_T3_pct",
  # Frozen K3 disease-effect summary.
# No new statistical testing is performed here.
k3_effect <- read_csv(
  F3_EFFECT,
  show_col_types = FALSE
)

k3_effect$median1[
    k3_effect$territory == "T3" &
      k3_effect$contrast == "AA_vs_FFA"
  ],

  "FFA_median_T3_pct",
  k3_effect$median2[
    k3_effect$territory == "T3" &
      k3_effect$contrast == "AA_vs_FFA"
  ],

  "AA_median_T2_pct",
  k3_effect$median1[
    k3_effect$territory == "T2" &
      k3_effect$contrast == "AA_vs_FFA"
  ],

  "FFA_median_T2_pct",
  k3_effect$median2[
    k3_effect$territory == "T2" &
      k3_effect$contrast == "AA_vs_FFA"
  ]
)

write_csv(
  figure_summary,
  file.path(
    SOURCE_DIR,
    "Fig03_key_statistics.csv"
  )
)


# ------------------------------------------------------------------------------
# 16. FORMAL CHECKS
# ------------------------------------------------------------------------------

formal_checks <- tibble(
  check = c(
    "Bulge_cells_5420",
    "Bulge_patients_16",
    "K3_three_territories",
    "T1_cells_1518",
    "T2_cells_1757",
    "T3_cells_2145",
    "T1_stem_high_19.5",
    "T2_stem_high_43.9",
    "T3_stem_high_76.6",
    "T2_global_FDR_lt_005",
    "T3_global_FDR_lt_005",
    "T2_AA_vs_FFA_FDR_lt_005",
    "T3_AA_vs_FFA_FDR_lt_005",
    "T2_FFA_vs_HC_FDR_lt_005",
    "T3_FFA_vs_HC_FDR_lt_005"
  ),

  pass = c(
    nrow(ledger) == 5420,

    n_distinct(
      ledger$patient_id
    ) == 16,

    setequal(
      unique(
        as.character(
          ledger$territory_K3
        )
      ),
      c(
        "T1",
        "T2",
        "T3"
      )
    ),

    sum(
      ledger$territory_K3 == "T1"
    ) == 1518,

    sum(
      ledger$territory_K3 == "T2"
    ) == 1757,

    sum(
      ledger$territory_K3 == "T3"
    ) == 2145,

    abs(
      E$pct_stem_high[
        E$territory == "T1"
      ] -
        19.5
    ) < 0.1,

    abs(
      E$pct_stem_high[
        E$territory == "T2"
      ] -
        43.9
    ) < 0.1,

    abs(
      E$pct_stem_high[
        E$territory == "T3"
      ] -
        76.6
    ) < 0.1,

    k3_global$FDR[
      k3_global$territory == "T2"
    ] < 0.05,

    k3_global$FDR[
      k3_global$territory == "T3"
    ] < 0.05,

    k3_pairwise$FDR[
      k3_pairwise$territory == "T2" &
        k3_pairwise$contrast == "AA_vs_FFA"
    ] < 0.05,

    k3_pairwise$FDR[
      k3_pairwise$territory == "T3" &
        k3_pairwise$contrast == "AA_vs_FFA"
    ] < 0.05,

    k3_pairwise$FDR[
      k3_pairwise$territory == "T2" &
        k3_pairwise$contrast == "FFA_vs_HC"
    ] < 0.05,

    k3_pairwise$FDR[
      k3_pairwise$territory == "T3" &
        k3_pairwise$contrast == "FFA_vs_HC"
    ] < 0.05
  )
)

write_csv(
  formal_checks,
  file.path(
    AUDIT_DIR,
    "Fig03_formal_checks.csv"
  )
)

if (
  !all(
    formal_checks$pass
  )
) {

  print(
    formal_checks
  )

  stop(
    "STATUS=FAIL - one or more Figure 3 formal checks failed."
  )
}


# ------------------------------------------------------------------------------
# 17. OUTPUT AUDIT
# ------------------------------------------------------------------------------

expected_stems <- c(
  "Fig03A_Bulge_abundance",
  "Fig03B_Bulge_functional_phenotype",
  "Fig03C_Bulge_K3_territories",
  "Fig03D_K3_disease_occupancy",
  "Fig03E_K3_functional_phenotype",
  "Fig03F_K3_stem_high_fraction",
  "Figure03_FINAL"
)

output_audit <- tibble(
  stem =
    expected_stems,

  svg =
    file.exists(
      file.path(
        SVG_DIR,
        paste0(
          expected_stems,
          ".svg"
        )
      )
    ),

  pdf =
    file.exists(
      file.path(
        PDF_DIR,
        paste0(
          expected_stems,
          ".pdf"
        )
      )
    ),

  png =
    file.exists(
      file.path(
        PNG_DIR,
        paste0(
          expected_stems,
          ".png"
        )
      )
    )
)

write_csv(
  output_audit,
  file.path(
    AUDIT_DIR,
    "Fig03_output_audit.csv"
  )
)

if (
  !all(
    output_audit$svg &
      output_audit$pdf &
      output_audit$png
  )
) {

  print(
    output_audit
  )

  stop(
    "STATUS=FAIL - incomplete Figure 3 output set."
  )
}


# ------------------------------------------------------------------------------
# 18. SESSION INFO
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file =
    file.path(
      AUDIT_DIR,
      "Fig03_sessionInfo.txt"
    )
)


# ------------------------------------------------------------------------------
# 19. FINAL REPORT
# ------------------------------------------------------------------------------

message(
  "\n============================================================"
)
message(
  "FIGURE 03 - FINAL REPORT"
)
message(
  "============================================================"
)

message(
  "Bulge/Isthmus cells: ",
  format(
    nrow(
      ledger
    ),
    big.mark = ","
  )
)

message(
  "Patients: ",
  n_distinct(
    ledger$patient_id
  )
)

message(
  "K3 territories:"
)

message(
  "  T1: ",
  sum(
    ledger$territory_K3 == "T1"
  ),
  " cells; stem-high = ",
  sprintf(
    "%.1f%%",
    E$pct_stem_high[
      E$territory == "T1"
    ]
  )
)

message(
  "  T2: ",
  sum(
    ledger$territory_K3 == "T2"
  ),
  " cells; stem-high = ",
  sprintf(
    "%.1f%%",
    E$pct_stem_high[
      E$territory == "T2"
    ]
  )
)

message(
  "  T3: ",
  sum(
    ledger$territory_K3 == "T3"
  ),
  " cells; stem-high = ",
  sprintf(
    "%.1f%%",
    E$pct_stem_high[
      E$territory == "T3"
    ]
  )
)

message(
  ""
)

message(
  "Biological summary:"
)

message(
  "  T1 = stem-low / cycling-enriched territory"
)

message(
  "  T2 = intermediate stem state with higher MHC-II activity"
)

message(
  "  T3 = stem-high / CD200-enriched regenerative territory"
)

message(
  "  AA preferentially occupies T3 relative to FFA."
)

message(
  "  FFA preferentially occupies T2 and is depleted from T3 relative to HC and AA."
)

message(
  ""
)

message(
  "Patient as statistical replicate: YES"
)

message(
  "Patient-balanced Bulge analyses: YES"
)

message(
  "New differential-expression testing: NO"
)

message(
  "Frozen K3 solution reused: YES"
)

message(
  "SVG + PDF + PNG outputs: YES"
)

message(
  "Output: ",
  OUT_ROOT
)

message(
  "STATUS=PASS"
)

