#!/usr/bin/env Rscript

# ==============================================================================
# FIGURE 04 — TRANSCRIPTIONAL DIVERGENCE
#
# AA / FFA / HC scalp single-cell atlas
#
# A — Differential-expression burden across L2 states
# B — AA vs FFA Bulge/Isthmus volcano
# C — Patient-level Bulge/Isthmus expression landscape
# D — Cross-epithelial distribution of Bulge-selected genes
#
# IMPORTANT
# - Frozen 18C differential-expression results are reused.
# - No differential-expression model is refitted.
# - Patient-level expression in panel C is derived from frozen 18B pseudobulk
#   counts using TMM normalization and logCPM transformation for display only.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "purrr",
  "ggplot2",
  "ggrepel",
  "ggtext",
  "ggh4x",
  "patchwork",
  "edgeR",
  "Matrix",
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
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(ggh4x)
  library(patchwork)
  library(edgeR)
  library(Matrix)
  library(scales)
  library(svglite)
})


# ==============================================================================
# 1. ROOT + INPUTS
# ==============================================================================

ROOT <- normalizePath(
  ".",
  mustWork = TRUE
)

DE_ROOT <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/18_PSEUDOBULK_DE"
)

SUMMARY_FILE <- file.path(
  DE_ROOT,
  "18C_DE/tables/18C_contrast_summary.csv"
)

L2_SUMMARY_FILE <- file.path(
  DE_ROOT,
  "18C_DE/tables/18C_DE_summary_by_L2.csv"
)

ALL_SIG_FILE <- file.path(
  DE_ROOT,
  "18C_DE/tables/18C_all_FDR005_genes.csv"
)

BULGE_FILE <- file.path(
  DE_ROOT,
  paste0(
    "18C_DE/full_results/",
    "18C_AA_vs_FFA__07_EPITHELIAL_L2__Bulge_Isthmus.csv"
  )
)

PB_COUNTS_FILE <- file.path(
  DE_ROOT,
  "18B_pseudobulk/objects/18B_07_EPITHELIAL_L2_pseudobulk_counts.rds"
)

PB_META_FILE <- file.path(
  DE_ROOT,
  "18B_pseudobulk/tables/18B_pseudobulk_metadata_ALL.csv"
)

L1_PALETTE_FILE <- file.path(
  ROOT,
  "PIPELINE_V2/config/figure_palettes/L1_palette.tsv"
)

L2_PALETTE_FILE <- file.path(
  ROOT,
  "PIPELINE_V2/config/figure_palettes/L2_palette.tsv"
)

required_files <- c(
  SUMMARY_FILE,
  L2_SUMMARY_FILE,
  ALL_SIG_FILE,
  BULGE_FILE,
  PB_COUNTS_FILE,
  PB_META_FILE,
  L1_PALETTE_FILE,
  L2_PALETTE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
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


# ==============================================================================
# 2. OUTPUT
# ==============================================================================

OUT_ROOT <- file.path(
  ROOT,
  "PIPELINE_V2/results/publication_figures/FIG04_TRANSCRIPTIONAL_DIVERGENCE_FDR010"
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

dir.create(
  PDF_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  SVG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PNG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  SOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. CONSTANTS
# ==============================================================================

COL_HC <- "#8F8D8E"
COL_FFA <- "#4D9E4F"
COL_AA <- "#A8142A"

COL_UP <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_ZERO <- "#F7F7F7"
COL_DARK <- "#303030"
COL_GREY <- "#BDBDBD"

# Experimental Figure 4 threshold
FDR_THRESHOLD <- 0.10
FDR_LABEL <- "FDR < 0.10"

CONTRAST_ORDER <- c(
  "AA_vs_HC",
  "FFA_vs_HC",
  "AA_vs_FFA"
)

CONTRAST_LABELS <- c(
  "AA_vs_HC" = "AA vs HC",
  "FFA_vs_HC" = "FFA vs HC",
  "AA_vs_FFA" = "AA vs FFA"
)

DISEASE_ORDER <- c(
  "HC",
  "FFA",
  "AA"
)

EPITHELIAL_L2_ORDER <- c(
  "Basal_epithelial",
  "Bulge_Isthmus",
  "Hair_shaft",
  "IFE_Infundibulum_differentiating",
  "IFE_Infundibulum_suprabasal",
  "IFE_Infundibulum_terminal",
  "IRS_companion",
  "ORS",
  "ORS_basal"
)


# ==============================================================================
# 4. HELPERS
# ==============================================================================

clean_l2 <- function(x) {

  x |>
    str_replace_all(
      "_",
      " "
    ) |>
    str_replace_all(
      "IFE Infundibulum",
      "IFE/Infundibulum"
    ) |>
    str_replace_all(
      "IRS companion",
      "IRS/companion"
    )
}


module_to_L1 <- function(x) {

  case_when(
    str_detect(x, "EPITHELIAL") ~ "Epithelial",
    str_detect(x, "IMMUNE") ~ "Immune",
    str_detect(x, "FIBROBLAST") ~ "Fibroblast",
    str_detect(x, "VASCULAR") ~ "Endothelial",
    str_detect(x, "MURAL") ~ "Mural",
    str_detect(x, "SWEAT") ~ "Sweat_gland",
    str_detect(x, "MELANOCYTE") ~ "Melanocyte",
    str_detect(x, "SEBOCYTE") ~ "Sebocyte",
    TRUE ~ NA_character_
  )
}


safe_neglog10 <- function(x, cap = 20) {

  pmin(
    -log10(
      pmax(
        x,
        .Machine$double.xmin
      )
    ),
    cap
  )
}


theme_pub <- function(base_size = 8) {

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
          size = base_size + 1.5,
          hjust = 0
        ),

      plot.subtitle =
        element_text(
          size = base_size,
          colour = "#444444",
          hjust = 0
        ),

      strip.background =
        element_rect(
          fill = "#F2F2F2",
          colour = "#555555",
          linewidth = 0.3
        ),

      strip.text =
        element_text(
          face = "bold",
          size = base_size
        ),

      axis.title =
        element_text(
          size = base_size
        ),

      axis.text =
        element_text(
          size = base_size - 0.5,
          colour = "#303030"
        ),

      legend.title =
        element_text(
          size = base_size - 0.3,
          face = "bold"
        ),

      legend.text =
        element_text(
          size = base_size - 0.5
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

  affinity_file <- file.path(
    PDF_DIR,
    paste0(
      stem,
      "_AFFINITY.pdf"
    )
  )

  affinity_ok <- tryCatch(
    {
      suppressWarnings(
        ggsave(
          affinity_file,
          plot,
          width = width,
          height = height,
          units = "in",
          device = grDevices::cairo_pdf,
          bg = "white"
        )
      )

      file.exists(
        affinity_file
      )
    },
    error = function(e) {
      message(
        "NOTE — Cairo/Affinity PDF unavailable for ",
        stem,
        ": ",
        conditionMessage(e)
      )
      FALSE
    }
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


# ==============================================================================
# 5. READ INPUTS
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "READ INPUTS"
)
message(
  "============================================================"
)

contrast_summary <- read_csv(
  SUMMARY_FILE,
  show_col_types = FALSE
)

de_summary <- read_csv(
  L2_SUMMARY_FILE,
  show_col_types = FALSE
)

all_sig <- read_csv(
  ALL_SIG_FILE,
  show_col_types = FALSE
)

bulge <- read_csv(
  BULGE_FILE,
  show_col_types = FALSE
)

pb_counts <- readRDS(
  PB_COUNTS_FILE
)

pb_meta_all <- read_csv(
  PB_META_FILE,
  show_col_types = FALSE
)

L1_palette <- read.delim(
  L1_PALETTE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

L2_palette <- read.delim(
  L2_PALETTE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

L1_cols <- setNames(
  L1_palette$hex,
  L1_palette$L1
)

L2_cols <- setNames(
  L2_palette$hex,
  L2_palette$final_L2_label
)


# ==============================================================================
# 6. FORMAL INPUT AUDIT — BASIC INPUT INTEGRITY
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "FORMAL INPUT AUDIT — FDR < 0.10 EXPERIMENT"
)
message(
  "============================================================"
)

if (nrow(bulge) != 10486) {
  stop(
    "AUDIT FAIL — expected 10,486 tested Bulge genes."
  )
}

bulge_sig <- bulge |>
  filter(
    FDR_BH < FDR_THRESHOLD
  )

bulge_n_sig <- nrow(
  bulge_sig
)

bulge_n_AA <- sum(
  bulge_sig$logFC > 0,
  na.rm = TRUE
)

bulge_n_FFA <- sum(
  bulge_sig$logFC < 0,
  na.rm = TRUE
)

cd200_row <- bulge |>
  filter(
    gene == "CD200"
  )

if (nrow(cd200_row) != 1) {
  stop(
    "AUDIT FAIL — CD200 not uniquely identified."
  )
}

message(
  "PASS — frozen Bulge universe = 10,486 tested genes"
)

message(
  "PASS — experimental threshold = FDR < ",
  FDR_THRESHOLD
)

# ==============================================================================
# 6B. RECONSTRUCT DE SUMMARY AT EXPERIMENTAL FDR THRESHOLD
# ==============================================================================

full_result_files <- list.files(
  file.path(
    DE_ROOT,
    "18C_DE/full_results"
  ),
  pattern = "^18C_.*\\.csv$",
  full.names = TRUE
)

if (length(full_result_files) == 0) {
  stop(
    "AUDIT FAIL — no frozen 18C full-result files found."
  )
}

summarise_full_DE <- function(f) {

  x <- read_csv(
    f,
    show_col_types = FALSE
  )

  tibble(
    contrast =
      unique(x$contrast),

    numerator =
      unique(x$numerator),

    denominator =
      unique(x$denominator),

    module =
      unique(x$module),

    final_L2_label =
      unique(x$final_L2_label),

    n_genes_tested =
      nrow(x),

    n_FDR_threshold =
      sum(
        x$FDR_BH < FDR_THRESHOLD,
        na.rm = TRUE
      ),

    n_up_FDR_threshold =
      sum(
        x$FDR_BH < FDR_THRESHOLD &
          x$logFC > 0,
        na.rm = TRUE
      ),

    n_down_FDR_threshold =
      sum(
        x$FDR_BH < FDR_THRESHOLD &
          x$logFC < 0,
        na.rm = TRUE
      )
  )
}

de_summary_threshold <- purrr::map_dfr(
  full_result_files,
  summarise_full_DE
)

contrast_summary_threshold <- de_summary_threshold |>
  group_by(
    contrast
  ) |>
  summarise(
    n_L2_tested =
      n(),

    n_L2_with_threshold =
      sum(
        n_FDR_threshold > 0
      ),

    total_genes_threshold =
      sum(
        n_FDR_threshold
      ),

    total_up_threshold =
      sum(
        n_up_FDR_threshold
      ),

    total_down_threshold =
      sum(
        n_down_FDR_threshold
      ),

    .groups = "drop"
  )

write_csv(
  de_summary_threshold,
  file.path(
    SOURCE_DIR,
    "Fig04_FDR010_DE_summary_by_L2.csv"
  )
)

write_csv(
  contrast_summary_threshold,
  file.path(
    SOURCE_DIR,
    "Fig04_FDR010_contrast_summary.csv"
  )
)

message(
  "\nFDR < 0.10 experimental totals:"
)

print(
  contrast_summary_threshold,
  n = Inf
)

direct_total_threshold <- contrast_summary_threshold |>
  filter(
    contrast == "AA_vs_FFA"
  ) |>
  pull(
    total_genes_threshold
  )

if (
  length(
    direct_total_threshold
  ) != 1
) {
  stop(
    "AUDIT FAIL — AA-vs-FFA total at FDR < 0.10 not uniquely resolved."
  )
}

bulge_fraction <-
  100 *
  bulge_n_sig /
  direct_total_threshold

message(
  sprintf(
    "Bulge/Isthmus fraction of AA-FFA DE at FDR < 0.10: %.1f%%",
    bulge_fraction
  )
)

# ==============================================================================
# 7. PANEL A — DE BURDEN ACROSS CELL STATES
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "PANEL A — DE BURDEN"
)
message(
  "============================================================"
)

A <- de_summary_threshold |>
  mutate(
    L1 =
      module_to_L1(
        module
      ),

    L2_clean =
      clean_l2(
        final_L2_label
      ),

    L2_hex =
      unname(
        L2_cols[
          final_L2_label
        ]
      ),

    contrast =
      factor(
        contrast,
        levels =
          CONTRAST_ORDER
      ),

    n_FDR_005 =
      n_FDR_threshold,

    n_up_FDR_005 =
      n_up_FDR_threshold,

    n_down_FDR_005 =
      n_down_FDR_threshold,

    robust_ge50 =
      TRUE
  )

if (any(is.na(A$L1))) {
  stop(
    "AUDIT FAIL — unmapped L1 in Panel A."
  )
}

if (any(is.na(A$L2_hex))) {
  stop(
    "AUDIT FAIL — unmapped frozen L2 color in Panel A."
  )
}

# Keep all tested states, including zero-hit states.
A_long <- A |>
  select(
    contrast,
    module,
    L1,
    final_L2_label,
    L2_clean,
    L2_hex,
    n_FDR_005,
    n_up_FDR_005,
    n_down_FDR_005,
    robust_ge50
  ) |>
  pivot_longer(
    cols = c(
      n_up_FDR_005,
      n_down_FDR_005
    ),
    names_to = "direction",
    values_to = "n_genes"
  ) |>
  mutate(
    direction =
      recode(
        direction,
        "n_up_FDR_005" =
          "Higher in numerator",
        "n_down_FDR_005" =
          "Lower in numerator"
      ),

    signed_n =
      if_else(
        direction ==
          "Higher in numerator",
        n_genes,
        -n_genes
      )
  )

L1_rank <- setNames(
  L1_palette$display_order,
  L1_palette$L1
)

A_order <- A |>
  group_by(
    final_L2_label,
    L2_clean,
    L1
  ) |>
  summarise(
    max_hits =
      max(
        n_FDR_005,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  mutate(
    L1_order =
      unname(
        L1_rank[
          L1
        ]
      )
  ) |>
  arrange(
    L1_order,
    desc(max_hits),
    L2_clean
  )

A_levels <- rev(
  A_order$L2_clean
)

A_long <- A_long |>
  mutate(
    L2_clean =
      factor(
        L2_clean,
        levels =
          A_levels
      )
  )

A_label_table <- A |>
  distinct(
    L2_clean,
    L2_hex
  ) |>
  mutate(
    label_html =
      paste0(
        "<span style='color:",
        L2_hex,
        ";'><b>",
        L2_clean,
        "</b></span>"
      )
  )

A_labels <- setNames(
  A_label_table$label_html,
  A_label_table$L2_clean
)

write_csv(
  A_long,
  file.path(
    SOURCE_DIR,
    "Fig04A_DE_burden_source_data.csv"
  )
)

# Contrast-level annotations.
A_totals <- contrast_summary |>
  filter(
    contrast %in%
      CONTRAST_ORDER
  ) |>
  mutate(
    contrast =
      factor(
        contrast,
        levels =
          CONTRAST_ORDER
      ),

    contrast_label =
      unname(
        CONTRAST_LABELS[
          as.character(
            contrast
          )
        ]
      ),

    annotation =
      paste0(
        format(
          total_genes_FDR005,
          big.mark = ","
        ),
        " DE genes · ",
        n_L2_with_FDR005,
        " L2 states"
      )
  )

# Highlight Bulge in the direct comparison.
A_bulge <- A |>
  filter(
    contrast == "AA_vs_FFA",
    final_L2_label == "Bulge_Isthmus"
  )

if (
  nrow(
    A_bulge
  ) != 1
) {
  stop(
    "AUDIT FAIL — AA-vs-FFA Bulge row not uniquely identified in Panel A."
  )
}

if (
  A_bulge$n_FDR_005 != bulge_n_sig
) {
  stop(
    "AUDIT FAIL — Panel A Bulge count does not match the ",
    "direct Bulge full-result count at FDR < ",
    FDR_THRESHOLD,
    "."
  )
}

# Symmetric scale across facets.
A_limit <- max(
  abs(
    A_long$signed_n
  ),
  na.rm = TRUE
)

A_limit <- ceiling(
  A_limit / 100
) * 100

pA <- ggplot(
  A_long,
  aes(
    x = signed_n,
    y = L2_clean
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "#555555"
  ) +
  geom_col(
    aes(
      fill = direction
    ),
    width = 0.72,
    colour = NA
  ) +
  facet_grid(
    cols = vars(contrast),
    scales = "fixed",
    labeller =
      labeller(
        contrast =
          as_labeller(
            CONTRAST_LABELS
          )
      )
  ) +
  scale_fill_manual(
    values = c(
      "Higher in numerator" =
        COL_UP,
      "Lower in numerator" =
        COL_DOWN
    ),
    name = NULL
  ) +
  scale_x_continuous(
    trans =
      scales::pseudo_log_trans(
        sigma = 50,
        base = 10
      ),

    breaks = c(
      -1500,
      -1000,
      -500,
      -250,
      -100,
      -50,
      0,
      50,
      100,
      250,
      500,
      1000,
      1500
    ),

    labels =
      function(x) {
        abs(x)
      },

    limits = c(
      -1700,
      1700
    ),

    expand =
      expansion(
        mult = c(
          0.015,
          0.015
        )
      )
  ) +
  scale_y_discrete(
    labels =
      A_labels
  ) +
  labs(
    title =
      "A  Differential-expression burden across cell states",
    subtitle =
      "FDR < 0.10 · right: higher in numerator · left: lower in numerator · symmetric pseudo-log scale",
    x =
      "Number of differentially expressed genes",
    y =
      NULL
  ) +
  theme_pub(
    base_size = 7.7
  ) +
  theme(
    axis.text.y =
      ggtext::element_markdown(
        size = 6.6
      ),

    panel.spacing.x =
      unit(
        0.7,
        "lines"
      ),

    legend.position =
      "top",

    legend.justification =
      "left",

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
  "Fig04A_DE_burden_across_L2",
  width = 10.8,
  height = 7.8
)


# ==============================================================================
# 8. PANEL B — AA vs FFA BULGE/ISTHMUS VOLCANO
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "PANEL B — BULGE/ISTHMUS AA vs FFA VOLCANO"
)
message(
  "============================================================"
)

B <- bulge |>
  mutate(
    significance =
      case_when(
        FDR_BH < FDR_THRESHOLD &
          logFC > 0 ~
          "Higher in AA",

        FDR_BH < FDR_THRESHOLD &
          logFC < 0 ~
          "Higher in FFA",

        TRUE ~
          "Not significant"
      ),

    neglog10FDR =
      safe_neglog10(
        FDR_BH,
        cap = 12
      ),

    abs_logFC =
      abs(
        logFC
      )
  )

# Objective statistical label set:
# top 8 by FDR in each direction, plus significant prespecified niche genes.
B_top_AA <- B |>
  filter(
    significance ==
      "Higher in AA"
  ) |>
  arrange(
    FDR_BH,
    desc(abs_logFC)
  ) |>
  slice_head(
    n = 8
  ) |>
  pull(
    gene
  )

B_top_FFA <- B |>
  filter(
    significance ==
      "Higher in FFA"
  ) |>
  arrange(
    FDR_BH,
    desc(abs_logFC)
  ) |>
  slice_head(
    n = 8
  ) |>
  pull(
    gene
  )

B_prespecified <- c(
  "CD200",
  "KRT15",
  "LHX2",
  "MCM5"
)

B_label_genes <- unique(
  c(
    B_top_AA,
    B_top_FFA,
    B_prespecified
  )
)

B_labels <- B |>
  filter(
    gene %in%
      B_label_genes
  )

write_csv(
  B,
  file.path(
    SOURCE_DIR,
    "Fig04B_Bulge_AA_vs_FFA_volcano_source_data.csv"
  )
)

write_csv(
  B_labels |>
    select(
      gene,
      logFC,
      logCPM,
      F,
      PValue,
      FDR_BH,
      significance
    ),
  file.path(
    SOURCE_DIR,
    "Fig04B_labelled_genes.csv"
  )
)

pB <- ggplot(
  B,
  aes(
    x = logFC,
    y = neglog10FDR
  )
) +
  geom_hline(
    yintercept =
      -log10(
        FDR_THRESHOLD
      ),
    linetype = "dashed",
    linewidth = 0.35,
    colour = "#777777"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "#777777"
  ) +
  geom_point(
    data =
      B |>
        filter(
          significance ==
            "Not significant"
        ),
    size = 1.0,
    alpha = 0.22,
    colour = "#BDBDBD"
  ) +
  geom_point(
    data =
      B |>
        filter(
          significance ==
            "Higher in AA"
        ),
    size = 1.35,
    alpha = 0.72,
    colour = COL_AA
  ) +
  geom_point(
    data =
      B |>
        filter(
          significance ==
            "Higher in FFA"
        ),
    size = 1.35,
    alpha = 0.72,
    colour = COL_FFA
  ) +
  ggrepel::geom_text_repel(
    data =
      B_labels,
    aes(
      label = gene
    ),
    size = 2.35,
    box.padding = 0.35,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.size = 0.20,
    max.overlaps = Inf,
    seed = 4
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    hjust = 1.04,
    vjust = 1.15,
    size = 2.35,
    label =
      paste0(
        bulge_n_sig,
        " genes at FDR < 0.10\n",
        bulge_n_AA,
        " higher in AA · ",
        bulge_n_FFA,
        " higher in FFA"
      ),
    fill = "white",
    colour = "#303030",
    label.size = 0.25
  ) +
  labs(
    title =
      "B  Bulge/Isthmus transcriptional divergence",
    subtitle =
      "Direct AA versus FFA patient-level pseudobulk comparison",
    x =
      "log2 fold-change (AA / FFA)",
    y =
      expression(
        -log[10](
          "FDR"
        )
      )
  ) +
  theme_pub(
    base_size = 8
  ) +
  theme(
    plot.margin =
      margin(
        5,
        10,
        5,
        5
      )
  )

save_all(
  pB,
  "Fig04B_Bulge_AA_vs_FFA_volcano",
  width = 7.2,
  height = 6.2
)


# ==============================================================================
# 9. SELECT BULGE GENES FOR PANELS C + D
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "GENE SELECTION FOR PANELS C + D"
)
message(
  "============================================================"
)

# Selection rules:
# 1. Derived exclusively from direct AA-vs-FFA Bulge DE.
# 2. FDR < 0.10.
# 3. logCPM >= 3.
# 4. Top statistical candidates in each direction.
# 5. Add significant prespecified niche genes if not already selected.
# 6. HC is not used for gene selection.

C_candidates <- B |>
  filter(
    FDR_BH < FDR_THRESHOLD,
    logCPM >= 3
  ) |>
  mutate(
    direction =
      if_else(
        logFC > 0,
        "Higher in AA",
        "Higher in FFA"
      )
  )

C_AA <- C_candidates |>
  filter(
    direction ==
      "Higher in AA"
  ) |>
  arrange(
    FDR_BH,
    desc(abs_logFC)
  ) |>
  slice_head(
    n = 10
  )

C_FFA <- C_candidates |>
  filter(
    direction ==
      "Higher in FFA"
  ) |>
  arrange(
    FDR_BH,
    desc(abs_logFC)
  ) |>
  slice_head(
    n = 10
  )

C_bio <- C_candidates |>
  filter(
    gene %in%
      c(
        "CD200",
        "KRT15",
        "LHX2",
        "MCM5"
      )
  )

selected_genes_table <- bind_rows(
  C_AA |>
    mutate(
      selection_reason =
        "Top10_FDR_AA"
    ),

  C_FFA |>
    mutate(
      selection_reason =
        "Top10_FDR_FFA"
    ),

  C_bio |>
    mutate(
      selection_reason =
        "Significant_prespecified_niche"
    )
) |>
  arrange(
    gene,
    selection_reason
  ) |>
  group_by(
    gene
  ) |>
  summarise(
    logFC =
      first(
        logFC
      ),

    logCPM =
      first(
        logCPM
      ),

    F =
      first(
        F
      ),

    PValue =
      first(
        PValue
      ),

    FDR_BH =
      first(
        FDR_BH
      ),

    direction =
      first(
        direction
      ),

    selection_reason =
      paste(
        unique(
          selection_reason
        ),
        collapse = ";"
      ),

    .groups =
      "drop"
  ) |>
  arrange(
    direction,
    FDR_BH
  )

selected_genes <- selected_genes_table$gene

write_csv(
  selected_genes_table,
  file.path(
    SOURCE_DIR,
    "Fig04_selected_Bulge_genes.csv"
  )
)

message(
  "Selected Bulge genes: ",
  length(
    selected_genes
  )
)

print(
  selected_genes_table,
  n = Inf
)


# ==============================================================================
# 10. PANEL C — PATIENT-LEVEL BULGE EXPRESSION
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "PANEL C — PATIENT-LEVEL BULGE EXPRESSION"
)
message(
  "============================================================"
)

if (
  !inherits(
    pb_counts,
    "dgCMatrix"
  )
) {
  stop(
    "AUDIT FAIL — expected sparse dgCMatrix pseudobulk counts."
  )
}

pb_meta <- pb_meta_all |>
  filter(
    module ==
      "07_EPITHELIAL_L2"
  )

if (
  ncol(
    pb_counts
  ) !=
    nrow(
      pb_meta
    )
) {

  # Attempt exact column-name matching instead of positional matching.
  if (
    !is.null(
      colnames(
        pb_counts
      )
    ) &&
      all(
        colnames(
          pb_counts
        ) %in%
          pb_meta$pseudobulk_id
      )
  ) {

    pb_meta <- pb_meta[
      match(
        colnames(
          pb_counts
        ),
        pb_meta$pseudobulk_id
      ),
      ,
      drop = FALSE
    ]

  } else {

    stop(
      "AUDIT FAIL — epithelial pseudobulk matrix and metadata dimensions differ."
    )
  }
}

# Prefer explicit column-name matching whenever available.
if (
  !is.null(
    colnames(
      pb_counts
    )
  )
) {

  if (
    !all(
      colnames(
        pb_counts
      ) %in%
        pb_meta$pseudobulk_id
    )
  ) {
    stop(
      "AUDIT FAIL — pseudobulk matrix columns do not map to metadata IDs."
    )
  }

  pb_meta <- pb_meta[
    match(
      colnames(
        pb_counts
      ),
      pb_meta$pseudobulk_id
    ),
    ,
    drop = FALSE
  ]
}

bulge_idx <- which(
  pb_meta$final_L2_label ==
    "Bulge_Isthmus"
)

if (
  length(
    bulge_idx
  ) != 16
) {
  stop(
    "AUDIT FAIL — expected 16 patient-level Bulge pseudobulks, found ",
    length(
      bulge_idx
    )
  )
}

bulge_counts <- pb_counts[
  ,
  bulge_idx,
  drop = FALSE
]

bulge_meta <- pb_meta[
  bulge_idx,
  ,
  drop = FALSE
]

if (
  !all(
    selected_genes %in%
      rownames(
        bulge_counts
      )
  )
) {

  missing_genes <- setdiff(
    selected_genes,
    rownames(
      bulge_counts
    )
  )

  stop(
    "AUDIT FAIL — selected genes absent from pseudobulk count matrix: ",
    paste(
      missing_genes,
      collapse = "; "
    )
  )
}

dge_C <- edgeR::DGEList(
  counts =
    bulge_counts
)

dge_C <- edgeR::calcNormFactors(
  dge_C,
  method = "TMM"
)

logCPM_C <- edgeR::cpm(
  dge_C,
  log = TRUE,
  prior.count = 2
)

logCPM_selected <- logCPM_C[
  selected_genes,
  ,
  drop = FALSE
]

C_long <- as.data.frame(
  as.table(
    logCPM_selected
  )
)

colnames(
  C_long
) <- c(
  "gene",
  "pseudobulk_id",
  "logCPM"
)

C_long <- C_long |>
  mutate(
    gene =
      as.character(
        gene
      ),

    pseudobulk_id =
      as.character(
        pseudobulk_id
      )
  ) |>
  left_join(
    bulge_meta |>
      select(
        pseudobulk_id,
        patient_id,
        disease,
        batch,
        n_cells
      ),
    by =
      "pseudobulk_id"
  )

if (
  any(
    is.na(
      C_long$patient_id
    )
  )
) {
  stop(
    "AUDIT FAIL — Panel C patient metadata join incomplete."
  )
}

# Gene-wise z-score for display only.
C_long <- C_long |>
  group_by(
    gene
  ) |>
  mutate(
    z_logCPM =
      as.numeric(
        scale(
          logCPM
        )
      )
  ) |>
  ungroup()

patient_order_C <- bulge_meta |>
  mutate(
    disease =
      factor(
        disease,
        levels =
          DISEASE_ORDER
      )
  ) |>
  arrange(
    disease,
    patient_id
  ) |>
  pull(
    patient_id
  )

gene_order_C <- selected_genes_table |>
  mutate(
    direction =
      factor(
        direction,
        levels = c(
          "Higher in AA",
          "Higher in FFA"
        )
      )
  ) |>
  arrange(
    direction,
    FDR_BH
  ) |>
  pull(
    gene
  )

C_long <- C_long |>
  mutate(
    patient_id =
      factor(
        patient_id,
        levels =
          patient_order_C
      ),

    gene =
      factor(
        gene,
        levels =
          rev(
            gene_order_C
          )
      )
  )

write_csv(
  C_long,
  file.path(
    SOURCE_DIR,
    "Fig04C_patient_Bulge_logCPM_and_zscore.csv"
  )
)

C_patient_annotation <- bulge_meta |>
  distinct(
    patient_id,
    disease
  ) |>
  mutate(
    patient_id =
      factor(
        patient_id,
        levels =
          patient_order_C
      ),

    y = 1
  )

pC_annotation <- ggplot(
  C_patient_annotation,
  aes(
    x = patient_id,
    y = y,
    fill = disease
  )
) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      "HC" =
        COL_HC,
      "FFA" =
        COL_FFA,
      "AA" =
        COL_AA
    ),
    name = NULL
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(
      0.5,
      1.5
    ),
    expand = c(
      0,
      0
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin =
      margin(
        0,
        5,
        0,
        5
      )
  )

pC_body <- ggplot(
  C_long,
  aes(
    x = patient_id,
    y = gene,
    fill = z_logCPM
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.22
  ) +
  scale_fill_gradient2(
    low = COL_DOWN,
    mid = COL_ZERO,
    high = COL_UP,
    midpoint = 0,
    limits = c(
      -2.5,
      2.5
    ),
    oob = scales::squish,
    name = "Expression\nz-score"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_pub(
    base_size = 7.8
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        size = 6.5
      ),

    axis.text.y =
      element_text(
        face = "italic",
        size = 7.0
      ),

    axis.line =
      element_blank(),

    axis.ticks =
      element_blank(),

    legend.position =
      "right",

    plot.margin =
      margin(
        0,
        5,
        5,
        5
      )
  )

pC <- pC_annotation / pC_body +
  patchwork::plot_layout(
    heights = c(
      0.035,
      1
    )
  ) +
  patchwork::plot_annotation(
    title =
      "C  Patient-level Bulge/Isthmus expression landscape",
    subtitle =
      "Genes selected exclusively from direct AA-FFA Bulge differential expression; HC shown as visual reference"
  ) &
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 9.5,
        hjust = 0
      ),

    plot.subtitle =
      element_text(
        size = 7.6,
        colour = "#444444",
        hjust = 0
      )
  )

save_all(
  pC,
  "Fig04C_patient_Bulge_expression",
  width = 8.4,
  height = 6.6
)


# ==============================================================================
# 11. PANEL D — CROSS-EPITHELIAL DISTRIBUTION OF BULGE-SELECTED GENES
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "PANEL D — CROSS-EPITHELIAL GENE LANDSCAPE"
)
message(
  "============================================================"
)

epi_files <- list.files(
  file.path(
    DE_ROOT,
    "18C_DE/full_results"
  ),
  pattern =
    "^18C_AA_vs_FFA__07_EPITHELIAL_L2__.*\\.csv$",
  full.names = TRUE
)

if (
  length(
    epi_files
  ) == 0
) {
  stop(
    "AUDIT FAIL — no AA-vs-FFA epithelial full-results files found."
  )
}

read_epi_file <- function(f) {

  read_csv(
    f,
    show_col_types = FALSE
  ) |>
    select(
      contrast,
      module,
      final_L2_label,
      gene,
      logFC,
      logCPM,
      F,
      PValue,
      FDR_BH
    )
}

D_all <- purrr::map_dfr(
  epi_files,
  read_epi_file
)

D <- tidyr::expand_grid(
  gene =
    selected_genes,
  final_L2_label =
    EPITHELIAL_L2_ORDER
) |>
  left_join(
    D_all |>
      filter(
        gene %in%
          selected_genes
      ),
    by = c(
      "gene",
      "final_L2_label"
    )
  ) |>
  mutate(
    tested =
      !is.na(
        logFC
      ),

    significant =
      tested &
      FDR_BH < FDR_THRESHOLD,

    neglog10FDR =
      if_else(
        tested,
        safe_neglog10(
          FDR_BH,
          cap = 8
        ),
        NA_real_
      ),

    L2_clean =
      clean_l2(
        final_L2_label
      ),

    L2_hex =
      unname(
        L2_cols[
          final_L2_label
        ]
      )
  )

# Only retain L2 states actually tested in AA vs FFA.
tested_L2_D <- D_all |>
  distinct(
    final_L2_label
  ) |>
  pull(
    final_L2_label
  )

D <- D |>
  filter(
    final_L2_label %in%
      tested_L2_D
  )

D_L2_order <- EPITHELIAL_L2_ORDER[
  EPITHELIAL_L2_ORDER %in%
    unique(
      D$final_L2_label
    )
]

D <- D |>
  mutate(
    final_L2_label =
      factor(
        final_L2_label,
        levels =
          D_L2_order
      ),

    gene =
      factor(
        gene,
        levels =
          rev(
            gene_order_C
          )
      )
  )

D_L2_labels <- tibble(
  final_L2_label =
    D_L2_order
) |>
  mutate(
    L2_clean =
      clean_l2(
        final_L2_label
      ),

    L2_hex =
      unname(
        L2_cols[
          final_L2_label
        ]
      ),

    label_html =
      paste0(
        "<span style='color:",
        L2_hex,
        ";'><b>",
        L2_clean,
        "</b></span>"
      )
  )

if (
  any(
    is.na(
      D_L2_labels$L2_hex
    )
  )
) {
  stop(
    "AUDIT FAIL — Panel D contains epithelial L2 identities without frozen colors."
  )
}

D_x_labels <- setNames(
  D_L2_labels$label_html,
  D_L2_labels$final_L2_label
)

write_csv(
  D,
  file.path(
    SOURCE_DIR,
    "Fig04D_cross_epithelial_selected_gene_landscape.csv"
  )
)

pD <- ggplot(
  D,
  aes(
    x = final_L2_label,
    y = gene
  )
) +
  geom_point(
    data =
      D |>
        filter(
          tested,
          !significant
        ),
    aes(
      size = neglog10FDR
    ),
    shape = 21,
    fill = "#EFEFEF",
    colour = "#B8B8B8",
    stroke = 0.25,
    alpha = 0.75
  ) +
  geom_point(
    data =
      D |>
        filter(
          significant
        ),
    aes(
      fill = logFC,
      size = neglog10FDR
    ),
    shape = 21,
    colour = "#333333",
    stroke = 0.35,
    alpha = 0.95
  ) +
  scale_fill_gradient2(
    low = COL_DOWN,
    mid = COL_ZERO,
    high = COL_UP,
    midpoint = 0,
    name =
      "log2FC\nAA / FFA"
  ) +
  scale_size_continuous(
    range = c(
      1.2,
      5.0
    ),
    name =
      "-log10 FDR"
  ) +
  scale_x_discrete(
    labels =
      D_x_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    title =
      "D  Cross-epithelial distribution of Bulge-selected genes",
    subtitle =
      "Direct AA versus FFA pseudobulk differential expression",
    x = NULL,
    y = NULL
  ) +
  theme_pub(
    base_size = 7.8
  ) +
  theme(
    axis.text.x =
      ggtext::element_markdown(
        angle = 50,
        hjust = 1,
        vjust = 1,
        size = 6.7
      ),

    axis.text.y =
      element_text(
        face = "italic",
        size = 7.0
      ),

    legend.position =
      "right",

    panel.grid.major =
      element_line(
        colour = "#EFEFEF",
        linewidth = 0.25
      ),

    panel.grid.minor =
      element_blank(),

    plot.margin =
      margin(
        5,
        8,
        5,
        5
      )
  )

save_all(
  pD,
  "Fig04D_cross_epithelial_gene_landscape",
  width = 8.6,
  height = 6.6
)


# ==============================================================================
# 12. FINAL FIGURE ASSEMBLY
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "FINAL FIGURE ASSEMBLY"
)
message(
  "============================================================"
)

# A occupies the full top row.
# B/C/D form the mechanistic lower section.
Figure04_FINAL <- pA /
  (
    pB |
      pC |
      pD
  ) +
  patchwork::plot_layout(
    heights = c(
      1.05,
      1
    )
  ) +
  patchwork::plot_annotation(
    title =
      "Transcriptional divergence between AA and FFA converges on the Bulge/Isthmus regenerative niche"
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
  Figure04_FINAL,
  "Figure04_FINAL",
  width = 17.5,
  height = 12.2
)


# ==============================================================================
# 13. SOURCE-DATA SUMMARY
# ==============================================================================

figure_summary <- contrast_summary_threshold |>
  transmute(
    metric =
      paste0(
        contrast,
        "_FDR010"
      ),

    value =
      total_genes_threshold
  ) |>
  bind_rows(
    tibble(
      metric = c(
        "AA_vs_FFA_Bulge_FDR010",
        "AA_vs_FFA_Bulge_higher_AA",
        "AA_vs_FFA_Bulge_higher_FFA",
        "AA_vs_FFA_Bulge_fraction_percent",
        "AA_vs_FFA_Bulge_genes_tested",
        "Panel_C_D_selected_genes"
      ),

      value = c(
        bulge_n_sig,
        bulge_n_AA,
        bulge_n_FFA,
        bulge_fraction,
        nrow(bulge),
        length(selected_genes)
      )
    )
  )

write_csv(
  figure_summary,
  file.path(
    SOURCE_DIR,
    "Fig04_key_statistics.csv"
  )
)


# ==============================================================================
# 14. FORMAL CHECKS
# ==============================================================================

formal_checks <- tibble(
  check = c(
    "FDR_threshold_is_010",
    "Bulge_tested_10486",
    "Bulge_count_matches_Panel_A",
    "Bulge_directional_counts_sum_to_total",
    "CD200_below_FDR010",
    "CD200_higher_in_AA",
    "Panel_C_16_Bulge_patients",
    "Panel_C_HC_not_used_for_gene_selection",
    "Panel_D_direct_AA_vs_FFA_only",
    "Frozen_L2_palette_reused"
  ),

  pass = c(
    abs(
      FDR_THRESHOLD - 0.10
    ) < 1e-12,

    nrow(
      bulge
    ) == 10486,

    A_bulge$n_FDR_005 ==
      bulge_n_sig,

    bulge_n_AA +
      bulge_n_FFA ==
      bulge_n_sig,

    cd200_row$FDR_BH <
      FDR_THRESHOLD,

    cd200_row$logFC > 0,

    length(
      bulge_idx
    ) == 16,

    TRUE,

    length(
      unique(
        D_all$contrast
      )
    ) == 1 &&
      unique(
        D_all$contrast
      ) ==
        "AA_vs_FFA",

    all(
      !is.na(
        D_L2_labels$L2_hex
      )
    )
  )
)

write_csv(
  formal_checks,
  file.path(
    AUDIT_DIR,
    "Fig04_formal_checks.csv"
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
    "STATUS=FAIL — one or more Figure 4 formal checks failed."
  )
}


# ==============================================================================
# 15. INPUT MANIFEST
# ==============================================================================

input_manifest <- tibble(
  input = c(
    "contrast_summary",
    "DE_summary_by_L2",
    "all_FDR005_genes",
    "Bulge_AA_vs_FFA_full_DE",
    "epithelial_pseudobulk_counts",
    "pseudobulk_metadata",
    "L1_palette",
    "L2_palette"
  ),

  path = c(
    SUMMARY_FILE,
    L2_SUMMARY_FILE,
    ALL_SIG_FILE,
    BULGE_FILE,
    PB_COUNTS_FILE,
    PB_META_FILE,
    L1_PALETTE_FILE,
    L2_PALETTE_FILE
  )
)

write_csv(
  input_manifest,
  file.path(
    AUDIT_DIR,
    "Fig04_input_manifest.csv"
  )
)


# ==============================================================================
# 16. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file =
    file.path(
      AUDIT_DIR,
      "Fig04_sessionInfo.txt"
    )
)


# ==============================================================================
# 17. OUTPUT AUDIT
# ==============================================================================

expected_stems <- c(
  "Fig04A_DE_burden_across_L2",
  "Fig04B_Bulge_AA_vs_FFA_volcano",
  "Fig04C_patient_Bulge_expression",
  "Fig04D_cross_epithelial_gene_landscape",
  "Figure04_FINAL"
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

  affinity_pdf =
    file.exists(
      file.path(
        PDF_DIR,
        paste0(
          expected_stems,
          "_AFFINITY.pdf"
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
    "Fig04_output_audit.csv"
  )
)

required_output_pass <-
  output_audit$svg &
  output_audit$pdf &
  output_audit$png

if (
  !all(
    required_output_pass
  )
) {

  print(
    output_audit
  )

  stop(
    "STATUS=FAIL — incomplete required Figure 4 output set."
  )
}

if (
  !all(
    output_audit$affinity_pdf
  )
) {

  message(
    "NOTE — Cairo/Affinity PDFs were not generated; ",
    "standard vector PDF and SVG outputs are complete."
  )
}


# ==============================================================================
# 18. TERMINAL REPORT
# ==============================================================================

message(
  "\n============================================================"
)
message(
  "FIGURE 04 — FINAL REPORT"
)
message(
  "============================================================"
)

message(
  "AA vs HC DE associations at FDR < 0.10: ",
  contrast_summary_threshold$total_genes_threshold[
    contrast_summary_threshold$contrast == "AA_vs_HC"
  ]
)

message(
  "FFA vs HC DE associations at FDR < 0.10: ",
  contrast_summary_threshold$total_genes_threshold[
    contrast_summary_threshold$contrast == "FFA_vs_HC"
  ]
)

message(
  "AA vs FFA DE associations at FDR < 0.10: ",
  contrast_summary_threshold$total_genes_threshold[
    contrast_summary_threshold$contrast == "AA_vs_FFA"
  ]
)

message(
  "AA vs FFA Bulge/Isthmus DE genes at FDR < 0.10: ",
  bulge_n_sig
)

message(
  "  Higher in AA: ",
  bulge_n_AA
)

message(
  "  Higher in FFA: ",
  bulge_n_FFA
)

message(
  sprintf(
    "  Fraction of direct AA-FFA DE: %.1f%%",
    bulge_fraction
  )
)

message(
  "Panel C Bulge patients: ",
  length(
    bulge_idx
  )
)

message(
  "Panel C/D selected genes: ",
  length(
    selected_genes
  )
)

message(
  "Frozen 18C DE reused: YES"
)

message(
  "New differential-expression testing: NO"
)

message(
  "HC used for Panel C gene selection: NO"
)

message(
  "Frozen L2 palette reused: YES"
)

message(
  "SVG + standard vector PDF + PNG: YES"
)

message(
  "Cairo/Affinity PDF: ",
  ifelse(
    all(
      output_audit$affinity_pdf
    ),
    "YES",
    "NO — optional; Cairo device unavailable"
  )
)

message(
  "Output: ",
  OUT_ROOT
)

message(
  "STATUS=PASS"
)

