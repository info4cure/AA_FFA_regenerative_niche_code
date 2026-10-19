#!/usr/bin/env Rscript

# ==============================================================================
# FIGS09 — ORTHOGONAL PATIENT-LEVEL GSVA VALIDATION OF HALLMARK GSEA
#
# a  Significant GSVA Hallmark landscape
# b  GSEA NES versus patient-level GSVA effect
# c  Quantitative concordance summary
#
# Frozen source: 20C_PATIENT_HALLMARK_GSVA
# No new statistical inference.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(svglite)
  library(grid)
})

# ==============================================================================
# 1. PATHS
# ==============================================================================

BASE <- file.path(
  "PIPELINE_V2/results/disease_biology",
  "20_PATHWAY_ANALYSIS",
  "20C_PATIENT_HALLMARK_GSVA",
  "tables"
)

OUT <- file.path(
  "PIPELINE_V2/results/publication_figures",
  "SUPPLEMENTARY",
  "FIGS09_GSVA_GSEA_VALIDATION"
)

SVG_DIR <- file.path(OUT, "SVG")
PDF_DIR <- file.path(OUT, "PDF")
PNG_DIR <- file.path(OUT, "PNG")
TAB_DIR <- file.path(OUT, "TABLES")

for (d in c(
  SVG_DIR,
  PDF_DIR,
  PNG_DIR,
  TAB_DIR
)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

GLOBAL_FILE <- file.path(
  BASE,
  "20C_global_GSVA_summary.csv"
)

CONCORD_FILE <- file.path(
  BASE,
  "20C_GSVA_GSEA_concordance_summary.csv"
)

PAIRS_FILE <- file.path(
  BASE,
  "20C_GSVA_GSEA_pairwise_concordance.csv"
)

DIFF_FILE <- file.path(
  BASE,
  "20C_Hallmark_GSVA_differential_scores.csv"
)

L2SUM_FILE <- file.path(
  BASE,
  "20C_L2_GSVA_summary.csv"
)

input_files <- c(
  GLOBAL_FILE,
  CONCORD_FILE,
  PAIRS_FILE,
  DIFF_FILE,
  L2SUM_FILE
)

if (!all(file.exists(input_files))) {
  stop(
    "Missing one or more frozen 20C source files."
  )
}

# ==============================================================================
# 2. LOAD FROZEN TABLES
# ==============================================================================

global <- read_csv(
  GLOBAL_FILE,
  show_col_types=FALSE
)

concord <- read_csv(
  CONCORD_FILE,
  show_col_types=FALSE
)

pairs <- read_csv(
  PAIRS_FILE,
  show_col_types=FALSE
)

diff <- read_csv(
  DIFF_FILE,
  show_col_types=FALSE
)

l2sum <- read_csv(
  L2SUM_FILE,
  show_col_types=FALSE
)

cat("\n============================================================\n")
cat("FIGS09 — GSVA / GSEA ORTHOGONAL VALIDATION\n")
cat("============================================================\n")
cat("Differential GSVA rows :", nrow(diff), "\n")
cat("GSEA-GSVA pairs        :", nrow(pairs), "\n")
cat("L2 summary rows        :", nrow(l2sum), "\n")

# ==============================================================================
# 3. FROZEN L2 PALETTE
# ==============================================================================

L2_PALETTE_FILE <- file.path(
  "PIPELINE_V2/config/figure_palettes",
  "L2_palette.tsv"
)

if (!file.exists(L2_PALETTE_FILE)) {
  stop(
    "Missing frozen L2 palette."
  )
}

L2_palette <- read.delim(
  L2_PALETTE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  filter(
    palette_status == "FROZEN"
  ) |>
  arrange(
    global_display_order
  )

L2_cols <- setNames(
  L2_palette$hex,
  L2_palette$final_L2_label
)

# ==============================================================================
# 4. FROZEN NUMERICAL AUDIT
# ==============================================================================

expected_global <- tibble(
  contrast=c(
    "AA_vs_FFA",
    "AA_vs_HC",
    "FFA_vs_HC"
  ),
  expected_L2=c(
    0,
    9,
    6
  ),
  expected_hits=c(
    0,
    33,
    35
  )
)

global_audit <- global |>
  left_join(
    expected_global,
    by="contrast"
  ) |>
  mutate(
    L2_match=
      L2_with_FDR005 ==
      expected_L2,

    hits_match=
      total_pathway_hits_FDR005 ==
      expected_hits
  )

expected_concord <- tibble(
  contrast=c(
    "AA_vs_FFA",
    "AA_vs_HC",
    "FFA_vs_HC"
  ),
  expected_rho=c(
    0.754,
    0.836,
    0.902
  ),
  expected_dir_GSEA_sig=c(
    0.966,
    0.984,
    1.000
  )
)

concord_audit <- concord |>
  left_join(
    expected_concord,
    by="contrast"
  ) |>
  mutate(
    rho_match=
      round(
        spearman_effect_NES,
        3
      ) ==
      expected_rho,

    dir_match=
      round(
        direction_concordance_GSEA_FDR005,
        3
      ) ==
      expected_dir_GSEA_sig
  )

cat("\n===== GLOBAL GSVA AUDIT =====\n")
print(
  global_audit,
  n=Inf
)

cat("\n===== GSEA-GSVA CONCORDANCE AUDIT =====\n")
print(
  concord_audit,
  n=Inf
)

if (
  !all(global_audit$L2_match) ||
  !all(global_audit$hits_match) ||
  !all(concord_audit$rho_match) ||
  !all(concord_audit$dir_match)
) {
  stop(
    "Frozen numerical audit FAILED."
  )
}

# ==============================================================================
# 5. DISPLAY HELPERS
# ==============================================================================

contrast_labels <- c(
  AA_vs_HC="AA vs HC",
  FFA_vs_HC="FFA vs HC",
  AA_vs_FFA="AA vs FFA"
)

contrast_order <- c(
  "AA vs HC",
  "FFA vs HC",
  "AA vs FFA"
)

clean_pathway <- function(x) {

  x |>
    str_remove("^HALLMARK_") |>
    str_replace_all("_", " ") |>
    str_to_title()
}

# ==============================================================================
# 6. PANEL A — DIRECTIONAL BUBBLE LANDSCAPE
# ==============================================================================

A <- diff |>
  filter(
    contrast %in%
      c(
        "AA_vs_HC",
        "FFA_vs_HC"
      ),
    FDR_BH < 0.05
  ) |>
  mutate(
    disease=case_when(
      contrast == "AA_vs_HC" ~ "AA",
      contrast == "FFA_vs_HC" ~ "FFA",
      TRUE ~ NA_character_
    ),

    direction=case_when(
      delta_GSVA > 0 ~ "Up",
      delta_GSVA < 0 ~ "Down",
      TRUE ~ "Zero"
    ),

    x_group=case_when(
      disease == "AA" &
        direction == "Down" ~ "AA down",

      disease == "AA" &
        direction == "Up" ~ "AA up",

      disease == "FFA" &
        direction == "Down" ~ "FFA down",

      disease == "FFA" &
        direction == "Up" ~ "FFA up",

      TRUE ~ NA_character_
    ),

    pathway_display=
      clean_pathway(pathway),

    sig_size=
      -log10(
        pmax(
          FDR_BH,
          1e-12
        )
      )
  ) |>
  filter(
    !is.na(x_group)
  )

if (nrow(A) != 68) {
  stop(
    "AUDIT FAIL - expected 68 significant GSVA associations; found ",
    nrow(A)
  )
}

if (
  sum(A$disease == "AA") != 33 ||
  sum(A$disease == "FFA") != 35
) {
  stop(
    "AUDIT FAIL - expected AA=33 and FFA=35 GSVA hits."
  )
}

missing_L2 <- setdiff(
  unique(A$final_L2_label),
  names(L2_cols)
)

if (length(missing_L2) > 0) {
  stop(
    "Missing frozen L2 colors for: ",
    paste(
      missing_L2,
      collapse="; "
    )
  )
}

A <- A |>
  mutate(
    L2_color=
      unname(
        L2_cols[
          final_L2_label
        ]
      )
  )

x_levels <- c(
  "AA down",
  "AA up",
  "FFA down",
  "FFA up"
)

A <- A |>
  mutate(
    x_group=factor(
      x_group,
      levels=x_levels
    )
  )

# Order pathways by:
# 1) number of significant disease/L2 associations
# 2) maximum significance
path_order <- A |>
  group_by(
    pathway_display
  ) |>
  summarise(
    n_hits=n(),
    max_sig=max(
      sig_size,
      na.rm=TRUE
    ),
    .groups="drop"
  ) |>
  arrange(
    n_hits,
    max_sig
  ) |>
  pull(
    pathway_display
  )

A <- A |>
  mutate(
    pathway_display=factor(
      pathway_display,
      levels=path_order
    )
  )

# Deterministic horizontal offsets by frozen L2 order.
A_L2_order <- c(
  # Non-follicular epithelium
  "Basal_epithelial",
  "IFE_Infundibulum_suprabasal",
  "IFE_Infundibulum_terminal",

  # Follicular epithelium
  "Bulge_Isthmus",
  "IRS_companion",
  "ORS",

  # Sebaceous lineage
  "Early sebocyte commitment",

  # Sweat gland
  "Clear secretory cell",
  "Dark/apocrine-like secretory cell",
  "Ductal/clear transitional cell",

  # Fibroblast
  "PTGDS/IGF1 fibroblast",

  # Immune
  "Macrophage/monocyte",

  # Mural / vascular
  "Contractile mural cell",
  "SLC10A1/CD83 endothelial state"
)

# Keep only states actually represented in Panel A.
A_L2_order <- A_L2_order[
  A_L2_order %in%
    unique(A$final_L2_label)
]

# Audit completeness: every Panel-A L2 must occur exactly once.
missing_ordered_L2 <- setdiff(
  unique(A$final_L2_label),
  A_L2_order
)

if (length(missing_ordered_L2) > 0) {
  stop(
    "AUDIT FAIL - biological L2 order missing: ",
    paste(
      missing_ordered_L2,
      collapse="; "
    )
  )
}

if (anyDuplicated(A_L2_order)) {
  stop(
    "AUDIT FAIL - duplicated L2 identity in biological order."
  )
}

offset_values <- seq(
  -0.28,
  0.28,
  length.out=length(A_L2_order)
)

L2_offsets <- setNames(
  offset_values,
  A_L2_order
)

A <- A |>
  mutate(
    x_base=
      as.numeric(
        x_group
      ),

    x_plot=
      x_base +
      unname(
        L2_offsets[
          final_L2_label
        ]
      )
  )

# x-axis labels use arrows to make direction immediately explicit.
x_breaks <- 1:4

x_labels <- c(
  "AA lower",
  "AA higher",
  "FFA lower",
  "FFA higher"
)

pA <- ggplot(
  A,
  aes(
    x=x_plot,
    y=pathway_display
  )
) +
  geom_vline(
    xintercept=2.5,
    linewidth=0.8,
    colour="#303030"
  ) +
  geom_vline(
    xintercept=c(
      1.5,
      3.5
    ),
    linewidth=0.35,
    linetype=2,
    colour="#808080"
  ) +
  geom_segment(
    aes(
      x=x_plot,
      xend=x_plot,

      y=as.numeric(pathway_display) +
        if_else(
          direction == "Up",
          -0.18,
          0.18
        ),

      yend=as.numeric(pathway_display) +
        if_else(
          direction == "Up",
          0.18,
          -0.18
        ),

      colour=final_L2_label,
      linewidth=sig_size
    ),

    arrow=arrow(
      type="closed",
      length=unit(
        0.075,
        "inches"
      )
    ),

    lineend="round",
    alpha=0.95
  ) +
  scale_colour_manual(
    values=L2_cols,
    breaks=A_L2_order,
    name="L2 state",
    guide=guide_legend(
      override.aes=list(
        linewidth=2.2,
        alpha=1
      ),
      keyheight=unit(
        0.42,
        "cm"
      ),
      keywidth=unit(
        0.65,
        "cm"
      ),
      ncol=2,
      byrow=TRUE
    )
  ) +
  scale_linewidth_continuous(
    range=c(
      0.45,
      1.35
    ),
    name=expression(-log[10]*" FDR"),
    guide=guide_legend(
      keyheight=unit(
        0.40,
        "cm"
      ),
      keywidth=unit(
        0.75,
        "cm"
      )
    )
  ) +
  scale_x_continuous(
    breaks=x_breaks,
    labels=x_labels,
    limits=c(
      0.55,
      4.45
    ),
    expand=c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    breaks=seq_along(
      levels(A$pathway_display)
    ),
    labels=levels(
      A$pathway_display
    ),
    limits=c(
      0.5,
      length(
        levels(A$pathway_display)
      ) + 0.5
    ),
    expand=c(
      0,
      0
    )
  ) +
  labs(
    title=
      "a  Patient-level GSVA resolves directional Hallmark activity across AA and FFA",
    subtitle=
      "Each arrow represents one L2-Hallmark association at FDR < 0.05; direction indicates lower or higher activity versus HC; AA vs FFA: 0 significant associations",
    x=NULL,
    y=NULL
  ) +
  theme_bw(
    base_size=8.5
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=8
    ),
    axis.text.x=element_text(
      face="bold",
      size=8.5
    ),
    axis.text.y=element_text(
      size=6.6
    ),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank(),
    legend.position="right"
  )

# ==============================================================================
# 7. PANEL B — GSEA NES vs GSVA EFFECT
# ==============================================================================

B <- pairs |>
  filter(
    is.finite(NES),
    is.finite(delta_GSVA)
  ) |>
  mutate(
    contrast_display=recode(
      contrast,
      !!!contrast_labels
    ),

    contrast_display=factor(
      contrast_display,
      levels=contrast_order
    ),

    GSEA_sig=
      !is.na(padj) &
      padj < 0.05
  )

if (nrow(B) != sum(concord$n_pairs)) {
  stop(
    "AUDIT FAIL - pair count differs from concordance summary."
  )
}

B_labels <- concord |>
  mutate(
    contrast_display=recode(
      contrast,
      !!!contrast_labels
    ),

    contrast_display=factor(
      contrast_display,
      levels=contrast_order
    ),

    label=paste0(
      "rho = ",
      sprintf(
        "%.3f",
        spearman_effect_NES
      ),
      "\n",
      "n = ",
      n_pairs
    )
  )

pB <- ggplot(
  B,
  aes(
    x=NES,
    y=delta_GSVA
  )
) +
  geom_hline(
    yintercept=0,
    linewidth=0.35,
    colour="#606060"
  ) +
  geom_vline(
    xintercept=0,
    linewidth=0.35,
    colour="#606060"
  ) +
  geom_point(
    aes(
      shape=GSEA_sig
    ),
    alpha=0.40,
    size=1.35,
    stroke=0.25
  ) +
  geom_smooth(
    method="lm",
    formula=y~x,
    se=FALSE,
    linewidth=0.75
  ) +
  geom_text(
    data=B_labels,
    aes(
      x=-Inf,
      y=Inf,
      label=label
    ),
    inherit.aes=FALSE,
    hjust=-0.08,
    vjust=1.12,
    size=2.8,
    fontface="bold"
  ) +
  facet_wrap(
    ~contrast_display,
    nrow=1
  ) +
  scale_shape_manual(
    values=c(
      `FALSE`=1,
      `TRUE`=16
    ),
    labels=c(
      "GSEA FDR >= 0.05",
      "GSEA FDR < 0.05"
    ),
    name=NULL
  ) +
  labs(
    title=
      "b  Patient-level GSVA effects reproduce the GSEA enrichment architecture",
    subtitle=
      "Each point represents an evaluable L2-Hallmark pair",
    x="GSEA normalized enrichment score (NES)",
    y=expression(Delta*"GSVA")
  ) +
  theme_bw(
    base_size=8.5
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=8
    ),
    strip.text=element_text(
      face="bold",
      size=8.5
    ),
    legend.position="bottom",
    panel.grid.minor=element_blank()
  )

# ==============================================================================
# 8. PANEL C — QUANTITATIVE VALIDATION SUMMARY
# ==============================================================================

C <- concord |>
  left_join(
    global |>
      select(
        contrast,
        L2_with_FDR005,
        total_pathway_hits_FDR005
      ),
    by="contrast"
  ) |>
  mutate(
    contrast_display=recode(
      contrast,
      !!!contrast_labels
    ),

    contrast_display=factor(
      contrast_display,
      levels=contrast_order
    ),

    GSEA_sig_direction_pct=
      100 *
      direction_concordance_GSEA_FDR005,

    GSVA_label=paste0(
      total_pathway_hits_FDR005,
      " GSVA hits / ",
      L2_with_FDR005,
      " L2"
    )
  )

C_long <- bind_rows(
  C |>
    transmute(
      contrast_display,
      metric="Spearman rho",
      value=spearman_effect_NES
    ),

  C |>
    transmute(
      contrast_display,
      metric="Directional concordance",
      value=
        direction_concordance_GSEA_FDR005
    )
) |>
  mutate(
    metric=factor(
      metric,
      levels=c(
        "Spearman rho",
        "Directional concordance"
      )
    )
  )

pC <- ggplot(
  C_long,
  aes(
    x=contrast_display,
    y=value,
    shape=metric
  )
) +
  geom_point(
    size=4.2,
    stroke=0.8
  ) +
  geom_text(
    data=C,
    aes(
      x=contrast_display,
      y=0.08,
      label=GSVA_label
    ),
    inherit.aes=FALSE,
    size=3.0,
    fontface="bold"
  ) +
  scale_shape_manual(
    values=c(
      "Spearman rho"=16,
      "Directional concordance"=17
    )
  ) +
  scale_y_continuous(
    limits=c(
      0,
      1.05
    ),
    breaks=c(
      0,
      0.25,
      0.50,
      0.75,
      1.00
    )
  ) +
  labs(
    title=
      "c  Orthogonal validation remains directionally concordant across contrasts",
    subtitle=
      "Directional concordance is calculated among Hallmark associations significant by GSEA",
    x=NULL,
    y="Concordance",
    shape=NULL
  ) +
  theme_bw(
    base_size=9
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=8
    ),
    axis.text.x=element_text(
      face="bold",
      size=8.5
    ),
    axis.text.y=element_text(
      size=8
    ),
    legend.position="bottom",
    legend.direction="horizontal",
    panel.grid.minor=element_blank()
  )

# ==============================================================================
# 9. PANEL-SPECIFIC AUDIT
# ==============================================================================

cat("\n===== PANEL A AUDIT =====\n")

A_audit <- A |>
  count(
    disease,
    direction,
    name="n_hits"
  ) |>
  arrange(
    disease,
    direction
  )

print(
  A_audit,
  n=Inf
)

cat("\n===== PANEL B AUDIT =====\n")

B_audit <- B |>
  count(
    contrast_display,
    name="n_pairs"
  )

print(
  B_audit,
  n=Inf
)

cat("\n===== PANEL C AUDIT =====\n")

print(
  C |>
    select(
      contrast,
      n_pairs,
      spearman_effect_NES,
      direction_concordance,
      n_GSEA_FDR005,
      direction_concordance_GSEA_FDR005,
      L2_with_FDR005,
      total_pathway_hits_FDR005
    ),
  n=Inf,
  width=Inf
)

# ==============================================================================
# 10. FINAL ASSEMBLY
# ==============================================================================

final <- wrap_plots(
  list(
    pA,
    pB,
    pC
  ),
  ncol=1,
  nrow=3,
  heights=c(
    1.45,
    1.00,
    0.72
  ),
  guides="keep"
)

# ==============================================================================
# 11. EXPORT
# ==============================================================================

stem <- "FigS09_GSVA_GSEA_validation"

svg_file <- file.path(
  SVG_DIR,
  paste0(
    stem,
    ".svg"
  )
)

pdf_file <- file.path(
  PDF_DIR,
  paste0(
    stem,
    ".pdf"
  )
)

png_file <- file.path(
  PNG_DIR,
  paste0(
    stem,
    ".png"
  )
)

ggsave(
  svg_file,
  final,
  width=11.5,
  height=16.0,
  units="in",
  device=svglite
)

ggsave(
  pdf_file,
  final,
  width=11.5,
  height=16.0,
  units="in"
)

ggsave(
  png_file,
  final,
  width=11.5,
  height=16.0,
  units="in",
  dpi=400
)

# ==============================================================================
# 12. SOURCE TABLES
# ==============================================================================

write_csv(
  A,
  file.path(
    TAB_DIR,
    "FigS09A_significant_GSVA_associations.csv"
  )
)

write_csv(
  B,
  file.path(
    TAB_DIR,
    "FigS09B_GSEA_GSVA_pairs.csv"
  )
)

write_csv(
  C,
  file.path(
    TAB_DIR,
    "FigS09C_validation_summary.csv"
  )
)

write_csv(
  global_audit,
  file.path(
    TAB_DIR,
    "FigS09_global_GSVA_audit.csv"
  )
)

write_csv(
  concord_audit,
  file.path(
    TAB_DIR,
    "FigS09_concordance_audit.csv"
  )
)

# ==============================================================================
# 13. FORMAL CHECKS
# ==============================================================================

expected_pair_counts <- c(
  "AA vs HC"=1169,
  "FFA vs HC"=886,
  "AA vs FFA"=861
)

observed_pair_counts <- B |>
  count(
    contrast_display
  ) |>
  mutate(
    expected=
      unname(
        expected_pair_counts[
          as.character(
            contrast_display
          )
        ]
      ),
    match=
      n == expected
  )

checks <- tibble(
  check=c(
    "3000 frozen GSVA differential rows",
    "2916 evaluable GSEA-GSVA pairs",
    "68 GSVA associations at FDR <0.05 in disease-vs-HC contrasts",
    "33 AA-vs-HC GSVA hits",
    "35 FFA-vs-HC GSVA hits",
    "0 AA-vs-FFA GSVA hits",
    "9 AA-vs-HC L2 identities with FDR<0.05",
    "6 FFA-vs-HC L2 identities with FDR<0.05",
    "AA-vs-HC GSEA-GSVA rho = 0.836",
    "FFA-vs-HC GSEA-GSVA rho = 0.902",
    "AA-vs-FFA GSEA-GSVA rho = 0.754",
    "AA-vs-HC GSEA-significant directional concordance = 98.4%",
    "FFA-vs-HC GSEA-significant directional concordance = 100%",
    "AA-vs-FFA GSEA-significant directional concordance = 96.6%",
    "Contrast-specific pair counts reproduce frozen summary",
    "All significant Panel A L2 states resolve to frozen colors",
    "SVG exists",
    "PDF exists",
    "PNG exists"
  ),

  pass=c(
    nrow(diff) == 3000,

    nrow(B) == 2916,

    nrow(A) == 68,

    sum(
      A$contrast ==
        "AA_vs_HC"
    ) == 33,

    sum(
      A$contrast ==
        "FFA_vs_HC"
    ) == 35,

    global |>
      filter(
        contrast ==
          "AA_vs_FFA"
      ) |>
      pull(
        total_pathway_hits_FDR005
      ) ==
      0,

    global |>
      filter(
        contrast ==
          "AA_vs_HC"
      ) |>
      pull(
        L2_with_FDR005
      ) ==
      9,

    global |>
      filter(
        contrast ==
          "FFA_vs_HC"
      ) |>
      pull(
        L2_with_FDR005
      ) ==
      6,

    concord |>
      filter(
        contrast ==
          "AA_vs_HC"
      ) |>
      pull(
        spearman_effect_NES
      ) |>
      round(3) ==
      0.836,

    concord |>
      filter(
        contrast ==
          "FFA_vs_HC"
      ) |>
      pull(
        spearman_effect_NES
      ) |>
      round(3) ==
      0.902,

    concord |>
      filter(
        contrast ==
          "AA_vs_FFA"
      ) |>
      pull(
        spearman_effect_NES
      ) |>
      round(3) ==
      0.754,

    concord |>
      filter(
        contrast ==
          "AA_vs_HC"
      ) |>
      pull(
        direction_concordance_GSEA_FDR005
      ) |>
      round(3) ==
      0.984,

    concord |>
      filter(
        contrast ==
          "FFA_vs_HC"
      ) |>
      pull(
        direction_concordance_GSEA_FDR005
      ) |>
      round(3) ==
      1.000,

    concord |>
      filter(
        contrast ==
          "AA_vs_FFA"
      ) |>
      pull(
        direction_concordance_GSEA_FDR005
      ) |>
      round(3) ==
      0.966,

    all(
      observed_pair_counts$match
    ),

    all(
      !is.na(
        A$L2_color
      )
    ),

    file.exists(
      svg_file
    ),

    file.exists(
      pdf_file
    ),

    file.exists(
      png_file
    )
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  checks,
  n=Inf
)

if (
  nrow(checks) !=
  length(checks$pass)
) {
  stop(
    "Malformed formal-check table."
  )
}

if (!all(checks$pass)) {
  stop(
    "FIGS09 formal checks FAILED."
  )
}

cat(
  "FORMAL_CHECKS=PASS\n"
)

# ==============================================================================
# 14. MANIFEST
# ==============================================================================

manifest <- tibble(
  figure=
    "Supplementary Figure S9",

  title=
    "Orthogonal patient-level GSVA validation of Hallmark GSEA architecture",

  panel_layout=
    "3 panels; 1 column x 3 rows",

  panel_A=
    "Directional significant GSVA bubble landscape",

  panel_B=
    "GSEA NES versus patient-level GSVA effect",

  panel_C=
    "Quantitative GSEA-GSVA validation summary",

  differential_rows=
    nrow(diff),

  evaluable_pairs=
    nrow(B),

  significant_GSVA_hits=
    nrow(A),

  svg=
    svg_file,

  pdf=
    pdf_file,

  png=
    png_file
)

write_csv(
  manifest,
  file.path(
    TAB_DIR,
    "FigS09_manifest.csv"
  )
)

# ==============================================================================
# 15. FINAL STATUS
# ==============================================================================

cat("\n============================================================\n")
cat("FIGS09 STATUS = PASS\n")
cat("============================================================\n")

cat(
  "GSVA hits AA vs HC   :",
  sum(
    A$contrast ==
      "AA_vs_HC"
  ),
  "\n"
)

cat(
  "GSVA hits FFA vs HC  :",
  sum(
    A$contrast ==
      "FFA_vs_HC"
  ),
  "\n"
)

cat(
  "GSVA hits AA vs FFA  :",
  global |>
    filter(
      contrast ==
        "AA_vs_FFA"
    ) |>
    pull(
      total_pathway_hits_FDR005
    ),
  "\n"
)

cat(
  "GSEA-GSVA pairs      :",
  nrow(B),
  "\n"
)

cat(
  "SVG                  :",
  svg_file,
  "\n"
)

cat(
  "PDF                  :",
  pdf_file,
  "\n"
)

cat(
  "PNG                  :",
  png_file,
  "\n"
)

cat("============================================================\n")

