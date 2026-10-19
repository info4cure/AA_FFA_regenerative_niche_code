#!/usr/bin/env Rscript

# ==============================================================================
# FigS07_Mesenchymal_Bulge_coupling.R
#
# Supplementary Figure S7
# Mesenchymal-vascular niche states covary with Bulge/Isthmus organization
# despite limited compositional remodeling
#
# READ-ONLY publication module.
# No new statistical inference.
# All statistics are read from frozen 22G3A outputs.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggh4x)
  library(svglite)
  library(ggtext)
})

# ==============================================================================
# 1. PATHS
# ==============================================================================

INPUT_RDS <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "22_MULTICOMPARTMENT_DISEASE_AXES/",
  "22G_BULGE_CENTERED_ECOSYSTEM/",
  "22G3_MESENCHYMAL_VASCULAR_NICHE/",
  "22G3A_PATIENT_MODELLING/objects/",
  "22G3A_Bulge_mesenchymal_vascular_niche_freeze.rds"
)

OUT_DIR <- paste0(
  "PIPELINE_V2/results/publication_figures/",
  "SUPPLEMENTARY/FIGS07_MESENCHYMAL_BULGE_COUPLING"
)

SVG_DIR <- file.path(OUT_DIR, "SVG")
PDF_DIR <- file.path(OUT_DIR, "PDF")
PNG_DIR <- file.path(OUT_DIR, "PNG")
TAB_DIR <- file.path(OUT_DIR, "TABLES")

dir.create(SVG_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PDF_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PNG_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TAB_DIR, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(INPUT_RDS)) {
  stop("Missing frozen input: ", INPUT_RDS)
}

# ==============================================================================
# 2. LOAD FROZEN OBJECT
# ==============================================================================


# ==============================================================================
# FROZEN L2 PALETTE
# Canonical atlas palette used by L2 UMAPs/publication figures
# ==============================================================================

L2_PALETTE_FILE <- file.path(
  "PIPELINE_V2",
  "config",
  "figure_palettes",
  "L2_palette.tsv"
)

if (!file.exists(L2_PALETTE_FILE)) {
  stop(
    "Missing frozen L2 palette: ",
    L2_PALETTE_FILE
  )
}

L2_palette_raw <- read.delim(
  L2_PALETTE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_L2_palette_cols <- c(
  "module",
  "L1",
  "final_L2_label",
  "hex",
  "L2_order_within_L1",
  "global_display_order",
  "palette_status"
)

if (!all(
  required_L2_palette_cols %in%
    colnames(L2_palette_raw)
)) {
  stop(
    "Frozen L2 palette missing required columns: ",
    paste(
      setdiff(
        required_L2_palette_cols,
        colnames(L2_palette_raw)
      ),
      collapse=", "
    )
  )
}

L2_palette_frozen <- L2_palette_raw |>
  filter(
    palette_status == "FROZEN"
  ) |>
  arrange(
    global_display_order
  )

if (anyDuplicated(
  L2_palette_frozen$final_L2_label
) > 0) {
  dup <- unique(
    L2_palette_frozen$final_L2_label[
      duplicated(
        L2_palette_frozen$final_L2_label
      )
    ]
  )

  stop(
    "Duplicated frozen L2 labels in palette: ",
    paste(dup, collapse="; ")
  )
}

L2_cols_frozen <- setNames(
  L2_palette_frozen$hex,
  L2_palette_frozen$final_L2_label
)

if (length(L2_cols_frozen) == 0) {
  stop(
    "Frozen L2 palette loaded but contains no colors."
  )
}

cat("\n===== FROZEN L2 PALETTE =====\n")
cat(
  "Frozen L2 colors loaded :",
  length(L2_cols_frozen),
  "\n"
)
cat("L2_PALETTE=PASS\n")


x <- readRDS(INPUT_RDS)

required_components <- c(
  "primary",
  "functional",
  "disease_effects",
  "model_support",
  "bulge_patient_data"
)

if (!all(required_components %in% names(x))) {
  stop(
    "Frozen object missing components: ",
    paste(
      setdiff(required_components, names(x)),
      collapse=", "
    )
  )
}

primary    <- x$primary
functional <- x$functional
disease    <- x$disease_effects

# ------------------------------------------------------------------
# Frozen L2 palette audit for every eligible S7 state
# ------------------------------------------------------------------

S7_L2 <- sort(
  unique(
    as.character(
      primary$L2
    )
  )
)

missing_S7_L2 <- setdiff(
  S7_L2,
  names(L2_cols_frozen)
)

cat("\n===== S7 L2 PALETTE AUDIT =====\n")
cat("Eligible frozen L2 states :", length(S7_L2), "\n")
cat("Matched frozen colors     :", length(S7_L2) - length(missing_S7_L2), "\n")

if (length(missing_S7_L2) > 0) {

  cat("\nUNMATCHED L2 LABELS:\n")
  print(missing_S7_L2)

  stop(
    "AUDIT FAIL - not all S7 L2 states match the frozen L2 palette."
  )
}

cat("L2_MATCH=PASS\n")

if (nrow(primary) != 23) {
  stop("Expected 23 primary niche states; found ", nrow(primary))
}

if (nrow(functional) != 92) {
  stop("Expected 92 functional relationships; found ", nrow(functional))
}

if (nrow(disease) != 23) {
  stop("Expected 23 disease-effect rows; found ", nrow(disease))
}

cat("\n============================================================\n")
cat("FIGS07 — MESENCHYMAL-VASCULAR / BULGE COUPLING\n")
cat("============================================================\n")
cat("Patients                 :", x$n_patients, "\n")
cat("Eligible niche states    :", nrow(primary), "\n")
cat("Functional relationships :", nrow(functional), "\n")
cat("Primary predictor        :", x$primary_predictor, "\n")
cat("Statistical unit         :", x$statistical_unit, "\n")

# ==============================================================================

# === FIGS07 FROZEN L1 PALETTE START ===

# Canonical frozen L1 palette used by Figures 1, 4 and 5.
L1_PALETTE_FILE <- file.path(
  "PIPELINE_V2",
  "config",
  "figure_palettes",
  "L1_palette.tsv"
)

if (!file.exists(L1_PALETTE_FILE)) {
  stop(
    "Missing frozen L1 palette: ",
    L1_PALETTE_FILE
  )
}

L1_palette_raw <- read.delim(
  L1_PALETTE_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_palette_cols <- c(
  "L1",
  "hex",
  "display_order"
)

if (!all(
  required_palette_cols %in%
    colnames(L1_palette_raw)
)) {
  stop(
    "Frozen L1 palette does not contain required columns: ",
    paste(required_palette_cols, collapse=", ")
  )
}

L1_palette_raw <- L1_palette_raw |>
  arrange(display_order)

L1_cols_frozen <- setNames(
  L1_palette_raw$hex,
  L1_palette_raw$L1
)

required_S7_L1 <- c(
  "Fibroblast",
  "Endothelial",
  "Mural"
)

missing_S7_L1 <- setdiff(
  required_S7_L1,
  names(L1_cols_frozen)
)

if (length(missing_S7_L1) > 0) {
  stop(
    "AUDIT FAIL - missing frozen L1 colors for: ",
    paste(missing_S7_L1, collapse=", ")
  )
}

if (any(is.na(
  L1_cols_frozen[required_S7_L1]
))) {
  stop(
    "AUDIT FAIL - NA values in required frozen L1 colors"
  )
}

cat("\n===== FROZEN L1 PALETTE =====\n")

print(
  L1_palette_raw |>
    filter(
      L1 %in% required_S7_L1
    ) |>
    select(
      L1,
      hex,
      display_order
    )
)

cat(
  "L1_PALETTE=PASS - canonical Figure 1 palette loaded\n"
)

# === FIGS07 FROZEN L1 PALETTE END ===





# 3. DISPLAY HELPERS
# ==============================================================================

short_label <- function(z) {
  z |>
    str_replace("^Fibroblast :: ", "") |>
    str_replace("^Vascular_EC :: ", "") |>
    str_replace("^Lymphatic_EC :: ", "") |>
    str_replace("^Pericyte_SMC :: ", "") |>
    str_replace("^Smooth_muscle :: ", "")
}


# Return exact frozen L2 color.
# No fallback color is allowed for publication figures.
get_L2_color <- function(z) {

  out <- unname(
    L2_cols_frozen[
      as.character(z)
    ]
  )

  if (any(is.na(out))) {

    missing <- unique(
      as.character(z)[
        is.na(out)
      ]
    )

    stop(
      "AUDIT FAIL - missing frozen L2 colors for: ",
      paste(
        missing,
        collapse="; "
      )
    )
  }

  out
}


make_L2_axis_labels <- function(z) {

  z <- as.character(z)

  cols <- unname(
    L2_cols_frozen[z]
  )

  if (any(is.na(cols))) {

    missing <- unique(
      z[is.na(cols)]
    )

    stop(
      "AUDIT FAIL - missing frozen L2 colors for axis labels: ",
      paste(
        missing,
        collapse="; "
      )
    )
  }

  paste0(
    "<span style='color:",
    cols,
    ";'><b>",
    z,
    "</b></span>"
  )
}

L1_labels <- c(
  Fibroblast    = "Fibroblast",
  Vascular_EC   = "Vascular EC",
  Lymphatic_EC  = "Lymphatic EC",
  Pericyte_SMC  = "Pericyte/SMC",
  Smooth_muscle = "Smooth muscle"
)

# Published L1 corresponding to each analytical niche compartment
niche_to_publication_L1 <- c(
  "Fibroblast"    = "Fibroblast",
  "Vascular EC"   = "Endothelial",
  "Lymphatic EC"  = "Endothelial",
  "Pericyte/SMC"  = "Mural",
  "Smooth muscle" = "Mural"
)

# Recover the exact frozen atlas colors used in Figures 1 and 5
display_L1_color <- c(
  "Fibroblast" =
    unname(L1_cols_frozen["Fibroblast"]),
  "Vascular EC" =
    unname(L1_cols_frozen["Endothelial"]),
  "Lymphatic EC" =
    unname(L1_cols_frozen["Endothelial"]),
  "Pericyte/SMC" =
    unname(L1_cols_frozen["Mural"]),
  "Smooth muscle" =
    unname(L1_cols_frozen["Mural"])
)

if (any(is.na(display_L1_color))) {
  stop("AUDIT FAIL - missing frozen colors for S7 L1 strips")
}

make_strip_background <- function(levels) {
  out <- lapply(
    levels,
    function(z) {
      element_rect(
        fill=unname(display_L1_color[z]),
        colour="#303030",
        linewidth=0.35
      )
    }
  )
  names(out) <- levels
  out
}

make_strip_text <- function(levels) {
  out <- lapply(
    levels,
    function(z) {
      element_text(
        colour="white",
        face="bold",
        size=7.0,
        angle=0
      )
    }
  )
  names(out) <- levels
  out
}


axis_labels <- c(
  stem_niche_mean = "Stem/niche",
  CD200_z_mean    = "CD200",
  cycling_mean    = "Cycling",
  MHCII_mean      = "MHC-II"
)

# ==============================================================================
# 4. PANEL A
# Mesenchymal-vascular differential abundance
# ==============================================================================

A_df <- disease |>
  mutate(
    label=short_label(feature),
    L1_display=recode(L1, !!!L1_labels),
    L2_color=get_L2_color(L2),
    label_html=paste0(
      "<span style='color:",
      L2_color,
      ";'><b>",
      label,
      "</b></span>"
    )
  ) |>
  select(
    feature,
    label,
    L1_display,
    KW_FDR,
    AA_FFA_FDR
  ) |>
  pivot_longer(
    cols=c(KW_FDR, AA_FFA_FDR),
    names_to="test",
    values_to="FDR"
  ) |>
  mutate(
    test=recode(
      test,
      KW_FDR="Global disease effect",
      AA_FFA_FDR="AA vs FFA"
    ),
    minuslog10FDR=-log10(pmax(FDR, 1e-300))
  )

A_order <- A_df |>
  group_by(label) |>
  summarise(
    best_FDR=min(FDR, na.rm=TRUE),
    .groups="drop"
  ) |>
  arrange(desc(best_FDR)) |>
  pull(label)

A_df <- A_df |>
  mutate(
    label=factor(label, levels=A_order)
  )

pA <- ggplot(
  A_df,
  aes(
    x=minuslog10FDR,
    y=label,
    shape=test
  )
) +
  geom_vline(
    xintercept=-log10(0.05),
    linetype="dashed",
    linewidth=0.35
  ) +
  geom_point(
    size=2.2,
    stroke=0.55
  ) +
  ggh4x::facet_grid2(
    L1_display ~ .,
    scales="free_y",
    space="free_y",
    strip=ggh4x::strip_themed(
      background_y=make_strip_background(
        levels(factor(A_df$L1_display))
      ),
      text_y=make_strip_text(
        levels(factor(A_df$L1_display))
      )
    )
  ) +
  scale_shape_manual(
    values=c(
      "Global disease effect"=16,
      "AA vs FFA"=1
    )
  ) +
  scale_y_discrete(
    labels=make_L2_axis_labels
  ) +
  labs(
    title="a  Mesenchymal-vascular composition",
    subtitle="No multiplicity-supported global disease or direct AA-FFA abundance differences",
    x=expression(-log[10]("FDR")),
    y=NULL,
    shape=NULL
  ) +
  theme_bw(base_size=8.5) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10
    ),
    plot.subtitle=element_text(
      size=7.3
    ),
    strip.text.y=element_text(
      angle=0,
      face="bold",
      size=7.3
    ),
    axis.text.y=ggtext::element_markdown(
      size=6.3
    ),
    legend.position="top",
    panel.grid.minor=element_blank()
  )

# ==============================================================================
# 5. PANEL B
# Disease-adjusted association with Bulge trajectory
# ==============================================================================

B_df <- primary |>
  mutate(
    label=short_label(feature),
    L1_display=recode(L1, !!!L1_labels),
    L2_color=get_L2_color(L2),
    label_html=paste0(
      "<span style='color:",
      L2_color,
      ";'><b>",
      label,
      "</b></span>"
    ),
    lo=M1_adjusted_beta - 1.96*M1_adjusted_SE,
    hi=M1_adjusted_beta + 1.96*M1_adjusted_SE,
    evidence=case_when(
      M1_adjusted_FDR < 0.05 ~ "FDR <0.05",
      M1_adjusted_FDR < 0.10 ~ "FDR 0.05-0.10",
      TRUE ~ "Not supported"
    )
  )

B_order <- B_df |>
  arrange(M1_adjusted_beta) |>
  pull(label)

B_df <- B_df |>
  mutate(
    label=factor(label, levels=unique(B_order))
  )

pB <- ggplot(
  B_df,
  aes(
    x=M1_adjusted_beta,
    y=label
  )
) +
  geom_vline(
    xintercept=0,
    linewidth=0.35
  ) +
  geom_errorbar(
    aes(
      xmin=lo,
      xmax=hi
    ),
    orientation="y",
    width=0,
    linewidth=0.4
  ) +
  geom_point(
    aes(shape=evidence),
    size=2.3,
    stroke=0.6
  ) +
  ggh4x::facet_grid2(
    L1_display ~ .,
    scales="free_y",
    space="free_y",
    strip=ggh4x::strip_themed(
      background_y=make_strip_background(
        levels(factor(B_df$L1_display))
      ),
      text_y=make_strip_text(
        levels(factor(B_df$L1_display))
      )
    )
  ) +
  scale_shape_manual(
    values=c(
      "FDR <0.05"=16,
      "FDR 0.05-0.10"=17,
      "Not supported"=1
    )
  ) +
  scale_y_discrete(
    labels=make_L2_axis_labels
  ) +
  labs(
    title="b  Coupling to Bulge trajectory position",
    subtitle="Disease-adjusted patient-level models",
    x=expression(beta~"(disease-adjusted)"),
    y=NULL,
    shape=NULL
  ) +
  theme_bw(base_size=8.5) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10
    ),
    plot.subtitle=element_text(
      size=7.3
    ),
    strip.text.y=element_text(
      angle=0,
      face="bold",
      size=7.3
    ),
    axis.text.y=ggtext::element_markdown(
      size=6.3
    ),
    legend.position="top",
    panel.grid.minor=element_blank()
  )

# ==============================================================================
# 6. PANEL C
# Three prioritized trajectory-associated niche states
# ==============================================================================

C_df <- primary |>
  filter(
    str_detect(
      feature,
      "COL11A1/DPEP1|PTGDS/IGF1|Epithelial-like EC"
    )
  ) |>
  mutate(
    plot_L1=case_when(
      L1 == "Fibroblast" ~ "Fibroblast",
      L1 == "Vascular_EC" ~ "Vascular EC",
      TRUE ~ NA_character_
    ),
    label=short_label(feature),
    L2_color=get_L2_color(L2),
    lo=M1_adjusted_beta - 1.96*M1_adjusted_SE,
    hi=M1_adjusted_beta + 1.96*M1_adjusted_SE,
    FDR_label=paste0(
      "FDR=",
      sprintf("%.3f", M1_adjusted_FDR)
    )
  )

if (nrow(C_df) != 3) {
  stop(
    "Expected exactly 3 prioritized trajectory states; found ",
    nrow(C_df)
  )
}

pC <- ggplot(
  C_df,
  aes(
    x=M1_adjusted_beta,
    y=reorder(label, M1_adjusted_beta)
  )
) +
  geom_vline(
    xintercept=0,
    linewidth=0.4
  ) +
  geom_errorbar(
    aes(
      xmin=lo,
      xmax=hi
    ),
    orientation="y",
    width=0.12,
    linewidth=0.55
  ) +
  geom_point(
    size=3
  ) +
  geom_text(
    aes(
      x=ifelse(
        M1_adjusted_beta >= 0,
        hi,
        lo
      ),
      label=FDR_label,
      hjust=ifelse(
        M1_adjusted_beta >= 0,
        -0.08,
        1.08
      )
    ),
    size=2.5
  ) +
  scale_x_continuous(
    expand=expansion(
      mult=c(0.30, 0.30)
    )
  ) +
  scale_colour_manual(
    values=setNames(
      C_df$L2_color,
      C_df$label
    ),
    guide="none"
  ) +
  labs(
    title="c  Prioritized niche states",
    subtitle="Association with patient-level Bulge median trajectory position",
    x=expression(beta~"(disease-adjusted)"),
    y=NULL
  ) +
  scale_y_discrete(
    labels=make_L2_axis_labels
  ) +
  theme_bw(base_size=9) +
  theme(
    axis.text.y=ggtext::element_markdown(
      size=7.2
    ),
    plot.title=element_text(
      face="bold",
      size=10
    ),
    plot.subtitle=element_text(
      size=7.5
    ),
    panel.grid.minor=element_blank()
  )

# ==============================================================================
# 7. PANEL D
# 92 niche-state × Bulge functional-feature relationships
# ==============================================================================

D_df <- functional |>
  mutate(
    label=short_label(feature),
    L2_color=get_L2_color(L2),
    label_html=paste0(
      "<span style='color:",
      L2_color,
      ";'><b>",
      label,
      "</b></span>"
    ),
    plot_L1=case_when(
      L1 == "Fibroblast" ~ "Fibroblast",
      L1 %in% c("Vascular_EC", "Lymphatic_EC") ~ "Vascular EC",
      L1 %in% c("Pericyte_SMC", "Smooth_muscle") ~ "Pericyte/SMC",
      TRUE ~ NA_character_
    ),
    axis=factor(
      recode(
        bulge_axis,
        !!!axis_labels
      ),
      levels=c(
        "Stem/niche",
        "CD200",
        "Cycling",
        "MHC-II"
      )
    ),
    evidence=case_when(
      adjusted_FDR < 0.05 ~ "FDR <0.05",
      adjusted_FDR < 0.10 ~ "FDR 0.05-0.10",
      TRUE ~ "NS"
    )
  )

D_order <- D_df |>
  group_by(label) |>
  summarise(
    best_FDR=min(adjusted_FDR, na.rm=TRUE),
    .groups="drop"
  ) |>
  arrange(desc(best_FDR)) |>
  pull(label)

D_df <- D_df |>
  mutate(
    label=factor(
      label,
      levels=D_order
    )
  )

pD <- ggplot(
  D_df,
  aes(
    x=axis,
    y=label,
    fill=adjusted_beta
  )
) +
  geom_tile(
    linewidth=0.25
  ) +
  geom_point(
    data=D_df |>
      filter(evidence != "NS"),
    aes(shape=evidence),
    size=2.2,
    stroke=0.65
  ) +
  scale_fill_gradient2(
    midpoint=0,
    name=expression(beta)
  ) +
  scale_shape_manual(
    values=c(
      "FDR <0.05"=8,
      "FDR 0.05-0.10"=4
    )
  ) +
  labs(
    title="d  Coupling to Bulge functional features",
    subtitle="Disease-adjusted models; symbols denote multiplicity-supported associations",
    x=NULL,
    y=NULL,
    shape=NULL
  ) +
  theme_bw(base_size=8.5) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10
    ),
    plot.subtitle=element_text(
      size=7.3
    ),
    axis.text.x=element_text(
      angle=0,
      hjust=0.5,
      vjust=0.5,
      size=7
    ),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid=element_blank(),
    legend.position="right"
  )


# ------------------------------------------------------------------------------
# Dedicated L2-colored row-label panel for D
# Vector-safe: no HTML/markdown rendering required.
# ------------------------------------------------------------------------------

D_label_df <- D_df |>
  distinct(
    label,
    L2,
    L2_color
  ) |>
  mutate(
    label=factor(
      as.character(label),
      levels=levels(D_df$label)
    )
  )

if (
  nrow(D_label_df) !=
  length(levels(D_df$label))
) {
  stop(
    "AUDIT FAIL - D label panel does not contain exactly one row per L2 state."
  )
}

if (any(is.na(D_label_df$L2_color))) {
  stop(
    "AUDIT FAIL - missing frozen L2 color in D label panel."
  )
}

pD_labels <- ggplot(
  D_label_df,
  aes(
    x=1,
    y=label,
    label=as.character(label),
    colour=L2_color
  )
) +
  geom_text(
    hjust=1,
    size=2.45,
    fontface="bold"
  ) +
  scale_colour_identity() +
  scale_x_continuous(
    limits=c(0,1),
    expand=c(0,0)
  ) +
  labs(
    x=NULL,
    y=NULL
  ) +
  theme_void() +
  theme(
    plot.margin=margin(
      t=0,
      r=2,
      b=0,
      l=0
    )
  )

# Combine vector-safe colored labels with beta heatmap
pD_final <- (
  pD_labels | pD
) +
  plot_layout(
    widths=c(0.57, 0.43)
  )


# ==============================================================================
# 8. ASSEMBLY
# ==============================================================================

final <- (
  pA | pB
) / (
  pC | pD_final
) +
  plot_layout(
    heights=c(1.12, 1),
    widths=c(1, 1)
  )

# ==============================================================================
# 9. EXPORT
# ==============================================================================

stem <- "FigS07_Mesenchymal_Bulge_coupling"

svg_file <- file.path(
  SVG_DIR,
  paste0(stem, ".svg")
)

pdf_file <- file.path(
  PDF_DIR,
  paste0(stem, ".pdf")
)

png_file <- file.path(
  PNG_DIR,
  paste0(stem, ".png")
)

ggsave(
  svg_file,
  final,
  width=13.2,
  height=9.4,
  units="in",
  device=svglite
)

ggsave(
  pdf_file,
  final,
  width=13.2,
  height=9.4,
  units="in"
)

ggsave(
  png_file,
  final,
  width=13.2,
  height=9.4,
  units="in",
  dpi=400
)

# ==============================================================================
# 10. EXPORT FROZEN SOURCE TABLES USED BY FIGURE
# ==============================================================================

write_csv(
  primary,
  file.path(
    TAB_DIR,
    "FigS07_primary_trajectory_coupling.csv"
  )
)

write_csv(
  functional,
  file.path(
    TAB_DIR,
    "FigS07_functional_coupling.csv"
  )
)

write_csv(
  disease,
  file.path(
    TAB_DIR,
    "FigS07_disease_effects.csv"
  )
)

# ==============================================================================
# 11. FORMAL CHECKS
# ==============================================================================

checks <- tibble(
  check=c(
    "16 patients",
    "23 eligible niche states",
    "92 functional relationships",
    "COL11A1/DPEP1 present",
    "PTGDS/IGF1 present",
    "Epithelial-like EC present",
    "APM-like smooth muscle / MHC-II present",
    "SVG exists",
    "PDF exists",
    "PNG exists"
  ),
  pass=c(
    x$n_patients == 16,
    nrow(primary) == 23,
    nrow(functional) == 92,
    any(
      str_detect(
        primary$feature,
        "COL11A1/DPEP1"
      )
    ),
    any(
      str_detect(
        primary$feature,
        "PTGDS/IGF1"
      )
    ),
    any(
      str_detect(
        primary$feature,
        "Epithelial-like EC"
      )
    ),
    any(
      str_detect(
        functional$feature,
        "APM-like smooth muscle"
      ) &
        functional$bulge_axis == "MHCII_mean"
    ),
    file.exists(svg_file),
    file.exists(pdf_file),
    file.exists(png_file)
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  checks,
  n=Inf
)

if (!all(checks$pass)) {
  stop("FIGS07 formal checks FAILED")
}

# ==============================================================================
# 12. FROZEN RESULT AUDIT
# ==============================================================================

cat("\n===== PRIMARY SUPPORTED FINDINGS =====\n")

primary_supported <- primary |>
  filter(
    !is.na(M1_adjusted_FDR),
    M1_adjusted_FDR < 0.10
  ) |>
  select(
    feature,
    L1,
    L2,
    M1_adjusted_beta,
    M1_adjusted_SE,
    M1_adjusted_p,
    M1_adjusted_FDR,
    adjusted_direction
  ) |>
  arrange(M1_adjusted_FDR)

print(
  primary_supported,
  n=Inf,
  width=Inf
)

cat("\n===== FUNCTIONAL SUPPORTED FINDINGS =====\n")

functional_supported <- functional |>
  filter(
    !is.na(adjusted_FDR),
    adjusted_FDR < 0.10
  ) |>
  select(
    feature,
    L1,
    L2,
    bulge_axis,
    rho,
    adjusted_beta,
    adjusted_SE,
    adjusted_p,
    adjusted_FDR
  ) |>
  arrange(adjusted_FDR)

print(
  functional_supported,
  n=Inf,
  width=Inf
)

# ==============================================================================
# 13. TARGETED MANUSCRIPT AUDIT
# ==============================================================================

cat("\n===== TARGETED PRIMARY FINDINGS =====\n")

target_primary <- primary |>
  filter(
    str_detect(
      feature,
      "COL11A1/DPEP1|PTGDS/IGF1|Epithelial-like EC"
    )
  ) |>
  select(
    feature,
    M1_adjusted_beta,
    M1_adjusted_SE,
    M1_adjusted_p,
    M1_adjusted_FDR
  ) |>
  arrange(M1_adjusted_FDR)

print(
  target_primary,
  n=Inf,
  width=Inf
)

cat("\n===== TARGETED FUNCTIONAL FINDING =====\n")

target_functional <- functional |>
  filter(
    str_detect(
      feature,
      "APM-like smooth muscle"
    ),
    bulge_axis == "MHCII_mean"
  ) |>
  select(
    feature,
    bulge_axis,
    rho,
    adjusted_beta,
    adjusted_SE,
    adjusted_p,
    adjusted_FDR
  )

print(
  target_functional,
  n=Inf,
  width=Inf
)

# ==============================================================================
# 14. EXACT RESULTS ↔ FROZEN OBJECT AUDIT
# ==============================================================================

get_primary_value <- function(pattern, variable) {

  z <- primary |>
    filter(
      str_detect(
        feature,
        fixed(pattern)
      )
    )

  if (nrow(z) != 1) {
    stop(
      "Expected exactly one primary row for: ",
      pattern,
      "; found ",
      nrow(z)
    )
  }

  z[[variable]][1]
}

get_functional_value <- function(
  pattern,
  axis,
  variable
) {

  z <- functional |>
    filter(
      str_detect(
        feature,
        fixed(pattern)
      ),
      bulge_axis == axis
    )

  if (nrow(z) != 1) {
    stop(
      "Expected exactly one functional row for: ",
      pattern,
      " / ",
      axis,
      "; found ",
      nrow(z)
    )
  }

  z[[variable]][1]
}

audit_values <- tibble(
  result=c(
    "COL11A1/DPEP1 -> Bulge trajectory",
    "Epithelial-like EC -> Bulge trajectory",
    "PTGDS/IGF1 -> Bulge trajectory",
    "APM-like smooth muscle -> Bulge MHC-II"
  ),

  expected_beta=c(
    2.737,
    5.964,
    -3.094,
    2.308
  ),

  expected_FDR=c(
    0.024,
    0.045,
    0.066,
    0.022
  ),

  observed_beta=c(
    get_primary_value(
      "COL11A1/DPEP1",
      "M1_adjusted_beta"
    ),
    get_primary_value(
      "Epithelial-like EC",
      "M1_adjusted_beta"
    ),
    get_primary_value(
      "PTGDS/IGF1",
      "M1_adjusted_beta"
    ),
    get_functional_value(
      "APM-like smooth muscle",
      "MHCII_mean",
      "adjusted_beta"
    )
  ),

  observed_FDR=c(
    get_primary_value(
      "COL11A1/DPEP1",
      "M1_adjusted_FDR"
    ),
    get_primary_value(
      "Epithelial-like EC",
      "M1_adjusted_FDR"
    ),
    get_primary_value(
      "PTGDS/IGF1",
      "M1_adjusted_FDR"
    ),
    get_functional_value(
      "APM-like smooth muscle",
      "MHCII_mean",
      "adjusted_FDR"
    )
  )
) |>
  mutate(
    beta_difference=observed_beta - expected_beta,
    FDR_difference=observed_FDR - expected_FDR,

    # Results values are rounded to three decimals.
    beta_match=
      round(observed_beta, 3) ==
      round(expected_beta, 3),

    FDR_match=
      round(observed_FDR, 3) ==
      round(expected_FDR, 3)
  )

cat("\n===== RESULTS ↔ FROZEN OBJECT AUDIT =====\n")

print(
  audit_values,
  n=Inf,
  width=Inf
)

write_csv(
  audit_values,
  file.path(
    TAB_DIR,
    "FigS07_manuscript_value_audit.csv"
  )
)

if (
  !all(audit_values$beta_match) ||
  !all(audit_values$FDR_match)
) {

  stop(
    paste0(
      "FIGS07 manuscript-value audit FAILED.\n",
      "At least one value in Results does not reproduce ",
      "the frozen 22G3A object after rounding to 3 decimals."
    )
  )
}

# ==============================================================================
# 15. COMPOSITIONAL AUDIT
# ==============================================================================

cat("\n===== COMPOSITIONAL DISEASE EFFECT AUDIT =====\n")

composition_supported <- disease |>
  filter(
    (!is.na(KW_FDR) & KW_FDR < 0.05) |
    (!is.na(AA_FFA_FDR) & AA_FFA_FDR < 0.05)
  ) |>
  select(
    feature,
    median_HC,
    median_FFA,
    median_AA,
    KW_p,
    KW_FDR,
    AA_FFA_p,
    AA_FFA_FDR
  ) |>
  arrange(
    pmin(
      KW_FDR,
      AA_FFA_FDR,
      na.rm=TRUE
    )
  )

print(
  composition_supported,
  n=Inf,
  width=Inf
)

if (nrow(composition_supported) != 0) {
  stop(
    paste0(
      "FIGS07 compositional audit FAILED: ",
      "at least one mesenchymal/vascular state has ",
      "FDR <0.05, inconsistent with the current Results text."
    )
  )
}

cat(
  "PASS: no mesenchymal/vascular state showed ",
  "FDR <0.05 for the global disease effect or ",
  "direct AA-vs-FFA abundance comparison.\n",
  sep=""
)

# ==============================================================================
# 16. SAVE FIGURE MANIFEST
# ==============================================================================

manifest <- tibble(
  figure="Supplementary Figure S7",
  module="22G3A",
  statistical_unit=x$statistical_unit,
  n_patients=x$n_patients,
  n_primary_states=nrow(primary),
  n_functional_relationships=nrow(functional),
  input_rds=INPUT_RDS,
  svg=svg_file,
  pdf=pdf_file,
  png=png_file
)

write_csv(
  manifest,
  file.path(
    TAB_DIR,
    "FigS07_manifest.csv"
  )
)

# ==============================================================================
# 17. FINAL STATUS
# ==============================================================================

cat("\n============================================================\n")
cat("FIGS07 STATUS = PASS\n")
cat("============================================================\n")
cat("Frozen module : 22G3A\n")
cat("Patients      :", x$n_patients, "\n")
cat("Niche states  :", nrow(primary), "\n")
cat("Relationships :", nrow(functional), "\n")
cat("SVG           :", svg_file, "\n")
cat("PDF           :", pdf_file, "\n")
cat("PNG           :", png_file, "\n")
cat("============================================================\n")

