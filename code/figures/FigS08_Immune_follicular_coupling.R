#!/usr/bin/env Rscript

# ==============================================================================
# FigS08_Immune_follicular_coupling.R
#
# Supplementary Figure S8
# Immune functional programmes covary with follicular remodeling in AA
#
# READ ONLY
# Frozen module: 22G4C3B
# Statistical unit: patient
# AA-only
# No new statistical inference
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(svglite)
})

# ==============================================================================
# 1. PATHS
# ==============================================================================

BASE <- file.path(
  "PIPELINE_V2/results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE",
  "22G4C3_IMMUNE_FUNCTION_FOLLICULAR_TOPOLOGY",
  "22G4C3B_AA_COUPLING"
)

INPUT_RDS <- file.path(
  BASE,
  "objects",
  "22G4C3B_AA_immune_function_follicular_coupling_freeze.rds"
)

OUT_DIR <- file.path(
  "PIPELINE_V2/results/publication_figures",
  "SUPPLEMENTARY",
  "FIGS08_IMMUNE_FOLLICULAR_COUPLING"
)

SVG_DIR <- file.path(OUT_DIR, "SVG")
PDF_DIR <- file.path(OUT_DIR, "PDF")
PNG_DIR <- file.path(OUT_DIR, "PNG")
TAB_DIR <- file.path(OUT_DIR, "TABLES")

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

if (!file.exists(INPUT_RDS)) {
  stop(
    "Missing frozen input: ",
    INPUT_RDS
  )
}

# ==============================================================================
# 2. LOAD FROZEN OBJECT
# ==============================================================================

x <- readRDS(INPUT_RDS)

required <- c(
  "primary_coupling",
  "GE5_sensitivity",
  "LOPO",
  "LOPO_summary",
  "integrated",
  "prioritized",
  "NK_coherence",
  "classification_rules",
  "interpretation_constraints"
)

missing_components <- setdiff(
  required,
  names(x)
)

if (length(missing_components) > 0) {
  stop(
    "Missing frozen components: ",
    paste(
      missing_components,
      collapse=", "
    )
  )
}

all_conn <- x$integrated
prior    <- x$prioritized
lopo     <- x$LOPO
nk_coh   <- x$NK_coherence

stopifnot(
  x$disease_universe == "AA only",
  x$statistical_unit == "patient",
  x$n_prespecified_connections == 63,
  nrow(all_conn) == 63
)

cat("\n============================================================\n")
cat("FIGS08 — IMMUNE FUNCTION × FOLLICULAR TOPOLOGY\n")
cat("============================================================\n")
cat("Disease universe          :", x$disease_universe, "\n")
cat("Statistical unit          :", x$statistical_unit, "\n")
cat("Prespecified connections  :", nrow(all_conn), "\n")
cat("Prioritized connections   :", nrow(prior), "\n")
cat("============================================================\n")

# ==============================================================================
# 3. FROZEN PALETTES
# ==============================================================================

L1_FILE <- file.path(
  "PIPELINE_V2/config/figure_palettes",
  "L1_palette.tsv"
)

L2_FILE <- file.path(
  "PIPELINE_V2/config/figure_palettes",
  "L2_palette.tsv"
)

stopifnot(
  file.exists(L1_FILE),
  file.exists(L2_FILE)
)

L1_pal <- read.delim(
  L1_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

L2_pal <- read.delim(
  L2_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  filter(
    palette_status == "FROZEN"
  )

L1_cols <- setNames(
  L1_pal$hex,
  L1_pal$L1
)

L2_cols <- setNames(
  L2_pal$hex,
  L2_pal$final_L2_label
)

# ==============================================================================
# 4. DISPLAY DEFINITIONS
# ==============================================================================

programme_labels <- c(
  CD8_CYTOTOXICITY="CD8 cytotoxicity",
  CD8_IFNA="CD8 IFN-alpha",
  CD8_IFNG="CD8 IFN-gamma",
  MAC_IFNA="Macrophage/monocyte IFN-alpha",
  MAC_TNFA_NFKB="Macrophage/monocyte TNF/NF-kB",
  NK_CYTOTOXICITY="NK cytotoxicity",
  NK_IFNG="NK IFN-gamma"
)

programme_order <- c(
  "CD8 cytotoxicity",
  "CD8 IFN-alpha",
  "CD8 IFN-gamma",
  "Macrophage/monocyte IFN-alpha",
  "Macrophage/monocyte TNF/NF-kB",
  "NK cytotoxicity",
  "NK IFN-gamma"
)

target_labels <- c(
  COCH_AQP4_dermal_papilla="Dermal papilla",
  COL11A1_DPEP1_dermal_sheath="Dermal sheath",
  Hair_shaft="Hair shaft",
  IRS_companion="IRS/companion",
  ORS="ORS",
  ORS_basal="ORS-basal",
  IFE_differentiating="IFE diff.",
  IFE_suprabasal="IFE suprabasal",
  IFE_terminal="IFE terminal"
)

target_order <- c(
  "Dermal papilla",
  "Dermal sheath",
  "Hair shaft",
  "IRS/companion",
  "ORS",
  "ORS-basal",
  "IFE diff.",
  "IFE suprabasal",
  "IFE terminal"
)

# Full frozen L2 identities corresponding to display targets
target_to_L2 <- c(
  "Dermal papilla"=
    "COCH/AQP4 dermal papilla",
  "Dermal sheath"=
    "COL11A1/DPEP1 dermal-sheath",
  "Hair shaft"=
    "Hair_shaft",
  "IRS/companion"=
    "IRS_companion",
  "ORS"=
    "ORS",
  "ORS-basal"=
    "ORS_basal",
  "IFE diff."=
    "IFE_Infundibulum_differentiating",
  "IFE suprabasal"=
    "IFE_Infundibulum_suprabasal",
  "IFE terminal"=
    "IFE_Infundibulum_terminal"
)

missing_target_L2 <- setdiff(
  unname(target_to_L2),
  names(L2_cols)
)

if (length(missing_target_L2) > 0) {
  stop(
    "AUDIT FAIL - missing target L2 palette entries: ",
    paste(missing_target_L2, collapse="; ")
  )
}

target_cols <- setNames(
  unname(
    L2_cols[
      unname(target_to_L2)
    ]
  ),
  names(target_to_L2)
)

if (any(is.na(target_cols))) {
  stop(
    "AUDIT FAIL - one or more target colors are NA."
  )
}

cat("
===== TARGET L2 PALETTE AUDIT =====
")
print(
  tibble(
    target=names(target_to_L2),
    frozen_L2=unname(target_to_L2),
    hex=unname(target_cols)
  ),
  n=Inf
)
cat("TARGET_L2_PALETTE=PASS
")



# ==============================================================================
# ANATOMICAL DISPLAY DOMAINS
# ==============================================================================

target_domain <- c(
  "Dermal papilla"="Lower mesenchyme",
  "Dermal sheath"="Lower mesenchyme",
  "Hair shaft"="Lower follicle",
  "IRS/companion"="Lower follicle",
  "ORS"="Lower follicle",
  "ORS-basal"="Lower follicle",
  "IFE diff."="Upper epithelium",
  "IFE suprabasal"="Upper epithelium",
  "IFE terminal"="Upper epithelium"
)

domain_cols <- c(
  "Lower mesenchyme"=
    unname(L1_cols["Fibroblast"]),
  "Lower follicle"=
    unname(L1_cols["Epithelial"]),
  "Upper epithelium"=
    unname(L1_cols["Epithelial"])
)

# Vector-safe L2 strip shared by panels a and d.
target_strip_df <- tibble(
  target=factor(
    target_order,
    levels=target_order
  ),
  y=1,
  L2=unname(
    target_to_L2[target_order]
  ),
  col=unname(
    target_cols[target_order]
  ),
  domain=unname(
    target_domain[target_order]
  )
)

if (nrow(target_strip_df) != 9) {
  stop(
    "AUDIT FAIL - target strip must contain exactly 9 states."
  )
}


# ==============================================================================
# 5. PANEL A — COMPLETE PRESPECIFIED AA COUPLING LANDSCAPE
# ==============================================================================

A <- all_conn |>
  mutate(
    programme=recode(
      immune_axis,
      !!!programme_labels
    ),
    target=recode(
      target_variable,
      !!!target_labels
    ),
    programme=factor(
      programme,
      levels=rev(programme_order)
    ),
    target=factor(
      target,
      levels=target_order
    ),
    tier=case_when(
      analysis_tier == "PRIMARY_CD8" ~ "Primary CD8",
      str_detect(analysis_tier, "^MAC") ~ "Macrophage/monocyte",
      str_detect(analysis_tier, "^NK") ~ "NK exploratory",
      TRUE ~ analysis_tier
    )
  )

pA <- ggplot(
  A,
  aes(
    x=target,
    y=programme,
    fill=rho
  )
) +
  geom_tile(
    colour="white",
    linewidth=0.35
  ) +
  geom_vline(
    xintercept=c(2.5,6.5),
    linewidth=0.65,
    colour="#303030"
  ) +
  geom_point(
    data=A |>
      filter(FDR_supported),
    shape=8,
    size=2.6
  ) +
  scale_fill_gradient2(
    midpoint=0,
    limits=c(-1,1),
    name=expression(rho)
  ) +
  labs(
    title="a  Prespecified immune-functional coupling landscape in AA",
    subtitle="Patient-level Spearman correlations; stars denote within-family FDR support (NK analyses exploratory)",
    x=NULL,
    y=NULL
  ) +
  theme_bw(base_size=8.5) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10
    ),
    plot.subtitle=element_text(
      size=7.4
    ),
    axis.text.x=element_text(
      angle=32,
      hjust=1,
      vjust=1,
      size=7.0
    ),
    axis.text.y=element_text(
      size=7.0,
      face="bold"
    ),
    panel.grid=element_blank(),
    legend.position="right"
  )



# ==============================================================================
# 6. PANEL B — PRIORITIZED CD8 / MACROPHAGE ASSOCIATIONS
# ==============================================================================

B <- prior |>
  filter(
    !str_starts(
      immune_axis,
      "NK_"
    )
  ) |>
  mutate(
    programme=recode(
      immune_axis,
      !!!programme_labels
    ),
    target=recode(
      target_variable,
      !!!target_labels
    ),
    connection_label=paste0(
      programme,
      " -> ",
      target
    ),
    evidence=case_when(
      FDR_supported ~
        "Within-family FDR supported",

      coupling_class == "PRIMARY_ROBUST" ~
        "Robust, not FDR-supported",

      coupling_class == "PRIMARY_SUGGESTIVE" ~
        "Suggestive",

      TRUE ~
        "Other"
    )
  ) |>
  arrange(
    rho
  ) |>
  mutate(
    connection_label=factor(
      connection_label,
      levels=connection_label
    )
  )

if (
  any(
    str_starts(
      B$immune_axis,
      "NK_"
    )
  )
) {
  stop(
    "AUDIT FAIL - NK programme present in Panel B."
  )
}

if (nrow(B) != 11) {
  stop(
    "AUDIT FAIL - Panel B expected 11 non-NK prioritized ",
    "connections; found ",
    nrow(B)
  )
}

pB <- ggplot(
  B,
  aes(
    x=rho,
    y=connection_label
  )
) +
  geom_vline(
    xintercept=0,
    linewidth=0.4,
    colour="#505050"
  ) +
  geom_segment(
    aes(
      x=0,
      xend=rho,
      y=connection_label,
      yend=connection_label
    ),
    linewidth=0.55,
    colour="#505050"
  ) +
  geom_point(
    aes(
      shape=evidence
    ),
    size=2.9,
    stroke=0.7
  ) +
  scale_shape_manual(
    values=c(
      "Within-family FDR supported"=8,
      "Robust, not FDR-supported"=16,
      "Suggestive"=17,
      "Other"=1
    ),
    drop=TRUE
  ) +
  scale_x_continuous(
    limits=c(-1,1),
    breaks=c(
      -1,
      -0.5,
      0,
      0.5,
      1
    )
  ) +
  labs(
    title=
      "b  Prioritized CD8 and macrophage/monocyte coupling",
    subtitle=
      "AA patients, n=6; patient is the statistical replicate",
    x=expression(
      "Spearman " * rho
    ),
    y=NULL,
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
    axis.text.y=element_text(
      size=7.7
    ),
    axis.text.x=element_text(
      size=8
    ),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.text=element_text(
      size=7.5
    ),
    panel.grid.minor=element_blank(),
    plot.margin=margin(
      6,
      8,
      6,
      6
    )
  )

# ==============================================================================
# 7. PANEL C — LOPO ROBUSTNESS
# ==============================================================================

C <- B |>
  transmute(
    connection_id,
    connection_label,
    rho,
    median_LOPO_rho,
    min_LOPO_rho,
    max_LOPO_rho,
    same_sign_fraction
  )

if (nrow(C) != 11) {
  stop(
    "AUDIT FAIL - Panel C expected 11 connections; found ",
    nrow(C)
  )
}

pC <- ggplot(
  C,
  aes(
    x=median_LOPO_rho,
    y=connection_label
  )
) +
  geom_vline(
    xintercept=0,
    linewidth=0.4,
    colour="#505050"
  ) +
  geom_errorbar(
    aes(
      xmin=min_LOPO_rho,
      xmax=max_LOPO_rho
    ),
    orientation="y",
    width=0.15,
    linewidth=0.6
  ) +
  geom_point(
    aes(
      size=same_sign_fraction
    )
  ) +
  geom_point(
    aes(
      x=rho
    ),
    shape=4,
    size=2.7,
    stroke=0.9
  ) +
  scale_x_continuous(
    limits=c(-1,1),
    breaks=c(
      -1,
      -0.5,
      0,
      0.5,
      1
    )
  ) +
  scale_size_continuous(
    range=c(
      1.8,
      3.5
    ),
    limits=c(
      0,
      1
    ),
    breaks=c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    name="LOPO same-sign fraction"
  ) +
  labs(
    title=
      "c  Leave-one-patient-out robustness",
    subtitle=
      "Circle, median LOPO rho; horizontal line, LOPO range; x, full-model rho.",
    x=expression(
      "Spearman " * rho
    ),
    y=NULL
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
    axis.text.y=element_text(
      size=7.7
    ),
    axis.text.x=element_text(
      size=8
    ),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.title=element_text(
      size=7.5
    ),
    legend.text=element_text(
      size=7.3
    ),
    panel.grid.minor=element_blank(),
    plot.margin=margin(
      6,
      8,
      6,
      6
    )
  )

# ==============================================================================
# 8. THREE-PANEL MEMBERSHIP AUDIT
# ==============================================================================

cat(
  "\n===== THREE-PANEL MEMBERSHIP AUDIT =====\n"
)

cat(
  "Panel A complete connections :",
  nrow(A),
  "\n"
)

cat(
  "Panel B non-NK prioritized   :",
  nrow(B),
  "\n"
)

cat(
  "Panel C non-NK LOPO          :",
  nrow(C),
  "\n"
)

if (nrow(A) != 63) {
  stop(
    "AUDIT FAIL - Panel A must contain 63 prespecified connections."
  )
}

if (nrow(B) != 11) {
  stop(
    "AUDIT FAIL - Panel B must contain 11 non-NK prioritized connections."
  )
}

if (nrow(C) != 11) {
  stop(
    "AUDIT FAIL - Panel C must contain the same 11 connections as Panel B."
  )
}

if (
  any(
    str_starts(
      B$immune_axis,
      "NK_"
    )
  )
) {
  stop(
    "AUDIT FAIL - NK programme present in Panel B."
  )
}

if (
  !all(
    c(
      "NK cytotoxicity",
      "NK IFN-gamma"
    ) %in%
      as.character(
        A$programme
      )
  )
) {
  stop(
    "AUDIT FAIL - Panel A does not retain both exploratory NK programmes."
  )
}

cat(
  "THREE_PANEL_MEMBERSHIP=PASS\n"
)

# ==============================================================================
# 9. FINAL ASSEMBLY — 1 COLUMN × 3 ROWS
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
    1.05,
    1.00,
    1.00
  ),
  guides="keep"
)


# ==============================================================================
# 10. EXPORT
# ==============================================================================

stem <- "FigS08_Immune_follicular_coupling"

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
  width=11.0,
  height=15.0,
  units="in",
  device=svglite
)

ggsave(
  pdf_file,
  final,
  width=11.0,
  height=15.0,
  units="in"
)

ggsave(
  png_file,
  final,
  width=11.0,
  height=15.0,
  units="in",
  dpi=400
)

# ==============================================================================
# 11. SOURCE TABLES
# ==============================================================================

write_csv(
  all_conn,
  file.path(
    TAB_DIR,
    "FigS08_all_63_connections.csv"
  )
)

write_csv(
  prior,
  file.path(
    TAB_DIR,
    "FigS08_prioritized_connections.csv"
  )
)

write_csv(
  x$LOPO_summary,
  file.path(
    TAB_DIR,
    "FigS08_LOPO_summary.csv"
  )
)

write_csv(
  nk_coh,
  file.path(
    TAB_DIR,
    "FigS08_NK_coherence.csv"
  )
)

# ==============================================================================
# 12. FROZEN RESULT AUDIT
# ==============================================================================

cat("\n===== PRIORITIZED CONNECTIONS =====\n")

print(
  prior |>
    select(
      connection_id,
      immune_L2,
      immune_axis,
      target_L2,
      n,
      rho,
      p,
      FDR_family,
      coupling_class,
      FDR_supported,
      same_sign_fraction
    ) |>
    arrange(
      desc(abs(rho))
    ),
  n=Inf,
  width=Inf
)

# ------------------------------------------------------------------------------
# Results-text targets
# ------------------------------------------------------------------------------

targets <- c(
  "CD8_CYTOTOXICITY__IFE_terminal",
  "CD8_CYTOTOXICITY__ORS_basal",
  "CD8_IFNG__IFE_terminal",
  "MAC_IFNA__IFE_terminal",
  "MAC_IFNA__ORS_basal",
  "MAC_TNFA_NFKB__IFE_differentiating",
  "NK_CYTOTOXICITY__COL11A1_DPEP1_dermal_sheath",
  "NK_CYTOTOXICITY__ORS_basal"
)

target_audit <- all_conn |>
  filter(
    connection_id %in% targets
  ) |>
  select(
    connection_id,
    immune_L2,
    immune_axis,
    target_L1,
    target_L2,
    n,
    rho,
    p,
    FDR_family,
    coupling_class,
    FDR_supported,
    median_LOPO_rho,
    min_LOPO_rho,
    max_LOPO_rho,
    same_sign_fraction
  ) |>
  arrange(
    match(
      connection_id,
      targets
    )
  )

cat("\n===== RESULTS-TEXT TARGET AUDIT =====\n")

print(
  target_audit,
  n=Inf,
  width=Inf
)

missing_targets <- setdiff(
  targets,
  target_audit$connection_id
)

if (length(missing_targets) > 0) {
  stop(
    "Missing Results-text targets: ",
    paste(
      missing_targets,
      collapse="; "
    )
  )
}

# ==============================================================================
# 13. EXACT NUMERICAL AUDIT AGAINST CURRENT RESULTS
# ==============================================================================

get_rho <- function(id) {

  z <- all_conn |>
    filter(
      connection_id == id
    )

  if (nrow(z) != 1) {
    stop(
      "Expected exactly one row for ",
      id,
      "; found ",
      nrow(z)
    )
  }

  z$rho[[1]]
}

get_n <- function(id) {

  z <- all_conn |>
    filter(
      connection_id == id
    )

  if (nrow(z) != 1) {
    stop(
      "Expected exactly one row for ",
      id
    )
  }

  z$n[[1]]
}

audit_values <- tibble(
  result=c(
    "CD8 cytotoxicity -> terminal IFE",
    "CD8 cytotoxicity -> ORS-basal",
    "CD8 IFNG -> terminal IFE",
    "Macrophage IFNA -> terminal IFE",
    "Macrophage IFNA -> ORS-basal",
    "Macrophage TNF/NFkB -> differentiating IFE",
    "NK cytotoxicity -> COL11A1/DPEP1",
    "NK cytotoxicity -> ORS-basal"
  ),

  connection_id=targets,

  expected_rho=c(
     0.886,
    -0.829,
     0.886,
     0.886,
    -0.829,
     0.886,
     1.000,
    -1.000
  ),

  observed_rho=vapply(
    targets,
    get_rho,
    numeric(1)
  ),

  observed_n=vapply(
    targets,
    get_n,
    integer(1)
  )
) |>
  mutate(
    rho_match=
      round(observed_rho, 3) ==
      round(expected_rho, 3),

    expected_n=case_when(
      str_starts(
        connection_id,
        "NK_"
      ) ~ 4L,
      TRUE ~ 6L
    ),

    n_match=
      observed_n ==
      expected_n
  )

cat("\n===== MANUSCRIPT VALUE AUDIT =====\n")

print(
  audit_values,
  n=Inf,
  width=Inf
)

write_csv(
  audit_values,
  file.path(
    TAB_DIR,
    "FigS08_manuscript_value_audit.csv"
  )
)

if (!all(audit_values$rho_match)) {
  stop(
    "FIGS08 rho audit FAILED: current Results values ",
    "do not reproduce the frozen object."
  )
}

if (!all(audit_values$n_match)) {
  stop(
    "FIGS08 n audit FAILED."
  )
}

# ==============================================================================
# 14. FDR AUDIT
# ==============================================================================

cat("\n===== FDR-SUPPORTED CONNECTIONS =====\n")

fdr_supported <- all_conn |>
  filter(
    FDR_supported
  ) |>
  select(
    connection_id,
    immune_L2,
    immune_axis,
    target_L2,
    n,
    rho,
    p,
    testing_family,
    FDR_family,
    coupling_class
  ) |>
  arrange(
    FDR_family
  )

print(
  fdr_supported,
  n=Inf,
  width=Inf
)

# Current Results state that the CD8/macrophage associations
# do not survive their prespecified testing-family correction.
non_NK_FDR_supported <- fdr_supported |>
  filter(
    !str_starts(
      connection_id,
      "NK_"
    )
  )

if (nrow(non_NK_FDR_supported) != 0) {
  stop(
    "FIGS08 FDR audit FAILED: a non-NK connection is ",
    "FDR-supported, inconsistent with current Results text."
  )
}

# ==============================================================================
# 15. NK EXPLORATORY AUDIT
# ==============================================================================

cat("\n===== NK EXPLORATORY AUDIT =====\n")

nk_targets <- all_conn |>
  filter(
    connection_id %in% c(
      "NK_CYTOTOXICITY__COL11A1_DPEP1_dermal_sheath",
      "NK_CYTOTOXICITY__ORS_basal"
    )
  ) |>
  select(
    connection_id,
    n,
    rho,
    p,
    testing_family,
    FDR_family,
    coupling_class,
    FDR_supported,
    same_sign_fraction
  )

print(
  nk_targets,
  n=Inf,
  width=Inf
)

if (!all(nk_targets$n == 4)) {
  stop(
    "FIGS08 NK audit FAILED: NK target analyses are not n=4."
  )
}

if (!all(nk_targets$FDR_supported)) {
  stop(
    paste0(
      "FIGS08 NK audit FAILED: both manuscript NK targets ",
      "are expected to survive correction within the ",
      "prespecified exploratory NK testing family."
    )
  )
}


# ==============================================================================
# TARGET PALETTE / DISPLAY AUDIT
# ==============================================================================

display_checks <- tibble(
  check=c(
    "9 anatomical display targets",
    "9 frozen target L2 colors resolved",
    "Panel A uses all 9 targets",
    "Panel A contains both NK programmes",
    "No NA target colors"
  ),
  pass=c(
    length(target_order) == 9,
    length(target_cols) == 9 &&
      all(!is.na(target_cols)),
    n_distinct(A$target) == 9,
    all(
      c(
        "NK cytotoxicity",
        "NK IFN-gamma"
      ) %in%
        as.character(A$programme)
    ),
    all(!is.na(target_cols))
  )
)

cat("\n===== DISPLAY / PALETTE CHECKS =====\n")
print(
  display_checks,
  n=Inf
)

if (!all(display_checks$pass)) {
  stop(
    "FIGS08 display/palette audit FAILED."
  )
}


# ==============================================================================
# 16. FORMAL FIGURE CHECKS
# ==============================================================================

checks <- tibble(
  check=c(
    "AA-only frozen universe",
    "Patient statistical unit",
    "63 prespecified connections",
    "15 frozen prioritized connections",
    "Panel A contains 63 connections",
    "Panel A contains all 9 anatomical targets",
    "Panel A retains both exploratory NK programmes",
    "Panel B contains 11 non-NK prioritized connections",
    "Panel C contains the same 11 connections as Panel B",
    "CD8 target rho values reproduce Results",
    "Macrophage target rho values reproduce Results",
    "NK target rho values reproduce Results",
    "NK target analyses n=4",
    "No non-NK FDR-supported connection",
    "Both manuscript NK targets FDR-supported",
    "SVG exists",
    "PDF exists",
    "PNG exists"
  ),

  pass=c(
    identical(
      x$disease_universe,
      "AA only"
    ),

    identical(
      x$statistical_unit,
      "patient"
    ),

    nrow(all_conn) == 63,

    nrow(prior) == 15,

    nrow(A) == 63,

    n_distinct(A$target) == 9,

    all(
      c(
        "NK cytotoxicity",
        "NK IFN-gamma"
      ) %in%
        as.character(A$programme)
    ),

    nrow(B) == 11,

    nrow(C) == 11 &&
      setequal(
        as.character(C$connection_id),
        as.character(B$connection_id)
      ),

    all(
      audit_values$rho_match[
        str_starts(
          audit_values$connection_id,
          "CD8_"
        )
      ]
    ),

    all(
      audit_values$rho_match[
        str_starts(
          audit_values$connection_id,
          "MAC_"
        )
      ]
    ),

    all(
      audit_values$rho_match[
        str_starts(
          audit_values$connection_id,
          "NK_"
        )
      ]
    ),

    all(
      audit_values$observed_n[
        str_starts(
          audit_values$connection_id,
          "NK_"
        )
      ] == 4
    ),

    nrow(non_NK_FDR_supported) == 0,

    all(
      nk_targets$FDR_supported
    ),

    file.exists(svg_file),

    file.exists(pdf_file),

    file.exists(png_file)
  )
)

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  checks,
  n=Inf
)

if (
  nrow(checks) != 18
) {
  stop(
    "FIGS08 formal-check audit malformed: expected 18 checks."
  )
}

if (!all(checks$pass)) {
  stop(
    "FIGS08 formal checks FAILED."
  )
}

cat(
  "FORMAL_CHECKS=PASS\n"
)


# ==============================================================================
# 17. MANIFEST
# ==============================================================================

manifest <- tibble(
  figure="Supplementary Figure S8",
  panel_layout="3 panels; 1 column x 3 rows",
  displayed_panels="A: complete AA coupling; B: prioritized non-NK; C: LOPO robustness",
  module=x$module,
  disease_universe=x$disease_universe,
  statistical_unit=x$statistical_unit,
  n_prespecified_connections=nrow(all_conn),
  n_prioritized_connections=nrow(prior),
  input_rds=INPUT_RDS,
  svg=svg_file,
  pdf=pdf_file,
  png=png_file
)

write_csv(
  manifest,
  file.path(
    TAB_DIR,
    "FigS08_manifest.csv"
  )
)

# ==============================================================================
# 18. FINAL STATUS
# ==============================================================================

cat("\n============================================================\n")
cat("FIGS08 STATUS = PASS\n")
cat("============================================================\n")
cat("Frozen module :", x$module, "\n")
cat("Universe      :", x$disease_universe, "\n")
cat("Unit          :", x$statistical_unit, "\n")
cat("Connections   :", nrow(all_conn), "\n")
cat("Prioritized   :", nrow(prior), "\n")
cat("SVG           :", svg_file, "\n")
cat("PDF           :", pdf_file, "\n")
cat("PNG           :", png_file, "\n")
cat("============================================================\n")

