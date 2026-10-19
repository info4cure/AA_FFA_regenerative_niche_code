
options(stringsAsFactors = FALSE, scipen = 999)

# ==============================================================================
# FIGURE 5 — SHARED MULTICOMPARTMENT PATHWAY ARCHITECTURE
#
# A — Shared Hallmark response
# B — Cell-state-resolved Hallmark landscape, grouped by biological programme
# C — High-confidence multicompartment pathway states
#
# Frozen outputs only.
# No statistical model is refitted.
# ==============================================================================

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

required <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "forcats",
  "patchwork",
  "scales",
  "svglite",
  "ggh4x",
  "ggtext"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(forcats)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ggh4x)
  library(ggtext)
})

HAS_GGREPEL <- requireNamespace(
  "ggrepel",
  quietly = TRUE
)

# ==============================================================================
# 1. ROOT + INPUTS
# ==============================================================================

ROOT <- normalizePath(".", mustWork = TRUE)

# ==============================================================================
# FROZEN L1 PALETTE — SAME SOURCE USED BY FIGURE 1
# ==============================================================================

L1_PALETTE_FILE <- file.path(
  ROOT,
  "PIPELINE_V2/config/figure_palettes/L1_palette.tsv"
)

if (!file.exists(L1_PALETTE_FILE)) {
  stop(
    "Missing frozen L1 palette: ",
    L1_PALETTE_FILE
  )
}

L1_palette_raw <- read.delim(
  L1_PALETTE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_palette_cols <- c(
  "L1",
  "hex",
  "display_order"
)

if (
  !all(
    required_palette_cols %in%
      colnames(L1_palette_raw)
  )
) {
  stop(
    "Frozen L1 palette does not contain required columns."
  )
}

L1_palette_raw <- L1_palette_raw |>
  arrange(
    display_order
  )

L1_cols_frozen <- setNames(
  L1_palette_raw$hex,
  L1_palette_raw$L1
)

L1_cols_B <- c(
  "Epithelial" =
    unname(L1_cols_frozen["Epithelial"]),

  "Fibroblast" =
    unname(L1_cols_frozen["Fibroblast"]),

  "Immune" =
    unname(L1_cols_frozen["Immune"]),

  "Endothelial" =
    unname(L1_cols_frozen["Endothelial"]),

  "Mural" =
    unname(L1_cols_frozen["Mural"]),

  "Sweat gland" =
    unname(L1_cols_frozen["Sweat_gland"]),

  "Melanocyte" =
    unname(L1_cols_frozen["Melanocyte"]),

  "Sebocyte" =
    unname(L1_cols_frozen["Sebocyte"])
)

if (any(is.na(L1_cols_B))) {
  stop(
    "AUDIT FAIL — one or more Figure 5 L1 colors could not be recovered."
  )
}

stopifnot(
  unname(
    L1_cols_B["Epithelial"]
  ) == "#8A034D"
)

message(
  "PASS — frozen Figure 1 L1 palette reused"
)

# ==============================================================================
# FROZEN L2 PALETTE — SAME SOURCE USED BY L2 UMAPs
# ==============================================================================

L2_PALETTE_FILE <- file.path(
  ROOT,
  "PIPELINE_V2/config/figure_palettes/L2_palette.tsv"
)

if (!file.exists(L2_PALETTE_FILE)) {
  stop(
    "Missing frozen L2 palette: ",
    L2_PALETTE_FILE
  )
}

L2_palette_raw <- read.delim(
  L2_PALETTE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
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

if (
  !all(
    required_L2_palette_cols %in%
      colnames(L2_palette_raw)
  )
) {
  stop(
    "Frozen L2 palette does not contain required columns."
  )
}

if (
  any(
    L2_palette_raw$palette_status !=
      "FROZEN"
  )
) {
  stop(
    "AUDIT FAIL — non-frozen entries detected in L2 palette."
  )
}

if (
  anyDuplicated(
    L2_palette_raw$final_L2_label
  )
) {
  stop(
    "AUDIT FAIL — duplicated final_L2_label in frozen L2 palette."
  )
}

L2_cols_frozen <- setNames(
  L2_palette_raw$hex,
  L2_palette_raw$final_L2_label
)

message(
  "PASS — frozen L2 UMAP palette loaded: ",
  nrow(L2_palette_raw),
  " identities"
)



PATH_GSEA_ALL <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS",
  "20A_HALLMARK_GSEA/tables",
  "20A_all_Hallmark_GSEA_results.csv"
)

PATH_GSEA_SIG <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS",
  "20A_HALLMARK_GSEA/tables",
  "20A_Hallmark_FDR005.csv"
)

PATH_INTEGRATED <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS",
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2B_integrated_prioritization/tables",
  "20F2B_integrated_state_prioritization.csv"
)

PATH_PRIMARY <- file.path(
  ROOT,
  "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS",
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2B_integrated_prioritization/tables",
  "20F2B_primary_publication_states.csv"
)

INPUTS <- c(
  PATH_GSEA_ALL,
  PATH_GSEA_SIG,
  PATH_INTEGRATED,
  PATH_PRIMARY
)

missing_inputs <- INPUTS[
  !file.exists(INPUTS)
]

if (length(missing_inputs) > 0) {
  stop(
    "\nMissing frozen input(s):\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

# ==============================================================================
# 2. OUTPUT — CENTRAL PUBLICATION FIGURE DIRECTORY
# ==============================================================================

PUBFIG_ROOT <- file.path(
  ROOT,
  "PIPELINE_V2/results/publication_figures"
)

OUT <- file.path(
  PUBFIG_ROOT,
  "FIG05_SHARED_PATHWAY_ARCHITECTURE"
)

SVG_DIR <- file.path(OUT, "SVG")
PDF_DIR <- file.path(OUT, "PDF")
PNG_DIR <- file.path(OUT, "PNG")
SOURCE_DIR <- file.path(OUT, "source_data")
AUDIT_DIR <- file.path(OUT, "audit")

for (d in c(
  PUBFIG_ROOT,
  OUT,
  SVG_DIR,
  PDF_DIR,
  PNG_DIR,
  SOURCE_DIR,
  AUDIT_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 3. CONSTANTS
# ==============================================================================

FDR_PRIMARY <- 0.05

COL_AA  <- "#D73027"
COL_FFA <- "#8FD175"

COL_UP   <- "#B2182B"
COL_DOWN <- "#3B6FB6"
COL_ZERO <- "#F7F7F7"

COL_DARK <- "#333333"
COL_GREY <- "#CFCFCF"

L1_ORDER <- c(
  "Epithelial",
  "Immune",
  "Fibroblast",
  "Endothelial",
  "Mural",
  "Melanocyte",
  "Sweat gland",
  "Sebocyte",
  "Other"
)

# ==============================================================================
# 4. HELPERS
# ==============================================================================

module_to_L1 <- function(x) {

  case_when(
    str_detect(x, "EPITHELIAL") ~ "Epithelial",
    str_detect(x, "IMMUNE") ~ "Immune",
    str_detect(x, "FIBROBLAST") ~ "Fibroblast",
    str_detect(x, "VASCULAR") ~ "Endothelial",
    str_detect(x, "MURAL") ~ "Mural",
    str_detect(x, "MELANOCYTE") ~ "Melanocyte",
    str_detect(x, "SWEAT") ~ "Sweat gland",
    str_detect(x, "SEBOCYTE") ~ "Sebocyte",
    TRUE ~ "Other"
  )
}

clean_pathway <- function(x) {

  x |>
    str_remove("^HALLMARK_") |>
    str_replace_all("_", " ") |>
    str_to_title()
}

clean_l2 <- function(x) {

  x |>
    str_replace_all("_", " ")
}

clean_axis <- function(x) {

  x |>
    str_replace_all("_", " ") |>
    str_to_title()
}

neglog10_cap <- function(
  x,
  cap = 12
) {

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

assert_n <- function(
  observed,
  expected,
  label
) {

  if (
    length(observed) != 1 ||
    is.na(observed) ||
    as.integer(observed) != as.integer(expected)
  ) {

    stop(
      sprintf(
        "AUDIT FAIL — %s: observed=%s expected=%s",
        label,
        paste(
          observed,
          collapse = ","
        ),
        expected
      )
    )
  }

  message(
    sprintf(
      "PASS — %s: %s",
      label,
      observed
    )
  )
}

theme_pub <- function(
  base_size = 9
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

      legend.background =
        element_rect(
          fill = "white",
          colour = NA
        ),

      legend.key =
        element_rect(
          fill = "white",
          colour = NA
        ),

      panel.grid.minor =
        element_blank(),

      panel.grid.major =
        element_line(
          linewidth = 0.18,
          colour = "#E9E9E9"
        ),

      panel.border =
        element_rect(
          colour = "#444444",
          linewidth = 0.4
        ),

      axis.text =
        element_text(
          colour = COL_DARK
        ),

      axis.title =
        element_text(
          face = "bold",
          colour = COL_DARK
        ),

      strip.background =
        element_rect(
          fill = "white",
          colour = "#444444",
          linewidth = 0.4
        ),

      strip.text =
        element_text(
          face = "bold",
          colour = COL_DARK
        ),

      legend.title =
        element_text(
          face = "bold"
        ),

      plot.title =
        element_text(
          face = "bold",
          size = rel(1.08)
        ),

      plot.margin =
        margin(
          6,
          7,
          6,
          6
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
      paste0(stem, ".svg")
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
      SVG_DIR,
      paste0(stem, "_AFFINITY.svg")
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::svg,
    bg = "white"
  )

  # Cairo PDF — vector master for Affinity/Illustrator interoperability
  ggsave(
    file.path(
      PDF_DIR,
      paste0(stem, "_AFFINITY.pdf")
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    file.path(
      PDF_DIR,
      paste0(stem, ".pdf")
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
      paste0(stem, ".png")
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
# 5. READ FROZEN INPUTS
# ==============================================================================

gsea_all <- read_csv(
  PATH_GSEA_ALL,
  show_col_types = FALSE
)

gsea_sig <- read_csv(
  PATH_GSEA_SIG,
  show_col_types = FALSE
)

integrated <- read_csv(
  PATH_INTEGRATED,
  show_col_types = FALSE
)

primary <- read_csv(
  PATH_PRIMARY,
  show_col_types = FALSE
)

# ==============================================================================
# 6. PREFLIGHT
# ==============================================================================

message("\n")
message("============================================================")
message("FIGURE 5 — PREFLIGHT")
message("============================================================")

assert_n(
  nrow(gsea_all),
  2916,
  "All Hallmark GSEA rows"
)

assert_n(
  nrow(gsea_sig),
  727,
  "Hallmark FDR<0.05 rows"
)

assert_n(
  sum(
    gsea_sig$contrast == "AA_vs_HC"
  ),
  313,
  "AA vs HC Hallmark hits"
)

assert_n(
  sum(
    gsea_sig$contrast == "FFA_vs_HC"
  ),
  268,
  "FFA vs HC Hallmark hits"
)

assert_n(
  sum(
    gsea_sig$contrast == "AA_vs_FFA"
  ),
  146,
  "AA vs FFA Hallmark hits"
)

assert_n(
  nrow(integrated),
  65,
  "Integrated pathway states"
)

assert_n(
  nrow(primary),
  22,
  "Frozen Tier 1/2 publication states"
)

assert_n(
  sum(
    integrated$integrated_evidence_tier == "TIER_1",
    na.rm = TRUE
  ),
  9,
  "Tier 1 states"
)

assert_n(
  sum(
    integrated$integrated_evidence_tier == "TIER_2",
    na.rm = TRUE
  ),
  13,
  "Tier 2 states"
)

# ==============================================================================
# 7. PANEL A — SHARED HALLMARK RESPONSE
# ==============================================================================

message("\n")
message("============================================================")
message("PANEL A — SHARED HALLMARK RESPONSE")
message("============================================================")

aa <- gsea_all |>
  filter(
    contrast == "AA_vs_HC"
  ) |>
  select(
    module,
    final_L2_label,
    pathway,
    NES_AA = NES,
    FDR_AA = padj
  )

ffa <- gsea_all |>
  filter(
    contrast == "FFA_vs_HC"
  ) |>
  select(
    module,
    final_L2_label,
    pathway,
    NES_FFA = NES,
    FDR_FFA = padj
  )

A <- inner_join(
  aa,
  ffa,
  by = c(
    "module",
    "final_L2_label",
    "pathway"
  )
) |>
  mutate(
    sig_AA =
      FDR_AA < FDR_PRIMARY,

    sig_FFA =
      FDR_FFA < FDR_PRIMARY,

    category =
      case_when(

        sig_AA &
          sig_FFA &
          sign(NES_AA) ==
            sign(NES_FFA) ~
          "Shared concordant",

        sig_AA &
          sig_FFA &
          sign(NES_AA) !=
            sign(NES_FFA) ~
          "Shared discordant",

        sig_AA ~
          "AA only",

        sig_FFA ~
          "FFA only",

        TRUE ~
          "Not significant"
      ),

    L1 =
      module_to_L1(
        module
      ),

    pathway_clean =
      clean_pathway(
        pathway
      ),

    L2_clean =
      clean_l2(
        final_L2_label
      ),

    label =
      paste0(
        L2_clean,
        " · ",
        pathway_clean
      )
  )

shared <- A |>
  filter(
    sig_AA,
    sig_FFA
  )

shared_n <- nrow(shared)

concordant_n <- sum(
  sign(
    shared$NES_AA
  ) ==
    sign(
      shared$NES_FFA
    )
)

concordance_pct <-
  100 *
  concordant_n /
  shared_n

assert_n(
  shared_n,
  137,
  "Shared Hallmark associations"
)

assert_n(
  concordant_n,
  128,
  "Shared directionally concordant associations"
)

discordant_n <- sum(
  A$category ==
    "Shared discordant"
)

assert_n(
  discordant_n,
  9,
  "Shared discordant associations"
)

message(
  sprintf(
    "PASS — shared directional concordance: %.3f%%",
    concordance_pct
  )
)

# Only the nine exceptions to the shared directional pattern are labelled.
labels_A <- A |>
  filter(
    category ==
      "Shared discordant"
  )

write_csv(
  A,
  file.path(
    SOURCE_DIR,
    "Fig05A_source_data.csv"
  )
)

write_csv(
  labels_A,
  file.path(
    SOURCE_DIR,
    "Fig05A_shared_discordant_labels.csv"
  )
)

CATEGORY_COLORS <- c(
  "Shared concordant" = "#222222",
  "Shared discordant" = "#7A3E9D",
  "AA only" = COL_AA,
  "FFA only" = COL_FFA
)

pA <- ggplot(
  A,
  aes(
    x = NES_AA,
    y = NES_FFA
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.30,
    colour = "#777777"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.30,
    colour = "#777777"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.45,
    colour = "#888888"
  ) +
  geom_point(
    data = A |>
      filter(
        category ==
          "Not significant"
      ),
    colour = COL_GREY,
    size = 1.0,
    alpha = 0.25
  ) +
  geom_point(
    data = A |>
      filter(
        category !=
          "Not significant"
      ),
    aes(
      colour = category
    ),
    size = 2.25,
    alpha = 0.88
  ) +
  scale_colour_manual(
    values = CATEGORY_COLORS,
    breaks = c(
      "Shared concordant",
      "Shared discordant",
      "AA only",
      "FFA only"
    ),
    name = NULL
  ) +
  annotate(
    "label",
    x = Inf,
    y = -Inf,
    hjust = 1.03,
    vjust = -0.15,
    size = 3.0,
    label = sprintf(
      "Shared FDR < 0.05: %d\nConcordant: %d/%d (%.1f%%)",
      shared_n,
      concordant_n,
      shared_n,
      concordance_pct
    ),
    label.size = 0.25,
    fill = "white"
  ) +
  labs(
    title =
      "A  Shared Hallmark response",

    x =
      "AA vs HC GSEA NES",

    y =
      "FFA vs HC GSEA NES"
  ) +
  coord_equal() +
  theme_pub(
    base_size = 9
  ) +
  theme(
    legend.position =
      "bottom"
  )

if (
  HAS_GGREPEL &&
    nrow(labels_A) > 0
) {

  pA <- pA +
    ggrepel::geom_text_repel(
      data =
        labels_A,

      aes(
        x = NES_AA,
        y = NES_FFA,
        label = label
      ),

      inherit.aes =
        FALSE,

      size =
        2.15,

      box.padding =
        0.30,

      point.padding =
        0.15,

      max.overlaps =
        Inf,

      min.segment.length =
        0,

      segment.size =
        0.20,

      seed =
        15
    )
}

save_all(
  pA,
  "Fig05A_shared_Hallmark_response",
  width = 7.4,
  height = 6.3
)

# ==============================================================================
# 8. PANEL B — BIOLOGICALLY GROUPED HALLMARK LANDSCAPE
# ==============================================================================

message("\n")
message("============================================================")
message("PANEL B — BIOLOGICALLY GROUPED HALLMARK LANDSCAPE")
message("============================================================")

# ------------------------------------------------------------------------------
# Prespecified display ontology for the 50 MSigDB Hallmark collections.
#
# This is DISPLAY ORGANIZATION ONLY.
# It does not alter GSEA,
# significance testing, effect estimates, or evidence tiers.
# ------------------------------------------------------------------------------

HALLMARK_GROUPS <- tibble::tribble(
  ~pathway, ~bio_group, ~group_order, ~pathway_order,

  # --------------------------------------------------------------------------
  # IMMUNE / INFLAMMATORY
  # --------------------------------------------------------------------------
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",       "Immune / inflammatory",          1L,  1L,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",       "Immune / inflammatory",          1L,  2L,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",         "Immune / inflammatory",          1L,  3L,
  "HALLMARK_INFLAMMATORY_RESPONSE",           "Immune / inflammatory",          1L,  4L,
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",         "Immune / inflammatory",          1L,  5L,
  "HALLMARK_IL2_STAT5_SIGNALING",             "Immune / inflammatory",          1L,  6L,
  "HALLMARK_COMPLEMENT",                      "Immune / inflammatory",          1L,  7L,
  "HALLMARK_ALLOGRAFT_REJECTION",             "Immune / inflammatory",          1L,  8L,

  # --------------------------------------------------------------------------
  # CELL CYCLE / PROLIFERATION
  # --------------------------------------------------------------------------
  "HALLMARK_E2F_TARGETS",                     "Cell cycle / proliferation",      2L,  1L,
  "HALLMARK_G2M_CHECKPOINT",                  "Cell cycle / proliferation",      2L,  2L,
  "HALLMARK_MITOTIC_SPINDLE",                 "Cell cycle / proliferation",      2L,  3L,
  "HALLMARK_MYC_TARGETS_V1",                  "Cell cycle / proliferation",      2L,  4L,
  "HALLMARK_MYC_TARGETS_V2",                  "Cell cycle / proliferation",      2L,  5L,

  # --------------------------------------------------------------------------
  # DNA DAMAGE / CELLULAR STRESS
  # --------------------------------------------------------------------------
  "HALLMARK_DNA_REPAIR",                      "DNA repair / cellular stress",    3L,  1L,
  "HALLMARK_P53_PATHWAY",                     "DNA repair / cellular stress",    3L,  2L,
  "HALLMARK_UV_RESPONSE_UP",                  "DNA repair / cellular stress",    3L,  3L,
  "HALLMARK_UV_RESPONSE_DN",                  "DNA repair / cellular stress",    3L,  4L,
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",       "DNA repair / cellular stress",    3L,  5L,
  "HALLMARK_APOPTOSIS",                       "DNA repair / cellular stress",    3L,  6L,
  "HALLMARK_HYPOXIA",                         "DNA repair / cellular stress",    3L,  7L,

  # --------------------------------------------------------------------------
  # METABOLISM / MITOCHONDRIAL
  # --------------------------------------------------------------------------
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",       "Metabolism / mitochondrial",      4L,  1L,
  "HALLMARK_GLYCOLYSIS",                      "Metabolism / mitochondrial",      4L,  2L,
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "Metabolism / mitochondrial",      4L,  3L,
  "HALLMARK_FATTY_ACID_METABOLISM",           "Metabolism / mitochondrial",      4L,  4L,
  "HALLMARK_ADIPOGENESIS",                    "Metabolism / mitochondrial",      4L,  5L,
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",         "Metabolism / mitochondrial",      4L,  6L,
  "HALLMARK_PEROXISOME",                      "Metabolism / mitochondrial",      4L,  7L,
  "HALLMARK_HEME_METABOLISM",                 "Metabolism / mitochondrial",      4L,  8L,
  "HALLMARK_BILE_ACID_METABOLISM",            "Metabolism / mitochondrial",      4L,  9L,
  "HALLMARK_XENOBIOTIC_METABOLISM",           "Metabolism / mitochondrial",      4L, 10L,

  # --------------------------------------------------------------------------
  # GROWTH / DEVELOPMENTAL SIGNALING
  # --------------------------------------------------------------------------
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",      "Growth / developmental",          5L,  1L,
  "HALLMARK_HEDGEHOG_SIGNALING",              "Growth / developmental",          5L,  2L,
  "HALLMARK_NOTCH_SIGNALING",                 "Growth / developmental",          5L,  3L,
  "HALLMARK_TGF_BETA_SIGNALING",              "Growth / developmental",          5L,  4L,
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",         "Growth / developmental",          5L,  5L,
  "HALLMARK_MTORC1_SIGNALING",                "Growth / developmental",          5L,  6L,
  "HALLMARK_KRAS_SIGNALING_UP",               "Growth / developmental",          5L,  7L,
  "HALLMARK_KRAS_SIGNALING_DN",               "Growth / developmental",          5L,  8L,

  # --------------------------------------------------------------------------
  # TISSUE REMODELING / EPITHELIAL ORGANIZATION
  # --------------------------------------------------------------------------
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","Tissue remodeling / epithelial", 6L,  1L,
  "HALLMARK_MYOGENESIS",                     "Tissue remodeling / epithelial",  6L,  2L,
  "HALLMARK_APICAL_JUNCTION",                 "Tissue remodeling / epithelial",  6L,  3L,
  "HALLMARK_APICAL_SURFACE",                  "Tissue remodeling / epithelial",  6L,  4L,
  "HALLMARK_ANGIOGENESIS",                    "Tissue remodeling / epithelial",  6L,  5L,
  "HALLMARK_COAGULATION",                     "Tissue remodeling / epithelial",  6L,  6L,

  # --------------------------------------------------------------------------
  # HORMONAL / SECRETORY
  # --------------------------------------------------------------------------
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",         "Hormonal / secretory",            7L,  1L,
  "HALLMARK_ESTROGEN_RESPONSE_LATE",          "Hormonal / secretory",            7L,  2L,
  "HALLMARK_ANDROGEN_RESPONSE",               "Hormonal / secretory",            7L,  3L,
  "HALLMARK_PROTEIN_SECRETION",               "Hormonal / secretory",            7L,  4L,

  # --------------------------------------------------------------------------
  # REMAINING SIGNALING / CELLULAR PROGRAMMES
  # --------------------------------------------------------------------------
  "HALLMARK_APICAL_SURFACE",                  "Tissue remodeling / epithelial",  6L,  4L,
  "HALLMARK_PANCREAS_BETA_CELLS",             "Hormonal / secretory",            7L,  5L,
  "HALLMARK_SPERMATOGENESIS",                 "Other cellular programmes",       8L,  1L
)

# Remove accidental duplicate definitions while retaining the first assignment.
HALLMARK_GROUPS <- HALLMARK_GROUPS |>
  distinct(
    pathway,
    .keep_all = TRUE
  )


# ==============================================================================
# COMPACT HALLMARK LABELS — DISPLAY ONLY
# ==============================================================================

HALLMARK_SHORT <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"       = "IFN-α",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"       = "IFN-γ",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"         = "TNF/NF-κB",
  "HALLMARK_INFLAMMATORY_RESPONSE"           = "Inflammatory response",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"          = "IL-6/JAK/STAT3",
  "HALLMARK_IL2_STAT5_SIGNALING"              = "IL-2/STAT5",
  "HALLMARK_COMPLEMENT"                       = "Complement",
  "HALLMARK_ALLOGRAFT_REJECTION"              = "Allograft rejection",

  "HALLMARK_E2F_TARGETS"                      = "E2F targets",
  "HALLMARK_G2M_CHECKPOINT"                   = "G2M checkpoint",
  "HALLMARK_MITOTIC_SPINDLE"                  = "Mitotic spindle",
  "HALLMARK_MYC_TARGETS_V1"                   = "MYC targets V1",
  "HALLMARK_MYC_TARGETS_V2"                   = "MYC targets V2",

  "HALLMARK_DNA_REPAIR"                       = "DNA repair",
  "HALLMARK_P53_PATHWAY"                      = "p53 pathway",
  "HALLMARK_UV_RESPONSE_UP"                   = "UV response ↑",
  "HALLMARK_UV_RESPONSE_DN"                   = "UV response ↓",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"        = "UPR",
  "HALLMARK_APOPTOSIS"                        = "Apoptosis",
  "HALLMARK_HYPOXIA"                          = "Hypoxia",

  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        = "OXPHOS",
  "HALLMARK_GLYCOLYSIS"                       = "Glycolysis",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  = "ROS",
  "HALLMARK_FATTY_ACID_METABOLISM"            = "Fatty acid metabolism",
  "HALLMARK_ADIPOGENESIS"                     = "Adipogenesis",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS"          = "Cholesterol",
  "HALLMARK_PEROXISOME"                       = "Peroxisome",
  "HALLMARK_HEME_METABOLISM"                  = "Heme metabolism",
  "HALLMARK_BILE_ACID_METABOLISM"             = "Bile acid metabolism",
  "HALLMARK_XENOBIOTIC_METABOLISM"            = "Xenobiotic metabolism",

  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"       = "WNT/β-catenin",
  "HALLMARK_HEDGEHOG_SIGNALING"               = "Hedgehog",
  "HALLMARK_NOTCH_SIGNALING"                  = "Notch",
  "HALLMARK_TGF_BETA_SIGNALING"               = "TGF-β",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"          = "PI3K/AKT/mTOR",
  "HALLMARK_MTORC1_SIGNALING"                 = "mTORC1",
  "HALLMARK_KRAS_SIGNALING_UP"                = "KRAS ↑",
  "HALLMARK_KRAS_SIGNALING_DN"                = "KRAS ↓",

  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_MYOGENESIS"                       = "Myogenesis",
  "HALLMARK_APICAL_JUNCTION"                  = "Apical junction",
  "HALLMARK_APICAL_SURFACE"                   = "Apical surface",
  "HALLMARK_ANGIOGENESIS"                     = "Angiogenesis",
  "HALLMARK_COAGULATION"                      = "Coagulation",

  "HALLMARK_ESTROGEN_RESPONSE_EARLY"          = "Estrogen early",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"           = "Estrogen late",
  "HALLMARK_ANDROGEN_RESPONSE"                = "Androgen response",
  "HALLMARK_PROTEIN_SECRETION"                = "Protein secretion",
  "HALLMARK_PANCREAS_BETA_CELLS"              = "β-cell programme",
  "HALLMARK_SPERMATOGENESIS"                  = "Spermatogenesis"
)

# ------------------------------------------------------------------------------
# Complete ontology audit against the actual frozen Hallmark universe.
# ------------------------------------------------------------------------------

all_hallmarks <- sort(
  unique(
    gsea_all$pathway
  )
)

missing_short_labels <- setdiff(
  all_hallmarks,
  names(HALLMARK_SHORT)
)

if (length(missing_short_labels) > 0) {
  stop(
    "AUDIT FAIL — missing short Hallmark labels: ",
    paste(
      missing_short_labels,
      collapse = "; "
    )
  )
}

message(
  "PASS — compact display labels mapped for all frozen Hallmarks"
)


unmapped_hallmarks <- setdiff(
  all_hallmarks,
  HALLMARK_GROUPS$pathway
)

extra_mapped_hallmarks <- setdiff(
  HALLMARK_GROUPS$pathway,
  all_hallmarks
)

if (
  length(
    unmapped_hallmarks
  ) > 0
) {

  cat(
    "\nUNMAPPED HALLMARKS:\n"
  )

  cat(
    paste(
      unmapped_hallmarks,
      collapse = "\n"
    ),
    "\n"
  )

  stop(
    "STATUS=FAIL — Hallmark display ontology is incomplete."
  )
}

message(
  "PASS — all ",
  length(
    all_hallmarks
  ),
  " frozen Hallmarks are mapped to biological display groups"
)

if (
  length(
    extra_mapped_hallmarks
  ) > 0
) {

  message(
    "INFO — mapped Hallmarks not present in this frozen universe: ",
    paste(
      extra_mapped_hallmarks,
      collapse = "; "
    )
  )
}

write_csv(
  HALLMARK_GROUPS,
  file.path(
    AUDIT_DIR,
    "Fig05B_Hallmark_biological_ontology.csv"
  )
)

# ------------------------------------------------------------------------------
# Significant disease-vs-HC Hallmark universe.
# No recurrence-based filtering.
# Every Hallmark with FDR<0.05 in AA-vs-HC or FFA-vs-HC is retained.
# ------------------------------------------------------------------------------

B <- gsea_sig |>
  filter(
    contrast %in%
      c(
        "AA_vs_HC",
        "FFA_vs_HC"
      )
  ) |>
  mutate(
    disease =
      case_when(
        contrast == "AA_vs_HC" ~ "AA",
        contrast == "FFA_vs_HC" ~ "FFA",
        TRUE ~ NA_character_
      ),

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
        L2_cols_frozen[
          final_L2_label
        ]
      ),

    L2_display =
      paste0(
        "<span style='color:",
        L2_hex,
        ";'><b>",
        L2_clean,
        "</b></span>"
      ),

    pathway_short =
      unname(
        HALLMARK_SHORT[
          pathway
        ]
      ),

    neglog10FDR =
      neglog10_cap(
        padj,
        cap = 12
      )
  ) |>
  left_join(
    HALLMARK_GROUPS,
    by = "pathway"
  )

if (
  any(is.na(B$bio_group)) ||
  any(is.na(B$pathway_short))
) {
  stop(
    "AUDIT FAIL — incomplete Panel B display annotation."
  )
}

if (
  any(
    is.na(
      B$L2_hex
    )
  )
) {

  missing_L2_palette <- B |>
    filter(
      is.na(
        L2_hex
      )
    ) |>
    distinct(
      final_L2_label
    ) |>
    pull(
      final_L2_label
    )

  stop(
    "AUDIT FAIL — Panel B L2 identities missing from frozen L2 palette: ",
    paste(
      missing_L2_palette,
      collapse = "; "
    )
  )
}

message(
  "PASS — all Panel B L2 identities mapped to frozen UMAP colors"
)


# ------------------------------------------------------------------------------
# Biological group order
# ------------------------------------------------------------------------------

GROUP_ORDER <- HALLMARK_GROUPS |>
  distinct(
    bio_group,
    group_order
  ) |>
  arrange(
    group_order
  ) |>
  pull(
    bio_group
  )

pathway_display_table <- HALLMARK_GROUPS |>
  filter(
    pathway %in%
      unique(
        B$pathway
      )
  ) |>
  mutate(
    pathway_short =
      unname(
        HALLMARK_SHORT[
          pathway
        ]
      )
  ) |>
  arrange(
    group_order,
    pathway_order
  ) |>
  mutate(
    x_index =
      row_number()
  )

pathway_order_B <-
  pathway_display_table$pathway_short

# ------------------------------------------------------------------------------
# Cell-state order
# ------------------------------------------------------------------------------

l2_order_B <- B |>
  count(
    L1,
    L2_clean,
    name = "n_hits"
  ) |>
  mutate(
    L1_rank =
      match(
        L1,
        L1_ORDER
      )
  ) |>
  arrange(
    L1_rank,
    desc(n_hits),
    L2_clean
  ) |>
  pull(
    L2_clean
  ) |>
  unique()

B <- B |>
  mutate(
    disease =
      factor(
        disease,
        levels = c(
          "AA",
          "FFA"
        )
      ),

    L1 =
      factor(
        L1,
        levels = L1_ORDER
      ),

    bio_group =
      factor(
        bio_group,
        levels = GROUP_ORDER
      ),

    pathway_short =
      factor(
        pathway_short,
        levels = pathway_order_B
      ),

    L2_clean =
      factor(
        L2_clean,
        levels = rev(
          l2_order_B
        )
      )
  )

write_csv(
  B,
  file.path(
    SOURCE_DIR,
    "Fig05B_source_data.csv"
  )
)

# ------------------------------------------------------------------------------
# Functional block geometry
# ------------------------------------------------------------------------------

group_boundaries <- pathway_display_table |>
  group_by(
    bio_group,
    group_order
  ) |>
  summarise(
    xmin =
      min(x_index) - 0.5,

    xmax =
      max(x_index) + 0.5,

    xmid =
      mean(
        range(x_index)
      ),

    .groups =
      "drop"
  ) |>
  arrange(
    group_order
  )

separator_x <- head(
  group_boundaries$xmax,
  -1
)

GROUP_SHORT <- c(
  "Immune / inflammatory" =
    "IMMUNE",

  "Cell cycle / proliferation" =
    "CELL CYCLE",

  "DNA repair / cellular stress" =
    "DNA / STRESS",

  "Metabolism / mitochondrial" =
    "METABOLISM",

  "Growth / developmental" =
    "GROWTH / DEV.",

  "Tissue remodeling / epithelial" =
    "REMODELING",

  "Hormonal / secretory" =
    "HORMONAL",

  "Other cellular programmes" =
    "OTHER"
)

group_boundaries <- group_boundaries |>
  mutate(
    group_label =
      unname(
        GROUP_SHORT[
          as.character(
            bio_group
          )
        ]
      )
  )

if (
  any(
    is.na(
      group_boundaries$group_label
    )
  )
) {
  stop(
    "AUDIT FAIL — missing functional-group display label."
  )
}

write_csv(
  group_boundaries,
  file.path(
    SOURCE_DIR,
    "Fig05B_functional_group_boundaries.csv"
  )
)

# ------------------------------------------------------------------------------
# L1 strip colors from the frozen Figure 1 palette
# ------------------------------------------------------------------------------

L1_present <- L1_ORDER[
  L1_ORDER %in%
    unique(
      as.character(
        B$L1
      )
    )
]

missing_L1_colors <- setdiff(
  L1_present,
  names(
    L1_cols_B
  )
)

if (
  length(
    missing_L1_colors
  ) > 0
) {
  stop(
    "AUDIT FAIL — missing frozen L1 colors for: ",
    paste(
      missing_L1_colors,
      collapse = "; "
    )
  )
}

# White text works for these saturated frozen atlas colors.
strip_background_y <- lapply(
  L1_present,
  function(x) {
    element_rect(
      fill =
        unname(
          L1_cols_B[x]
        ),
      colour =
        "#303030",
      linewidth =
        0.35
    )
  }
)

strip_text_y <- lapply(
  L1_present,
  function(x) {
    element_text(
      colour =
        "white",
      face =
        "bold",
      size =
        7.0,
      angle =
        0
    )
  }
)

names(strip_background_y) <- L1_present
names(strip_text_y) <- L1_present

# ------------------------------------------------------------------------------
# Main dot matrix
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Colored L2 axis labels — frozen L2 UMAP palette
# ------------------------------------------------------------------------------

L2_label_map <- B |>
  distinct(
    L2_clean,
    L2_display
  )

if (
  anyDuplicated(
    L2_label_map$L2_clean
  )
) {
  stop(
    "AUDIT FAIL — non-unique colored display labels for Panel B L2 identities."
  )
}

L2_display_labels <- setNames(
  L2_label_map$L2_display,
  L2_label_map$L2_clean
)

pB_body <- ggplot(
  B,
  aes(
    x = pathway_short,
    y = L2_clean
  )
)

if (
  length(
    separator_x
  ) > 0
) {

  pB_body <- pB_body +
    geom_vline(
      xintercept =
        separator_x,
      colour =
        "#8A8A8A",
      linewidth =
        0.42
    )
}

pB_body <- pB_body +
  geom_point(
    aes(
      fill = NES,
      size = neglog10FDR
    ),
    shape = 21,
    colour = "#3A3A3A",
    stroke = 0.22,
    alpha = 0.96
  ) +
  ggh4x::facet_grid2(
    rows = vars(L1),
    cols = vars(disease),
    scales = "free_y",
    space = "free_y",
    independent = "none",
    axes = "all",
    remove_labels = "x",
    strip = ggh4x::strip_themed(
      background_x = list(
        element_rect(
          fill = "#F2F2F2",
          colour = "#303030",
          linewidth = 0.35
        )
      ),

      text_x = list(
        element_text(
          colour = "#202020",
          face = "bold",
          size = 8.2
        )
      ),

      background_y =
        strip_background_y,

      text_y =
        strip_text_y
    )
  ) +
  scale_y_discrete(
    labels =
      L2_display_labels
  ) +
  scale_fill_gradient2(
    low = COL_DOWN,
    mid = COL_ZERO,
    high = COL_UP,
    midpoint = 0,
    name = "GSEA NES"
  ) +
  scale_size_continuous(
    name = "-log10 FDR",
    range = c(
      1.6,
      4.9
    )
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
        angle = 50,
        hjust = 1,
        vjust = 1,
        size = 6.7
      ),

    axis.text.y =
      ggtext::element_markdown(
        size = 7.0,
        colour = "#303030"
      ),

    panel.spacing.x =
      unit(
        0.12,
        "lines"
      ),

    panel.spacing.y =
      unit(
        0.10,
        "lines"
      ),

    legend.position =
      "right",

    plot.margin =
      margin(
        0,
        7,
        5,
        6
      )
  )

# ------------------------------------------------------------------------------
# Functional-group header
#
# Separate plot aligned to the Hallmark columns. This avoids overloading the
# disease/L1 facet strips and makes the biological hierarchy explicit.
# ------------------------------------------------------------------------------

# Alternating neutral backgrounds only. Functional groups are intentionally
# NOT assigned saturated colors because color is already used for NES and L1.
header_fill <- rep(
  c(
    "#F2F2F2",
    "#E5E5E5"
  ),
  length.out =
    nrow(
      group_boundaries
    )
)

B_header <- group_boundaries |>
  mutate(
    header_fill =
      header_fill,

    block_width =
      xmax - xmin,

    narrow_block =
      block_width <= 4
  )

# Header is duplicated over the AA and FFA column blocks so that each disease
# matrix has exactly the same biological-group annotation.

B_header_plot <- bind_rows(
  B_header |>
    mutate(
      disease = "AA"
    ),

  B_header |>
    mutate(
      disease = "FFA"
    )
) |>
  mutate(
    disease =
      factor(
        disease,
        levels = c(
          "AA",
          "FFA"
        )
      )
  )

pB_header <- ggplot(
  B_header_plot
) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0,
      ymax = 1,
      fill = header_fill
    ),
    colour = "#555555",
    linewidth = 0.30
  ) +
  geom_text(
    data =
      B_header_plot |>
        filter(
          !narrow_block
        ),

    aes(
      x = xmid,
      y = 0.5,
      label = group_label
    ),

    size = 2.45,
    fontface = "bold",
    colour = "#202020"
  ) +

  geom_text(
    data =
      B_header_plot |>
        filter(
          narrow_block
        ),

    aes(
      x = xmid,
      y = 0.5,
      label = group_label
    ),

    angle = 90,
    size = 2.15,
    fontface = "bold",
    colour = "#202020"
  ) +
  facet_grid(
    cols = vars(disease)
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(
      0.5,
      length(
        pathway_order_B
      ) + 0.5
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_void() +
  theme(
    strip.text =
      element_blank(),

    strip.background =
      element_blank(),

    panel.spacing.x =
      unit(
        0.12,
        "lines"
      ),

    plot.margin =
      margin(
        3,
        7,
        0,
        6
      )
  )

# ------------------------------------------------------------------------------
# Panel B assembly
# ------------------------------------------------------------------------------

pB <- pB_header / pB_body +
  patchwork::plot_layout(
    heights = c(
      0.13,
      1.00
    )
  ) +
  patchwork::plot_annotation(
    title =
      "B  Cell-state-resolved Hallmark landscape",
    subtitle =
      "Hallmarks grouped by prespecified biological programme"
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
        size = 8.0,
        colour = "#444444",
        hjust = 0
      )
  )

save_all(
  pB,
  "Fig05B_biologically_grouped_Hallmark_landscape",
  width = 15.5,
  height = 9.2
)


# ==============================================================================
# 9. PANEL C — HIGH-CONFIDENCE MULTICOMPARTMENT PATHWAY STATES
# ==============================================================================

message("\n")
message("============================================================")
message("PANEL C — HIGH-CONFIDENCE MULTICOMPARTMENT PATHWAY STATES")
message("============================================================")

C0 <- primary |>
  filter(
    integrated_evidence_tier %in%
      c(
        "TIER_1",
        "TIER_2"
      )
  )

assert_n(
  nrow(C0),
  22,
  "Panel C Tier 1/2 states"
)

C_AA <- C0 |>
  transmute(
    module,
    final_L2_label,
    axis,

    disease =
      "AA",

    median_NES =
      AA_median_database_NES,

    primary_status =
      AA_primary_status,

    integrated_evidence_tier,

    architecture,

    publication_role
  )

C_FFA <- C0 |>
  transmute(
    module,
    final_L2_label,
    axis,

    disease =
      "FFA",

    median_NES =
      FFA_median_database_NES,

    primary_status =
      FFA_primary_status,

    integrated_evidence_tier,

    architecture,

    publication_role
  )

C <- bind_rows(
  C_AA,
  C_FFA
) |>
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
        L2_cols_frozen[
          final_L2_label
        ]
      ),

    L2_display =
      paste0(
        "<span style='color:",
        L2_hex,
        ";'><b>",
        L2_clean,
        "</b></span>"
      ),

    axis_clean =
      clean_axis(
        axis
      ),

    disease =
      factor(
        disease,
        levels = c(
          "AA",
          "FFA"
        )
      ),

    tier =
      case_when(
        integrated_evidence_tier ==
          "TIER_1" ~ "Tier 1",

        integrated_evidence_tier ==
          "TIER_2" ~ "Tier 2",

        TRUE ~ NA_character_
      ),

    tier =
      factor(
        tier,
        levels = c(
          "Tier 2",
          "Tier 1"
        )
      ),

    median_NES =
      as.numeric(
        median_NES
      )
  )

if (
  any(
    is.na(
      C$L2_hex
    )
  )
) {

  missing_C_L2 <- C |>
    filter(
      is.na(
        L2_hex
      )
    ) |>
    distinct(
      final_L2_label
    ) |>
    pull(
      final_L2_label
    )

  stop(
    "AUDIT FAIL — Panel C L2 identities missing from frozen L2 palette: ",
    paste(
      missing_C_L2,
      collapse = "; "
    )
  )
}

message(
  "PASS — all Panel C L2 identities mapped to frozen UMAP colors"
)

C_plot <- C |>
  filter(
    !is.na(
      primary_status
    ),
    primary_status != "NONE",
    !is.na(
      median_NES
    )
  )

# ------------------------------------------------------------------------------
# Biological axis order:
# inflammation -> follicular/developmental biology -> epithelial organization ->
# proliferation/repair -> fibrosis -> metabolism -> vascular remodeling.
# ------------------------------------------------------------------------------

AXIS_RAW_ORDER <- c(
  "IFN_JAK_STAT_ANTIGEN",
  "T_CELL_CYTOTOXIC",
  "WNT_HEDGEHOG_FOLLICLE",
  "EPITHELIAL_DIFFERENTIATION_JUNCTIONAL",
  "CELL_CYCLE_DNA_REPAIR",
  "ECM_COLLAGEN_FIBROSIS",
  "METABOLISM_OXPHOS_ROS",
  "VASCULAR_ANGIOGENESIS"
)

axis_inventory <- C0 |>
  distinct(
    axis
  )

unmapped_axes <- setdiff(
  axis_inventory$axis,
  AXIS_RAW_ORDER
)

if (
  length(
    unmapped_axes
  ) > 0
) {
  stop(
    "AUDIT FAIL — unmapped Panel C biological axes: ",
    paste(
      unmapped_axes,
      collapse = "; "
    )
  )
}

axis_order_C <- AXIS_RAW_ORDER[
  AXIS_RAW_ORDER %in%
    axis_inventory$axis
]

axis_clean_order_C <- clean_axis(
  axis_order_C
)

# ------------------------------------------------------------------------------
# Cell-state order:
# biological L1 order; Tier 1 before Tier 2 within each compartment.
# ------------------------------------------------------------------------------

l2_order_C <- C0 |>
  mutate(
    L1 =
      module_to_L1(
        module
      ),

    L2_clean =
      clean_l2(
        final_L2_label
      ),

    L1_rank =
      match(
        L1,
        L1_ORDER
      ),

    tier_rank_display =
      case_when(
        integrated_evidence_tier ==
          "TIER_1" ~ 1L,

        integrated_evidence_tier ==
          "TIER_2" ~ 2L,

        TRUE ~ 3L
      )
  ) |>
  group_by(
    L1,
    L2_clean
  ) |>
  summarise(
    L1_rank =
      min(
        L1_rank,
        na.rm = TRUE
      ),

    tier_rank_display =
      min(
        tier_rank_display,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) |>
  arrange(
    L1_rank,
    tier_rank_display,
    L2_clean
  ) |>
  pull(
    L2_clean
  ) |>
  unique()

C_plot <- C_plot |>
  mutate(
    L1 =
      factor(
        L1,
        levels =
          L1_ORDER
      ),

    L2_clean =
      factor(
        L2_clean,
        levels =
          rev(
            l2_order_C
          )
      ),

    axis_clean =
      factor(
        axis_clean,
        levels =
          axis_clean_order_C
      )
  )

# ------------------------------------------------------------------------------
# Colored Panel C L2 labels — frozen L2 UMAP palette
# ------------------------------------------------------------------------------

C_L2_label_map <- C_plot |>
  distinct(
    L2_clean,
    L2_display
  )

if (
  anyDuplicated(
    C_L2_label_map$L2_clean
  )
) {
  stop(
    "AUDIT FAIL — non-unique colored Panel C L2 display labels."
  )
}

C_L2_display_labels <- setNames(
  C_L2_label_map$L2_display,
  as.character(
    C_L2_label_map$L2_clean
  )
)

write_csv(
  C,
  file.path(
    SOURCE_DIR,
    "Fig05C_source_data_all_Tier1_Tier2.csv"
  )
)

write_csv(
  C_plot,
  file.path(
    SOURCE_DIR,
    "Fig05C_source_data_displayed.csv"
  )
)

write_csv(
  C0,
  file.path(
    AUDIT_DIR,
    "Fig05C_frozen_22_state_inventory.csv"
  )
)

# ------------------------------------------------------------------------------
# Visual encoding:
#
# Fill = direction and magnitude of cross-database consensus NES.
# Size = integrated evidence tier.
#
# Hallmark and GSVA evidence are NOT represented with additional symbols here.
# They contribute to the frozen integrated classification and remain available
# in the source-data tables.
# ------------------------------------------------------------------------------

pC <- ggplot(
  C_plot,
  aes(
    x = axis_clean,
    y = L2_clean
  )
) +
  geom_point(
    aes(
      fill =
        median_NES,

      size =
        tier
    ),

    shape =
      21,

    colour =
      "#303030",

    stroke =
      0.35,

    alpha =
      0.97
  ) +
  facet_grid(
    rows =
      vars(L1),

    cols =
      vars(disease),

    scales =
      "free_y",

    space =
      "free_y",

    drop =
      TRUE
  ) +
  scale_y_discrete(
    labels =
      C_L2_display_labels
  ) +
  scale_fill_gradient2(
    low =
      COL_DOWN,

    mid =
      COL_ZERO,

    high =
      COL_UP,

    midpoint =
      0,

    name =
      "Cross-database NES"
  ) +
  scale_size_manual(
    values =
      c(
        "Tier 2" = 3.4,
        "Tier 1" = 5.6
      ),

    name =
      "Integrated evidence"
  ) +
  labs(
    title =
      "C  High-confidence multicompartment pathway states",

    x =
      NULL,

    y =
      NULL
  ) +
  theme_pub(
    base_size =
      8.2
  ) +
  theme(
    axis.text.x =
      element_text(
        angle =
          38,

        hjust =
          1,

        vjust =
          1,

        size =
          7.0
      ),

    axis.text.y =
      ggtext::element_markdown(
        size = 7.2,
        colour = "#303030"
      ),

    strip.text.x =
      element_text(
        size =
          8.3
      ),

    strip.text.y =
      element_text(
        angle =
          0,

        size =
          7.0
      ),

    panel.spacing.x =
      unit(
        0.12,
        "lines"
      ),

    panel.spacing.y =
      unit(
        0.10,
        "lines"
      ),

    legend.position =
      "right"
  )

save_all(
  pC,
  "Fig05C_high_confidence_multicompartment_states",
  width = 12.7,
  height = 7.2
)

# ==============================================================================
# 10. FINAL FIGURE ASSEMBLY
# ==============================================================================

message("\n")
message("============================================================")
message("FINAL FIGURE ASSEMBLY")
message("============================================================")

# Global concordance and high-confidence integrated architecture on top.
# Detailed Hallmark landscape occupies the complete lower row.

top <- pA + pC +
  plot_layout(
    widths =
      c(
        0.80,
        1.35
      )
  )

final_fig <- top / pB +
  plot_layout(
    heights =
      c(
        0.82,
        1.35
      )
  )

save_all(
  final_fig,
  "Figure05_FINAL",
  width = 16.5,
  height = 14.0
)

# ==============================================================================
# 11. FORMAL AUDIT
# ==============================================================================

formal_checks <- tibble(
  check = c(
    "all_GSEA_rows",
    "GSEA_FDR005_rows",
    "AA_vs_HC_hits",
    "FFA_vs_HC_hits",
    "AA_vs_FFA_hits",
    "shared_disease_vs_HC_hits",
    "shared_concordant",
    "shared_discordant",
    "shared_concordance_pct",
    "Hallmark_universe_size",
    "Hallmark_unmapped",
    "integrated_states",
    "Tier1_states",
    "Tier2_states",
    "PanelC_primary_states",
    "PanelC_unmapped_axes"
  ),

  observed = c(
    nrow(
      gsea_all
    ),

    nrow(
      gsea_sig
    ),

    sum(
      gsea_sig$contrast ==
        "AA_vs_HC"
    ),

    sum(
      gsea_sig$contrast ==
        "FFA_vs_HC"
    ),

    sum(
      gsea_sig$contrast ==
        "AA_vs_FFA"
    ),

    shared_n,

    concordant_n,

    discordant_n,

    round(
      concordance_pct,
      3
    ),

    length(
      all_hallmarks
    ),

    length(
      unmapped_hallmarks
    ),

    nrow(
      integrated
    ),

    sum(
      integrated$integrated_evidence_tier ==
        "TIER_1",
      na.rm = TRUE
    ),

    sum(
      integrated$integrated_evidence_tier ==
        "TIER_2",
      na.rm = TRUE
    ),

    nrow(
      C0
    ),

    length(
      unmapped_axes
    )
  ),

  expected = c(
    2916,
    727,
    313,
    268,
    146,
    137,
    128,
    9,
    round(
      100 * 128 / 137,
      3
    ),
    length(
      all_hallmarks
    ),
    0,
    65,
    9,
    13,
    22,
    0
  )
)

formal_checks <- formal_checks |>
  mutate(
    status =
      if_else(
        as.character(observed) ==
          as.character(expected),
        "PASS",
        "CHECK"
      )
  )

write_csv(
  formal_checks,
  file.path(
    AUDIT_DIR,
    "Fig05_formal_checks.csv"
  )
)

# Hard-stop on any formal check failure.
if (
  any(
    formal_checks$status !=
      "PASS"
  )
) {

  print(
    formal_checks |>
      filter(
        status != "PASS"
      )
  )

  stop(
    "AUDIT FAIL — one or more Figure 5 formal checks failed."
  )
}

# ==============================================================================
# 12. INPUT MANIFEST
# ==============================================================================

input_manifest <- tibble(
  analytical_block = c(
    "20A_HALLMARK_GSEA",
    "20A_HALLMARK_GSEA",
    "20F2B_INTEGRATED_PRIORITIZATION",
    "20F2B_INTEGRATED_PRIORITIZATION"
  ),

  file = c(
    PATH_GSEA_ALL,
    PATH_GSEA_SIG,
    PATH_INTEGRATED,
    PATH_PRIMARY
  ),

  role = c(
    "Complete frozen Hallmark GSEA universe",
    "Frozen Hallmark associations at FDR < 0.05",
    "Integrated 65-state cross-database universe",
    "Frozen 22 Tier 1/Tier 2 publication states"
  )
)

write_csv(
  input_manifest,
  file.path(
    AUDIT_DIR,
    "Fig05_input_manifest.csv"
  )
)

# ==============================================================================
# 13. SESSION INFO
# ==============================================================================

capture.output(
  sessionInfo(),
  file =
    file.path(
      AUDIT_DIR,
      "Fig05_sessionInfo.txt"
    )
)

# ==============================================================================
# 14. OUTPUT AUDIT
# ==============================================================================

expected_stems <- c(
  "Fig05A_shared_Hallmark_response",
  "Fig05B_biologically_grouped_Hallmark_landscape",
  "Fig05C_high_confidence_multicompartment_states",
  "Figure05_FINAL"
)

expected_files <- unlist(
  lapply(
    expected_stems,
    function(stem) {

      c(
        file.path(
          SVG_DIR,
          paste0(
            stem,
            ".svg"
          )
        ),

        file.path(
          PDF_DIR,
          paste0(
            stem,
            ".pdf"
          )
        ),

        file.path(
          PNG_DIR,
          paste0(
            stem,
            ".png"
          )
        )
      )
    }
  )
)

missing_outputs <- expected_files[
  !file.exists(
    expected_files
  )
]

if (
  length(
    missing_outputs
  ) > 0
) {

  stop(
    "\nOUTPUT AUDIT FAIL — missing files:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}

output_manifest <- tibble(
  file =
    expected_files,

  exists =
    file.exists(
      expected_files
    ),

  size_bytes =
    file.info(
      expected_files
    )$size
)

if (
  any(
    is.na(
      output_manifest$size_bytes
    )
  ) ||
    any(
      output_manifest$size_bytes <= 0
    )
) {

  stop(
    "OUTPUT AUDIT FAIL — one or more output files are empty."
  )
}

write_csv(
  output_manifest,
  file.path(
    AUDIT_DIR,
    "Fig05_output_manifest.csv"
  )
)

# ==============================================================================
# 15. FINAL REPORT
# ==============================================================================

report_lines <- c(
  "================================================================================",
  "FIGURE 5 — SHARED MULTICOMPARTMENT PATHWAY ARCHITECTURE",
  "================================================================================",
  "",
  paste0(
    "ROOT: ",
    ROOT
  ),
  paste0(
    "OUTPUT: ",
    OUT
  ),
  "",
  "PANEL A — SHARED HALLMARK RESPONSE",
  paste0(
    "Shared FDR<0.05 associations: ",
    shared_n
  ),
  paste0(
    "Directionally concordant: ",
    concordant_n,
    "/",
    shared_n,
    " (",
    sprintf(
      "%.1f",
      concordance_pct
    ),
    "%)"
  ),
  paste0(
    "Directionally discordant: ",
    discordant_n
  ),
  "Only the shared-discordant associations are labelled.",
  "",
  "PANEL B — BIOLOGICALLY GROUPED HALLMARK LANDSCAPE",
  paste0(
    "Frozen Hallmark universe: ",
    length(
      all_hallmarks
    )
  ),
  "All significant AA-vs-HC and FFA-vs-HC Hallmark associations retained.",
  "No recurrence-based pathway filtering.",
  "Columns organized by prespecified biological programme.",
  "Biological display groups:",
  paste0(
    "  - ",
    GROUP_ORDER
  ),
  "",
  "PANEL C — HIGH-CONFIDENCE MULTICOMPARTMENT PATHWAY STATES",
  paste0(
    "Frozen Tier 1/2 states: ",
    nrow(
      C0
    )
  ),
  paste0(
    "Tier 1: ",
    sum(
      C0$integrated_evidence_tier ==
        "TIER_1",
      na.rm = TRUE
    )
  ),
  paste0(
    "Tier 2: ",
    sum(
      C0$integrated_evidence_tier ==
        "TIER_2",
      na.rm = TRUE
    )
  ),
  "Biological axes ordered from inflammation through follicular biology and tissue remodeling.",
  "",
  "STATISTICAL POLICY",
  "No statistical model refitted.",
  "Frozen 20A and 20F2B outputs reused.",
  "FDR < 0.05 retained as primary Hallmark threshold.",
  "Panel B grouping is display organization only.",
  "Panel C uses frozen integrated evidence tiers.",
  "",
  "EXPORTS",
  "SVG: PASS",
  "PDF: PASS",
  "PNG 400 dpi: PASS",
  "",
  "STATUS=PASS"
)

writeLines(
  report_lines,
  con =
    file.path(
      AUDIT_DIR,
      "Fig05_FINAL_REPORT.txt"
    )
)

# ==============================================================================
# 16. TERMINAL REPORT
# ==============================================================================

message("\n")
message("============================================================")
message("FIGURE 5 COMPLETE")
message("============================================================")

message(
  "Output: ",
  OUT
)

message("")

message(
  "A — shared Hallmarks: ",
  shared_n,
  "; concordant: ",
  concordant_n,
  "/",
  shared_n,
  " (",
  sprintf(
    "%.1f",
    concordance_pct
  ),
  "%)"
)

message(
  "B — Hallmark universe mapped: ",
  length(
    all_hallmarks
  ),
  "; unmapped: ",
  length(
    unmapped_hallmarks
  )
)

message(
  "C — high-confidence states: ",
  nrow(
    C0
  ),
  " (Tier 1=",
  sum(
    C0$integrated_evidence_tier ==
      "TIER_1",
    na.rm = TRUE
  ),
  "; Tier 2=",
  sum(
    C0$integrated_evidence_tier ==
      "TIER_2",
    na.rm = TRUE
  ),
  ")"
)

message("")
message("SVG/PDF/PNG: PASS")
message("SOURCE DATA: PASS")
message("AUDIT: PASS")
message("")
message("STATUS=PASS")

