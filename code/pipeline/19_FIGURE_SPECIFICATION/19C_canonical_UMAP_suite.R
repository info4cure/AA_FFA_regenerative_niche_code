suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 19C CANONICAL UMAP FIGURE SUITE\n")
cat("============================================================\n\n")

# ============================================================
# 0. PATHS
# ============================================================

GLOBAL_RDS <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "global_atlas_v2_L1_frozen.rds"
)

LEDGER_RDS <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17A_frozen_analysis_atlas/tables/",
  "17A_master_cell_ledger_FINAL.rds"
)

PAL_DIR <- "PIPELINE_V2/config/figure_palettes"

OUT <- paste0(
  "PIPELINE_V2/results/figure_specification/",
  "19C_UMAP_suite"
)

FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD FROZEN SOURCES
# ============================================================

cat("STEP 1/9 — loading frozen sources\n")

obj <- readRDS(GLOBAL_RDS)
ledger <- readRDS(LEDGER_RDS)

disease_pal <- read_tsv(
  file.path(PAL_DIR, "disease_palette.tsv"),
  show_col_types = FALSE
)

L1_pal <- read_tsv(
  file.path(PAL_DIR, "L1_palette.tsv"),
  show_col_types = FALSE
)

L2_pal <- read_tsv(
  file.path(PAL_DIR, "L2_palette.tsv"),
  show_col_types = FALSE
)

stopifnot(inherits(obj, "Seurat"))
stopifnot(ncol(obj) == 82893)
stopifnot(nrow(ledger) == 82893)
stopifnot(nrow(disease_pal) == 3)
stopifnot(nrow(L1_pal) == 9)
stopifnot(nrow(L2_pal) == 59)

cat("Global cells:", ncol(obj), "\n")
cat("Ledger cells:", nrow(ledger), "\n")
cat("Frozen L2 colours:", nrow(L2_pal), "\n")

# ============================================================
# 2. CANONICAL HARMONY UMAP
# ============================================================

cat("\nSTEP 2/9 — recovering canonical Harmony UMAP\n")

stopifnot("umap_harmony" %in% Reductions(obj))

emb <- Embeddings(obj, "umap_harmony")

stopifnot(nrow(emb) == 82893)
stopifnot(ncol(emb) >= 2)
stopifnot(!anyDuplicated(rownames(emb)))

umap <- tibble(
  cell = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2]
)

cat(
  "Canonical embedding:",
  nrow(umap),
  "cells\n"
)

# ============================================================
# 3. NORMALISE LEDGER KEYS / JOIN
# ============================================================

cat("\nSTEP 3/9 — joining final analysis ledger\n")

# Find cell column safely.
cell_candidates <- c(
  "cell",
  "barcode",
  "cell_id",
  "Cell",
  "CELL"
)

cell_col <- cell_candidates[
  cell_candidates %in% colnames(ledger)
][1]

if (is.na(cell_col)) {
  stop(
    "No cell identifier column found in 17A final ledger.\n",
    "Available columns:\n",
    paste(colnames(ledger), collapse = ", ")
  )
}

ledger2 <- ledger %>%
  rename(cell = all_of(cell_col))

stopifnot(n_distinct(ledger2$cell) == 82893)

# Disease is canonical in the Seurat object.
md <- obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

stopifnot("disease" %in% colnames(md))

base <- umap %>%
  left_join(
    md %>% select(cell, disease),
    by = "cell"
  ) %>%
  left_join(
    ledger2,
    by = "cell",
    suffix = c("", "_ledger")
  )

stopifnot(nrow(base) == 82893)
stopifnot(!anyDuplicated(base$cell))
stopifnot(all(!is.na(base$disease)))

cat("Joined cells:", nrow(base), "\n")

# ============================================================
# 4. RESOLVE FINAL L1 / L2 COLUMNS
# ============================================================

cat("\nSTEP 4/9 — resolving final L1/L2 analysis universes\n")

required_ledger <- c(
  "analysis_eligible_global",
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label"
)

missing_required <- setdiff(
  required_ledger,
  colnames(base)
)

if (length(missing_required) > 0) {
  stop(
    "Missing required final-ledger columns: ",
    paste(missing_required, collapse = ", ")
  )
}

# Disease universe = all original singlets.
disease_df <- base

# Final L1 universe = globally eligible cells.
L1_df <- base %>%
  filter(analysis_eligible_global)

# Frozen L2 universe only.
L2_df <- base %>%
  filter(analysis_eligible_L2)

stopifnot(nrow(disease_df) == 82893)
stopifnot(nrow(L1_df) == 82680)
stopifnot(nrow(L2_df) == 82494)

stopifnot(
  all(L1_df$final_compartment %in% L1_pal$L1)
)

stopifnot(
  all(L2_df$final_L2_label %in% L2_pal$final_L2_label)
)

stopifnot(
  n_distinct(L2_df$final_L2_label) == 59
)

cat("Disease universe:", nrow(disease_df), "\n")
cat("L1 universe:", nrow(L1_df), "\n")
cat("L2 universe:", nrow(L2_df), "\n")

# ============================================================
# 5. FROZEN PALETTE VECTORS
# ============================================================

cat("\nSTEP 5/9 — binding frozen colour specification\n")

disease_cols <- setNames(
  disease_pal$hex,
  disease_pal$disease
)

L1_cols <- setNames(
  L1_pal$hex,
  L1_pal$L1
)

L2_cols <- setNames(
  L2_pal$hex,
  L2_pal$final_L2_label
)

disease_levels <- disease_pal %>%
  arrange(display_order) %>%
  pull(disease)

L1_levels <- L1_pal %>%
  arrange(display_order) %>%
  pull(L1)

L2_levels <- L2_pal %>%
  arrange(global_display_order) %>%
  pull(final_L2_label)

disease_df$disease <- factor(
  disease_df$disease,
  levels = disease_levels
)

L1_df$final_compartment <- factor(
  L1_df$final_compartment,
  levels = L1_levels
)

L2_df$final_L2_label <- factor(
  L2_df$final_L2_label,
  levels = L2_levels
)

# ============================================================
# 6. PLOT HELPERS
# ============================================================

cat("\nSTEP 6/9 — defining canonical UMAP graphics\n")

theme_umap <- function() {
  theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0,
        margin = margin(b = 8)
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(4.0, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

make_umap <- function(
  dat,
  colour_var,
  colours,
  title,
  subtitle = NULL,
  point_size = 0.32,
  point_alpha = 0.88,
  stroke_size = 0.08,
  legend_ncol = 1
) {

  dat$.plot_colour <- dat[[colour_var]]

  # Draw larger groups first so smaller populations are not buried.
  ord <- table(dat$.plot_colour)

  dat$.draw_n <- unname(
    ord[as.character(dat$.plot_colour)]
  )

  dat <- dat %>%
    arrange(desc(.draw_n))

  ggplot(
    dat,
    aes(
      x = UMAP_1,
      y = UMAP_2
    )
  ) +
    geom_point(
      aes(fill = .plot_colour),
      shape = 21,
      colour = "#303030",
      stroke = stroke_size,
      size = point_size,
      alpha = point_alpha
    ) +
    scale_fill_manual(
      values = colours,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 3.2,
          alpha = 1,
          colour = "grey25",
          stroke = 0.25
        ),
        ncol = legend_ncol
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    coord_equal() +
    theme_umap()
}

save_plot <- function(
  p,
  stem,
  width,
  height
) {

  pdf_file <- file.path(
    FIG,
    paste0(stem, ".pdf")
  )

  png_file <- file.path(
    FIG,
    paste0(stem, ".png")
  )

  ggsave(
    pdf_file,
    p,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )

  ggsave(
    png_file,
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

# ============================================================
# 7. GLOBAL UMAP SUITE
# ============================================================

cat("\nSTEP 7/9 — generating global UMAP suite\n")

# ------------------------------------------------------------
# 7A Disease
# ------------------------------------------------------------

p_disease <- make_umap(
  disease_df,
  "disease",
  disease_cols,
  "Global atlas by disease",
  "Harmony UMAP | 82,893 cells",
  point_size = 0.65,
  point_alpha = 0.88,
  stroke_size = 0.045
)

save_plot(
  p_disease,
  "19C_01_global_UMAP_disease",
  width = 7.2,
  height = 5.8
)

cat("Generated: disease UMAP\n")

# ------------------------------------------------------------
# 7B Final L1
# ------------------------------------------------------------

p_L1 <- make_umap(
  L1_df,
  "final_compartment",
  L1_cols,
  "Global atlas — final L1 compartments",
  "Harmony UMAP | 82,680 analysis-eligible cells",
  point_size = 0.65,
  point_alpha = 0.90,
  stroke_size = 0.045
)

save_plot(
  p_L1,
  "19C_02_global_UMAP_L1",
  width = 8.2,
  height = 6.2
)

cat("Generated: L1 UMAP\n")

# ------------------------------------------------------------
# 7C Frozen 59-L2 global atlas
# ------------------------------------------------------------

p_L2 <- make_umap(
  L2_df,
  "final_L2_label",
  L2_cols,
  "Global atlas — frozen L2 identities",
  "Harmony UMAP | 82,494 cells | 59 frozen L2 identities",
  point_size = 0.52,
  point_alpha = 0.90,
  stroke_size = 0.035,
  legend_ncol = 2
)

save_plot(
  p_L2,
  "19C_03_global_UMAP_L2_59",
  width = 12.2,
  height = 8.4
)

cat("Generated: 59-L2 global UMAP\n")

# ============================================================
# 8. COMPARTMENT-SPECIFIC L2 UMAPs
# ============================================================

cat("\nSTEP 8/9 — generating compartment-specific L2 UMAPs\n")

parent_order <- L1_pal %>%
  filter(L1 != "Neural_Schwann") %>%
  arrange(display_order) %>%
  pull(L1)

figure_inventory <- tibble(
  figure = c(
    "19C_01_global_UMAP_disease",
    "19C_02_global_UMAP_L1",
    "19C_03_global_UMAP_L2_59"
  ),
  universe = c(
    "all_original_singlets",
    "final_L1_analysis_eligible",
    "frozen_L2_analysis_eligible"
  ),
  n_cells = c(
    nrow(disease_df),
    nrow(L1_df),
    nrow(L2_df)
  ),
  n_categories = c(
    n_distinct(disease_df$disease),
    n_distinct(L1_df$final_compartment),
    n_distinct(L2_df$final_L2_label)
  )
)

for (i in seq_along(parent_order)) {

  parent <- parent_order[i]

  allowed_L2 <- L2_pal %>%
    filter(L1 == parent) %>%
    arrange(L2_order_within_L1)

  dd <- L2_df %>%
    filter(final_compartment == parent)

  stopifnot(nrow(dd) > 0)

  dd$final_L2_label <- factor(
    as.character(dd$final_L2_label),
    levels = allowed_L2$final_L2_label
  )

  cols <- setNames(
    allowed_L2$hex,
    allowed_L2$final_L2_label
  )

  title_parent <- gsub("_", " ", parent)

  p <- make_umap(
    dd,
    "final_L2_label",
    cols,
    paste0(title_parent, " — frozen L2 identities"),
    paste0(
      "Harmony UMAP | ",
      format(nrow(dd), big.mark = ","),
      " cells | ",
      nrow(allowed_L2),
      " L2 identities"
    ),
    point_size = 0.90,
    point_alpha = 0.92,
    stroke_size = 0.12
  )

  stem <- sprintf(
    "19C_%02d_L2_UMAP_%s",
    i + 3,
    parent
  )

  save_plot(
    p,
    stem,
    width = 7.8,
    height = 6.2
  )

  figure_inventory <- bind_rows(
    figure_inventory,
    tibble(
      figure = stem,
      universe = paste0(parent, "_frozen_L2"),
      n_cells = nrow(dd),
      n_categories = nrow(allowed_L2)
    )
  )

  cat(
    sprintf(
      "%-15s : %6d cells | %2d L2\n",
      parent,
      nrow(dd),
      nrow(allowed_L2)
    )
  )
}

# ============================================================
# 9. AUDIT / EXPORT
# ============================================================

cat("\nSTEP 9/9 — formal figure audit\n")

write_csv(
  figure_inventory,
  file.path(
    TAB,
    "19C_figure_inventory.csv"
  )
)

disease_counts <- disease_df %>%
  count(disease, name = "n_cells")

L1_counts <- L1_df %>%
  count(
    final_compartment,
    name = "n_cells",
    .drop = FALSE
  )

L2_counts <- L2_df %>%
  count(
    final_compartment,
    final_L2_label,
    name = "n_cells",
    .drop = FALSE
  ) %>%
  filter(n_cells > 0)

write_csv(
  disease_counts,
  file.path(TAB, "19C_disease_counts.csv")
)

write_csv(
  L1_counts,
  file.path(TAB, "19C_L1_counts.csv")
)

write_csv(
  L2_counts,
  file.path(TAB, "19C_L2_counts.csv")
)

# Coordinate conservation.
coord_check_L1 <- L1_df %>%
  select(cell, UMAP_1, UMAP_2) %>%
  inner_join(
    umap,
    by = "cell",
    suffix = c("_plot", "_source")
  )

coord_check_L2 <- L2_df %>%
  select(cell, UMAP_1,
UMAP_2) %>%
  inner_join(
    umap,
    by = "cell",
    suffix = c("_plot", "_source")
  )

coordinates_identical_L1 <-
  nrow(coord_check_L1) == nrow(L1_df) &&
  isTRUE(all.equal(
    coord_check_L1$UMAP_1_plot,
    coord_check_L1$UMAP_1_source,
    tolerance = 0
  )) &&
  isTRUE(all.equal(
    coord_check_L1$UMAP_2_plot,
    coord_check_L1$UMAP_2_source,
    tolerance = 0
  ))

coordinates_identical_L2 <-
  nrow(coord_check_L2) == nrow(L2_df) &&
  isTRUE(all.equal(
    coord_check_L2$UMAP_1_plot,
    coord_check_L2$UMAP_1_source,
    tolerance = 0
  )) &&
  isTRUE(all.equal(
    coord_check_L2$UMAP_2_plot,
    coord_check_L2$UMAP_2_source,
    tolerance = 0
  ))

# ------------------------------------------------------------
# Palette coverage
# ------------------------------------------------------------

disease_observed <- unique(
  as.character(disease_df$disease)
)

L1_observed <- unique(
  as.character(L1_df$final_compartment)
)

L2_observed <- unique(
  as.character(L2_df$final_L2_label)
)

disease_palette_complete <-
  all(disease_observed %in% names(disease_cols))

L1_palette_complete <-
  all(L1_observed %in% names(L1_cols))

L2_palette_complete <-
  all(L2_observed %in% names(L2_cols))

# Exact category conservation.
disease_categories_exact <-
  setequal(
    disease_observed,
    disease_pal$disease
  )

L1_categories_exact <-
  setequal(
    L1_observed,
    L1_pal$L1
  )

L2_categories_exact <-
  setequal(
    L2_observed,
    L2_pal$final_L2_label
  )

# ------------------------------------------------------------
# Cell-universe nesting
# ------------------------------------------------------------

L1_subset_of_disease <-
  all(L1_df$cell %in% disease_df$cell)

L2_subset_of_L1 <-
  all(L2_df$cell %in% L1_df$cell)

expected_L1_minus_L2 <-
  nrow(L1_df) - nrow(L2_df)

# Expected:
# 112 epithelial returns
# 73 neural/schwann returns
# 1 immune return
# = 186 L1-eligible but not frozen-L2-eligible cells.
L1_minus_L2_is_186 <-
  expected_L1_minus_L2 == 186

# Global exclusions:
# 82,893 - 82,680 = 213.
global_exclusions_are_213 <-
  nrow(disease_df) - nrow(L1_df) == 213

# ------------------------------------------------------------
# L2 -> L1 mapping consistency
# ------------------------------------------------------------

observed_L2_mapping <- L2_df %>%
  distinct(
    final_compartment,
    final_L2_label
  ) %>%
  mutate(
    final_compartment =
      as.character(final_compartment),
    final_L2_label =
      as.character(final_L2_label)
  )

frozen_L2_mapping <- L2_pal %>%
  transmute(
    final_compartment = as.character(L1),
    final_L2_label = as.character(final_L2_label)
  )

L2_parent_mapping_exact <-
  nrow(observed_L2_mapping) == 59 &&
  nrow(frozen_L2_mapping) == 59 &&
  nrow(
    anti_join(
      observed_L2_mapping,
      frozen_L2_mapping,
      by = c(
        "final_compartment",
        "final_L2_label"
      )
    )
  ) == 0 &&
  nrow(
    anti_join(
      frozen_L2_mapping,
      observed_L2_mapping,
      by = c(
        "final_compartment",
        "final_L2_label"
      )
    )
  ) == 0

# ------------------------------------------------------------
# Figure-file audit
# ------------------------------------------------------------

expected_stems <- c(
  "19C_01_global_UMAP_disease",
  "19C_02_global_UMAP_L1",
  "19C_03_global_UMAP_L2_59",
  sprintf(
    "19C_%02d_L2_UMAP_%s",
    seq_along(parent_order) + 3,
    parent_order
  )
)

expected_pdf <- file.path(
  FIG,
  paste0(expected_stems, ".pdf")
)

expected_png <- file.path(
  FIG,
  paste0(expected_stems, ".png")
)

all_PDF_written <- all(file.exists(expected_pdf))
all_PNG_written <- all(file.exists(expected_png))

n_expected_figures <- length(expected_stems)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "global_source_cells_82893",
    "disease_universe_82893",
    "L1_universe_82680",
    "L2_universe_82494",
    "global_exclusions_213",
    "L1_minus_L2_186",
    "disease_categories_3",
    "L1_categories_9",
    "L2_categories_59",
    "disease_palette_complete",
    "L1_palette_complete",
    "L2_palette_complete",
    "disease_categories_exact",
    "L1_categories_exact",
    "L2_categories_exact",
    "L2_parent_mapping_exact",
    "L1_subset_of_disease",
    "L2_subset_of_L1",
    "L1_coordinates_identical_to_source",
    "L2_coordinates_identical_to_source",
    "canonical_reduction_umap_harmony",
    "no_UMAP_recalculation",
    "no_Harmony_recalculation",
    "no_annotation_modified",
    "expected_11_figures",
    "all_PDF_written",
    "all_PNG_written"
  ),
  pass = c(
    nrow(umap) == 82893,
    nrow(disease_df) == 82893,
    nrow(L1_df) == 82680,
    nrow(L2_df) == 82494,
    global_exclusions_are_213,
    L1_minus_L2_is_186,
    n_distinct(disease_df$disease) == 3,
    n_distinct(L1_df$final_compartment) == 9,
    n_distinct(L2_df$final_L2_label) == 59,
    disease_palette_complete,
    L1_palette_complete,
    L2_palette_complete,
    disease_categories_exact,
    L1_categories_exact,
    L2_categories_exact,
    L2_parent_mapping_exact,
    L1_subset_of_disease,
    L2_subset_of_L1,
    coordinates_identical_L1,
    coordinates_identical_L2,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    n_expected_figures == 11,
    all_PDF_written,
    all_PNG_written
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "19C_formal_checks.csv"
  )
)

cat("\n===== DISEASE COUNTS =====\n")
print(
  as.data.frame(disease_counts),
  row.names = FALSE
)

cat("\n===== FINAL L1 COUNTS =====\n")
print(
  as.data.frame(L1_counts),
  row.names = FALSE
)

cat("\n===== FROZEN L2 FAMILY COUNTS =====\n")
print(
  as.data.frame(
    L2_pal %>%
      count(
        L1,
        name = "n_L2"
      ) %>%
      arrange(
        match(
          L1,
          L1_pal$L1
        )
      )
  ),
  row.names = FALSE
)

cat("\n===== FIGURE INVENTORY =====\n")
print(
  as.data.frame(figure_inventory),
  row.names = FALSE
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  as.data.frame(checks),
  row.names = FALSE
)

if (!all(checks$pass)) {

  failed <- checks %>%
    filter(!pass) %>%
    pull(check)

  stop(
    paste0(
      "19C FAILED formal checks: ",
      paste(failed, collapse = ", ")
    )
  )
}

cat("\n============================================================\n")
cat("19C COMPLETE\n")
cat("============================================================\n")
cat("Canonical geometry: umap_harmony\n")
cat("UMAP recalculated: NO\n")
cat("Harmony recalculated: NO\n")
cat("Disease universe: 82,893 cells\n")
cat("Final L1 universe: 82,680 cells\n")
cat("Frozen L2 universe: 82,494 cells\n")
cat("Disease colours: frozen 19B palette\n")
cat("L1 colours: frozen 19B palette\n")
cat("L2 colours: frozen 19B palette\n")
cat("Global UMAPs generated: 3\n")
cat("Compartment-specific L2 UMAPs generated: 8\n")
cat("Total figure designs: 11\n")
cat("Formats: PDF + PNG 600 dpi\n")
cat("Annotation modified: NO\n")
cat("Palette modified: NO\n")
cat("READY FOR VISUAL FIGURE QC\n")
cat("PASS\n")
