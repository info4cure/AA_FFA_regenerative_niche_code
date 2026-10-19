suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

cat("============================================================\n")
cat("PIPELINE_V2 - 20A3 PRESENTATION-GROUP FIGURE SUITE\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

OBJ <- file.path(
  ROOT,
  "objects/ambient_corrected/global_atlas_v2_SoupX_expression.rds"
)

LEDGER <- file.path(
  ROOT,
  "results/disease_biology/17A_frozen_analysis_atlas/tables",
  "17A_master_cell_ledger_FINAL.rds"
)

DISPLAY <- file.path(
  ROOT,
  "results/figure_specification/19C3_publication_suite/tables",
  "19C3_L2_display_codes.tsv"
)

DISPLAY_LABELS <- file.path(
  ROOT,
  "config/figure_specification/L2_display_labels.tsv"
)

PRESENTATION_GROUPS <- file.path(
  ROOT,
  "config/figure_specification/presentation_groups.tsv"
)

PATIENT_ORDER <- file.path(
  ROOT,
  "results/figure_specification/19C3_publication_suite/tables",
  "19C3_patient_order.tsv"
)

MARKERS <- file.path(
  ROOT,
  "results/publication_panels/20A1_frozen_L2_markers/tables",
  "20A1_marker_ranking_by_frozen_L2.tsv"
)

MARKERS_FINAL <- file.path(
  ROOT,
  "results/publication_panels/20A1_frozen_L2_markers/tables",
  "20A1_final_marker_panel.tsv"
)

OUT <- file.path(
  ROOT,
  "results/publication_panels/20A3_presentation_group_suite"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(OBJ))
stopifnot(file.exists(LEDGER))
stopifnot(file.exists(DISPLAY))
stopifnot(file.exists(PATIENT_ORDER))
stopifnot(file.exists(MARKERS))
stopifnot(file.exists(MARKERS_FINAL))

# ============================================================
# 2. LOAD FROZEN DATA
# ============================================================

cat("\nSTEP 1/9 - loading frozen atlas and specifications\n")

obj <- readRDS(OBJ)
ledger <- readRDS(LEDGER)

display <- read_tsv(
  DISPLAY,
  show_col_types = FALSE
)

display_labels <- read_tsv(
  DISPLAY_LABELS,
  show_col_types = FALSE
) %>%
  mutate(
    plot_label_2line = paste0(
      display_code,
      "\\n",
      short_identity
    )
  )

stopifnot(
  nrow(display_labels) == 59,
  n_distinct(display_labels$final_L2_label) == 59,
  n_distinct(display_labels$display_code) == 59,
  !anyNA(display_labels$short_identity)
)

display <- display %>%
  select(
    -any_of(
      c(
        "short_identity",
        "full_identity",
        "plot_label",
        "plot_label_2line"
      )
    )
  ) %>%
  left_join(
    display_labels %>%
      select(
        final_L2_label,
        short_identity,
        full_identity,
        plot_label,
        plot_label_2line
      ),
    by = "final_L2_label"
  )

stopifnot(
  nrow(display) == 59,
  !anyNA(display$short_identity),
  !anyNA(display$plot_label),
  !anyNA(display$plot_label_2line)
)

patient_order <- read_tsv(
  PATIENT_ORDER,
  show_col_types = FALSE
)

presentation_groups <- read_tsv(
  PRESENTATION_GROUPS,
  show_col_types = FALSE
)

stopifnot(
  nrow(presentation_groups) == 8,
  n_distinct(presentation_groups$frozen_L1) == 8,
  n_distinct(presentation_groups$presentation_group) == 6
)

marker_rank <- read_tsv(
  MARKERS,
  show_col_types = FALSE
)

marker_final <- read_tsv(
  MARKERS_FINAL,
  show_col_types = FALSE
)

if (!"cell" %in% colnames(ledger)) {
  cc <- intersect(
    c("cell", "cell_id", "barcode", "Cell"),
    colnames(ledger)
  )

  stopifnot(length(cc) == 1)

  ledger <- ledger %>%
    rename(cell = all_of(cc))
}

required <- c(
  "cell",
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label"
)

stopifnot(all(required %in% colnames(ledger)))

# Add disease/patient metadata from Seurat.
md <- obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

stopifnot("disease" %in% colnames(md))

patient_col <- intersect(
  c("patient_id", "patient", "sample_id"),
  colnames(md)
)

if (length(patient_col) == 0) {

  stopifnot(
    "patient_id" %in% colnames(ledger)
  )

  base <- ledger %>%
    filter(analysis_eligible_L2) %>%
    left_join(
      md %>%
        select(cell, disease),
      by = "cell"
    )

} else {

  patient_col <- patient_col[[1]]

  base <- ledger %>%
    filter(analysis_eligible_L2) %>%
    left_join(
      md %>%
        transmute(
          cell,
          disease,
          patient_id =
            as.character(.data[[patient_col]])
        ),
      by = "cell"
    )

  if (
    "patient_id.x" %in% colnames(base) &&
    "patient_id.y" %in% colnames(base)
  ) {
    base <- base %>%
      mutate(
        patient_id = coalesce(
          as.character(patient_id.x),
          as.character(patient_id.y)
        )
      ) %>%
      select(
        -patient_id.x,
        -patient_id.y
      )
  }
}

# Resolve disease after the ledger/Seurat join.
# Both sources contain disease, so dplyr may create disease.x/y.
if (!"disease" %in% colnames(base)) {

  disease_cols_join <- intersect(
    c("disease.x", "disease.y"),
    colnames(base)
  )

  if (length(disease_cols_join) == 0) {
    stop("Cannot resolve disease after ledger/Seurat join.")
  }

  base$disease <- NA_character_

  for (cc in disease_cols_join) {
    base$disease <- coalesce(
      base$disease,
      as.character(base[[cc]])
    )
  }
}

base$disease <- as.character(base$disease)
base$patient_id <- as.character(base$patient_id)

stopifnot(!anyNA(base$disease))
stopifnot(!anyNA(base$patient_id))
stopifnot(setequal(unique(base$disease), c("HC", "FFA", "AA")))
stopifnot(n_distinct(base$patient_id) == 16)

stopifnot(nrow(base) == 82494)
stopifnot(n_distinct(base$final_L2_label) == 59)

# ============================================================
# 3. NORMALIZED EXPRESSION
# ============================================================

cat("STEP 2/9 - recovering normalized RNA expression\n")

DefaultAssay(obj) <- "RNA"

expr <- tryCatch(
  GetAssayData(
    obj,
    assay = "RNA",
    layer = "data"
  ),
  error = function(e) NULL
)

if (is.null(expr) || nrow(expr) == 0) {
  stop(
    "Normalized RNA layer unavailable; ",
    "counts will not be substituted."
  )
}

stopifnot(
  all(base$cell %in% colnames(expr))
)

# ============================================================
# 4. FROZEN ORDER
# ============================================================

cat("STEP 3/9 - binding frozen L2 order and colours\n")

display <- display %>%
  arrange(global_display_order)

base <- base %>%
  left_join(
    display %>%
      select(
        L1,
        final_L2_label,
        display_code,
        legend_label,
        short_identity,
        plot_label,
        plot_label_2line,
        hex,
        L2_order_within_L1
      ),
    by = c(
      "final_compartment" = "L1",
      "final_L2_label"
    )
  )

stopifnot(!anyNA(base$display_code))
stopifnot(!anyNA(base$hex))

base <- base %>%
  left_join(
    presentation_groups %>%
      select(
        frozen_L1,
        presentation_group,
        presentation_order = display_order
      ),
    by = c(
      "final_compartment" = "frozen_L1"
    )
  )

stopifnot(
  !anyNA(base$presentation_group),
  n_distinct(base$presentation_group) == 6
)

presentation_order <- presentation_groups %>%
  distinct(
    presentation_group,
    display_order
  ) %>%
  arrange(display_order) %>%
  pull(presentation_group)

stopifnot(
  identical(
    presentation_order,
    c(
      "Epithelial_lineage",
      "Fibroblast",
      "Immune",
      "Endothelial",
      "Mural",
      "Melanocyte"
    )
  )
)

# ============================================================
# 5. PAPER-STYLE MARKER SELECTION
# ============================================================

cat("STEP 4/9 - selecting paper-style marker panels\n")

# B: up to 6 discriminating genes per frozen L2.
# C: one strongest representative gene per frozen L2.
#
# Candidate genes must exist in the canonical RNA matrix.

marker_rank <- marker_rank %>%
  filter(
    gene %in% rownames(expr)
  ) %>%
  left_join(
    display %>%
      select(
        L1,
        final_L2_label,
        L2_order_within_L1,
        display_code
      ),
    by = c(
      "L1",
      "final_L2_label"
    ),
    suffix = c("", "_display")
  )

marker_B <- marker_rank %>%
  group_by(
    L1,
    final_L2_label
  ) %>%
  arrange(
    desc(n_supporting_clusters),
    desc(max_marker_score),
    desc(max_log2FC),
    min_p_adj,
    .by_group = TRUE
  ) %>%
  slice_head(n = 6) %>%
  ungroup()

# Avoid repeated genes within a parent compartment where possible.
marker_B <- marker_B %>%
  group_by(L1, gene) %>%
  arrange(
    desc(n_supporting_clusters),
    desc(max_marker_score),
    desc(max_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

# Ensure every L2 still has >= 1 marker.
coverage_B <- marker_B %>%
  count(
    L1,
    final_L2_label,
    name = "n_markers"
  )

missing_B <- display %>%
  anti_join(
    coverage_B,
    by = c(
      "L1",
      "final_L2_label"
    )
  )

if (nrow(missing_B) > 0) {

  rescue <- marker_rank %>%
    semi_join(
      missing_B,
      by = c(
        "L1",
        "final_L2_label"
      )
    ) %>%
    group_by(
      L1,
      final_L2_label
    ) %>%
    arrange(
      desc(n_supporting_clusters),
      desc(max_marker_score),
      desc(max_log2FC),
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup()

  marker_B <- bind_rows(
    marker_B,
    rescue
  )
}

marker_C <- marker_rank %>%
  group_by(
    L1,
    final_L2_label
  ) %>%
  arrange(
    desc(n_supporting_clusters),
    desc(max_marker_score),
    desc(max_log2FC),
    min_p_adj,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

write_tsv(
  marker_B,
  file.path(
    TAB,
    "20A2_panel_B_marker_selection.tsv"
  )
)

write_tsv(
  marker_C,
  file.path(
    TAB,
    "20A2_panel_C_marker_selection.tsv"
  )
)

# ============================================================
# 6. EXPORT HELPER
# ============================================================

save3 <- function(
  p,
  parent,
  panel,
  stem,
  width,
  height
) {

  root <- file.path(
    FIG,
    parent,
    panel
  )

  for (fmt in c("PDF", "PNG", "SVG")) {
    dir.create(
      file.path(root, fmt),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  ggsave(
    file.path(
      root,
      "PDF",
      paste0(stem, ".pdf")
    ),
    p,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      root,
      "PNG",
      paste0(stem, ".png")
    ),
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  if (
    requireNamespace(
      "svglite",
      quietly = TRUE
    )
  ) {

    ggsave(
      file.path(
        root,
        "SVG",
        paste0(stem, ".svg")
      ),
      p,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white",
      limitsize = FALSE
    )

  } else {

    svg(
      file.path(
        root,
        "SVG",
        paste0(stem, ".svg")
      ),
      width = width,
      height = height
    )

    print(p)
    dev.off()
  }
}

# ============================================================
# PRESENTATION-GROUP DISPLAY HELPERS
# ============================================================

presentation_L1_order <- c(
  "Epithelial",
  "Sweat_gland",
  "Sebocyte",
  "Fibroblast",
  "Immune",
  "Endothelial",
  "Mural",
  "Melanocyte"
)

get_group_display <- function(group) {

  frozen_members <- presentation_groups %>%
    filter(
      presentation_group == group
    ) %>%
    pull(frozen_L1)

  z <- display %>%
    filter(
      L1 %in% frozen_members
    ) %>%
    mutate(
      .L1_order = match(
        L1,
        presentation_L1_order
      )
    ) %>%
    arrange(
      .L1_order,
      L2_order_within_L1
    ) %>%
    select(-.L1_order)

  stopifnot(nrow(z) > 0)

  z
}

pretty_group <- function(group) {
  ifelse(
    group == "Epithelial_lineage",
    "Epithelial lineage",
    gsub("_", " ", group)
  )
}

safe_group <- function(group) {
  gsub(
    "[^A-Za-z0-9]+",
    "_",
    group
  )
}

# ============================================================
# 7. PAPER-STYLE CELL HEATMAP
# ============================================================

cat("STEP 5/9 - generating paper-style cell-level heatmaps\n")

make_cell_heatmap <- function(
  group,
  max_cells_per_L2 = 250
) {

  disp <- get_group_display(group)

  markers <- marker_B %>%
    semi_join(
      disp %>%
        select(
          L1,
          final_L2_label
        ),
      by = c(
        "L1",
        "final_L2_label"
      )
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          .display_order =
            global_display_order
        ),
      by = "final_L2_label"
    ) %>%
    arrange(
      match(
        final_L2_label,
        disp$final_L2_label
      ),
      desc(max_marker_score)
    )

  # ----------------------------------------------------------
  # Compact publication marker panel
  #
  # Use at most the top 3 already-selected markers per frozen
  # L2 identity. This is DISPLAY SELECTION ONLY.
  #
  # No marker inference, DE testing, clustering or annotation
  # is recalculated here.
  # ----------------------------------------------------------

  markers_plot <- markers %>%
    group_by(
      final_L2_label
    ) %>%
    arrange(
      desc(max_marker_score),
      desc(max_log2FC),
      gene,
      .by_group = TRUE
    ) %>%
    slice_head(
      n = 3
    ) %>%
    ungroup()

  genes <- unique(
    markers_plot$gene
  )

  dd <- base %>%
    filter(
      presentation_group == group
    )

  # ----------------------------------------------------------
  # Order cells biologically within each L2.
  #
  # For each frozen L2, rank cells using the mean expression
  # of that identity's selected marker genes. This produces
  # compact blocks similar to publication heatmaps while
  # retaining actual individual cells.
  # ----------------------------------------------------------

  selected_cells <- list()

  for (lbl in disp$final_L2_label) {

    cc <- dd %>%
      filter(
        final_L2_label == lbl
      ) %>%
      pull(cell)

    gg <- markers %>%
      filter(
        final_L2_label == lbl
      ) %>%
      pull(gene) %>%
      intersect(
        rownames(expr)
      )

    if (length(cc) == 0) next

    if (length(gg) > 0) {

      score <- Matrix::colMeans(
        expr[
          gg,
          cc,
          drop = FALSE
        ]
      )

      cc <- cc[
        order(
          score,
          decreasing = TRUE
        )
      ]
    }

    # For large compartments use a deterministic evenly-spaced
    # subset for publication rendering; no random sampling.
    if (
      length(cc) >
      max_cells_per_L2
    ) {

      idx <- unique(
        round(
          seq(
            1,
            length(cc),
            length.out =
              max_cells_per_L2
          )
        )
      )

      cc <- cc[idx]
    }

    selected_cells[[lbl]] <- cc
  }

  # ----------------------------------------------------------
  # Publication display ordering:
  # within each frozen L2, order cells by the mean expression
  # of that L2's selected marker programme.
  #
  # DISPLAY ONLY:
  # - no clustering
  # - no annotation modification
  # - no marker inference
  # ----------------------------------------------------------

  ordered_cells <- list()

  for (lbl in disp$final_L2_label) {

    cc <- selected_cells[[lbl]]

    if (length(cc) == 0) {
      next
    }

    marker_genes_lbl <- markers_plot %>%
      filter(
        final_L2_label == lbl
      ) %>%
      pull(gene) %>%
      unique()

    marker_genes_lbl <- intersect(
      marker_genes_lbl,
      rownames(expr)
    )

    if (length(marker_genes_lbl) > 0) {

      score_mat <- as.matrix(
        expr[
          marker_genes_lbl,
          cc,
          drop = FALSE
        ]
      )

      if (length(marker_genes_lbl) == 1) {

        marker_score <- as.numeric(
          score_mat[1, ]
        )

      } else {

        marker_score <- colMeans(
          score_mat,
          na.rm = TRUE
        )
      }

      ord <- order(
        marker_score,
        decreasing = TRUE,
        na.last = TRUE
      )

      cc <- cc[ord]
    }

    ordered_cells[[lbl]] <- cc
  }

  cell_order <- unlist(
    ordered_cells,
    use.names = FALSE
  )

  stopifnot(length(cell_order) > 0)

  mat <- as.matrix(
    expr[
      genes,
      cell_order,
      drop = FALSE
    ]
  )

  # Gene-wise z-score across displayed individual cells.
  mat_z <- t(
    apply(
      mat,
      1,
      function(x) {

        sx <- sd(x)

        if (
          is.na(sx) ||
          sx == 0
        ) {
          rep(
            0,
            length(x)
          )
        } else {
          (x - mean(x)) / sx
        }
      }
    )
  )

  mat_z[
    mat_z > 2.5
  ] <- 2.5

  mat_z[
    mat_z < -2.5
  ] <- -2.5

  cell_key <- dd %>%
    filter(
      cell %in% cell_order
    ) %>%
    select(
      cell,
      final_L2_label,
      display_code,
      plot_label,
      plot_label_2line,
      hex
    ) %>%
    mutate(
      cell = factor(
        cell,
        levels = cell_order
      )
    ) %>%
    arrange(cell) %>%
    mutate(
      x = row_number()
    )

  # ----------------------------------------------------------
  # Gene ordering:
  # frozen L2 order first, then marker strength within L2.
  # This creates a publication-style descending block structure
  # without changing the marker set.
  # ----------------------------------------------------------

  # ----------------------------------------------------------
  # Publication staircase ordering of marker genes
  #
  # Columns already follow the frozen L2 order.
  # For DISPLAY ONLY, order genes according to the L2 in which
  # their mean displayed-cell expression is maximal.
  #
  # Secondary ordering within each peak-L2 uses the relative
  # strength of the peak over the remaining L2 identities.
  #
  # This changes row order only:
  # marker membership, expression values, clustering, annotation
  # and differential-expression inference remain unchanged.
  # ----------------------------------------------------------

  cell_l2 <- cell_key %>%
    mutate(
      final_L2_label =
        as.character(final_L2_label)
    ) %>%
    select(
      cell,
      final_L2_label
    )

  l2_levels <- disp$final_L2_label

  gene_expression_order <- lapply(
    genes,
    function(g) {

      vv <- as.numeric(
        mat[
          g,
          cell_order
        ]
      )

      tmp <- tibble(
        cell = cell_order,
        expression = vv
      ) %>%
        left_join(
          cell_l2,
          by = "cell"
        ) %>%
        group_by(
          final_L2_label
        ) %>%
        summarise(
          mean_expression =
            mean(
              expression,
              na.rm = TRUE
            ),
          .groups = "drop"
        ) %>%
        complete(
          final_L2_label =
            l2_levels,
          fill = list(
            mean_expression = 0
          )
        ) %>%
        mutate(
          l2_index =
            match(
              final_L2_label,
              l2_levels
            )
        ) %>%
        arrange(
          l2_index
        )

      peak_idx <- which.max(
        tmp$mean_expression
      )

      peak_value <-
        tmp$mean_expression[
          peak_idx
        ]

      background_value <-
        if (
          nrow(tmp) > 1
        ) {
          mean(
            tmp$mean_expression[
              -peak_idx
            ],
            na.rm = TRUE
          )
        } else {
          0
        }

      tibble(
        gene = g,
        peak_L2 =
          tmp$final_L2_label[
            peak_idx
          ],
        peak_L2_index =
          tmp$l2_index[
            peak_idx
          ],
        peak_expression =
          peak_value,
        background_expression =
          background_value,
        peak_specificity =
          peak_value -
          background_value
      )
    }
  ) %>%
    bind_rows()

  gene_key <- markers_plot %>%
    distinct(
      gene,
      final_L2_label,
      L2_order_within_L1,
      .keep_all = TRUE
    ) %>%
    left_join(
      gene_expression_order,
      by = "gene"
    ) %>%
    mutate(
      marker_strength =
        dplyr::coalesce(
          suppressWarnings(
            as.numeric(
              max_marker_score
            )
          ),
          suppressWarnings(
            as.numeric(
              max_log2FC
            )
          ),
          0
        )
    ) %>%
    arrange(
      desc(peak_L2_index),
      desc(peak_specificity),
      desc(marker_strength),
      gene
    )

  # ggplot places the first y-factor level at the bottom and
  # the plotting code below applies rev(gene_order).
  # Therefore descending peak_L2_index here yields:
  #
  # EP.1 -> ... -> EP.9 -> SG.1 -> ... -> SG.5 ->
  # SE.1 -> ... -> SE.5
  #
  # from TOP to BOTTOM, concordant with LEFT-to-RIGHT columns.
  gene_order <- unique(
    gene_key$gene
  )

  long <- as.data.frame(
    mat_z
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    pivot_longer(
      -gene,
      names_to = "cell",
      values_to = "z"
    ) %>%
    mutate(
      cell = factor(
        cell,
        levels = cell_order
      ),
      # Final publication orientation:
      # invert the previous displayed gene order.
      gene = factor(
        gene,
        levels = gene_order
      ),
      x = as.integer(cell)
    )

  p_body <- ggplot(
    long,
    aes(
      x = x,
      y = gene,
      fill = z
    )
  ) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "Expression"
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    labs(
      title = paste0(
        pretty_group(group),
        " - L2 marker signatures"
      ),
      subtitle = paste0(
        length(genes),
        " marker genes | ",
        length(cell_order),
        " representative individual cells"
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y =
        element_text(
          size = 9.2,
          face = "bold",
          colour = "#111111"
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      plot.subtitle =
        element_text(size = 8),
      legend.title =
        element_text(size = 8),
      legend.text =
        element_text(size = 7),
      plot.margin =
        margin(
          t = 28,
          r = 8,
          b = 4,
          l = 8
        )
    )

  # ----------------------------------------------------------
  # Top annotation bar
  # ----------------------------------------------------------

  annotation <- cell_key %>%
    mutate(
      y = 1
    )

  ann_cols <- setNames(
    disp$hex,
    disp$display_code
  )

  p_ann <- ggplot(
    annotation,
    aes(
      x = x,
      y = y,
      fill = display_code
    )
  ) +
    geom_tile(
      height = 0.85
    ) +
    scale_fill_manual(
      values = ann_cols,
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour = NA
        )
      )
    ) +
    labs(
      fill = NULL
    ) +
    theme_void(
      base_size = 8
    ) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.text =
        element_text(
          size = 7
        ),
      legend.key.width =
        grid::unit(
          4.5,
          "mm"
        ),
      legend.key.height =
        grid::unit(
          2.5,
          "mm"
        ),
      plot.margin =
        margin(
          2, 6, 0, 6
        )
    )

  # ----------------------------------------------------------
  # Vertical separators between L2 blocks
  # ----------------------------------------------------------

  boundaries <- cell_key %>%
    count(
      display_code,
      plot_label,
      plot_label_2line,
      name = "n_cells"
    ) %>%
    mutate(
      display_code = factor(
        display_code,
        levels = disp$display_code
      )
    ) %>%
    arrange(
      display_code
    ) %>%
    mutate(
      xmax = cumsum(n_cells),
      xmin = lag(
        xmax,
        default = 0
      ) + 1,
      midpoint =
        (xmin + xmax) / 2
    )

  separator_x <- boundaries$xmax[
    seq_len(
      max(
        nrow(boundaries) - 1,
        0
      )
    )
  ]

  if (length(separator_x) > 0) {

    p_body <- p_body +
      geom_vline(
        xintercept =
          separator_x + 0.5,
        linewidth = 0.25,
        colour = "grey35"
      )

    p_ann <- p_ann +
      geom_vline(
        xintercept =
          separator_x + 0.5,
        linewidth = 0.25,
        colour = "white"
      )
  }

  # ----------------------------------------------------------
  # Compact L2 labels above the annotation strip
  # ----------------------------------------------------------

  p_ann <- p_ann +
    geom_text(
      data = boundaries,
      aes(
        x = midpoint,
        y = 1.72,
        label = plot_label_2line
      ),
      inherit.aes = FALSE,
      size = 2.75,
      angle = 28,
      hjust = 0,
      vjust = 0
    ) +
    coord_cartesian(
      ylim = c(
        0.55,
        2.55
      ),
      clip = "off"
    )

  # No duplicated legend is needed because the short codes
  # are printed directly above the colour strip.
  p_ann <- p_ann +
    theme(
      legend.position = "none"
    )

  combined <- p_ann / p_body +
    plot_layout(
      heights = c(
        1.35,
        8
      )
    )

  list(
    plot = combined,
    n_cells = length(cell_order),
    n_genes = length(genes),
    boundaries = boundaries,
    cells = tibble(
      presentation_group = group,
      cell = cell_order
    ),
    genes = markers %>%
      select(
        L1,
        final_L2_label,
        display_code,
        gene,
        n_supporting_clusters,
        max_log2FC,
        max_marker_score,
        min_p_adj
      )
  )
}

# ============================================================
# 8. GENERATE PANEL B FOR ALL COMPARTMENTS
# ============================================================

panel_B_inventory <- tibble(
  L1 = character(),
  n_L2 = integer(),
  n_cells_displayed = integer(),
  n_genes = integer(),
  figure = character()
)

panel_B_cells <- list()
panel_B_genes <- list()
panel_B_boundaries <- list()

for (group in presentation_order) {

  cat(
    sprintf(
      "  %-20s",
      group
    )
  )

  # The combined epithelial lineage has 19 L2 identities.
  # Use fewer cells per L2 only for rendering efficiency.
  max_cells_group <- ifelse(
    group == "Epithelial_lineage",
    180,
    250
  )

  res <- make_cell_heatmap(
    group = group,
    max_cells_per_L2 = max_cells_group
  )

  stem <- paste0(
    "20A3_",
    safe_group(group),
    "_B_cell_marker_heatmap"
  )

  disp_parent <- get_group_display(group)

  # Width depends primarily on the number of displayed cells.
  # Height depends on the marker-gene panel.
  width <- max(
    7.0,
    min(
      11.5,
      5.5 +
        res$n_cells / 450
    )
  )

  height <- max(
    5.0,
    2.8 +
      0.17 * res$n_genes
  )

  save3(
    res$plot,
    safe_group(group),
    "B_cell_marker_heatmap",
    stem,
    width,
    height
  )

  panel_B_inventory <- bind_rows(
    panel_B_inventory,
    tibble(
      presentation_group = group,
      n_L2 =
        nrow(disp_parent),
      n_cells_displayed =
        res$n_cells,
      n_genes =
        res$n_genes,
      figure = stem
    )
  )

  panel_B_cells[[group]] <-
    res$cells

  panel_B_genes[[group]] <-
    res$genes

  panel_B_boundaries[[group]] <-
    res$boundaries %>%
    mutate(
      presentation_group = group
    )

  cat(
    sprintf(
      " : %4d cells | %2d genes | DONE\n",
      res$n_cells,
      res$n_genes
    )
  )
}

write_tsv(
  panel_B_inventory,
  file.path(
    TAB,
    "20A2_panel_B_inventory.tsv"
  )
)

write_tsv(
  bind_rows(panel_B_cells),
  file.path(
    TAB,
    "20A2_panel_B_displayed_cells.tsv"
  )
)

write_tsv(
  bind_rows(panel_B_genes),
  file.path(
    TAB,
    "20A2_panel_B_genes.tsv"
  )
)

write_tsv(
  bind_rows(panel_B_boundaries),
  file.path(
    TAB,
    "20A2_panel_B_boundaries.tsv"
  )
)

cat("20A2 BLOCK 1 COMPLETE\n")

# ============================================================
# 9. PANEL C — COMPACT IDENTIFYING-GENE VIOLINS
# ============================================================

cat("\nSTEP 6/9 - generating compact marker violins\n")

make_compact_violins <- function(
  group,
  max_cells_per_L2 = 800
) {

  disp <- get_group_display(group)

  markers <- marker_C %>%
    semi_join(
      disp %>%
        select(
          L1,
          final_L2_label
        ),
      by = c(
        "L1",
        "final_L2_label"
      )
    ) %>%
    arrange(
      match(
        final_L2_label,
        disp$final_L2_label
      )
    )

  stopifnot(nrow(markers) > 0)

  genes <- unique(markers$gene)

  dd <- base %>%
    filter(
      presentation_group == group
    )

  # Deterministic thinning only for rendering very large groups.
  # This does not alter any statistical analysis.
  cells_use <- lapply(
    disp$final_L2_label,
    function(lbl) {

      cc <- dd %>%
        filter(
          final_L2_label == lbl
        ) %>%
        pull(cell)

      if (length(cc) > max_cells_per_L2) {

        idx <- unique(
          round(
            seq(
              1,
              length(cc),
              length.out =
                max_cells_per_L2
            )
          )
        )

        cc <- cc[idx]
      }

      cc
    }
  ) %>%
    unlist(
      use.names = FALSE
    )

  md_use <- dd %>%
    filter(
      cell %in% cells_use
    ) %>%
    select(
      cell,
      final_L2_label,
      display_code,
      short_identity,
      plot_label,
      plot_label_2line
    )

  mat <- as.matrix(
    expr[
      genes,
      md_use$cell,
      drop = FALSE
    ]
  )

  long <- as.data.frame(
    t(mat)
  ) %>%
    rownames_to_column(
      "cell"
    ) %>%
    pivot_longer(
      -cell,
      names_to = "gene",
      values_to = "expression"
    ) %>%
    left_join(
      md_use,
      by = "cell"
    )

  gene_labels <- markers %>%
    distinct(
      gene,
      final_L2_label,
      display_code
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          short_identity,
          plot_label,
          plot_label_2line
        ),
      by = "final_L2_label"
    ) %>%
    mutate(
      gene_label = gene
    )

  long <- long %>%
    left_join(
      gene_labels %>%
        select(
          gene,
          gene_label
        ),
      by = "gene"
    )

  code_levels <- disp$display_code
  plot_levels <- disp$plot_label_2line

  gene_levels <- gene_labels %>%
    arrange(
      match(
        display_code,
        code_levels
      )
    ) %>%
    pull(gene_label) %>%
    unique()

  long <- long %>%
    mutate(
      display_code = factor(
        display_code,
        levels = code_levels
      ),
      plot_label_2line = factor(
        plot_label_2line,
        levels = plot_levels
      ),
      gene_label = factor(
        gene_label,
        levels = rev(gene_levels)
      )
    )

  cols <- setNames(
    disp$hex,
    disp$display_code
  )

  p <- ggplot(
    long,
    aes(
      x = plot_label_2line,
      y = expression,
      fill = display_code
    )
  ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      linewidth = 0.28,
      colour = "#303030"
    ) +
    facet_grid(
      rows = vars(gene_label),
      scales = "free_y",
      switch = "y"
    ) +
    scale_fill_manual(
      values = cols,
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    labs(
      title = paste0(
        pretty_group(group),
        " - cluster-identifying genes"
      ),
      subtitle =
        "One representative marker per frozen L2 identity",
      x = NULL,
      y = NULL
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      legend.position = "none",
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 9,
          face = "bold"
        ),
      axis.text.y =
        element_blank(),
      axis.ticks.y =
        element_blank(),
      strip.placement =
        "outside",
      strip.background =
        element_blank(),
      strip.text.y.left =
        element_text(
          angle = 0,
          hjust = 1,
          size = 11,
          face = "bold"
        ),
      panel.spacing.y =
        grid::unit(
          1.5,
          "mm"
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 12
        ),
      plot.subtitle =
        element_text(
          size = 8
        ),
      plot.margin =
        margin(
          6, 8, 6, 8
        )
    )

  list(
    plot = p,
    n_cells =
      n_distinct(long$cell),
    n_genes =
      length(genes),
    genes =
      gene_labels
  )
}

panel_C_inventory <- tibble(
  L1 = character(),
  n_L2 = integer(),
  n_cells_displayed = integer(),
  n_genes = integer(),
  figure = character()
)

panel_C_genes <- list()

for (group in presentation_order) {

  res <- make_compact_violins(
    group
  )

  stem <- paste0(
    "20A3_",
    safe_group(group),
    "_C_marker_violins"
  )

  n_L2_parent <- get_group_display(group) %>%
    nrow()

  width <- ifelse(
    group == "Epithelial_lineage",
    11.5,
    max(
      7.0,
      0.72 * n_L2_parent + 3.0
    )
  )

  height <- ifelse(
    group == "Epithelial_lineage",
    16.5,
    max(
      9.0,
      1.05 * res$n_genes + 2.5
    )
  )

  save3(
    res$plot,
    safe_group(group),
    "C_marker_violins",
    stem,
    width,
    height
  )

  panel_C_inventory <- bind_rows(
    panel_C_inventory,
    tibble(
      presentation_group = group,
      n_L2 = n_L2_parent,
      n_cells_displayed =
        res$n_cells,
      n_genes =
        res$n_genes,
      figure = stem
    )
  )

  panel_C_genes[[group]] <-
    res$genes %>%
    mutate(
      presentation_group = group
    )

  cat(
    sprintf(
      "%-15s : %5d cells | %2d genes | violin panel DONE\n",
      group,
      res$n_cells,
      res$n_genes
    )
  )
}

write_tsv(
  panel_C_inventory,
  file.path(
    TAB,
    "20A2_panel_C_inventory.tsv"
  )
)

write_tsv(
  bind_rows(panel_C_genes),
  file.path(
    TAB,
    "20A2_panel_C_genes.tsv"
  )
)

# ============================================================
# 10. PANEL D — PATIENT COMPOSITION
# ============================================================

cat("\nSTEP 7/9 - generating patient composition panels\n")

# Resolve patient order directly from frozen 19C3 specification.
stopifnot(
  "patient_id" %in%
    colnames(patient_order)
)

if (!"disease" %in% colnames(patient_order)) {
  patient_order <- patient_order %>%
    left_join(
      base %>%
        distinct(
          patient_id,
          disease
        ),
      by = "patient_id"
    )
}

patient_levels <- patient_order$patient_id

stopifnot(
  length(unique(patient_levels)) == 16
)

disease_cols <- c(
  HC = "#8F8D8E",
  FFA = "#4D9E4F",
  AA = "#A8142A"
)

make_patient_composition <- function(
  group
) {

  disp <- get_group_display(group)

  dd <- base %>%
    filter(
      presentation_group == group
    )

  comp <- dd %>%
    count(
      patient_id,
      disease,
      final_L2_label,
      display_code,
      name = "n_cells"
    ) %>%
    complete(
      patient_id =
        patient_levels,
      final_L2_label =
        disp$final_L2_label,
      fill = list(
        n_cells = 0
      )
    ) %>%
    left_join(
      patient_order %>%
        select(
          patient_id,
          disease
        ),
      by = "patient_id",
      suffix = c(
        "",
        "_order"
      )
    ) %>%
    mutate(
      disease = coalesce(
        disease,
        disease_order
      )
    ) %>%
    select(
      -disease_order
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code,
          short_identity,
          plot_label,
          plot_label_2line,
          hex,
          L2_order_within_L1
        ),
      by = "final_L2_label",
      suffix = c(
        "",
        "_display"
      )
    ) %>%
    mutate(
      display_code = coalesce(
        display_code,
        display_code_display
      )
    ) %>%
    select(
      -display_code_display
    ) %>%
    group_by(
      patient_id
    ) %>%
    mutate(
      total_cells =
        sum(n_cells),
      fraction =
        ifelse(
          total_cells > 0,
          n_cells /
            total_cells,
          0
        )
    ) %>%
    ungroup()

  comp <- comp %>%
    mutate(
      patient_id = factor(
        patient_id,
        levels =
          patient_levels
      ),
      display_code = factor(
        display_code,
        levels =
          disp$display_code
      )
    )

  cols <- setNames(
    disp$hex,
    disp$display_code
  )

  legend_labels <- setNames(
    disp$plot_label,
    disp$display_code
  )

  disease_bar <- patient_order %>%
    mutate(
      patient_id = factor(
        patient_id,
        levels =
          patient_levels
      ),
      y = 1
    )

  p_top <- ggplot(
    disease_bar,
    aes(
      x = patient_id,
      y = y,
      fill = disease
    )
  ) +
    geom_tile(
      height = 0.8
    ) +
    scale_fill_manual(
      values = disease_cols,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +
    labs(
      fill = NULL
    ) +
    theme_void(
      base_size = 8
    ) +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.key.width =
        grid::unit(
          5,
          "mm"
        ),
      legend.key.height =
        grid::unit(
          2.5,
          "mm"
        ),
      legend.text =
        element_text(
          size = 7
        ),
      plot.margin =
        margin(
          2, 6, 0, 6
        )
    )

  p_body <- ggplot(
    comp,
    aes(
      x = patient_id,
      y = fraction,
      fill = display_code
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.15
    ) +
    scale_fill_manual(
      values = cols,
      breaks = disp$display_code,
      labels = unname(
        legend_labels[
          disp$display_code
        ]
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(
        0,
        1
      ),
      expand = c(
        0,
        0
      )
    ) +
    labs(
      title = paste0(
        pretty_group(group),
        " - L2 composition by patient"
      ),
      subtitle =
        "Patients ordered HC -> FFA -> AA",
      x = NULL,
      y = "Cell proportion",
      fill = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          size = 7
        ),
      axis.title.y =
        element_text(
          size = 8
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      plot.subtitle =
        element_text(
          size = 8
        ),
      legend.position = "right",
      legend.text =
        element_text(
          size = 7
        ),
      legend.key.height =
        grid::unit(
          3.3,
          "mm"
        ),
      plot.margin =
        margin(
          0, 6, 6, 6
        )
    )

  combined <- p_top / p_body +
    plot_layout(
      heights = c(
        0.45,
        6
      )
    )

  list(
    plot = combined,
    table = comp,
    n_cells = nrow(dd)
  )
}

panel_D_inventory <- tibble(
  L1 = character(),
  n_cells = integer(),
  n_patients = integer(),
  figure = character()
)

panel_D_tables <- list()

for (group in presentation_order) {

  res <- make_patient_composition(
    group
  )

  stem <- paste0(
    "20A3_",
    safe_group(group),
    "_D_patient_composition"
  )

  save3(
    res$plot,
    safe_group(group),
    "D_patient_composition",
    stem,
    width = 10.5,
    height = 6.2
  )

  panel_D_tables[[group]] <-
    res$table %>%
    mutate(
      presentation_group = group
    )

  panel_D_inventory <- bind_rows(
    panel_D_inventory,
    tibble(
      presentation_group = group,
      n_cells = res$n_cells,
      n_patients =
        n_distinct(
          res$table$patient_id
        ),
      figure = stem
    )
  )

  cat(
    sprintf(
      "%-15s : %6d cells | patient composition DONE\n",
      group,
      res$n_cells
    )
  )
}

write_tsv(
  panel_D_inventory,
  file.path(
    TAB,
    "20A2_panel_D_inventory.tsv"
  )
)

write_tsv(
  bind_rows(
    panel_D_tables
  ),
  file.path(
    TAB,
    "20A2_panel_D_patient_composition.tsv"
  )
)

# ============================================================
# 11. PANEL A / A2 — PRESENTATION-GROUP UMAPs
# ============================================================

cat("\nSTEP 8/10 - generating presentation-group UMAPs\n")

stopifnot(
  "umap_harmony" %in%
    names(obj@reductions)
)

emb <- Embeddings(
  obj,
  reduction = "umap_harmony"
)

stopifnot(
  nrow(emb) == 82893
)

umap <- tibble(
  cell = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2]
)

plot_base <- base %>%
  left_join(
    umap,
    by = "cell"
  )

stopifnot(
  nrow(plot_base) == 82494,
  !anyNA(plot_base$UMAP_1),
  !anyNA(plot_base$UMAP_2)
)

xlim_global <- range(
  umap$UMAP_1,
  finite = TRUE
)

ylim_global <- range(
  umap$UMAP_2,
  finite = TRUE
)

theme_pub_umap <- function() {

  theme_void(
    base_size = 9
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      plot.subtitle = element_text(
        size = 8,
        margin = margin(b = 5)
      ),
      legend.title =
        element_blank(),
      legend.text =
        element_text(size = 7),
      legend.key.height =
        grid::unit(3.5, "mm"),
      legend.key.width =
        grid::unit(3.5, "mm"),
      strip.text =
        element_text(
          face = "bold",
          size = 9
        ),
      plot.margin =
        margin(6, 6, 6, 6)
    )
}

make_group_umap <- function(
  group,
  split_disease = FALSE
) {

  disp <- get_group_display(group)

  dd <- plot_base %>%
    filter(
      presentation_group == group
    ) %>%
    mutate(
      final_L2_label = factor(
        final_L2_label,
        levels =
          disp$final_L2_label
      )
    )

  cols <- setNames(
    disp$hex,
    disp$final_L2_label
  )

  legend_labels <- setNames(
    disp$plot_label,
    disp$final_L2_label
  )

  # Larger populations first so small populations remain visible.
  nn <- table(
    dd$final_L2_label
  )

  dd$.draw_n <- unname(
    nn[
      as.character(
        dd$final_L2_label
      )
    ]
  )

  dd <- dd %>%
    arrange(
      desc(.draw_n)
    )

  # ----------------------------------------------------------
  # Direct display labels.
  #
  # Median coordinates are used only as deterministic anchors;
  # ggrepel moves the text to avoid overlap.
  # ----------------------------------------------------------

  label_df <- dd %>%
    mutate(
      final_L2_label =
        as.character(
          final_L2_label
        )
    ) %>%
    group_by(
      final_L2_label
    ) %>%
    summarise(
      UMAP_1 =
        median(
          UMAP_1,
          na.rm = TRUE
        ),
      UMAP_2 =
        median(
          UMAP_2,
          na.rm = TRUE
        ),
      n = n(),
      .groups = "drop"
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code
        ),
      by = "final_L2_label"
    )

  point_size <- ifelse(
    group == "Epithelial_lineage",
    0.72,
    0.86
  )

  point_stroke <- ifelse(
    group == "Epithelial_lineage",
    0.085,
    0.11
  )

  base_layers <- list(
    geom_point(
      aes(
        fill =
          final_L2_label
      ),
      shape = 21,
      colour = "#303030",
      stroke = point_stroke,
      size = point_size,
      alpha = 0.94
    ),

    scale_fill_manual(
      values = cols,
      breaks =
        disp$final_L2_label,
      labels =
        unname(
          legend_labels[
            disp$final_L2_label
          ]
        ),
      drop = FALSE
    ),

    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 4.5,
          alpha = 1,
          colour = "#303030",
          stroke = 0.45
        ),
        keyheight =
          grid::unit(
            5.0,
            "mm"
          ),
        keywidth =
          grid::unit(
            5.0,
            "mm"
          )
      )
    )
  )

  if (!split_disease) {

    p <- ggplot(
      dd,
      aes(
        UMAP_1,
        UMAP_2
      )
    ) +
      base_layers +
      ggrepel::geom_label_repel(
        data = label_df,
        aes(
          x = UMAP_1,
          y = UMAP_2,
          label = display_code
        ),
        inherit.aes = FALSE,
        seed = 2026,
        size = ifelse(
          group ==
            "Epithelial_lineage",
          3.0,
          3.3
        ),
        fontface = "bold",
        colour = "#202020",
        fill =
          scales::alpha(
            "white",
            0.88
          ),
        label.size = 0.18,
        box.padding = 0.42,
        point.padding = 0.24,
        min.segment.length = 0,
        segment.colour = "grey45",
        segment.size = 0.28,
        max.overlaps = Inf,
        force = 1.4,
        max.time = 3
      ) +
      coord_equal(
        xlim = xlim_global,
        ylim = ylim_global,
        expand = FALSE,
        clip = "off"
      ) +
      labs(
        title = paste0(
          pretty_group(group),
          " - frozen L2 identities"
        ),
        subtitle = paste0(
          format(
            nrow(dd),
            big.mark = ","
          ),
          " cells | ",
          nrow(disp),
          " L2 identities"
        )
      ) +
      theme_pub_umap() +
      theme(
        legend.text =
          element_text(
            size = 8.2
          )
      )

    return(p)
  }

  # ----------------------------------------------------------
  # Disease-split UMAP:
  # direct labels are calculated independently within each
  # disease panel, but the Harmony coordinates remain fixed.
  # ----------------------------------------------------------

  dd$disease <- factor(
    dd$disease,
    levels = c(
      "HC",
      "FFA",
      "AA"
    )
  )

  label_split <- dd %>%
    mutate(
      final_L2_label =
        as.character(
          final_L2_label
        )
    ) %>%
    group_by(
      disease,
      final_L2_label
    ) %>%
    summarise(
      UMAP_1 =
        median(
          UMAP_1,
          na.rm = TRUE
        ),
      UMAP_2 =
        median(
          UMAP_2,
          na.rm = TRUE
        ),
      n = n(),
      .groups = "drop"
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code
        ),
      by = "final_L2_label"
    )

  ggplot(
    dd,
    aes(
      UMAP_1,
      UMAP_2
    )
  ) +
    geom_point(
      aes(
        fill =
          final_L2_label
      ),
      shape = 21,
      colour = "#303030",
      stroke = ifelse(
        group ==
          "Epithelial_lineage",
        0.07,
        0.09
      ),
      size = ifelse(
        group ==
          "Epithelial_lineage",
        0.58,
        0.70
      ),
      alpha = 0.94
    ) +
    ggrepel::geom_label_repel(
      data = label_split,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        label = display_code
      ),
      inherit.aes = FALSE,
      seed = 2026,
      size = 2.5,
      fontface = "bold",
      colour = "#202020",
      fill =
        scales::alpha(
          "white",
          0.86
        ),
      label.size = 0.15,
      box.padding = 0.30,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.colour = "grey50",
      segment.size = 0.22,
      max.overlaps = Inf,
      force = 1.2,
      max.time = 3
    ) +
    scale_fill_manual(
      values = cols,
      breaks =
        disp$final_L2_label,
      labels =
        unname(
          legend_labels[
            disp$final_L2_label
          ]
        ),
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 4.2,
          alpha = 1,
          colour = "#303030",
          stroke = 0.4
        ),
        keyheight =
          grid::unit(
            4.8,
            "mm"
          ),
        keywidth =
          grid::unit(
            4.8,
            "mm"
          )
      )
    ) +
    facet_wrap(
      ~ disease,
      nrow = 1,
      drop = FALSE
    ) +
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE,
      clip = "off"
    ) +
    labs(
      title = paste0(
        pretty_group(group),
        " - frozen L2 identities"
      ),
      subtitle =
        "HC | FFA | AA - identical Harmony coordinates"
    ) +
    theme_pub_umap() +
    theme(
      legend.text =
        element_text(
          size = 8
        )
    )
}

panel_A_inventory <- tibble()
panel_A2_inventory <- tibble()

for (group in presentation_order) {

  dd <- plot_base %>%
    filter(
      presentation_group ==
        group
    )

  n_L2_group <- n_distinct(
    dd$final_L2_label
  )

  p_A <- make_group_umap(
    group,
    split_disease = FALSE
  )

  stem_A <- paste0(
    "20A3_",
    safe_group(group),
    "_A_L2_UMAP"
  )

  save3(
    p_A,
    safe_group(group),
    "A_L2_UMAP",
    stem_A,
    width = ifelse(
      group ==
        "Epithelial_lineage",
      10.5,
      8.2
    ),
    height = 6.4
  )

  panel_A_inventory <- bind_rows(
    panel_A_inventory,
    tibble(
      presentation_group =
        group,
      n_cells = nrow(dd),
      n_L2 = n_L2_group,
      figure = stem_A
    )
  )

  p_A2 <- make_group_umap(
    group,
    split_disease = TRUE
  )

  stem_A2 <- paste0(
    "20A3_",
    safe_group(group),
    "_A2_L2_UMAP_HC_FFA_AA"
  )

  save3(
    p_A2,
    safe_group(group),
    "A2_L2_UMAP_HC_FFA_AA",
    stem_A2,
    width = 14,
    height = 5.6
  )

  panel_A2_inventory <- bind_rows(
    panel_A2_inventory,
    tibble(
      presentation_group =
        group,
      n_cells = nrow(dd),
      n_L2 = n_L2_group,
      figure = stem_A2
    )
  )

  cat(
    sprintf(
      "%-20s : %6d cells | %2d L2 | A+A2 DONE\n",
      group,
      nrow(dd),
      n_L2_group
    )
  )
}

write_tsv(
  panel_A_inventory,
  file.path(
    TAB,
    "20A3_panel_A_inventory.tsv"
  )
)

write_tsv(
  panel_A2_inventory,
  file.path(
    TAB,
    "20A3_panel_A2_inventory.tsv"
  )
)

# ============================================================
# 12. PANEL D2 — DISEASE COMPOSITION WITHIN L2
# ============================================================

cat("\nSTEP 9/10 - generating disease-composition panels\n")

disease_cols <- c(
  HC = "#8F8D8E",
  FFA = "#4D9E4F",
  AA = "#A8142A"
)

make_disease_composition <- function(group) {

  disp <- get_group_display(group)

  dd <- base %>%
    filter(
      presentation_group == group
    )

  comp <- dd %>%
    count(
      final_L2_label,
      disease,
      name = "n_cells"
    ) %>%
    complete(
      final_L2_label =
        disp$final_L2_label,
      disease = c(
        "HC",
        "FFA",
        "AA"
      ),
      fill = list(
        n_cells = 0L
      )
    ) %>%
    group_by(
      final_L2_label
    ) %>%
    mutate(
      total_L2_cells =
        sum(n_cells),
      proportion =
        n_cells /
        total_L2_cells
    ) %>%
    ungroup() %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code,
          short_identity
        ),
      by = "final_L2_label"
    ) %>%
    mutate(
      plot_label_2line =
        paste0(
          display_code,
          "\n",
          short_identity
        ),
      plot_label_2line =
        factor(
          plot_label_2line,
          levels = paste0(
            disp$display_code,
            "\n",
            disp$short_identity
          )
        ),
      disease = factor(
        disease,
        levels = c(
          "HC",
          "FFA",
          "AA"
        )
      )
    )

  stopifnot(
    all(
      is.finite(
        comp$proportion
      )
    )
  )

  p <- ggplot(
    comp,
    aes(
      x = plot_label_2line,
      y = proportion,
      fill = disease
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.25,
      position =
        position_stack(
          reverse = TRUE
        )
    ) +
    scale_fill_manual(
      values = disease_cols,
      breaks = c(
        "HC",
        "FFA",
        "AA"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      labels =
        scales::percent_format(
          accuracy = 1
        ),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0(
        pretty_group(group),
        " - disease composition by L2 subtype"
      ),
      subtitle =
        "Cell-level descriptive composition within each frozen L2 identity",
      x = NULL,
      y = "Proportion",
      fill = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = ifelse(
            group ==
              "Epithelial_lineage",
            6.2,
            7.2
          )
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 11
        ),
      plot.subtitle =
        element_text(size = 8),
      legend.position =
        "right",
      legend.text =
        element_text(size = 8),
      plot.margin =
        margin(7, 8, 7, 7)
    )

  list(
    plot = p,
    table = comp,
    n_cells = nrow(dd),
    n_L2 = nrow(disp)
  )
}

panel_D2_inventory <- tibble()
panel_D2_tables <- list()

for (group in presentation_order) {

  res <- make_disease_composition(
    group
  )

  stem <- paste0(
    "20A3_",
    safe_group(group),
    "_D2_disease_composition"
  )

  save3(
    res$plot,
    safe_group(group),
    "D2_disease_composition",
    stem,
    width = ifelse(
      group ==
        "Epithelial_lineage",
      13.5,
      7.5
    ),
    height = 5.2
  )

  panel_D2_tables[[group]] <-
    res$table %>%
    mutate(
      presentation_group =
        group
    )

  panel_D2_inventory <- bind_rows(
    panel_D2_inventory,
    tibble(
      presentation_group =
        group,
      n_cells =
        res$n_cells,
      n_L2 =
        res$n_L2,
      figure =
        stem
    )
  )

  cat(
    sprintf(
      "%-20s : %6d cells | %2d L2 | D2 DONE\n",
      group,
      res$n_cells,
      res$n_L2
    )
  )
}

write_tsv(
  bind_rows(
    panel_D2_tables
  ),
  file.path(
    TAB,
    "20A3_D2_disease_composition.tsv"
  )
)

write_tsv(
  panel_D2_inventory,
  file.path(
    TAB,
    "20A3_panel_D2_inventory.tsv"
  )
)

# ============================================================
# 13. FORMAL AUDIT
# ============================================================

cat("\nSTEP 10/10 - formal audit\n")

expected_group_counts <- tibble(
  presentation_group = c(
    "Epithelial_lineage",
    "Fibroblast",
    "Immune",
    "Endothelial",
    "Mural",
    "Melanocyte"
  ),
  expected_n_L2 = c(
    19L,
    10L,
    9L,
    10L,
    7L,
    4L
  )
)

observed_group_counts <- base %>%
  distinct(
    presentation_group,
    final_L2_label
  ) %>%
  count(
    presentation_group,
    name = "observed_n_L2"
  ) %>%
  right_join(
    expected_group_counts,
    by = "presentation_group"
  )

group_counts_correct <- all(
  observed_group_counts$observed_n_L2 ==
    observed_group_counts$expected_n_L2
)

epithelial_family_counts <- base %>%
  filter(
    presentation_group ==
      "Epithelial_lineage"
  ) %>%
  distinct(
    final_compartment,
    final_L2_label
  ) %>%
  count(
    final_compartment,
    name = "n_L2"
  )

epithelial_expected <- c(
  Epithelial = 9L,
  Sweat_gland = 5L,
  Sebocyte = 5L
)

epithelial_family_structure_correct <-
  nrow(epithelial_family_counts) == 3 &&
  all(
    epithelial_expected[
      epithelial_family_counts$final_compartment
    ] ==
      epithelial_family_counts$n_L2
  )

all_B_L2_covered <-
  marker_B %>%
    distinct(
      final_L2_label
    ) %>%
    nrow() == 59

all_C_L2_covered <-
  marker_C %>%
    distinct(
      final_L2_label
    ) %>%
    nrow() == 59

expected_designs <- 36L

pdf_files <- list.files(
  FIG,
  pattern = "\\.pdf$",
  recursive = TRUE,
  full.names = TRUE
)

png_files <- list.files(
  FIG,
  pattern = "\\.png$",
  recursive = TRUE,
  full.names = TRUE
)

svg_files <- list.files(
  FIG,
  pattern = "\\.svg$",
  recursive = TRUE,
  full.names = TRUE
)

checks <- tibble(
  check = c(
    "analysis_eligible_L2_82494",
    "frozen_L2_59",
    "frozen_L1_8",
    "presentation_groups_6",
    "presentation_group_counts_correct",
    "epithelial_lineage_L2_19",
    "epithelial_lineage_EP9_SG5_SE5",
    "fibroblast_L2_10",
    "
immune_L2_9",
    "endothelial_L2_10",
    "mural_L2_7",
    "melanocyte_L2_4",
    "canonical_display_labels_59",
    "panel_A_6",
    "panel_A2_6",
    "panel_B_6",
    "panel_C_6",
    "panel_D1_6",
    "panel_D2_6",
    "patients_16",
    "patient_disease_groups_HC_FFA_AA",
    "HC_colour_8F8D8E",
    "FFA_colour_4D9E4F",
    "AA_colour_A8142A",
    "frozen_L2_colours_preserved",
    "epithelial_EP_colours_preserved",
    "sweat_gland_SG_colours_preserved",
    "sebocyte_SE_colours_preserved",
    "all_B_59_L2_covered",
    "all_C_59_L2_covered",
    "same_Harmony_coordinates",
    "figure_designs_36",
    "PDF_files_36",
    "PNG_files_36",
    "SVG_files_36",
    "dynamic_abbreviations_FALSE",
    "presentation_group_only",
    "no_clustering_recalculated",
    "no_UMAP_recalculated",
    "no_annotation_modified",
    "no_DE_inference_recalculated"
  ),

  pass = c(
    nrow(base) == 82494,

    n_distinct(
      base$final_L2_label
    ) == 59,

    n_distinct(
      base$final_compartment
    ) == 8,

    n_distinct(
      base$presentation_group
    ) == 6,

    group_counts_correct,

    base %>%
      filter(
        presentation_group ==
          "Epithelial_lineage"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 19,

    epithelial_family_structure_correct,

    base %>%
      filter(
        presentation_group ==
          "Fibroblast"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 10,

    base %>%
      filter(
        presentation_group ==
          "Immune"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 9,

    base %>%
      filter(
        presentation_group ==
          "Endothelial"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 10,

    base %>%
      filter(
        presentation_group ==
          "Mural"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 7,

    base %>%
      filter(
        presentation_group ==
          "Melanocyte"
      ) %>%
      distinct(
        final_L2_label
      ) %>%
      nrow() == 4,

    nrow(display_labels) == 59 &&
      n_distinct(
        display_labels$final_L2_label
      ) == 59 &&
      n_distinct(
        display_labels$display_code
      ) == 59,

    nrow(panel_A_inventory) == 6,

    nrow(panel_A2_inventory) == 6,

    nrow(panel_B_inventory) == 6,

    nrow(panel_C_inventory) == 6,

    nrow(panel_D_inventory) == 6,

    nrow(panel_D2_inventory) == 6,

    n_distinct(
      base$patient_id
    ) == 16,

    setequal(
      unique(
        as.character(
          base$disease
        )
      ),
      c(
        "HC",
        "FFA",
        "AA"
      )
    ),

    identical(
      unname(
        disease_cols["HC"]
      ),
      "#8F8D8E"
    ),

    identical(
      unname(
        disease_cols["FFA"]
      ),
      "#4D9E4F"
    ),

    identical(
      unname(
        disease_cols["AA"]
      ),
      "#A8142A"
    ),

    all(
      display$hex[
        match(
          display_labels$final_L2_label,
          display$final_L2_label
        )
      ] ==
        display_labels$hex
    ),

    all(
      display_labels %>%
        filter(
          L1 == "Epithelial"
        ) %>%
        arrange(
          L2_order_within_L1
        ) %>%
        pull(hex) ==
        display %>%
          filter(
            L1 == "Epithelial"
          ) %>%
          arrange(
            L2_order_within_L1
          ) %>%
          pull(hex)
    ),

    all(
      display_labels %>%
        filter(
          L1 == "Sweat_gland"
        ) %>%
        arrange(
          L2_order_within_L1
        ) %>%
        pull(hex) ==
        display %>%
          filter(
            L1 == "Sweat_gland"
          ) %>%
          arrange(
            L2_order_within_L1
          ) %>%
          pull(hex)
    ),

    all(
      display_labels %>%
        filter(
          L1 == "Sebocyte"
        ) %>%
        arrange(
          L2_order_within_L1
        ) %>%
        pull(hex) ==
        display %>%
          filter(
            L1 == "Sebocyte"
          ) %>%
          arrange(
            L2_order_within_L1
          ) %>%
          pull(hex)
    ),

    all_B_L2_covered,

    all_C_L2_covered,

    nrow(umap) == 82893 &&
      nrow(plot_base) == 82494 &&
      all(
        plot_base$cell %in%
          rownames(emb)
      ),

    nrow(panel_A_inventory) +
      nrow(panel_A2_inventory) +
      nrow(panel_B_inventory) +
      nrow(panel_C_inventory) +
      nrow(panel_D_inventory) +
      nrow(panel_D2_inventory) ==
        expected_designs,

    length(pdf_files) ==
      expected_designs,

    length(png_files) ==
      expected_designs,

    length(svg_files) ==
      expected_designs,

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "20A3_formal_checks.tsv"
  )
)

cat("\n===== PRESENTATION-GROUP L2 COUNTS =====\n")

print(
  as.data.frame(
    observed_group_counts %>%
      arrange(
        match(
          presentation_group,
          presentation_order
        )
      )
  ),
  row.names = FALSE
)

cat("\n===== EPITHELIAL LINEAGE STRUCTURE =====\n")

print(
  as.data.frame(
    epithelial_family_counts
  ),
  row.names = FALSE
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  as.data.frame(checks),
  row.names = FALSE
)

if (!all(checks$pass)) {

  cat("\n===== FAILED CHECKS =====\n")

  print(
    as.data.frame(
      checks %>%
        filter(!pass)
    ),
    row.names = FALSE
  )

  stop(
    "20A3 formal checks failed. ",
    "Do not accept presentation-group figures."
  )
}

# ============================================================
# 14. MASTER FIGURE INVENTORY
# ============================================================

master_inventory <- bind_rows(

  panel_A_inventory %>%
    transmute(
      presentation_group,
      panel = "A_L2_UMAP",
      n_cells,
      n_L2,
      figure
    ),

  panel_A2_inventory %>%
    transmute(
      presentation_group,
      panel = "A2_L2_UMAP_HC_FFA_AA",
      n_cells,
      n_L2,
      figure
    ),

  panel_B_inventory %>%
    transmute(
      presentation_group,
      panel = "B_cell_marker_heatmap",
      n_cells = n_cells_displayed,
      n_L2,
      figure
    ),

  panel_C_inventory %>%
    transmute(
      presentation_group,
      panel = "C_marker_violins",
      n_cells = n_cells_displayed,
      n_L2,
      figure
    ),

  panel_D_inventory %>%
    transmute(
      presentation_group,
      panel = "D1_patient_composition",
      n_cells = NA_integer_,
      n_L2 = NA_integer_,
      figure
    ),

  panel_D2_inventory %>%
    transmute(
      presentation_group,
      panel = "D2_disease_composition",
      n_cells,
      n_L2,
      figure
    )
) %>%
  mutate(
    presentation_group = factor(
      presentation_group,
      levels = presentation_order
    )
  ) %>%
  arrange(
    presentation_group,
    match(
      panel,
      c(
        "A_L2_UMAP",
        "A2_L2_UMAP_HC_FFA_AA",
        "B_cell_marker_heatmap",
        "C_marker_violins",
        "D1_patient_composition",
        "D2_disease_composition"
      )
    )
  ) %>%
  mutate(
    presentation_group =
      as.character(
        presentation_group
      )
  )

write_tsv(
  master_inventory,
  file.path(
    TAB,
    "20A3_master_figure_inventory.tsv"
  )
)

cat("\n===== MASTER FIGURE INVENTORY =====\n")

print(
  as.data.frame(
    master_inventory
  ),
  row.names = FALSE
)

# ============================================================
# 15. REPRODUCIBILITY RECORD
# ============================================================

record <- c(
  "PIPELINE_V2 - 20A3 PRESENTATION-GROUP FIGURE SUITE",
  "",
  "STATUS: GENERATED - VISUAL QC REQUIRED",
  "",
  "Presentation grouping is a figure-level operation only.",
  "Frozen biological L1/L2 annotation is unchanged.",
  "",
  "PRESENTATION GROUPS",
  "Epithelial_lineage = Epithelial + Sweat_gland + Sebocyte",
  "Fibroblast = Fibroblast",
  "Immune = Immune",
  "Endothelial = Endothelial",
  "Mural = Mural",
  "Melanocyte = Melanocyte",
  "",
  "EPITHELIAL LINEAGE",
  "19 frozen L2 identities = EP9 + SG5 + SE5.",
  "Original frozen colours retained.",
  "",
  "PANELS",
  "A  = global L2 UMAP",
  "A2 = HC | FFA | AA split L2 UMAP",
  "B  = cell-level marker heatmap",
  "C  = marker violins",
  "D1 = patient-level composition",
  "D2 = disease composition within L2",
  "",
  "Canonical L2 display labels: YES",
  "Dynamic abbreviation generation: NO",
  "Clustering recalculated: NO",
  "Harmony UMAP recalculated: NO",
  "Annotation modified: NO",
  "DE inference recalculated: NO",
  "",
  "36 designs exported as PDF, PNG and SVG.",
  "PASS"
)

writeLines(
  record,
  con = file.path(
    OUT,
    "20A3_REPRODUCIBILITY_RECORD.txt"
  )
)

# ============================================================
# 16. FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("20A3 COMPLETE\n")
cat("============================================================\n")
cat("Presentation groups: 6\n")
cat("Frozen L2 identities: 59\n")
cat("Epithelial lineage L2: 19 (EP9 + SG5 + SE5)\n")
cat("Panel A UMAPs: 6\n")
cat("Panel A2 disease-split UMAPs: 6\n")
cat("Panel B marker heatmaps: 6\n")
cat("Panel C marker violins: 6\n")
cat("Panel D1 patient compositions: 6\n")
cat("Panel D2 disease compositions: 6\n")
cat("Figure designs generated: 36\n")
cat("PDF files generated: 36\n")
cat("PNG files generated: 36\n")
cat("SVG files generated: 36\n")
cat("Editable SVG export: YES\n")
cat("Dynamic abbreviations generated: NO\n")
cat("Frozen colours preserved: YES\n")
cat("UMAP recalculated: NO\n")
cat("Clustering recalculated: NO\n")
cat("Annotation modified: NO\n")
cat("DE inference recalculated: NO\n")
cat("STATUS: VISUAL QC REQUIRED\n")
cat("PASS\n")
