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
cat("PIPELINE_V2 - 20A2 PAPER-STYLE CELLTYPE PANELS\n")
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
  "results/publication_panels/20A2_paper_style_celltype"
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
# 7. PAPER-STYLE CELL HEATMAP
# ============================================================

cat("STEP 5/9 - generating paper-style cell-level heatmaps\n")

make_cell_heatmap <- function(
  parent,
  max_cells_per_L2 = 250
) {

  disp <- display %>%
    filter(L1 == parent) %>%
    arrange(L2_order_within_L1)

  markers <- marker_B %>%
    filter(L1 == parent) %>%
    arrange(
      L2_order_within_L1,
      desc(max_marker_score)
    )

  genes <- unique(markers$gene)

  dd <- base %>%
    filter(
      final_compartment == parent
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

  cell_order <- unlist(
    selected_cells,
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

  gene_key <- markers %>%
    distinct(
      gene,
      final_L2_label,
      L2_order_within_L1
    ) %>%
    arrange(
      L2_order_within_L1,
      desc(
        match(
          gene,
          genes
        )
      )
    )

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
      gene = factor(
        gene,
        levels = rev(gene_order)
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
        gsub("_", " ", parent),
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
        element_text(size = 7),
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
        margin(0, 6, 4, 6)
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
      size = 2.45,
      angle = 28,
      hjust = 0,
      vjust = 0.5
    ) +
    coord_cartesian(
      ylim = c(
        0.55,
        2.2
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
        0.55,
        8
      )
    )

  list(
    plot = combined,
    n_cells = length(cell_order),
    n_genes = length(genes),
    boundaries = boundaries,
    cells = tibble(
      parent = parent,
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

for (parent in unique(display$L1)) {

  if (parent == "Neural_Schwann") {
    next
  }

  cat(
    sprintf(
      "  %-15s",
      parent
    )
  )

  res <- make_cell_heatmap(
    parent = parent,
    max_cells_per_L2 = 250
  )

  stem <- paste0(
    "20A2_",
    parent,
    "_B_cell_marker_heatmap"
  )

  disp_parent <- display %>%
    filter(
      L1 == parent
    )

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
    parent,
    "B_cell_marker_heatmap",
    stem,
    width,
    height
  )

  panel_B_inventory <- bind_rows(
    panel_B_inventory,
    tibble(
      L1 = parent,
      n_L2 =
        nrow(disp_parent),
      n_cells_displayed =
        res$n_cells,
      n_genes =
        res$n_genes,
      figure = stem
    )
  )

  panel_B_cells[[parent]] <-
    res$cells

  panel_B_genes[[parent]] <-
    res$genes

  panel_B_boundaries[[parent]] <-
    res$boundaries %>%
    mutate(
      L1 = parent
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
  parent,
  max_cells_per_L2 = 800
) {

  disp <- display %>%
    filter(L1 == parent) %>%
    arrange(L2_order_within_L1)

  markers <- marker_C %>%
    filter(L1 == parent) %>%
    arrange(L2_order_within_L1)

  stopifnot(nrow(markers) > 0)

  genes <- unique(markers$gene)

  dd <- base %>%
    filter(
      final_compartment == parent
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
      linewidth = 0.18,
      colour = "grey25"
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
        gsub("_", " ", parent),
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
          size = 7
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
          size = 7,
          face = "italic"
        ),
      panel.spacing.y =
        grid::unit(
          0.8,
          "mm"
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

for (parent in unique(display$L1)) {

  if (parent == "Neural_Schwann") {
    next
  }

  res <- make_compact_violins(
    parent
  )

  stem <- paste0(
    "20A2_",
    parent,
    "_C_marker_violins"
  )

  n_L2_parent <- display %>%
    filter(
      L1 == parent
    ) %>%
    nrow()

  width <- max(
    6.5,
    0.62 * n_L2_parent + 3.0
  )

  height <- max(
    5.5,
    0.55 * res$n_genes + 2.2
  )

  save3(
    res$plot,
    parent,
    "C_marker_violins",
    stem,
    width,
    height
  )

  panel_C_inventory <- bind_rows(
    panel_C_inventory,
    tibble(
      L1 = parent,
      n_L2 = n_L2_parent,
      n_cells_displayed =
        res$n_cells,
      n_genes =
        res$n_genes,
      figure = stem
    )
  )

  panel_C_genes[[parent]] <-
    res$genes %>%
    mutate(
      L1 = parent
    )

  cat(
    sprintf(
      "%-15s : %5d cells | %2d genes | violin panel DONE\n",
      parent,
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
  parent
) {

  disp <- display %>%
    filter(
      L1 == parent
    ) %>%
    arrange(
      L2_order_within_L1
    )

  dd <- base %>%
    filter(
      final_compartment == parent
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
        gsub("_", " ", parent),
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

for (parent in unique(display$L1)) {

  if (parent == "Neural_Schwann") {
    next
  }

  res <- make_patient_composition(
    parent
  )

  stem <- paste0(
    "20A2_",
    parent,
    "_D_patient_composition"
  )

  save3(
    res$plot,
    parent,
    "D_patient_composition",
    stem,
    width = 10.5,
    height = 6.2
  )

  panel_D_tables[[parent]] <-
    res$table %>%
    mutate(
      L1 = parent
    )

  panel_D_inventory <- bind_rows(
    panel_D_inventory,
    tibble(
      L1 = parent,
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
      parent,
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
# 11. PANEL A — REGISTER EXISTING 19C3 UMAPS
# ============================================================

cat("\nSTEP 8/9 - registering existing panel A UMAPs\n")

UMAP_ROOT <- file.path(
  ROOT,
  "results/figure_specification/19C3_publication_suite"
)

panel_A_registry <- display %>%
  distinct(
    L1
  ) %>%
  filter(
    L1 != "Neural_Schwann"
  ) %>%
  mutate(
    panel = "A_L2_UMAP",
    source =
      "19C3 publication suite",
    status =
      "REUSE_EXISTING_FROZEN_UMAP"
  )

write_tsv(
  panel_A_registry,
  file.path(
    TAB,
    "20A2_panel_A_UMAP_registry.tsv"
  )
)

# ============================================================
# 12. FORMAL AUDIT
# ============================================================

cat("\nSTEP 9/9 - formal audit\n")

n_parents <- display %>%
  filter(
    L1 != "Neural_Schwann"
  ) %>%
  distinct(
    L1
  ) %>%
  nrow()

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

# 8 B heatmaps + 8 C violin panels + 8 D compositions.
expected_new_designs <- 24L

all_B_L2_covered <- marker_B %>%
  distinct(
    L1,
    final_L2_label
  ) %>%
  nrow() == 59

all_C_L2_covered <- marker_C %>%
  distinct(
    L1,
    final_L2_label
  ) %>%
  nrow() == 59

patient_order_diseases <- patient_order %>%
  distinct(
    disease
  ) %>%
  pull(
    disease
  )

checks <- tibble(
  check = c(
    "L2_universe_82494",
    "frozen_L2_59",
    "parent_compartments_8",
    "panel_A_registered_8",
    "panel_B_heatmaps_8",
    "panel_C_violins_8",
    "panel_D_compositions_8",
    "panel_B_all_59_L2_covered",
    "panel_C_all_59_L2_covered",
    "patients_16",
    "patient_disease_groups_HC_FFA_AA",
    "HC_colour_8F8D8E",
    "FFA_colour_4D9E4F",
    "AA_colour_A8142A",
    "new_figure_designs_24",
    "PDF_files_24",
    "PNG_files_24",
    "SVG_files_24",
    "no_clustering_recalculated",
    "no_annotation_modified",
    "no_DE_inference_recalculated"
  ),
  pass = c(
    nrow(base) == 82494,
    n_distinct(
      base
$final_L2_label
    ) == 59,

    n_distinct(
      base$final_compartment
    ) == 8,

    nrow(panel_A_registry) == 8,

    nrow(panel_B_inventory) == 8,

    nrow(panel_C_inventory) == 8,

    nrow(panel_D_inventory) == 8,

    all_B_L2_covered,

    all_C_L2_covered,

    n_distinct(
      base$patient_id
    ) == 16,

    setequal(
      patient_order_diseases,
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

    nrow(panel_B_inventory) +
      nrow(panel_C_inventory) +
      nrow(panel_D_inventory) ==
      expected_new_designs,

    length(pdf_files) ==
      expected_new_designs,

    length(png_files) ==
      expected_new_designs,

    length(svg_files) ==
      expected_new_designs,

    TRUE,
    TRUE,
    TRUE
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "20A2_formal_checks.tsv"
  )
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
    "20A2 formal checks failed. ",
    "Do not accept publication panels until resolved."
  )
}

# ============================================================
# 13. MASTER FIGURE INVENTORY
# ============================================================

master_inventory <- bind_rows(

  panel_A_registry %>%
    transmute(
      L1,
      panel = "A_L2_UMAP",
      figure =
        "REUSE_FROM_19C3",
      status =
        "EXISTING"
    ),

  panel_B_inventory %>%
    transmute(
      L1,
      panel =
        "B_cell_marker_heatmap",
      figure,
      status =
        "GENERATED_20A2"
    ),

  panel_C_inventory %>%
    transmute(
      L1,
      panel =
        "C_marker_violins",
      figure,
      status =
        "GENERATED_20A2"
    ),

  panel_D_inventory %>%
    transmute(
      L1,
      panel =
        "D_patient_composition",
      figure,
      status =
        "GENERATED_20A2"
    )
)

write_tsv(
  master_inventory,
  file.path(
    TAB,
    "20A2_master_panel_inventory.tsv"
  )
)

cat("\n===== MASTER PANEL INVENTORY =====\n")

print(
  as.data.frame(
    master_inventory
  ),
  row.names = FALSE
)

# ============================================================
# 14. REPRODUCIBILITY RECORD
# ============================================================

record <- c(
  "PIPELINE_V2 - 20A2 PAPER-STYLE CELLTYPE PANELS",
  "",
  "STATUS: GENERATED - VISUAL QC REQUIRED",
  "",
  "PURPOSE",
  "Reproduce the cell-type characterization structure used",
  "in the reference single-cell publication figure.",
  "",
  "PANEL A",
  "Frozen L2 UMAPs are reused from 19C3.",
  "No UMAP coordinates were recalculated.",
  "",
  "PANEL B",
  "Cell-level marker heatmaps.",
  "Individual cells are columns and marker genes are rows.",
  "Cells are grouped by frozen L2 identity.",
  "Within each L2, cells are ordered by the mean expression",
  "of that identity's selected marker genes.",
  "For rendering, at most 250 cells per L2 are selected",
  "deterministically at evenly spaced ranks.",
  "This thinning is graphical only.",
  "Gene-wise expression is scaled across displayed cells.",
  "A frozen L2 colour strip is shown above the heatmap.",
  "",
  "PANEL C",
  "Compact violin expression plots.",
  "One representative cluster-identifying marker is shown",
  "for each frozen L2 identity.",
  "At most 800 cells per L2 are retained deterministically",
  "for rendering only.",
  "No statistical inference is performed in this panel.",
  "",
  "PANEL D",
  "Patient-level 100% stacked L2 composition.",
  "One bar represents one biological patient.",
  "Patients follow the frozen HC -> FFA -> AA order.",
  "",
  "COLOUR SPECIFICATION",
  "Frozen 19B/19C3 L2 colours reused.",
  "HC  = #8F8D8E",
  "FFA = #4D9E4F",
  "AA  = #A8142A",
  "",
  "OUTPUTS",
  "8 panel-B cell-level heatmaps.",
  "8 panel-C compact violin panels.",
  "8 panel-D patient composition panels.",
  "24 new figure designs in PDF, PNG and editable SVG.",
  "8 panel-A UMAPs registered from 19C3.",
  "",
  "Clustering recalculated: NO",
  "UMAP recalculated: NO",
  "Annotation modified: NO",
  "DE inference recalculated: NO",
  "PASS"
)

writeLines(
  record,
  con = file.path(
    OUT,
    "20A2_REPRODUCIBILITY_RECORD.txt"
  )
)

# ============================================================
# 15. FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("20A2 COMPLETE\n")
cat("============================================================\n")
cat("Panel A UMAPs registered: 8\n")
cat("Panel B cell-level heatmaps generated: 8\n")
cat("Panel C marker violin panels generated: 8\n")
cat("Panel D patient compositions generated: 8\n")
cat("New figure designs generated: 24\n")
cat("PDF files generated: 24\n")
cat("PNG files generated: 24\n")
cat("SVG files generated: 24\n")
cat("Editable SVG export: YES\n")
cat("Frozen L2 identities represented: 59\n")
cat("Patient order: HC -> FFA -> AA\n")
cat("UMAP recalculated: NO\n")
cat("Clustering recalculated: NO\n")
cat("Annotation modified: NO\n")
cat("DE inference recalculated: NO\n")
cat("STATUS: VISUAL QC REQUIRED\n")
cat("PASS\n")
