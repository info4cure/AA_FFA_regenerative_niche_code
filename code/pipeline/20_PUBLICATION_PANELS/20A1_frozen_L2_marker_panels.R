suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
})

cat("============================================================\n")
cat("PIPELINE_V2 - 20A1 FROZEN L2 MARKER PANELS\n")
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

OUT <- file.path(
  ROOT,
  "results/publication_panels/20A1_frozen_L2_markers"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. SOURCE SPECIFICATION
# ============================================================

cat("\nSTEP 1/10 - defining frozen source specification\n")

spec <- tribble(
  ~L1, ~module, ~marker_file, ~resolution_filter,

  "Epithelial",
  "07_EPITHELIAL_L2",
  "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07D_markers/07D_markers_res0.5_all.csv",
  NA_real_,

  "Immune",
  "08_IMMUNE_L2",
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08H_final_L2_adjudication/tables/08H_markers_all.csv",
  NA_real_,

  "Fibroblast",
  "09_FIBROBLAST_L2",
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/09H5A_postreturn_marker_evidence/tables/09H5A_all_postreturn_markers.csv",
  NA_real_,

  "Endothelial",
  "10_VASCULAR_L2",
  "PIPELINE_V2/results/annotation/10_VASCULAR_L2/10B5_final_identity_adjudication/tables/10B5_r06_all_positive_markers.csv",
  NA_real_,

  "Mural",
  "11_MURAL_L2",
  "PIPELINE_V2/results/annotation/11_MURAL_L2/11B2_identity_evidence/tables/11B2_all_candidate_markers.csv",
  0.6,

  "Sweat_gland",
  "12_SWEAT_GLAND_L2",
  "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12C_adjudication/tables/12C_r04_all_markers.csv",
  NA_real_,

  "Melanocyte",
  "13_MELANOCYTE_L2",
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/13B_recluster/tables/13B_all_candidate_markers.csv",
  0.4,

  "Sebocyte",
  "14_SEBOCYTE_L2",
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14B_recluster/tables/14B_all_candidate_markers.csv",
  0.4
)

mapping_files <- c(
  Epithelial =
    "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07I_frozen/tables/07I_epithelial_L2_frozen_map.csv",

  Immune =
    "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08I_L2_freeze/tables/08I_cluster_to_L2_mapping.csv",

  Fibroblast =
    "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/09H5B_final_L2_freeze/tables/09H5B_cluster_L2_summary.csv",

  Endothelial =
    "PIPELINE_V2/results/annotation/10_VASCULAR_L2/10B5_final_identity_adjudication/tables/10B5_final_mapping_compact.csv",

  Mural =
    "PIPELINE_V2/results/annotation/11_MURAL_L2/11B5_freeze/tables/11B5_frozen_mural_mapping.csv",

  Sweat_gland =
    "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/12D_freeze/tables/12D_frozen_sweat_gland_mapping.csv",

  Melanocyte =
    "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/13D_freeze/tables/13D_frozen_melanocyte_mapping.csv",

  Sebocyte =
    "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14D_freeze/tables/14D_frozen_sebocyte_mapping.csv"
)

stopifnot(file.exists(OBJ))
stopifnot(file.exists(LEDGER))
stopifnot(file.exists(DISPLAY))
stopifnot(all(file.exists(spec$marker_file)))
stopifnot(all(file.exists(mapping_files)))

# ============================================================
# 3. LOAD CANONICAL DATA
# ============================================================

cat("STEP 2/10 - loading canonical expression and frozen ledger\n")

obj <- readRDS(OBJ)
ledger <- readRDS(LEDGER)
display <- read_tsv(
  DISPLAY,
  show_col_types = FALSE
)

stopifnot(inherits(obj, "Seurat"))

if (!"cell" %in% colnames(ledger)) {
  candidate <- intersect(
    c("cell", "cell_id", "barcode", "Cell"),
    colnames(ledger)
  )

  if (length(candidate) != 1) {
    stop("Cannot uniquely resolve cell identifier in final ledger.")
  }

  ledger <- ledger %>%
    rename(cell = all_of(candidate))
}

required_ledger <- c(
  "cell",
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label"
)

stopifnot(
  all(required_ledger %in% colnames(ledger))
)

ledger_L2 <- ledger %>%
  filter(analysis_eligible_L2) %>%
  select(
    cell,
    final_compartment,
    final_L2_label
  )

stopifnot(nrow(ledger_L2) == 82494)
stopifnot(!anyDuplicated(ledger_L2$cell))
stopifnot(n_distinct(ledger_L2$final_L2_label) == 59)

stopifnot(
  all(ledger_L2$cell %in% colnames(obj))
)

# ============================================================
# 4. EXPRESSION MATRIX
# ============================================================

cat("STEP 3/10 - recovering normalized expression matrix\n")

DefaultAssay(obj) <- "RNA"

get_data_matrix <- function(x) {

  out <- tryCatch(
    GetAssayData(
      x,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) NULL
  )

  if (is.null(out) || nrow(out) == 0) {
    out <- tryCatch(
      GetAssayData(
        x,
        assay = "RNA",
        slot = "data"
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out) || nrow(out) == 0) {
    stop(
      "No normalized RNA data layer found. ",
      "20A1 will not silently substitute counts."
    )
  }

  out
}

expr <- get_data_matrix(obj)

stopifnot(
  all(ledger_L2$cell %in% colnames(expr))
)

cat(
  "Expression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "cells\n"
)

# ============================================================
# 5. DISPLAY SPECIFICATION
# ============================================================

cat("STEP 4/10 - validating frozen display specification\n")

required_display <- c(
  "module",
  "L1",
  "final_L2_label",
  "hex",
  "L2_order_within_L1",
  "display_code",
  "legend_label"
)

stopifnot(
  all(required_display %in% colnames(display))
)

stopifnot(nrow(display) == 59)
stopifnot(!anyDuplicated(display$final_L2_label))
stopifnot(!anyDuplicated(display$display_code))

# ============================================================
# 6. MAPPING NORMALIZATION
# ============================================================

cat("STEP 5/10 - recovering cluster-to-frozen-L2 mappings\n")

find_col <- function(df, candidates, what, source_name) {

  hit <- candidates[
    candidates %in% colnames(df)
  ]

  if (length(hit) == 0) {
    stop(
      "Cannot resolve ", what,
      " in ", source_name,
      ". Columns: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  hit[[1]]
}

read_mapping <- function(L1, path) {

  x <- read_csv(
    path,
    show_col_types = FALSE
  )

  cluster_col <- find_col(
    x,
    c(
      "cluster",
      "immune_L2_cluster",
      "fibroblast_cluster_frozen",
      "vascular_cluster_frozen",
      "mural_cluster_frozen",
      "sweat_cluster_frozen",
      "melanocyte_cluster_frozen",
      "sebocyte_cluster_frozen"
    ),
    "cluster column",
    path
  )

  label_col <- find_col(
    x,
    c(
      "L2_identity",
      "immune_L2",
      "fibroblast_L2_label",
      "vascular_L2_label",
      "mural_L2_label",
      "sweat_L2_label",
      "melanocyte_L2_label",
      "sebocyte_L2_label",
      "final_L2_label"
    ),
    "L2 label column",
    path
  )

  x %>%
    transmute(
      L1 = L1,
      cluster = as.character(
        .data[[cluster_col]]
      ),
      final_L2_label = as.character(
        .data[[label_col]]
      )
    ) %>%
    filter(
      !is.na(cluster),
      !is.na(final_L2_label),
      final_L2_label != ""
    ) %>%
    distinct()
}

mapping_all <- bind_rows(
  lapply(
    names(mapping_files),
    function(nm) {
      read_mapping(
        nm,
        mapping_files[[nm]]
      )
    }
  )
)

mapping_all <- mapping_all %>%
  left_join(
    display %>%
      select(
        L1,
        final_L2_label,
        display_code,
        legend_label,
        hex,
        L2_order_within_L1
      ),
    by = c(
      "L1",
      "final_L2_label"
    )
  )

if (anyNA(mapping_all$display_code)) {
  bad <- mapping_all %>%
    filter(is.na(display_code))

  print(bad)

  stop(
    "Frozen cluster mapping does not match 19C3 display identities."
  )
}

write_tsv(
  mapping_all,
  file.path(
    TAB,
    "20A1_cluster_to_frozen_L2_mapping.tsv"
  )
)

cat("\n===== CLUSTER -> FROZEN L2 SUMMARY =====\n")

print(
  as.data.frame(
    mapping_all %>%
      count(
        L1,
        final_L2_label,
        display_code,
        name = "n_source_clusters"
      ) %>%
      arrange(
        match(L1, spec$L1),
        display_code
      )
  ),
  row.names = FALSE
)

# ============================================================
# 7. MARKER CANDIDATE RECOVERY
# ============================================================

cat("\nSTEP 6/10 - selecting candidate marker genes\n")

read_markers <- function(
  L1,
  marker_file,
  resolution_filter
) {

  m <- read_csv(
    marker_file,
    show_col_types = FALSE
  )

  required <- c(
    "cluster",
    "gene",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj"
  )

  stopifnot(
    all(required %in% colnames(m))
  )

  if (
    "resolution" %in% colnames(m) &&
    !is.na(resolution_filter)
  ) {
    m <- m %>%
      filter(
        abs(
          as.numeric(resolution) -
            resolution_filter
        ) < 1e-8
      )
  }

  parent_L1 <- L1

  map <- mapping_all %>%
    filter(.data$L1 == parent_L1) %>%
    select(
      cluster,
      final_L2_label,
      display_code
    )

  m %>%
    mutate(
      cluster = as.character(cluster),
      gene = as.character(gene)
    ) %>%
    inner_join(
      map,
      by = "cluster"
    ) %>%
    filter(
      is.finite(avg_log2FC),
      avg_log2FC > 0,
      pct.1 >= 0.10,
      p_val_adj <= 0.05
    ) %>%
    mutate(
      L1 = parent_L1,
      specificity_delta = pct.1 - pct.2,
      marker_score =
        avg_log2FC *
        pmax(
          specificity_delta,
          0.05
        )
    )
}

marker_candidates <- bind_rows(
  lapply(
    seq_len(nrow(spec)),
    function(i) {

      read_markers(
        spec$L1[i],
        spec$marker_file[i],
        spec$resolution_filter[i]
      ) %>%
        mutate(
          module = spec$module[i]
        )
    }
  )
)

stopifnot(nrow(marker_candidates) > 0)

write_tsv(
  marker_candidates,
  file.path(
    TAB,
    "20A1_marker_candidates_cluster_evidence.tsv"
  )
)

# ============================================================
# 8. SELECT TOP GENES PER FROZEN L2
# ============================================================

cat("STEP 7/10 - selecting reproducible marker panel\n")

# Multiple source clusters may converge onto one frozen L2.
# We use cluster-level marker evidence ONLY for candidate selection.
# Expression shown in final panels is recomputed on frozen L2 cells.

marker_ranked <- marker_candidates %>%
  group_by(
    L1,
    final_L2_label,
    display_code,
    gene
  ) %>%
  summarise(
    n_supporting_clusters =
      n_distinct(cluster),
    max_log2FC =
      max(avg_log2FC, na.rm = TRUE),
    max_pct =
      max(pct.1, na.rm = TRUE),
    max_specificity =
      max(
        specificity_delta,
        na.rm = TRUE
      ),
    max_marker_score =
      max(
        marker_score,
        na.rm = TRUE
      ),
    min_p_adj =
      min(
        p_val_adj,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    L1,
    final_L2_label,
    desc(n_supporting_clusters),
    desc(max_marker_score),
    desc(max_log2FC),
    gene
  )

# Only genes represented in the canonical normalized RNA matrix
# are eligible for publication marker panels.
marker_ranked <- marker_ranked %>%
  filter(gene %in% rownames(expr))

# Prefer genes supported by multiple source clusters when an L2
# arose through cluster consolidation, but do not require it.
marker_selected <- marker_ranked %>%
  group_by(
    L1,
    final_L2_label,
    display_code
  ) %>%
  arrange(
    desc(n_supporting_clusters),
    desc(max_marker_score),
    desc(max_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 5) %>%
  ungroup()

# Prevent the same gene from dominating many identities.
# First assignment follows frozen L2 display order and marker rank.
marker_selected <- marker_selected %>%
  left_join(
    display %>%
      select(
        L1,
        final_L2_label,
        L2_order_within_L1
      ),
    by = c(
      "L1",
      "final_L2_label"
    )
  ) %>%
  group_by(L1) %>%
  arrange(
    L2_order_within_L1,
    desc(n_supporting_clusters),
    desc(max_marker_score),
    .by_group = TRUE
  ) %>%
  mutate(
    selection_rank =
      row_number()
  ) %>%
  ungroup()

write_tsv(
  marker_ranked,
  file.path(
    TAB,
    "20A1_marker_ranking_by_frozen_L2.tsv"
  )
)

write_tsv(
  marker_selected,
  file.path(
    TAB,
    "20A1_marker_panel_selected.tsv"
  )
)

cat(
  "Selected marker rows:",
  nrow(marker_selected),
  "\n"
)

cat("20A1 BLOCK 1 COMPLETE\n")

# ============================================================
# 9. FROZEN-L2 EXPRESSION SUMMARIES
# ============================================================

cat("\nSTEP 8/10 - computing expression on frozen L2 identities\n")

# Restrict expression matrix to the analysis-eligible frozen-L2 cells.
cell_order <- ledger_L2$cell

expr_L2 <- expr[
  ,
  cell_order,
  drop = FALSE
]

cell_meta <- ledger_L2 %>%
  mutate(
    cell_index = seq_len(n())
  )

# Only genes selected for the publication marker panels are needed.
selected_genes <- intersect(
  unique(marker_selected$gene),
  rownames(expr_L2)
)

missing_selected_genes <- setdiff(
  unique(marker_selected$gene),
  rownames(expr_L2)
)

if (length(missing_selected_genes) > 0) {
  cat(
    "Selected genes absent from normalized matrix:",
    length(missing_selected_genes),
    "\n"
  )

  write_tsv(
    tibble(gene = missing_selected_genes),
    file.path(
      TAB,
      "20A1_selected_genes_missing_from_expression.tsv"
    )
  )
}

stopifnot(length(selected_genes) > 0)

marker_selected_present <- marker_selected %>%
  filter(gene %in% selected_genes)

# ------------------------------------------------------------
# Sparse descriptive aggregation
# ------------------------------------------------------------

summarise_frozen_expression <- function(
  L1_name,
  genes
) {

  md <- cell_meta %>%
    filter(
      final_compartment == L1_name
    )

  stopifnot(nrow(md) > 0)

  disp <- display %>%
    filter(
      L1 == L1_name
    ) %>%
    arrange(
      L2_order_within_L1
    )

  stopifnot(nrow(disp) > 0)

  genes <- intersect(
    genes,
    rownames(expr_L2)
  )

  stopifnot(length(genes) > 0)

  mat <- expr_L2[
    genes,
    md$cell,
    drop = FALSE
  ]

  out <- lapply(
    disp$final_L2_label,
    function(lbl) {

      idx <- which(
        md$final_L2_label == lbl
      )

      if (length(idx) == 0) {
        return(
          tibble(
            final_L2_label = lbl,
            gene = genes,
            n_cells = 0L,
            mean_expression = NA_real_,
            pct_expressing = NA_real_
          )
        )
      }

      mm <- mat[
        ,
        idx,
        drop = FALSE
      ]

      mean_expr <- Matrix::rowMeans(mm)

      pct_expr <- 100 *
        Matrix::rowMeans(
          mm > 0
        )

      tibble(
        final_L2_label = lbl,
        gene = genes,
        n_cells = length(idx),
        mean_expression =
          as.numeric(mean_expr),
        pct_expressing =
          as.numeric(pct_expr)
      )
    }
  ) %>%
    bind_rows()

  out %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          display_code,
          legend_label,
          hex,
          L2_order_within_L1
        ),
      by = "final_L2_label"
    ) %>%
    mutate(
      L1 = L1_name
    )
}

expression_summary <- bind_rows(
  lapply(
    spec$L1,
    function(parent) {

      genes <- marker_selected_present %>%
        filter(
          L1 == parent
        ) %>%
        pull(gene) %>%
        unique()

      summarise_frozen_expression(
        parent,
        genes
      )
    }
  )
)

stopifnot(nrow(expression_summary) > 0)
stopifnot(!anyNA(expression_summary$display_code))

# ------------------------------------------------------------
# Gene-wise z score across frozen identities
# ------------------------------------------------------------

expression_summary <- expression_summary %>%
  group_by(
    L1,
    gene
  ) %>%
  mutate(
    expression_z = {
      x <- mean_expression

      if (
        all(is.na(x)) ||
        sd(x, na.rm = TRUE) == 0 ||
        is.na(sd(x, na.rm = TRUE))
      ) {
        rep(0, length(x))
      } else {
        as.numeric(
          scale(x)
        )
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    expression_z_clipped =
      pmax(
        pmin(
          expression_z,
          2.5
        ),
        -2.5
      )
  )

write_tsv(
  expression_summary,
  file.path(
    TAB,
    "20A1_frozen_L2_expression_summary.tsv"
  )
)

# ============================================================
# 10. PLOT HELPERS
# ============================================================

cat("STEP 9/10 - generating marker heatmaps and DotPlots\n")

theme_marker <- function(base_size = 9) {

  theme_bw(
    base_size = base_size
  ) +
    theme(
      panel.grid.major =
        element_line(
          linewidth = 0.20,
          colour = "grey90"
        ),
      panel.grid.minor =
        element_blank(),
      axis.title =
        element_blank(),
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 8
        ),
      axis.text.y =
        element_text(
          size = 8
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 12
        ),
      plot.subtitle =
        element_text(
          size = 9,
          margin = margin(
            b = 8
          )
        ),
      legend.title =
        element_text(
          size = 8
        ),
      legend.text =
        element_text(
          size = 8
        ),
      strip.text =
        element_text(
          face = "bold"
        ),
      plot.margin =
        margin(
          8, 10, 8, 8
        )
    )
}

save_three_formats <- function(
  p,
  directory,
  stem,
  width,
  height
) {

  pdf_dir <- file.path(
    directory,
    "PDF"
  )

  png_dir <- file.path(
    directory,
    "PNG"
  )

  svg_dir <- file.path(
    directory,
    "SVG"
  )

  dir.create(
    pdf_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    png_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    svg_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  ggsave(
    file.path(
      pdf_dir,
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
      png_dir,
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
        svg_dir,
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
      filename = file.path(
        svg_dir,
        paste0(stem, ".svg")
      ),
      width = width,
      height = height,
      onefile = TRUE
    )

    print(p)
    dev.off()
  }
}

# ============================================================
# 11. BUILD EACH COMPARTMENT
# ============================================================

figure_inventory <- tibble(
  L1 = character(),
  panel = character(),
  figure = character(),
  n_L2 = integer(),
  n_genes = integer(),
  n_cells = integer()
)

selected_panel_final <- list()

for (parent in spec$L1) {

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "Compartment:",
    parent,
    "\n"
  )

  disp <- display %>%
    filter(
      L1 == parent
    ) %>%
    arrange(
      L2_order_within_L1
    )

  dd <- expression_summary %>%
    filter(
      L1 == parent
    )

  selected <- marker_selected_present %>%
    filter(
      L1 == parent
    ) %>%
    left_join(
      disp %>%
        select(
          final_L2_label,
          L2_order_within_L1
        ),
      by = "final_L2_label",
      suffix = c("", "_display")
    )

  # ----------------------------------------------------------
  # Resolve duplicated marker genes across L2 identities.
  # A gene is assigned to the identity for which marker
  # evidence is strongest. This keeps the publication panel
  # compact and avoids repeated y-axis genes.
  # ----------------------------------------------------------

  selected_unique <- selected %>%
    arrange(
      desc(n_supporting_clusters),
      desc(max_marker_score),
      desc(max_log2FC),
      min_p_adj,
      L2_order_within_L1
    ) %>%
    group_by(gene) %>%
    slice_head(n = 1) %>%
    ungroup()

  # Guarantee representation of every frozen L2.
  represented <- unique(
    selected_unique$final_L2_label
  )

  missing_L2 <- setdiff(
    disp$final_L2_label,
    represented
  )

  if (length(missing_L2) > 0) {

    rescue <- marker_ranked %>%
      filter(
        L1 == parent,
        final_L2_label %in% missing_L2,
        gene %in% rownames(expr_L2)
      ) %>%
      anti_join(
        selected_unique %>%
          select(gene),
        by = "gene"
      ) %>%
      group_by(
        final_L2_label
      ) %>%
      arrange(
        desc(n_supporting_clusters),
        desc(max_marker_score),
        desc(max_log2FC),
        .by_group = TRUE
      ) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      left_join(
        disp %>%
          select(
            final_L2_label,
            L2_order_within_L1
          ),
        by = "final_L2_label"
      )

    selected_unique <- bind_rows(
      selected_unique,
      rescue
    )
  }

  # Limit to max 4 markers per identity after deduplication.
  selected_unique <- selected_unique %>%
    group_by(
      final_L2_label
    ) %>%
    arrange(
      desc(n_supporting_clusters),
      desc(max_marker_score),
      desc(max_log2FC),
      .by_group = TRUE
    ) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    arrange(
      L2_order_within_L1,
      desc(max_marker_score)
    )

  genes_plot <- unique(
    selected_unique$gene
  )

  stopifnot(length(genes_plot) > 0)

  selected_panel_final[[parent]] <-
    selected_unique %>%
    mutate(
      L1 = parent
    )

  plot_dat <- dd %>%
    filter(
      gene %in% genes_plot
    )

  # Gene order follows the frozen identity from which each
  # marker was selected.
  gene_key <- selected_unique %>%
    select(
      gene,
      final_L2_label,
      L2_order_within_L1
    ) %>%
    distinct() %>%
    mutate(
      marker_order =
        row_number()
    )

  gene_levels <- gene_key %>%
    arrange(
      L2_order_within_L1,
      marker_order
    ) %>%
    pull(gene) %>%
    unique()

  plot_dat <- plot_dat %>%
    mutate(
      display_code = factor(
        display_code,
        levels =
          disp$display_code
      ),
      gene = factor(
        gene,
        levels =
          rev(gene_levels)
      )
    )

  # ----------------------------------------------------------
  # B. MARKER HEATMAP
  # ----------------------------------------------------------

  p_heat <- ggplot(
    plot_dat,
    aes(
      x = display_code,
      y = gene,
      fill = expression_z_clipped
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.25
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "Scaled\nexpression"
    ) +
    labs(
      title = paste0(
        gsub("_", " ", parent),
        " - frozen L2 marker heatmap"
      ),
      subtitle = paste0(
        length(genes_plot),
        " discriminating genes | ",
        nrow(disp),
        " frozen L2 identities"
      )
    ) +
    theme_marker()

  # ----------------------------------------------------------
  # C. MARKER DOTPLOT
  # ----------------------------------------------------------

  p_dot <- ggplot(
    plot_dat,
    aes(
      x = display_code,
      y = gene
    )
  ) +
    geom_point(
      aes(
        size = pct_expressing,
        colour = expression_z_clipped
      ),
      alpha = 0.92
    ) +
    scale_size_continuous(
      range = c(0.6, 6.0),
      limits = c(0, 100),
      breaks = c(
        25,
        50,
        75,
        100
      ),
      name = "% expressing"
    ) +
    scale_colour_gradient2(
      low = "#2166AC",
      mid = "grey92",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = squish,
      name = "Scaled\nexpression"
    ) +
    labs(
      title = paste0(
        gsub("_", " ", parent),
        " - frozen L2 marker DotPlot"
      ),
      subtitle = paste0(
        "Dot size = % expressing | colour = gene-wise scaled mean expression"
      )
    ) +
    theme_marker()

  # Dynamic figure dimensions.
  fig_width <- max(
    6.8,
    0.65 * nrow(disp) + 4.0
  )

  fig_height <- max(
    5.5,
    0.22 * length(genes_plot) + 2.8
  )

  parent_dir <- file.path(
    FIG,
    parent
  )

  save_three_formats(
    p_heat,
    file.path(
      parent_dir,
      "B_marker_heatmap"
    ),
    paste0(
      "20A1_",
      parent,
      "_marker_heatmap"
    ),
    width = fig_width,
    height = fig_height
  )

  save_three_formats(
    p_dot,
    file.path(
      parent_dir,
      "C_marker_dotplot"
    ),
    paste0(
      "20A1_",
      parent,
      "_marker_dotplot"
    ),
    width = fig_width + 0.8,
    height = fig_height
  )

  n_cells_parent <- ledger_L2 %>%
    filter(
      final_compartment == parent
    ) %>%
    nrow()

  figure_inventory <- bind_rows(
    figure_inventory,
    tibble(
      L1 = parent,
      panel = c(
        "B_marker_heatmap",
        "C_marker_dotplot"
      ),
      figure = c(
        paste0(
          "20A1_",
          parent,
          "_marker_heatmap"
        ),
        paste0(
          "20A1_",
          parent,
          "_marker_dotplot"
        )
      ),
      n_L2 = nrow(disp),
      n_genes = length(genes_plot),
      n_cells = n_cells_parent
    )
  )

  cat(
    sprintf(
      "%-15s : %6d cells | %2d L2 | %2d genes | heatmap + DotPlot\n",
      parent,
      n_cells_parent,
      nrow(disp),
      length(genes_plot)
    )
  )
}

# ============================================================
# 12. SAVE FINAL PANEL SELECTION
# ============================================================

final_marker_panel <- bind_rows(
  selected_panel_final
)

write_tsv(
  final_marker_panel,
  file.path(
    TAB,
    "20A1_final_marker_panel.tsv"
  )
)

write_tsv(
  figure_inventory,
  file.path(
    TAB,
    "20A1_figure_inventory.tsv"
  )
)

# ============================================================
# 13. FORMAL CHECKS
# ============================================================

cat("\nSTEP 10/10 - formal checks\n")

expected_L2_counts <- tibble(
  L1 = c(
    "Epithelial",
    "Fibroblast",
    "Immune",
    "Endothelial",
    "Mural",
    "Sweat_gland",
    "Melanocyte",
    "Sebocyte"
  ),
  expected_n = c(
    9L,
    10L,
    9L,
    10L,
    7L,
    5L,
    4L,
    5L
  )
)

observed_L2_counts <- display %>%
  count(
    L1,
    name = "observed_n"
  ) %>%
  inner_join(
    expected_L2_counts,
    by = "L1"
  )

pdf_files <- list.files(
  FIG,
  pattern = "\\.pdf$",
  recursive = TRUE,
  full.names = TRUE
)

png_files <- list.files(
  FIG,
  pattern = "\\.png$",
  recursive
 = TRUE,
  full.names = TRUE
)

svg_files <- list.files(
  FIG,
  pattern = "\\.svg$",
  recursive = TRUE,
  full.names = TRUE
)

# Every frozen L2 must be represented in the display specification.
L2_count_ok <-
  nrow(observed_L2_counts) == 8 &&
  all(
    observed_L2_counts$observed_n ==
      observed_L2_counts$expected_n
  )

# Every frozen L2 should have at least one selected marker
# in the final publication marker panel.
marker_L2_coverage <- final_marker_panel %>%
  distinct(
    L1,
    final_L2_label
  ) %>%
  count(
    L1,
    name = "n_L2_with_marker"
  ) %>%
  right_join(
    expected_L2_counts,
    by = "L1"
  ) %>%
  mutate(
    n_L2_with_marker =
      replace_na(
        n_L2_with_marker,
        0L
      ),
    complete =
      n_L2_with_marker ==
        expected_n
  )

write_tsv(
  marker_L2_coverage,
  file.path(
    TAB,
    "20A1_marker_L2_coverage.tsv"
  )
)

all_L2_marker_coverage <-
  nrow(marker_L2_coverage) == 8 &&
  all(marker_L2_coverage$complete)

# Expression summary must cover every frozen L2 identity.
expression_L2_coverage <-
  expression_summary %>%
  distinct(
    L1,
    final_L2_label
  ) %>%
  nrow() == 59

# All plotted genes must exist in the normalized expression matrix.
all_final_genes_present <-
  all(
    final_marker_panel$gene %in%
      rownames(expr_L2)
  )

# There should be two figure designs per compartment.
expected_figure_designs <- 16L

checks <- tibble(
  check = c(
    "analysis_eligible_L2_82494",
    "frozen_L2_identities_59",
    "eight_parent_compartments",
    "expected_L2_counts_by_parent",
    "cluster_mapping_matches_display",
    "all_frozen_L2_have_marker",
    "expression_summary_covers_59_L2",
    "all_final_marker_genes_in_expression",
    "marker_selection_from_existing_evidence",
    "expression_recomputed_on_frozen_L2",
    "no_marker_inference_recalculated",
    "no_clustering_recalculated",
    "no_annotation_modified",
    "figure_designs_16",
    "PDF_files_16",
    "PNG_files_16",
    "SVG_files_16"
  ),
  pass = c(
    nrow(ledger_L2) == 82494,
    n_distinct(
      ledger_L2$final_L2_label
    ) == 59,
    n_distinct(
      ledger_L2$final_compartment
    ) == 8,
    L2_count_ok,
    !anyNA(
      mapping_all$display_code
    ),
    all_L2_marker_coverage,
    expression_L2_coverage,
    all_final_genes_present,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    nrow(figure_inventory) ==
      expected_figure_designs,
    length(pdf_files) ==
      expected_figure_designs,
    length(png_files) ==
      expected_figure_designs,
    length(svg_files) ==
      expected_figure_designs
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "20A1_formal_checks.tsv"
  )
)

write_tsv(
  observed_L2_counts,
  file.path(
    TAB,
    "20A1_observed_L2_counts.tsv"
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

  cat("\n===== MARKER L2 COVERAGE =====\n")

  print(
    as.data.frame(
      marker_L2_coverage
    ),
    row.names = FALSE
  )

  stop(
    "20A1 formal checks failed. ",
    "Inspect tables before accepting marker panels."
  )
}

# ============================================================
# 14. REPRODUCIBILITY RECORD
# ============================================================

record <- c(
  "PIPELINE_V2 - 20A1 FROZEN L2 MARKER PANELS",
  "",
  "STATUS: GENERATED - VISUAL QC REQUIRED",
  "",
  "Analysis universe:",
  "  82,494 analysis-eligible cells with frozen L2 identities.",
  "  59 frozen L2 identities across 8 parent compartments.",
  "",
  "Marker selection:",
  "  Existing annotation-stage positive-marker evidence was reused.",
  "  Cluster-level marker statistics were used only to select candidate genes.",
  "  No marker differential-expression test was rerun.",
  "",
  "Expression visualization:",
  "  Mean normalized RNA expression was recomputed descriptively",
  "  after grouping cells by the final frozen L2 identities.",
  "  Percent expressing was computed directly on the same frozen-L2 cells.",
  "  Heatmap colour represents gene-wise scaled mean expression.",
  "  DotPlot size represents percent expressing.",
  "  DotPlot colour represents gene-wise scaled mean expression.",
  "",
  "Identity handling:",
  "  Frozen cluster-to-L2 mappings were recovered from annotation outputs.",
  "  Multiple source clusters may map to one final frozen L2 identity.",
  "  Biological annotations were not modified.",
  "",
  "Display:",
  "  19C3 short display codes and frozen colours were reused.",
  "  Display codes do not replace biological identity labels.",
  "",
  "Outputs:",
  "  8 marker heatmaps.",
  "  8 marker DotPlots.",
  "  PDF + PNG + editable SVG.",
  "",
  "Clustering recalculated: NO",
  "Annotation modified: NO",
  "Marker inference recalculated: NO",
  "PASS"
)

writeLines(
  record,
  con = file.path(
    OUT,
    "20A1_REPRODUCIBILITY_RECORD.txt"
  )
)

# ============================================================
# 15. FINAL SUMMARY
# ============================================================

cat("\n===== FIGURE INVENTORY =====\n")

print(
  as.data.frame(
    figure_inventory
  ),
  row.names = FALSE
)

cat("\n===== MARKER COVERAGE =====\n")

print(
  as.data.frame(
    marker_L2_coverage
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("20A1 COMPLETE\n")
cat("============================================================\n")
cat("Frozen L2 cells: 82,494\n")
cat("Frozen L2 identities: 59\n")
cat("Parent compartments: 8\n")
cat("Marker heatmaps generated: 8\n")
cat("Marker DotPlots generated: 8\n")
cat("Figure designs generated: 16\n")
cat("PDF files generated: 16\n")
cat("PNG files generated: 16\n")
cat("SVG files generated: 16\n")
cat("Editable SVG export: YES\n")
cat("Cluster-level marker evidence reused: YES\n")
cat("Expression summarized on frozen L2 cells: YES\n")
cat("Marker inference recalculated: NO\n")
cat("Clustering recalculated: NO\n")
cat("Annotation modified: NO\n")
cat("STATUS: VISUAL QC REQUIRED\n")
cat("PASS\n")
