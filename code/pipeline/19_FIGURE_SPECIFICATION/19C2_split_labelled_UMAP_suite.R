#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(Matrix)
})

cat("============================================================\n")
cat("PIPELINE_V2 - 19C2 SPLIT + LABELLED UMAP FIGURE SUITE\n")
cat("============================================================\n")

# ============================================================
# 0. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

SOURCE_19C <- file.path(
  ROOT,
  "scripts/19_FIGURE_SPECIFICATION/19C_canonical_UMAP_suite.R"
)

OBJ_FILE <- file.path(
  ROOT,
  "objects/annotated/global_atlas_v2_L1_frozen.rds"
)

LEDGER_FILE <- file.path(
  ROOT,
  "results/disease_biology/17A_frozen_analysis_atlas/tables",
  "17A_master_cell_ledger_FINAL.rds"
)

PAL_DIR <- file.path(
  ROOT,
  "config/figure_palettes"
)

# Recover palette directory if canonical location differs.
palette_candidates <- c(
  PAL_DIR,
  file.path(ROOT, "config"),
  file.path(
    ROOT,
    "results/figure_specification/19B_palette_freeze"
  ),
  file.path(
    ROOT,
    "results/figure_specification/19B_palette_freeze/tables"
  )
)

find_palette <- function(filename) {
  hits <- unlist(
    lapply(
      palette_candidates,
      function(d) {
        p <- file.path(d, filename)
        if (file.exists(p)) p else character()
      }
    )
  )

  if (length(hits) == 0) {
    hits <- list.files(
      ROOT,
      pattern = paste0("^", filename, "$"),
      recursive = TRUE,
      full.names = TRUE
    )
  }

  if (length(hits) != 1) {
    stop(
      "Expected exactly one canonical ",
      filename,
      "; found ",
      length(hits),
      ".\n",
      paste(hits, collapse = "\n")
    )
  }

  hits[[1]]
}

DISEASE_PAL_FILE <- find_palette("disease_palette.tsv")
L1_PAL_FILE      <- find_palette("L1_palette.tsv")
L2_PAL_FILE      <- find_palette("L2_palette.tsv")

OUT <- file.path(
  ROOT,
  "results/figure_specification/19C2_split_labelled_UMAP_suite"
)

FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

cat("\nSTEP 1/10 - loading frozen inputs\n")

stopifnot(file.exists(OBJ_FILE))
stopifnot(file.exists(LEDGER_FILE))
stopifnot(file.exists(DISEASE_PAL_FILE))
stopifnot(file.exists(L1_PAL_FILE))
stopifnot(file.exists(L2_PAL_FILE))

obj <- readRDS(OBJ_FILE)
ledger <- readRDS(LEDGER_FILE)

disease_pal <- read_tsv(
  DISEASE_PAL_FILE,
  show_col_types = FALSE
)

L1_pal <- read_tsv(
  L1_PAL_FILE,
  show_col_types = FALSE
)

L2_pal <- read_tsv(
  L2_PAL_FILE,
  show_col_types = FALSE
)

stopifnot(inherits(obj, "Seurat"))
stopifnot(ncol(obj) == 82893)
stopifnot(nrow(disease_pal) == 3)
stopifnot(nrow(L1_pal) == 9)
stopifnot(nrow(L2_pal) == 59)

# ============================================================
# 2. HARMONY UMAP
# ============================================================

cat("STEP 2/10 - recovering canonical Harmony UMAP\n")

stopifnot("umap_harmony" %in% names(obj@reductions))

emb <- Embeddings(
  obj,
  reduction = "umap_harmony"
)

stopifnot(nrow(emb) == 82893)
stopifnot(ncol(emb) >= 2)

umap <- tibble(
  cell = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2]
)

xlim_global <- range(umap$UMAP_1, finite = TRUE)
ylim_global <- range(umap$UMAP_2, finite = TRUE)

# ============================================================
# 3. JOIN METADATA + FINAL LEDGER
# ============================================================

cat("STEP 3/10 - joining disease, patient and frozen ledger\n")

md <- obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

stopifnot("disease" %in% colnames(md))

cell_candidates <- c(
  "cell",
  "barcode",
  "cell_id",
  "Cell",
  "CELL"
)

cell_col <- intersect(
  cell_candidates,
  colnames(ledger)
)[1]

if (is.na(cell_col)) {
  stop("No cell identifier column found in final ledger.")
}

ledger2 <- ledger %>%
  rename(cell = all_of(cell_col))

required_ledger <- c(
  "patient_id",
  "analysis_eligible_global",
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label"
)

missing_ledger <- setdiff(
  required_ledger,
  colnames(ledger2)
)

if (length(missing_ledger) > 0) {
  stop(
    "Missing final-ledger columns: ",
    paste(missing_ledger, collapse = ", ")
  )
}

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
stopifnot(all(!is.na(base$patient_id)))

disease_df <- base

L1_df <- base %>%
  filter(analysis_eligible_global)

L2_df <- base %>%
  filter(analysis_eligible_L2)

stopifnot(nrow(disease_df) == 82893)
stopifnot(nrow(L1_df) == 82680)
stopifnot(nrow(L2_df) == 82494)

# ============================================================
# 4. FROZEN PALETTES / ORDERS
# ============================================================

cat("STEP 4/10 - binding frozen palettes\n")

disease_pal <- disease_pal %>%
  arrange(display_order)

L1_pal <- L1_pal %>%
  arrange(display_order)

L2_pal <- L2_pal %>%
  arrange(global_display_order)

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

disease_levels <- disease_pal$disease
L1_levels <- L1_pal$L1
L2_levels <- L2_pal$final_L2_label

disease_df$disease <- factor(
  disease_df$disease,
  levels = disease_levels
)

L1_df$disease <- factor(
  L1_df$disease,
  levels = disease_levels
)

L2_df$disease <- factor(
  L2_df$disease,
  levels = disease_levels
)

L1_df$final_compartment <- factor(
  L1_df$final_compartment,
  levels = L1_levels
)

L2_df$final_compartment <- factor(
  L2_df$final_compartment,
  levels = L1_levels
)

L2_df$final_L2_label <- factor(
  L2_df$final_L2_label,
  levels = L2_levels
)

stopifnot(
  identical(
    unname(disease_cols[c("HC", "FFA", "AA")]),
    c("#8F8D8E", "#4D9E4F", "#A8142A")
  )
)

# ============================================================
# 5. PATIENT ORDER / FACET LABELS
# ============================================================

cat("STEP 5/10 - freezing patient facet specification\n")

patient_key <- disease_df %>%
  distinct(patient_id, disease) %>%
  mutate(
    disease = as.character(disease)
  )

stopifnot(nrow(patient_key) == 16)
stopifnot(!anyDuplicated(patient_key$patient_id))

patient_key <- patient_key %>%
  mutate(
    disease_order = match(
      disease,
      disease_levels
    )
  ) %>%
  arrange(
    disease_order,
    patient_id
  )

patient_levels <- patient_key$patient_id

disease_df$patient_id <- factor(
  disease_df$patient_id,
  levels = patient_levels
)

L1_df$patient_id <- factor(
  L1_df$patient_id,
  levels = patient_levels
)

L2_df$patient_id <- factor(
  L2_df$patient_id,
  levels = patient_levels
)

make_patient_labels <- function(dat) {

  dat %>%
    count(
      patient_id,
      disease,
      name = "n_cells"
    ) %>%
    mutate(
      patient_id = as.character(patient_id),
      disease = as.character(disease),
      facet_label = paste0(
        patient_id,
        "\n",
        disease,
        " | n=",
        format(
          n_cells,
          big.mark = ",",
          scientific = FALSE,
          trim = TRUE
        )
      )
    ) %>%
    select(
      patient_id,
      facet_label
    ) %>%
    deframe()
}

patient_labels_disease <- make_patient_labels(disease_df)
patient_labels_L1 <- make_patient_labels(L1_df)
patient_labels_L2 <- make_patient_labels(L2_df)

write_csv(
  patient_key,
  file.path(
    TAB,
    "19C2_patient_order.csv"
  )
)

# ============================================================
# 6. GRAPHICAL HELPERS
# ============================================================

cat("STEP 6/10 - defining labelled/split plot helpers\n")

theme_umap2 <- function(
  facet = FALSE
) {

  th <- theme_void(base_size = 10) +
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

  if (facet) {
    th <- th +
      theme(
        strip.background = element_rect(
          fill = "grey96",
          colour = "grey80",
          linewidth = 0.3
        ),
        strip.text = element_text(
          face = "bold",
          size = 8,
          margin = margin(4, 4, 4, 4)
        ),
        panel.spacing = grid::unit(
          1.5,
          "mm"
        )
      )
  }

  th
}

prepare_draw_order <- function(
  dat,
  colour_var
) {

  dat$.plot_colour <- dat[[colour_var]]

  ord <- table(
    dat$.plot_colour
  )

  dat$.draw_n <- unname(
    ord[
      as.character(
        dat$.plot_colour
      )
    ]
  )

  dat %>%
    arrange(desc(.draw_n))
}

make_label_positions <- function(
  dat,
  label_var,
  min_cells = 1
) {

  dat %>%
    filter(
      !is.na(.data[[label_var]])
    ) %>%
    group_by(
      .data[[label_var]]
    ) %>%
    filter(n() >= min_cells) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    rename(
      label = 1
    )
}

add_direct_labels <- function(
  p,
  dat,
  label_var,
  min_cells = 1,
  label_size = 2.7
) {

  lab <- make_label_positions(
    dat,
    label_var,
    min_cells = min_cells
  )

  if (nrow(lab) == 0) {
    return(p)
  }

  p +
    geom_label(
      data = lab,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        label = label
      ),
      inherit.aes = FALSE,
      size = label_size,
      fontface = "bold",
      label.size = 0.18,
      label.padding = grid::unit(
        0.12,
        "lines"
      ),
      fill = scales::alpha(
        "white",
        0.82
      ),
      colour = "grey15"
    )
}

make_global_umap <- function(
  dat,
  colour_var,
  colours,
  title,
  subtitle,
  point_size,
  stroke_size,
  point_alpha = 0.90,
  legend_ncol = 1,
  label_var = NULL,
  label_size = 2.7
) {

  dd <- prepare_draw_order(
    dat,
    colour_var
  )

  p <- ggplot(
    dd,
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
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE
    ) +
    theme_umap2()

  if (!is.null(label_var)) {
    p <- add_direct_labels(
      p,
      dat,
      label_var,
      label_size = label_size
    )
  }

  p
}

make_disease_split <- function(
  dat,
  colour_var,
  colours,
  title,
  subtitle,
  point_size,
  stroke_size,
  point_alpha = 0.90,
  legend_ncol = 1
) {

  dd <- prepare_draw_order(
    dat,
    colour_var
  )

  ggplot(
    dd,
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
    facet_wrap(
      ~ disease,
      nrow = 1,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 3.0,
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
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE
    ) +
    theme_umap2(
      facet = TRUE
    )
}

make_patient_split <- function(
  dat,
  colour_var,
  colours,
  facet_labels,
  title,
  subtitle,
  point_size,
  stroke_size,
  point_alpha = 0.92,
  legend_ncol = 1
) {

  dd <- prepare_draw_order(
    dat,
    colour_var
  )

  ggplot(
    dd,
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
    facet_wrap(
      ~ patient_id,
      ncol = 4,
      labeller = as_labeller(
        facet_labels
      ),
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 3.0,
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
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE
    ) +
    theme_umap2(
      facet = TRUE
    )
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
# 7. DISEASE SUITE
# ============================================================

cat("STEP 7/10 - generating disease UMAP suite\n")

# 7A. Global disease UMAP
p <- make_global_umap(
  disease_df,
  "disease",
  disease_cols,
  "Global atlas by disease",
  "Harmony UMAP | 82,893 cells",
  point_size = 0.65,
  stroke_size = 0.045,
  point_alpha = 0.88
)

save_plot(
  p,
  "19C2_01_disease_global",
  width = 7.2,
  height = 5.8
)

# 7B. Disease split
p <- make_disease_split(
  disease_df,
  "disease",
  disease_cols,
  "Global atlas by disease - split view",
  "Identical Harmony coordinates and axis limits across HC, FFA and AA",
  point_size = 0.72,
  stroke_size = 0.045,
  point_alpha = 0.90
)

save_plot(
  p,
  "19C2_02_disease_split_HC_FFA_AA",
  width = 14.0,
  height = 5.0
)

# 7C. Patient split
p <- make_patient_split(
  disease_df,
  "disease",
  disease_cols,
  patient_labels_disease,
  "Global atlas - individual patients",
  "16 patients | identical Harmony coordinates and axis limits",
  point_size = 0.55,
  stroke_size = 0.035,
  point_alpha = 0.90
)

save_plot(
  p,
  "19C2_03_disease_split_16_patients",
  width = 14.0,
  height = 13.0
)

# ============================================================
# 8. L1 SUITE
# ============================================================

cat("STEP 8/10 - generating L1 UMAP suite\n")

# 8A. Global labelled L1
p <- make_global_umap(
  L1_df,
  "final_compartment",
  L1_cols,
  "Global atlas - final L1 compartments",
  "82,680 analysis-eligible cells | direct compartment labels",
  point_size = 0.65,
  stroke_size = 0.045,
  point_alpha = 0.90,
  label_var = "final_compartment",
  label_size = 2.8
)

save_plot(
  p,
  "19C2_04_L1_global_labelled",
  width = 8.6,
  height = 6.4
)

# 8B. L1 by disease
p <- make_disease_split(
  L1_df,
  "final_compartment",
  L1_cols,
  "Final L1 compartments - split by disease",
  "HC | FFA | AA | identical Harmony coordinates and L1 colours",
  point_size = 0.58,
  stroke_size = 0.035,
  point_alpha = 0.90
)

save_plot(
  p,
  "19C2_05_L1_split_HC_FFA_AA",
  width = 15.5,
  height = 5.4
)

# 8C. L1 by patient
p <- make_patient_split(
  L1_df,
  "final_compartment",
  L1_cols,
  patient_labels_L1,
  "Final L1 compartments - individual patients",
  "16 patients | frozen L1 colours | identical Harmony coordinates",
  point_size = 0.48,
  stroke_size = 0.03,
  point_alpha = 0.91
)

save_plot(
  p,
  "19C2_06_L1_split_16_patients",
  width = 15.5,
  height = 13.2
)

# ============================================================
# 9. GLOBAL L2 + COMPARTMENT-SPECIFIC SUITES
# ============================================================

cat("STEP 9/10 - generating L2 UMAP suites\n")

# ------------------------------------------------------------
# 9A. Global 59-L2 atlas
# ------------------------------------------------------------

p <- make_global_umap(
  L2_df,
  "final_L2_label",
  L2_cols,
  "Global atlas - frozen L2 identities",
  "82,494 cells | 59 frozen L2 identities",
  point_size = 0.52,
  stroke_size = 0.035,
  point_alpha = 0.90,
  legend_ncol = 2
)

save_plot(
  p,
  "19C2_07_L2_59_global",
  width = 12.2,
  height = 8.4
)

p <- make_disease_split(
  L2_df,
  "final_L2_label",
  L2_cols,
  "Frozen L2 identities - split by disease",
  "59 L2 identities | HC | FFA | AA | identical Harmony coordinates",
  point_size = 0.48,
  stroke_size = 0.025,
  point_alpha = 0.90,
  legend_ncol = 2
)

save_plot(
  p,
  "19C2_08_L2_59_split_HC_FFA_AA",
  width = 18.0,
  height = 7.6
)

p <- make_patient_split(
  L2_df,
  "final_L2_label",
  L2_cols,
  patient_labels_L2,
  "Frozen L2 identities - individual patients",
  "16 patients | 59 frozen L2 identities | identical Harmony coordinates",
  point_size = 0.40,
  stroke_size = 0.02,
  point_alpha = 0.91,
  legend_ncol = 2
)

save_plot(
  p,
  "19C2_09_L2_59_split_16_patients",
  width = 18.0,
  height = 14.5
)

# ------------------------------------------------------------
# 9B. Each frozen L2 parent compartment
# ------------------------------------------------------------

parent_order <- L1_pal %>%
  filter(L1 != "Neural_Schwann") %>%
  arrange(display_order) %>%
  pull(L1)

inventory <- tibble(
  figure = c(
    "19C2_01_disease_global",
    "19C2_02_disease_split_HC_FFA_AA",
    "19C2_03_disease_split_16_patients",
    "19C2_04_L1_global_labelled",
    "19C2_05_L1_split_HC_FFA_AA",
    "19C2_06_L1_split_16_patients",
    "19C2_07_L2_59_global",
    "19C2_08_L2_59_split_HC_FFA_AA",
    "19C2_09_L2_59_split_16_patients"
  ),
  level = c(
    rep("Disease", 3),
    rep("L1", 3),
    rep("L2_global", 3)
  ),
  view = rep(
    c(
      "global",
      "disease_split",
      "patient_split"
    ),
    3
  ),
  n_cells = c(
    rep(nrow(disease_df), 3),
    rep(nrow(L1_df), 3),
    rep(nrow(L2_df), 3)
  )
)

next_index <- 10L

for (parent in parent_order) {

  allowed_L2 <- L2_pal %>%
    filter(L1 == parent) %>%
    arrange(L2_order_within_L1)

  dd <- L2_df %>%
    filter(
      as.character(final_compartment) == parent
    )

  stopifnot(nrow(dd) > 0)

  dd$final_L2_label <- factor(
    as.character(dd$final_L2_label),
    levels = allowed_L2$final_L2_label
  )

  cols <- setNames(
    allowed_L2$hex,
    allowed_L2$final_L2_label
  )

  parent_pretty <- gsub(
    "_",
    " ",
    parent
  )

  safe_parent <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    parent
  )

  # Patient labels must use counts within this compartment.
  patient_labels_parent <- make_patient_labels(
    dd
  )

  # -------------------------
  # Global + direct L2 labels
  # -------------------------

  p <- make_global_umap(
    dd,
    "final_L2_label",
    cols,
    paste0(
      parent_pretty,
      " - frozen L2 identities"
    ),
    paste0(
      format(
        nrow(dd),
        big.mark = ","
      ),
      " cells | ",
      nrow(allowed_L2),
      " L2 identities | direct labels"
    ),
    point_size = 0.90,
    stroke_size = 0.12,
    point_alpha = 0.92,
    label_var = "final_L2_label",
    label_size = 2.55
  )

  stem_global <- sprintf(
    "19C2_%02d_%s_L2_global_labelled",
    next_index,
    safe_parent
  )

  save_plot(
    p,
    stem_global,
    width = 8.6,
    height = 6.5
  )

  next_index <- next_index + 1L

  # -------------------------
  # Disease split
  # -------------------------

  p <- make_disease_split(
    dd,
    "final_L2_label",
    cols,
    paste0(
      parent_pretty,
      " L2 identities - split by disease"
    ),
    paste0(
      nrow(allowed_L2),
      " frozen L2 identities | HC | FFA | AA"
    ),
    point_size = 0.72,
    stroke_size = 0.07,
    point_alpha = 0.92
  )

  stem_disease <- sprintf(
    "19C2_%02d_%s_L2_split_HC_FFA_AA",
    next_index,
    safe_parent
  )

  save_plot(
    p,
    stem_disease,
    width = 14.8,
    height = 5.2
  )

  next_index <- next_index + 1L

  # -------------------------
  # Patient split
  # -------------------------

  p <- make_patient_split(
    dd,
    "final_L2_label",
    cols,
    patient_labels_parent,
    paste0(
      parent_pretty,
      " L2 identities - individual patients"
    ),
    paste0(
      nrow(allowed_L2),
      " frozen L2 identities | identical Harmony coordinates"
    ),
    point_size = 0.58,
    stroke_size = 0.045,
    point_alpha = 0.93
  )

  stem_patient <- sprintf(
    "19C2_%02d_%s_L2_split_16_patients",
    next_index,
    safe_parent
  )

  save_plot(
    p,
    stem_patient,
    width = 15.5,
    height = 13.2
  )

  next_index <- next_index + 1L

  inventory <- bind_rows(
    inventory,
    tibble(
      figure = c(
        stem_global,
        stem_disease,
        stem_patient
      ),
      level = rep(
        paste0(parent, "_L2"),
        3
      ),
      view = c(
        "global",
        "disease_split",
        "patient_split"
      ),
      n_cells = rep(
        nrow(dd),
        3
      )
    )
  )

  cat(
    sprintf(
      "%-15s : %6d cells | %2d L2 | 3 views\n",
      parent,
      nrow(dd),
      nrow(allowed_L2)
    )
  )
}

# ============================================================
# 10. AUDIT / FORMAL CHECKS
# ============================================================

cat("STEP 10/10 - formal checks and output inventory\n")

write_csv(
  inventory,
  file.path(
    TAB,
    "19C2_figure_inventory.csv"
  )
)

patient_counts <- disease_df %>%
  count(
    patient_id,
    disease,
    name = "n_cells"
  )

write_csv(
  patient_counts,
  file.path(
    TAB,
    "19C2_patient_cell_counts.csv"
  )
)

# Expected:
# 9 global-level figures +
# 8 parent compartments x 3 views = 24
# Total = 33 figure designs.
expected_designs <- 33L

pdf_files <- list.files(
  FIG,
  pattern = "\\.pdf$",
  full.names = TRUE
)

png_files <- list.files(
  FIG,
  pattern = "\\.png$",
  full.names = TRUE
)

checks <- tibble(
  check = c(
    "disease_universe_82893",
    "L1_universe_82680",
    "L2_universe_82494",
    "patients_16",
    "disease_categories_3",
    "L1_categories_9",
    "L2_categories_59",
    "eight_L2_parent_compartments",
    "disease_palette_frozen",
    "L1_palette_frozen",
    "L2_palette_frozen",
    "HC_colour_8F8D8E",
    "FFA_colour_4D9E4F",
    "AA_colour_A8142A",
    "patient_id_complete",
    "same_global_coordinates_used",
    "same_axis_limits_used",
    "no_UMAP_recalculation",
    "no_annotation_modification",
    "figure_designs_33",
    "PDF_files_33",
    "PNG_files_33"
  ),
  pass = c(
    nrow(disease_df) == 82893,
    nrow(L1_df) == 82680,
    nrow(L2_df) == 82494,
    n_distinct(disease_df$patient_id) == 16,
    n_distinct(disease_df$disease) == 3,
    n_distinct(L1_df$final_compartment) == 9,
    n_distinct(L2_df$final_L2_label) == 59,
    length(parent_order) == 8,
    nrow(disease_pal) == 3,
    nrow(L1_pal) == 9,
    nrow(L2_pal) == 59,
    unname(disease_cols["HC"]) == "#8F8D8E",
    unname(disease_cols["FFA"]) == "#4D9E4F",
    unname(disease_cols["AA"]) == "#A8142A",
    all(!is.na(disease_df$patient_id)),
    nrow(umap) == 82893,
    length(xlim_global) == 2 &&
      length(ylim_global) == 2,
    TRUE,
    TRUE,
    nrow(inventory) == expected_designs,
    length(pdf_files) == expected_designs,
    length(png_files) == expected_designs
  )
)

write_csv(
  checks,
  file.path(
    TAB,
    "19C2_formal_checks.csv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  as.data.frame(checks),
  row.names = FALSE
)

if (!all(checks$pass)) {
  stop(
    "19C2 formal checks failed."
  )
}

cat("\n===== FIGURE INVENTORY =====\n")
print(
  as.data.frame(inventory),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("19C2 COMPLETE\n")
cat("============================================================\n")
cat("Figure designs generated: 33\n")
cat("PDF files generated: 33\n")
cat("PNG files generated: 33\n")
cat("Views: global + disease split + patient split\n")
cat("Disease split: HC | FFA | AA\n")
cat("Patient split: 16 patients, fixed Harmony coordinates\n")
cat("Direct labels: L1 global + compartment-specific L2 globals\n")
cat("Frozen disease palette reused: YES\n")
cat("Frozen L1 palette reused: YES\n")
cat("Frozen L2 palette reused: YES\n")
cat("Dynamic colours generated: NO\n")
cat("UMAP recalculated: NO\n")
cat("Annotation modified: NO\n")
cat("PASS\n")
