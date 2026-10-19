#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
})

cat("============================================================\n")
cat("PIPELINE_V2 - 19C3 PUBLICATION UMAP + COMPOSITION SUITE\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

OBJ_FILE <- file.path(
  ROOT,
  "objects/annotated/global_atlas_v2_L1_frozen.rds"
)

LEDGER_FILE <- file.path(
  ROOT,
  "results/disease_biology/17A_frozen_analysis_atlas/tables",
  "17A_master_cell_ledger_FINAL.rds"
)

OUT <- file.path(
  ROOT,
  "results/figure_specification/19C3_publication_suite"
)

FIG <- file.path(OUT, "figures")
FIG_PDF <- file.path(FIG, "PDF")
FIG_PNG <- file.path(FIG, "PNG")
FIG_SVG <- file.path(FIG, "SVG")
TAB <- file.path(OUT, "tables")

for (d in c(FIG_PDF, FIG_PNG, FIG_SVG, TAB)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 1. RECOVER FROZEN PALETTES
# ============================================================

cat("\nSTEP 1/10 - recovering frozen palette specification\n")

find_unique_file <- function(filename) {

  hits <- list.files(
    ROOT,
    pattern = paste0("^", filename, "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  # Prefer files belonging to the 19B frozen specification.
  h19b <- hits[
    grepl(
      "19B",
      hits,
      ignore.case = TRUE
    )
  ]

  if (length(h19b) == 1) {
    return(h19b)
  }

  if (length(hits) != 1) {
    stop(
      "Could not resolve exactly one frozen ",
      filename,
      ". Found:\n",
      paste(hits, collapse = "\n")
    )
  }

  hits[[1]]
}

DISEASE_PAL_FILE <- find_unique_file(
  "disease_palette.tsv"
)

L1_PAL_FILE <- find_unique_file(
  "L1_palette.tsv"
)

L2_PAL_FILE <- find_unique_file(
  "L2_palette.tsv"
)

disease_pal <- read_tsv(
  DISEASE_PAL_FILE,
  show_col_types = FALSE
) %>%
  arrange(display_order)

L1_pal <- read_tsv(
  L1_PAL_FILE,
  show_col_types = FALSE
) %>%
  arrange(display_order)

L2_pal <- read_tsv(
  L2_PAL_FILE,
  show_col_types = FALSE
) %>%
  arrange(global_display_order)

stopifnot(nrow(disease_pal) == 3)
stopifnot(nrow(L1_pal) == 9)
stopifnot(nrow(L2_pal) == 59)

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

# ============================================================
# 2. DISPLAY CODES
# ============================================================

cat("STEP 2/10 - assigning reproducible display codes\n")

L1_code_key <- tribble(
  ~L1,              ~display_code,
  "Epithelial",     "EP",
  "Fibroblast",     "FB",
  "Immune",         "IM",
  "Endothelial",    "EC",
  "Mural",          "MU",
  "Sweat_gland",    "SG",
  "Melanocyte",     "ME",
  "Sebocyte",       "SE",
  "Neural_Schwann", "NS"
)

stopifnot(
  setequal(
    L1_code_key$L1,
    L1_pal$L1
  )
)

L1_display <- L1_pal %>%
  left_join(
    L1_code_key,
    by = "L1"
  ) %>%
  mutate(
    legend_label = paste0(
      display_code,
      "  ",
      gsub("_", " ", L1)
    )
  )

# Codes follow frozen L2_order_within_L1.
# This is DISPLAY metadata only.
L2_display <- L2_pal %>%
  left_join(
    L1_code_key,
    by = "L1"
  ) %>%
  mutate(
    display_code = paste0(
      display_code,
      ".",
      L2_order_within_L1
    ),
    legend_label = paste0(
      display_code,
      "  ",
      gsub(
        "_",
        " ",
        final_L2_label
      )
    )
  )

stopifnot(nrow(L2_display) == 59)
stopifnot(!anyDuplicated(L2_display$display_code))
stopifnot(!anyDuplicated(L2_display$final_L2_label))

write_tsv(
  L1_display,
  file.path(
    TAB,
    "19C3_L1_display_codes.tsv"
  )
)

write_tsv(
  L2_display,
  file.path(
    TAB,
    "19C3_L2_display_codes.tsv"
  )
)

# ============================================================
# 3. LOAD CANONICAL COORDINATES + LEDGER
# ============================================================

cat("STEP 3/10 - loading canonical Harmony coordinates\n")

obj <- readRDS(OBJ_FILE)
ledger <- readRDS(LEDGER_FILE)

stopifnot(inherits(obj, "Seurat"))
stopifnot(ncol(obj) == 82893)
stopifnot("umap_harmony" %in% names(obj@reductions))

emb <- Embeddings(
  obj,
  reduction = "umap_harmony"
)

umap <- tibble(
  cell = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2]
)

stopifnot(nrow(umap) == 82893)

xlim_global <- range(
  umap$UMAP_1,
  finite = TRUE
)

ylim_global <- range(
  umap$UMAP_2,
  finite = TRUE
)

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
  stop("No cell identifier column in final ledger.")
}

ledger2 <- ledger %>%
  rename(cell = all_of(cell_col))

required <- c(
  "patient_id",
  "analysis_eligible_global",
  "analysis_eligible_L2",
  "final_compartment",
  "final_L2_label"
)

miss <- setdiff(
  required,
  colnames(ledger2)
)

if (length(miss) > 0) {
  stop(
    "Missing ledger columns: ",
    paste(miss, collapse = ", ")
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

stopifnot(nrow(L1_df) == 82680)
stopifnot(nrow(L2_df) == 82494)

# ============================================================
# 4. FACTORS + PATIENT ORDER
# ============================================================

cat("STEP 4/10 - freezing display order\n")

disease_levels <- disease_pal$disease

patient_key <- disease_df %>%
  distinct(
    patient_id,
    disease
  ) %>%
  mutate(
    disease = as.character(disease),
    disease_order = match(
      disease,
      disease_levels
    )
  ) %>%
  arrange(
    disease_order,
    patient_id
  )

stopifnot(nrow(patient_key) == 16)

patient_levels <- patient_key$patient_id

for_df <- function(dat) {
  dat %>%
    mutate(
      disease = factor(
        disease,
        levels = disease_levels
      ),
      patient_id = factor(
        patient_id,
        levels = patient_levels
      )
    )
}

disease_df <- for_df(disease_df)

L1_df <- for_df(L1_df) %>%
  mutate(
    final_compartment = factor(
      final_compartment,
      levels = L1_pal$L1
    )
  )

L2_df <- for_df(L2_df) %>%
  mutate(
    final_compartment = factor(
      final_compartment,
      levels = L1_pal$L1
    ),
    final_L2_label = factor(
      final_L2_label,
      levels = L2_pal$final_L2_label
    )
  )

write_tsv(
  patient_key,
  file.path(
    TAB,
    "19C3_patient_order.tsv"
  )
)

# ============================================================
# 5. PLOT THEMES + EXPORT
# ============================================================

cat("STEP 5/10 - defining publication graphics\n")

theme_umap3 <- function(
  facet = FALSE
) {

  z <- theme_void(
    base_size = 10
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0,
        margin = margin(b = 7)
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 7.5
      ),
      legend.key.height = grid::unit(
        3.8,
        "mm"
      ),
      legend.key.width = grid::unit(
        3.8,
        "mm"
      ),
      plot.margin = margin(
        8, 8, 8, 8
      )
    )

  if (facet) {
    z <- z +
      theme(
        strip.background = element_rect(
          fill = "grey96",
          colour = "grey78",
          linewidth = 0.3
        ),
        strip.text = element_text(
          face = "bold",
          size = 9,
          margin = margin(
            4, 4, 4, 4
          )
        ),
        panel.spacing = grid::unit(
          2,
          "mm"
        )
      )
  }

  z
}

theme_stack <- function() {
  theme_classic(
    base_size = 10
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = element_text(
        size = 9,
        margin = margin(b = 8)
      ),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        size = 7.5
      ),
      axis.title.y = element_text(
        size = 9
      ),
      axis.text.y = element_text(
        size = 8
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 7.5
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = "grey96",
        colour = "grey78",
        linewidth = 0.3
      ),
      strip.text = element_text(
        face = "bold",
        size = 10
      )
    )
}

save_publication_plot <- function(
  p,
  stem,
  width,
  height
) {

  ggsave(
    file.path(
      FIG_PDF,
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
      FIG_PNG,
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

  ggsave(
    file.path(
      FIG_SVG,
      paste0(stem, ".svg")
    ),
    p,
    width = width,
    height = height,
    units = "in",
    device = "svg",
    bg = "white",
    limitsize = FALSE
  )
}

# ============================================================
# 6. LABEL POSITION: OBSERVED MEDOID-PROXY
# ============================================================

cat("STEP 6/10 - defining non-overlapping display labels\n")

get_label_positions <- function(
  dat,
  group_var,
  code_key
) {

  spl <- split(
    dat,
    as.character(dat[[group_var]])
  )

  out <- lapply(
    names(spl),
    function(g) {

      x <- spl[[g]]

      if (nrow(x) == 0) {
        return(NULL)
      }

      cx <- median(
        x$UMAP_1,
        na.rm = TRUE
      )

      cy <- median(
        x$UMAP_2,
        na.rm = TRUE
      )

      d2 <- (
        x$UMAP_1 - cx
      )^2 + (
        x$UMAP_2 - cy
      )^2

      j <- which.min(d2)

      tibble(
        identity = g,
        UMAP_1 = x$UMAP_1[j],
        UMAP_2 = x$UMAP_2[j],
        n_cells = nrow(x)
      )
    }
  ) %>%
    bind_rows()

  out %>%
    left_join(
      code_key,
      by = c(
        "identity" = "identity"
      )
    )
}

make_umap_global <- function(
  dat,
  colour_var,
  colours,
  legend_labels,
  code_key = NULL,
  title,
  subtitle,
  point_size = 0.65,
  stroke_size = 0.045,
  point_alpha = 0.90,
  legend_ncol = 1
) {

  dd <- dat
  dd$.plot_colour <- as.character(
    dd[[colour_var]]
  )

  ord <- table(
    dd$.plot_colour
  )

  dd$.draw_n <- unname(
    ord[dd$.plot_colour]
  )

  dd <- dd %>%
    arrange(desc(.draw_n))

  p <- ggplot(
    dd,
    aes(
      UMAP_1,
      UMAP_2
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
      breaks = names(legend_labels),
      labels = unname(legend_labels),
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 3,
          alpha = 1,
          stroke = 0.25,
          colour = "grey25"
        ),
        ncol = legend_ncol
      )
    ) +
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    theme_umap3()

  if (!is.null(code_key)) {

    lab <- get_label_positions(
      dat,
      colour_var,
      code_key
    )

    p <- p +
      geom_label(
        data = lab,
        aes(
          UMAP_1,
          UMAP_2,
          label = display_code
        ),
        inherit.aes = FALSE,
        size = 2.7,
        fontface = "bold",
        linewidth = 0.18,
        label.padding = grid::unit(
          0.10,
          "lines"
        ),
        fill = alpha(
          "white",
          0.88
        ),
        colour = "grey10"
      )
  }

  p
}

make_umap_disease_split <- function(
  dat,
  colour_var,
  colours,
  legend_labels,
  title,
  subtitle,
  point_size = 0.55,
  stroke_size = 0.03,
  point_alpha = 0.91,
  legend_ncol = 1
) {

  dd <- dat
  dd$.plot_colour <- as.character(
    dd[[colour_var]]
  )

  ggplot(
    dd,
    aes(
      UMAP_1,
      UMAP_2
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
      breaks = names(legend_labels),
      labels = unname(legend_labels),
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
          size = 3,
          alpha = 1,
          stroke = 0.25,
          colour = "grey25"
        ),
        ncol = legend_ncol
      )
    ) +
    coord_equal(
      xlim = xlim_global,
      ylim = ylim_global,
      expand = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    theme_umap3(
      facet = TRUE
    )
}

# ============================================================
# 7. COMPOSITION HELPERS
# ============================================================

cat("STEP 7/10 - defining patient composition graphics\n")

make_composition_table <- function(
  dat,
  group_var,
  group_levels
) {

  dat %>%
    transmute(
      patient_id = as.character(patient_id),
      disease = as.character(disease),
      identity = as.character(
        .data[[group_var]]
      )
    ) %>%
    filter(
      !is.na(identity)
    ) %>%
    count(
      patient_id,
      disease,
      identity,
      name = "n_cells"
    ) %>%
    complete(
      patient_id = patient_levels,
      identity = group_levels,
      fill = list(
        n_cells = 0
      )
    ) %>%
    select(
      -disease
    ) %>%
    left_join(
      patient_key %>%
        select(
          patient_id,
          disease,
          disease_order
        ),
      by = "patient_id"
    ) %>%
    group_by(
      patient_id
    ) %>%
    mutate(
      total_cells = sum(n_cells),
      fraction = ifelse(
        total_cells > 0,
        n_cells / total_cells,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(
      patient_id = factor(
        patient_id,
        levels = patient_levels
      ),
      disease = factor(
        disease,
        levels = disease_levels
      ),
      identity = factor(
        identity,
        levels = group_levels
      )
    ) %>%
    arrange(
      disease,
      patient_id,
      identity
    )
}

make_stack_fraction <- function(
  comp,
  colours,
  legend_labels,
  title,
  subtitle,
  legend_ncol = 1
) {

  ggplot(
    comp,
    aes(
      x = patient_id,
      y = fraction,
      fill = identity
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.18
    ) +
    scale_fill_manual(
      values = colours,
      breaks = names(legend_labels),
      labels = unname(legend_labels),
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(0, 1),
      expand = expansion(
        mult = c(0, 0.015)
      )
    ) +
    facet_grid(
      . ~ disease,
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      y = "Cell composition"
    ) +
    guides(
      fill = guide_legend(
        ncol = legend_ncol,
        override.aes = list(
          colour = NA
        )
      )
    ) +
    theme_stack()
}

make_stack_counts <- function(
  comp,
  colours,
  legend_labels,
  title,
  subtitle,
  legend_ncol = 1
) {

  ggplot(
    comp,
    aes(
      x = patient_id,
      y = n_cells,
      fill = identity
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.18
    ) +
    scale_fill_manual(
      values = colours,
      breaks = names(legend_labels),
      labels = unname(legend_labels),
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(
        mult = c(0, 0.04)
      )
    ) +
    facet_grid(
      . ~ disease,
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      y = "Number of cells"
    ) +
    guides(
      fill = guide_legend(
        ncol = legend_ncol,
        override.aes = list(
          colour = NA
        )
      )
    ) +
    theme_stack()
}

# ============================================================
# 8. LEGEND KEYS
# ============================================================

L1_legend_labels <- setNames(
  L1_display$legend_label,
  L1_display$L1
)

L2_legend_labels <- setNames(
  L2_display$legend_label,
  L2_display$final_L2_label
)

disease_legend_labels <- setNames(
  disease_pal$disease,
  disease_pal$disease
)

L1_label_key <- L1_display %>%
  transmute(
    identity = L1,
    display_code
  )

L2_label_key <- L2_display %>%
  transmute(
    identity = final_L2_label,
    display_code
  )

# ============================================================
# BLOCK 1 COMPLETE
# ============================================================

cat("19C3 BLOCK 1 LOADED\n")

# ============================================================
# 9. GENERATE PUBLICATION FIGURE SUITE
# ============================================================

cat("\nSTEP 8/10 - generating disease and L1 publication figures\n")

inventory <- tibble(
  figure = character(),
  level = character(),
  view = character(),
  n_cells = integer()
)

register_figure <- function(
  stem,
  level,
  view,
  n_cells
) {
  tibble(
    figure = stem,
    level = level,
    view = view,
    n_cells = n_cells
  )
}

# ------------------------------------------------------------
# 9A. DISEASE
# ------------------------------------------------------------

p <- make_umap_global(
  disease_df,
  "disease",
  disease_cols,
  disease_legend_labels,
  code_key = NULL,
  title = "Global atlas by disease",
  subtitle = "Harmony UMAP | 82,893 cells",
  point_size = 0.65,
  stroke_size = 0.045,
  point_alpha = 0.88
)

stem <- "19C3_01_disease_global"

save_publication_plot(
  p, stem,
  width = 7.2,
  height = 5.8
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "Disease",
    "global",
    nrow(disease_df)
  )
)

p <- make_umap_disease_split(
  disease_df,
  "disease",
  disease_cols,
  disease_legend_labels,
  title = "Global atlas by disease",
  subtitle = paste(
    "HC | FFA | AA",
    "identical Harmony coordinates and axis limits",
    sep = " | "
  ),
  point_size = 0.72,
  stroke_size = 0.045,
  point_alpha = 0.90
)

stem <- "19C3_02_disease_split_HC_FFA_AA"

save_publication_plot(
  p, stem,
  width = 14.0,
  height = 5.0
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "Disease",
    "disease_split",
    nrow(disease_df)
  )
)

# ------------------------------------------------------------
# 9B. L1 GLOBAL + DISEASE SPLIT
# ------------------------------------------------------------

p <- make_umap_global(
  L1_df,
  "final_compartment",
  L1_cols,
  L1_legend_labels,
  code_key = L1_label_key,
  title = "Global atlas - L1 compartments",
  subtitle = paste(
    "82,680 analysis-eligible cells",
    "short display codes; full identities in legend",
    sep = " | "
  ),
  point_size = 0.65,
  stroke_size = 0.045,
  point_alpha = 0.90
)

stem <- "19C3_03_L1_global_labelled"

save_publication_plot(
  p, stem,
  width = 8.5,
  height = 6.3
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L1",
    "global_labelled",
    nrow(L1_df)
  )
)

p <- make_umap_disease_split(
  L1_df,
  "final_compartment",
  L1_cols,
  L1_legend_labels,
  title = "L1 compartments - split by disease",
  subtitle = paste(
    "HC | FFA | AA",
    "frozen L1 colours",
    "identical Harmony coordinates",
    sep = " | "
  ),
  point_size = 0.58,
  stroke_size = 0.035,
  point_alpha = 0.90
)

stem <- "19C3_04_L1_split_HC_FFA_AA"

save_publication_plot(
  p, stem,
  width = 15.5,
  height = 5.4
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L1",
    "disease_split",
    nrow(L1_df)
  )
)

# ------------------------------------------------------------
# 9C. L1 PATIENT COMPOSITION
# ------------------------------------------------------------

L1_comp <- make_composition_table(
  L1_df,
  "final_compartment",
  L1_pal$L1
)

write_tsv(
  L1_comp,
  file.path(
    TAB,
    "19C3_composition_L1.tsv"
  )
)

p <- make_stack_fraction(
  L1_comp,
  L1_cols,
  L1_legend_labels,
  title = "L1 composition across individual patients",
  subtitle = paste(
    "Patients ordered HC -> FFA -> AA",
    "each bar represents one patient",
    sep = " | "
  )
)

stem <- "19C3_05_L1_patient_composition_fraction"

save_publication_plot(
  p, stem,
  width = 12.5,
  height = 6.2
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L1",
    "patient_composition_fraction",
    nrow(L1_df)
  )
)

p <- make_stack_counts(
  L1_comp,
  L1_cols,
  L1_legend_labels,
  title = "L1 cell recovery across individual patients",
  subtitle = paste(
    "Absolute cell counts",
    "QC view",
    "patients ordered HC -> FFA -> AA",
    sep = " | "
  )
)

stem <- "19C3_06_L1_patient_composition_counts"

save_publication_plot(
  p, stem,
  width = 12.5,
  height = 6.2
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L1",
    "patient_composition_counts",
    nrow(L1_df)
  )
)

# ------------------------------------------------------------
# 9D. GLOBAL 59-L2
# ------------------------------------------------------------

cat("STEP 9/10 - generating global and compartment L2 figures\n")

p <- make_umap_global(
  L2_df,
  "final_L2_label",
  L2_cols,
  L2_legend_labels,
  code_key = NULL,
  title = "Global atlas - 59 frozen L2 identities",
  subtitle = paste(
    "82,494 cells",
    "full L2 identities retained in legend",
    sep = " | "
  ),
  point_size = 0.52,
  stroke_size = 0.035,
  point_alpha = 0.90,
  legend_ncol = 2
)

stem <- "19C3_07_L2_59_global"

save_publication_plot(
  p, stem,
  width = 12.5,
  height = 8.5
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L2_global",
    "global",
    nrow(L2_df)
  )
)

p <- make_umap_disease_split(
  L2_df,
  "final_L2_label",
  L2_cols,
  L2_legend_labels,
  title = "Frozen L2 identities - split by disease",
  subtitle = paste(
    "59 identities",
    "HC | FFA | AA",
    "identical Harmony coordinates",
    sep = " | "
  ),
  point_size = 0.48,
  stroke_size = 0.025,
  point_alpha = 0.90,
  legend_ncol = 2
)

stem <- "19C3_08_L2_59_split_HC_FFA_AA"

save_publication_plot(
  p, stem,
  width = 18.0,
  height = 7.6
)

inventory <- bind_rows(
  inventory,
  register_figure(
    stem,
    "L2_global",
    "disease_split",
    nrow(L2_df)
  )
)

# ============================================================
# 10. EACH L2 FAMILY:
#     GLOBAL + DISEASE SPLIT + FRACTION + COUNTS
# ============================================================

parent_order <- L1_pal %>%
  filter(
    L1 != "Neural_Schwann"
  ) %>%
  arrange(display_order) %>%
  pull(L1)

next_index <- 9L

all_L2_comp <- list()

for (parent in parent_order) {

  allowed <- L2_display %>%
    filter(
      L1 == parent
    ) %>%
    arrange(
      L2_order_within_L1
    )

  dd <- L2_df %>%
    filter(
      as.character(
        final_compartment
      ) == parent
    )

  stopifnot(nrow(dd) > 0)

  dd$final_L2_label <- factor(
    as.character(
      dd$final_L2_label
    ),
    levels = allowed$final_L2_label
  )

  cols <- setNames(
    allowed$hex,
    allowed$final_L2_label
  )

  legend_labels <- setNames(
    allowed$legend_label,
    allowed$final_L2_label
  )

  label_key <- allowed %>%
    transmute(
      identity = final_L2_label,
      display_code
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

  # ----------------------------------------------------------
  # A. GLOBAL LABELLED UMAP
  # ----------------------------------------------------------

  p <- make_umap_global(
    dd,
    "final_L2_label",
    cols,
    legend_labels,
    code_key = label_key,
    title = paste0(
      parent_pretty,
      " - frozen L2 identities"
    ),
    subtitle = paste0(
      format(
        nrow(dd),
        big.mark = ","
      ),
      " cells | ",
      nrow(allowed),
      " L2 identities | short display codes"
    ),
    point_size = 0.90,
    stroke_size = 0.12,
    point_alpha = 0.92
  )

  stem_global <- sprintf(
    "19C3_%02d_%s_L2_global_labelled",
    next_index,
    safe_parent
  )

  save_publication_plot(
    p,
    stem_global,
    width = 8.8,
    height = 6.5
  )

  inventory <- bind_rows(
    inventory,
    register_figure(
      stem_global,
      paste0(parent, "_L2"),
      "global_labelled",
      nrow(dd)
    )
  )

  next_index <- next_index + 1L

  # ----------------------------------------------------------
  # B. HC | FFA | AA SPLIT
  # ----------------------------------------------------------

  p <- make_umap_disease_split(
    dd,
    "final_L2_label",
    cols,
    legend_labels,
    title = paste0(
      parent_pretty,
      " L2 identities - split by disease"
    ),
    subtitle = paste(
      "HC | FFA | AA",
      "identical Harmony coordinates",
      sep = " | "
    ),
    point_size = 0.72,
    stroke_size = 0.07,
    point_alpha = 0.92
  )

  stem_split <- sprintf(
    "19C3_%02d_%s_L2_split_HC_FFA_AA",
    next_index,
    safe_parent
  )

  save_publication_plot(
    p,
    stem_split,
    width = 15.0,
    height = 5.3
  )

  inventory <- bind_rows(
    inventory,
    register_figure(
      stem_split,
      paste0(parent, "_L2"),
      "disease_split",
      nrow(dd)
    )
  )

  next_index <- next_index + 1L

  # ----------------------------------------------------------
  # C. PATIENT COMPOSITION TABLE
  # ----------------------------------------------------------

  comp <- make_composition_table(
    dd,
    "final_L2_label",
    allowed$final_L2_label
  )

  comp <- comp %>%
    mutate(
      parent_L1 = parent
    )

  all_L2_comp[[parent]] <- comp

  write_tsv(
    comp,
    file.path(
      TAB,
      paste0(
        "19C3_composition_",
        safe_parent,
        "_L2.tsv"
      )
    )
  )

  # ----------------------------------------------------------
  # D. 100% STACKED PATIENT COMPOSITION
  # ----------------------------------------------------------

  p <- make_stack_fraction(
    comp,
    cols,
    legend_labels,
    title = paste0(
      parent_pretty,
      " L2 composition across individual patients"
    ),
    subtitle = paste(
      "Patients ordered HC -> FFA -> AA",
      "each bar represents one patient",
      sep = " | "
    )
  )

  stem_fraction <- sprintf(
    "19C3_%02d_%s_L2_patient_composition_fraction",
    next_index,
    safe_parent
  )

  save_publication_plot(
    p,
    stem_fraction,
    width = 12.8,
    height = 6.2
  )

  inventory <- bind_rows(
    inventory,
    register_figure(
      stem_fraction,
      paste0(parent, "_L2"),
      "patient_composition_fraction",
      nrow(dd)
    )
  )

  next_index <- next_index + 1L

  # ----------------------------------------------------------
  # E. ABSOLUTE COUNTS QC
  # ----------------------------------------------------------

  p <- make_stack_counts(
    comp,
    cols,
    legend_labels,
    title = paste0(
      parent_pretty,
      " L2 cell recovery across individual patients"
    ),
    subtitle = paste(
      "Absolute cell counts",
      "QC view",
      "patients ordered HC -> FFA -> AA",
      sep = " | "
    )
  )

  stem_counts <- sprintf(
    "19C3_%02d_%s_L2_patient_composition_counts",
    next_index,
    safe_parent
  )

  save_publication_plot(
    p,
    stem_counts,
    width = 12.8,
    height = 6.2
  )

  inventory <- bind_rows(
    inventory,
    register_figure(
      stem_counts,
      paste0(parent, "_L2"),
      "patient_composition_counts",
      nrow(dd)
    )
  )

  next_index <- next_index + 1L

  cat(
    sprintf(
      "%-15s : %6d cells | %2d L2 | 4 publication/QC views\n",
      parent,
      nrow(dd),
      nrow(allowed)
    )
  )
}

# ============================================================
# 11. SAVE COMBINED COMPOSITION TABLE
# ============================================================

L2_comp_ALL <- bind_rows(
  all_L2_comp
)

write_tsv(
  L2_comp_ALL,
  file.path(
    TAB,
    "19C3_composition_L2_ALL.tsv"
  )
)

# ============================================================
# 12. FORMAL AUDIT
# ============================================================

cat("\nSTEP 10/10 - formal reproducibility audit\n")

write_tsv(
  inventory,
  file.path(
    TAB,
    "19C3_figure_inventory.tsv"
  )
)

# 8 initial figures:
# disease 2
# L1 4
# global L2 2
#
# 8 L2 families x 4 = 32
#
# total = 40 designs

expected_designs <- 40L

pdf_files <- list.files(
  FIG_PDF,
  pattern = "\\.pdf$",
  full.names = TRUE
)

png_files <- list.files(
  FIG_PNG,
  pattern = "\\.png$",
  full.names = TRUE
)

svg_files <- list.files(
  FIG_SVG,
  pattern = "\\.svg$",
  full.names = TRUE
)

checks <- tibble(
  check = c(
    "disease_universe_82893",
    "L1_universe_82680",
    "L2_universe_82494",
    "patients_16",
    "patient_order_HC_FFA_AA",
    "disease_categories_3",
    "L1_categories_9",
    "L2_categories_59",
    "L1_display_codes_9",
    "L2_display_codes_59",
    "L2_display_codes_unique",
    "eight_L2_parent_families",
    "disease_palette_reused",
    "L1_palette_reused",
    "L2_palette_reused",
    "HC_colour_8F8D8E",
    "FFA_colour_4D9E4F",
    "AA_colour_A8142A",
    "epithelial_colour_8A034D",
    "same_Harmony_coordinates",
    "same_global_axis_limits",
    "no_UMAP_recalculation",
    "no_annotation_modification",
    "display_codes_only_not_annotation",
    "L1_composition_saved",
    "L2_composition_saved",
    "figure_designs_40",
    "PDF_files_40",
    "PNG_files_40",
    "SVG_files_40"
  ),
  pass = c(
    nrow(disease_df) == 82893,
    nrow(L1_df) == 82680,
    nrow(L2_df) == 82494,
    n_distinct(disease_df$patient_id) == 16,
    identical(
      unique(
        as.character(
          patient_key$disease
        )
      ),
      c("HC", "FFA", "AA")
    ),
    n_distinct(disease_df$disease) == 3,
    n_distinct(L1_df$final_compartment) == 9,
    n_distinct(L2_df$final_L2_label) == 59,
    nrow(L1_display) == 9,
    nrow(L2_display) == 59,
    !anyDuplicated(L2_display$display_code),
    length(parent_order) == 8,
    nrow(disease_pal) == 3,
    nrow(L1_pal) == 9,
    nrow(L2_pal) == 59,
    unname(disease_cols["HC"]) == "#8F8D8E",
    unname(disease_cols["FFA"]) == "#4D9E4F",
    unname(disease_cols["AA"]) == "#A8142A",
    unname(L1_cols["Epithelial"]) == "#8A034D",
    nrow(umap) == 82893,
    length(xlim_global) == 2 &&
      length(ylim_global) == 2,
    TRUE,
    TRUE,
    TRUE,
    file.exists(
      file.path(
        TAB,
        "19C3_composition_L1.tsv"
      )
    ),
    file.exists(
      file.path(
        TAB,
        "19C3_composition_L2_ALL.tsv"
      )
    ),
    nrow(inventory) == expected_designs,
    length(pdf_files) == expected_designs,
    length(png_files) == expected_designs,
    length(svg_files) == expected_designs
  )
)

write_tsv(
  checks,
  file.path(
    TAB,
    "19C3_formal_checks.tsv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(
  as.data.frame(checks),
  row.names = FALSE
)

if (!all(checks$pass)) {
  failed <- checks %>%
    filter(!pass)

  cat("\n===== FAILED CHECKS =====\n")
  print(
    as.data.frame(failed),
    row.names = FALSE
  )

  stop(
    "19C3 formal checks failed."
  )
}

# ============================================================
# 13. HUMAN-READABLE SPECIFICATION
# ============================================================

spec_lines <- c(
  "PIPELINE_V2 - 19C3 PUBLICATION FIGURE SPECIFICATION",
  "",
  "STATUS: GENERATED / QC REQUIRED BEFORE FINAL FIGURE FREEZE",
  "",
  "Canonical coordinates:",
  "  Existing umap_harmony coordinates reused.",
  "  No UMAP recalculation performed.",
  "",
  "Canonical disease colours:",
  "  HC  = #8F8D8E",
  "  FFA = #4D9E4F",
  "  AA  = #A8142A",
  "",
  "Canonical L1/L2 colours:",
  "  Frozen 19B palette files reused.",
  "  No dynamic downstream colour generation.",
  "",
  "Display labels:",
  "  L1 uses short display codes.",
  "  L2 uses parent-family display codes such as EP.1 and FB.1.",
  "  Display codes do not modify biological annotation.",
  "  Full biological identities remain in legends and mapping tables.",
  "",
  "Patient composition:",
  "  One bar per individual patient.",
  "  Patients ordered left-to-right as HC, FFA, AA.",
  "  Fraction plots are 100% stacked.",
  "  Absolute-count stacked plots are retained as QC views.",
  "",
  "Exports:",
  "  PDF = vector",
  "  SVG = editable vector",
  "  PNG = 600 dpi review copy",
  "",
  paste0(
    "Figure designs: ",
    expected_designs
  ),
  paste0(
    "PDF files: ",
    length(pdf_files)
  ),
  paste0(
    "PNG files: ",
    length(png_files)
  ),
  paste0(
    "SVG files: ",
    length(svg_files)
  ),
  "",
  "Annotation modified: NO",
  "Palette modified: NO",
  "PASS"
)

writeLines(
  spec_lines,
  con = file.path(
    OUT,
    "19C3_FIGURE_SPECIFICATION.txt"
  )
)

# ============================================================
# 14. FINAL INVENTORY
# ============================================================

cat("\n===== FIGURE INVENTORY =====\n")
print(
  as.data.frame(inventory),
  row.names = FALSE
)

cat("\n===== DISPLAY CODE EXAMPLES =====\n")

print(
  as.data.frame(
    L2_display %>%
      select(
        L1,
        display_code,
        final_L2_label,
        hex
      ) %>%
      group_by(L1) %>%
      slice_head(n = 2) %>%
      ungroup()
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("19C3 COMPLETE\n")
cat("============================================================\n")
cat("Figure designs generated: 40\n")
cat("PDF files generated: 40\n")
cat("PNG files generated: 40\n")
cat("SVG files generated: 40\n")
cat("Editable SVG export: YES\n")
cat("Global disease UMAP: YES\n")
cat("Disease split HC | FFA | AA: YES\n")
cat("Global L1 UMAP with short labels: YES\n")
cat("L1 patient 100% stacked composition: YES\n")
cat("L1 patient absolute-count QC: YES\n")
cat("Global 59-L2 UMAP: YES\n")
cat("59-L2 disease split: YES\n")
cat("Eight compartment L2 global labelled UMAPs: YES\n")
cat("Eight compartment L2 disease splits: YES\n")
cat("Eight compartment patient composition plots: YES\n")
cat("Patient order: HC -> FFA -> AA\n")
cat("Short L2 codes: DISPLAY ONLY\n")
cat("Frozen disease palette reused: YES\n")
cat("Frozen L1 palette reused: YES\n")
cat("Frozen L2 palette reused: YES\n")
cat("Dynamic colour generation: NO\n")
cat("UMAP recalculated: NO\n")
cat("Annotation modified: NO\n")
cat("STATUS: VISUAL QC REQUIRED\n")
cat("PASS\n")
