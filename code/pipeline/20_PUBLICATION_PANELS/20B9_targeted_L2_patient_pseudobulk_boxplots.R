suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(stringr)
  library(tibble)
})

set.seed(20261008)

cat("============================================================\n")
cat("20B9 — TARGETED L2 PATIENT PSEUDOBULK BOXPLOTS\n")
cat("============================================================\n")

# ============================================================
# PATHS
# ============================================================

PB_ROOT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18B_pseudobulk"
)

META_FILE <- file.path(
  PB_ROOT,
  "tables",
  "18B_pseudobulk_metadata_ALL.csv"
)

DE_ROOT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18C_DE/full_results"
)

OUT <- paste0(
  "PIPELINE_V2/results/publication_panels/",
  "20B9_targeted_L2_patient_pseudobulk_boxplots"
)

SVG_DIR <- file.path(
  OUT,
  "figures",
  "SVG"
)

PDF_DIR <- file.path(
  OUT,
  "figures",
  "PDF"
)

PNG_DIR <- file.path(
  OUT,
  "figures",
  "PNG"
)

TAB_DIR <- file.path(
  OUT,
  "tables"
)

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

# ============================================================
# FROZEN PARAMETERS
# ============================================================

MIN_CELLS <- 20
FDR_DISPLAY <- 0.10

disease_order <- c(
  "HC",
  "FFA",
  "AA"
)

contrast_order <- c(
  "FFA_vs_HC",
  "AA_vs_HC",
  "AA_vs_FFA"
)

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

# ============================================================
# FROZEN 20B9 TARGETS
# ============================================================

targets <- tribble(
  ~panel, ~axis_group, ~module, ~L2, ~L2_file, ~gene,

  1L,
  "Follicular niche",
  "07_EPITHELIAL_L2",
  "Bulge_Isthmus",
  "Bulge_Isthmus",
  "CD200",

  2L,
  "Follicular niche",
  "07_EPITHELIAL_L2",
  "Bulge_Isthmus",
  "Bulge_Isthmus",
  "HLA-DRB1",

  3L,
  "Fibroblast remodeling",
  "09_FIBROBLAST_L2",
  "PTGDS/IGF1 fibroblast",
  "PTGDS_IGF1_fibroblast",
  "CTHRC1",

  4L,
  "Fibroblast remodeling",
  "09_FIBROBLAST_L2",
  "PTGDS/IGF1 fibroblast",
  "PTGDS_IGF1_fibroblast",
  "NOX4",

  5L,
  "Fibroblast remodeling",
  "09_FIBROBLAST_L2",
  "COL11A1/DPEP1 dermal-sheath",
  "COL11A1_DPEP1_dermal-sheath",
  "CXCL9",

  6L,
  "Immune",
  "08_IMMUNE_L2",
  "CD8 T",
  "CD8_T",
  "IFITM1",

  7L,
  "Immune",
  "08_IMMUNE_L2",
  "CD8 T",
  "CD8_T",
  "GZMB",

  8L,
  "Immune",
  "08_IMMUNE_L2",
  "Macrophage/monocyte",
  "Macrophage_monocyte",
  "CXCL9",

  9L,
  "Vascular",
  "10_VASCULAR_L2",
  "ACKR1+ venular EC",
  "ACKR1_venular_EC",
  "ABCA1",

  10L,
  "Vascular",
  "10_VASCULAR_L2",
  "ACKR1+ venular EC",
  "ACKR1_venular_EC",
  "CD36"
)

# ============================================================
# MODULE -> 18B OBJECT
# ============================================================

module_objects <- c(
  "07_EPITHELIAL_L2" =
    "18B_07_EPITHELIAL_L2_pseudobulk_counts.rds",

  "08_IMMUNE_L2" =
    "18B_08_IMMUNE_L2_pseudobulk_counts.rds",

  "09_FIBROBLAST_L2" =
    "18B_09_FIBROBLAST_L2_pseudobulk_counts.rds",

  "10_VASCULAR_L2" =
    "18B_10_VASCULAR_L2_pseudobulk_counts.rds"
)

# ============================================================
# INPUT METADATA
# ============================================================

stopifnot(
  file.exists(META_FILE)
)

pb_meta <- read_csv(
  META_FILE,
  show_col_types=FALSE
)

# ============================================================
# FDR SYMBOL
# ============================================================

fdr_symbol <- function(x) {

  case_when(
    is.na(x)  ~ "",
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.10  ~ "\u2020",
    TRUE      ~ ""
  )
}

# ============================================================
# LOAD ALL THREE FROZEN 18C RESULTS FOR TARGET
# ============================================================

load_target_DE <- function(
  module,
  L2_file,
  gene
) {

  out <- list()

  for (cc in contrast_order) {

    f <- file.path(
      DE_ROOT,
      paste0(
        "18C_",
        cc,
        "__",
        module,
        "__",
        L2_file,
        ".csv"
      )
    )

    if (!file.exists(f)) {
      next
    }

    z <- read_csv(
      f,
      show_col_types=FALSE
    ) %>%
      filter(
        .data$gene == .env$gene
      )

    if (nrow(z) == 0) {
      next
    }

    out[[cc]] <- z
  }

  bind_rows(out)
}

# ============================================================
# STORAGE
# ============================================================

plot_list <- list()
expression_list <- list()
stats_list <- list()
support_list <- list()

# Cache module matrices
matrix_cache <- list()

# ============================================================
# MAIN TARGET LOOP
# ============================================================

for (ii in seq_len(nrow(targets))) {

  tt <- targets[ii,]

  cat("\n------------------------------------------------------------\n")
  cat(
    sprintf(
      "[%02d/%02d] %s | %s | %s\n",
      ii,
      nrow(targets),
      tt$axis_group,
      tt$L2,
      tt$gene
    )
  )
  cat("------------------------------------------------------------\n")

  module <- tt$module
  l2 <- tt$L2
  gene_i <- tt$gene

  # ----------------------------------------------------------
  # Load module matrix once
  # ----------------------------------------------------------

  if (is.null(
    matrix_cache[[module]]
  )) {

    obj_file <- file.path(
      PB_ROOT,
      "objects",
      module_objects[[module]]
    )

    stopifnot(
      file.exists(obj_file)
    )

    matrix_cache[[module]] <- readRDS(
      obj_file
    )
  }

  mat <- matrix_cache[[module]]

  # ----------------------------------------------------------
  # Align matrix metadata
  # ----------------------------------------------------------

  md_module <- pb_meta %>%
    filter(
      .data$module == .env$module
    )

  idx <- match(
    colnames(mat),
    md_module$pseudobulk_id
  )

  if (anyNA(idx)) {
    stop(
      "Metadata mismatch: ",
      module
    )
  }

  md_module <- md_module[
    idx,
    ,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(mat),
      md_module$pseudobulk_id
    )
  )

  # ----------------------------------------------------------
  # Exact 18C display eligibility:
  # target L2 and n_cells >=20
  # ----------------------------------------------------------

  keep_cols <- which(
    md_module$final_L2_label == l2 &
    md_module$n_cells >= MIN_CELLS &
    md_module$disease %in% disease_order
  )

  if (length(keep_cols) == 0) {
    stop(
      "No eligible pseudobulks: ",
      l2
    )
  }

  md <- md_module[
    keep_cols,
    ,
    drop=FALSE
  ]

  counts <- mat[
    ,
    keep_cols,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(counts),
      md$pseudobulk_id
    )
  )

  # one L2 pseudobulk per patient expected
  stopifnot(
    !anyDuplicated(
      md$patient_id
    )
  )

  # ----------------------------------------------------------
  # Patient ordering
  # ----------------------------------------------------------

  md <- md %>%
    mutate(
      disease=factor(
        disease,
        levels=disease_order
      )
    ) %>%
    arrange(
      disease,
      patient_id
    )

  counts <- counts[
    ,
    match(
      md$pseudobulk_id,
      colnames(counts)
    ),
    drop=FALSE
  ]

  # ----------------------------------------------------------
  # Common three-group TMM normalization for DISPLAY ONLY
  # No new inferential testing.
  # ----------------------------------------------------------

  y_display <- DGEList(
    counts=counts
  )

  y_display <- calcNormFactors(
    y_display,
    method="TMM"
  )

  logCPM_display <- cpm(
    y_display,
    log=TRUE,
    prior.count=2
  )

  # ----------------------------------------------------------
  # Gene presence
  # ----------------------------------------------------------

  if (!gene_i %in% rownames(
    logCPM_display
  )) {

    stop(
      "Gene absent from matrix: ",
      gene_i,
      " | ",
      l2
    )
  }

  # ----------------------------------------------------------
  # Expression table
  # ----------------------------------------------------------

  dd <- tibble(
    patient_id=md$patient_id,
    disease=as.character(
      md$disease
    ),
    batch=md$batch,
    n_cells=md$n_cells,
    module=module,
    L2=l2,
    gene=gene_i,
    logCPM=as.numeric(
      logCPM_display[
        gene_i,
        md$pseudobulk_id
      ]
    )
  ) %>%
    mutate(
      disease=factor(
        disease,
        levels=disease_order
      ),
      panel=tt$panel,
      axis_group=tt$axis_group
    )

  expression_list[[ii]] <- dd

  # ----------------------------------------------------------
  # Frozen 18C statistics
  # ----------------------------------------------------------

  ss_all <- load_target_DE(
    module=module,
    L2_file=tt$L2_file,
    gene=gene_i
  )

  if (nrow(ss_all) == 0) {
    stop(
      "No 18C statistics found: ",
      l2,
      " | ",
      gene_i
    )
  }

  ss_all <- ss_all %>%
    mutate(
      panel=tt$panel,
      axis_group=tt$axis_group,
      display_significant=
        FDR_BH < FDR_DISPLAY,
      symbol=fdr_symbol(
        FDR_BH
      )
    )

  stats_list[[ii]] <- ss_all

  ss <- ss_all %>%
    filter(
      FDR_BH < FDR_DISPLAY
    )

  # ----------------------------------------------------------
  # Support table
  # ----------------------------------------------------------

  support <- dd %>%
    mutate(
      disease=as.character(
        disease
      )
    ) %>%
    group_by(
      module,
      L2,
      gene,
      disease
    ) %>%
    summarise(
      n_patients=n_distinct(
        patient_id
      ),
      total_cells=sum(
        n_cells
      ),
      median_cells=median(
        n_cells
      ),
      min_cells=min(
        n_cells
      ),
      max_cells=max(
        n_cells
      ),
      .groups="drop"
    ) %>%
    mutate(
      panel=tt$panel,
      axis_group=tt$axis_group
    )

  support_list[[ii]] <- support

  # ----------------------------------------------------------
  # Bracket map
  # ----------------------------------------------------------

  bracket_map <- tribble(
    ~contrast,    ~x1, ~x2, ~level,
    "FFA_vs_HC",    1,   2,      1,
    "AA_vs_HC",     1,   3,      2,
    "AA_vs_FFA",    2,   3,      3
  )

  ann <- bracket_map %>%
    inner_join(
      ss %>%
        select(
          contrast,
          logFC,
          FDR_BH,
          symbol
        ),
      by="contrast"
    ) %>%
    arrange(
      level
    )

  # ----------------------------------------------------------
  # Y geometry
  # ----------------------------------------------------------

  yr <- range(
    dd$logCPM,
    na.rm=TRUE
  )

  span <- diff(yr)

  if (
    !is.finite(span) ||
    span == 0
  ) {
    span <- 1
  }

  ymax <- max(
    dd$logCPM,
    na.rm=TRUE
  )

  if (nrow(ann) > 0) {

    ann <- ann %>%
      mutate(
        display_level=row_number(),
        y=ymax +
          span * (
            0.12 +
            0.13 *
            (
              display_level - 1
            )
          )
      )
  }

  # ----------------------------------------------------------
  # Display labels
  # ----------------------------------------------------------

  l2_display <- l2 %>%
    str_replace_all(
      "_",
      " "
    )

  panel_title <- gene_i

  panel_subtitle <- l2_display

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot(
    dd,
    aes(
      x=disease,
      y=logCPM,
      fill=disease
    )
  ) +

    geom_boxplot(
      width=0.56,
      outlier.shape=NA,
      alpha=0.48,
      linewidth=0.55,
      colour="#333333"
    ) +

    geom_jitter(
      aes(
        colour=disease
      ),
      width=0.08,
      height=0,
      size=2.75,
      alpha=0.96,
      stroke=0.3,
      show.legend=FALSE
    ) +

    scale_fill_manual(
      values=DISEASE_COLORS,
      drop=FALSE
    ) +

    scale_colour_manual(
      values=DISEASE_COLORS,
      drop=FALSE
    ) +

    labs(
      title=panel_title,
      subtitle=panel_subtitle,
      x=NULL,
      y="TMM-normalized logCPM"
    ) +

    theme_classic(
      base_size=10.5
    ) +

    theme(
      legend.position="none",

      plot.title=element_text(
        face="bold",
        size=11.5,
        hjust=0.5
      ),

      plot.subtitle=element_text(
        size=8.5,
        hjust=0.5,
        colour="#4A4A4A"
      ),

      axis.text.x=element_text(
        face="bold",
        size=9
      ),

      axis.text.y=element_text(
        size=8
      ),

      axis.title.y=element_text(
        size=8.5
      ),

      axis.line=element_line(
        linewidth=0.45,
        colour="#333333"
      ),

      plot.margin=margin(
        7,
        9,
        7,
        9
      )
    )

  # ----------------------------------------------------------
  # Significant frozen 18C brackets only
  # ----------------------------------------------------------

  if (nrow(ann) > 0) {

    for (jj in seq_len(
      nrow(ann)
    )) {

      a <- ann[jj,]

      p <- p +

        annotate(
          "segment",
          x=a$x1,
          xend=a$x2,
          y=a$y,
          yend=a$y,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "segment",
          x=a$x1,
          xend=a$x1,
          y=a$y,
          yend=a$y -
            span * 0.025,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "segment",
         
          "segment",
          x=a$x2,
          xend=a$x2,
          y=a$y,
          yend=a$y -
            span * 0.025,
          linewidth=0.42,
          colour="#333333"
        ) +

        annotate(
          "text",
          x=(a$x1 + a$x2)/2,
          y=a$y + span * 0.025,
          label=a$symbol,
          size=3.7,
          fontface="bold",
          colour="#222222"
        )
    }
  }

  # ----------------------------------------------------------
  # Y-axis limits
  # ----------------------------------------------------------

  if (nrow(ann) > 0) {

    upper <- max(
      ann$y,
      na.rm=TRUE
    ) + span * 0.13

  } else {

    upper <- ymax +
      span * 0.10
  }

  lower <- min(
    dd$logCPM,
    na.rm=TRUE
  ) - span * 0.08

  p <- p +
    coord_cartesian(
      ylim=c(
        lower,
        upper
      ),
      clip="off"
    )

  plot_list[[ii]] <- p
}

# ============================================================
# CONSOLIDATE TABLES
# ============================================================

expression_tbl <- bind_rows(
  expression_list
)

stats_tbl <- bind_rows(
  stats_list
)

support_tbl <- bind_rows(
  support_list
)

# ============================================================
# RESTORE TARGET ORDER
# ============================================================

expression_tbl <- expression_tbl %>%
  mutate(
    panel=as.integer(
      panel
    )
  ) %>%
  arrange(
    panel,
    disease,
    patient_id
  )

stats_tbl <- stats_tbl %>%
  mutate(
    panel=as.integer(
      panel
    )
  ) %>%
  arrange(
    panel,
    match(
      contrast,
      contrast_order
    )
  )

support_tbl <- support_tbl %>%
  mutate(
    panel=as.integer(
      panel
    )
  ) %>%
  arrange(
    panel,
    match(
      disease,
      disease_order
    )
  )

# ============================================================
# SAVE AUDIT TABLES
# ============================================================

write_csv(
  targets,
  file.path(
    TAB_DIR,
    "20B9_frozen_target_panel.csv"
  )
)

write_csv(
  expression_tbl,
  file.path(
    TAB_DIR,
    "20B9_patient_L2_logCPM.csv"
  )
)

write_csv(
  stats_tbl,
  file.path(
    TAB_DIR,
    "20B9_frozen_18C_statistics.csv"
  )
)

write_csv(
  support_tbl,
  file.path(
    TAB_DIR,
    "20B9_patient_support.csv"
  )
)

# ============================================================
# STATISTICAL DISPLAY LEDGER
# ============================================================

display_ledger <- stats_tbl %>%
  select(
    panel,
    axis_group,
    module,
    final_L2_label,
    gene,
    contrast,
    numerator,
    denominator,
    n_numerator,
    n_denominator,
    batch_adjusted,
    robust_ge50,
    logFC,
    logCPM,
    PValue,
    FDR_BH,
    display_significant,
    symbol
  ) %>%
  arrange(
    panel,
    match(
      contrast,
      contrast_order
    )
  )

write_csv(
  display_ledger,
  file.path(
    TAB_DIR,
    "20B9_statistical_display_ledger.csv"
  )
)

# ============================================================
# SUPPORT QC
# ============================================================

support_wide <- support_tbl %>%
  select(
    panel,
    axis_group,
    module,
    L2,
    gene,
    disease,
    n_patients,
    total_cells,
    median_cells
  ) %>%
  pivot_wider(
    names_from=disease,
    values_from=c(
      n_patients,
      total_cells,
      median_cells
    )
  ) %>%
  mutate(
    min_patients_across_groups=pmin(
      n_patients_HC,
      n_patients_FFA,
      n_patients_AA
    ),

    all_groups_ge3=
      min_patients_across_groups >= 3,

    all_groups_ge4=
      min_patients_across_groups >= 4
  ) %>%
  arrange(
    panel
  )

write_csv(
  support_wide,
  file.path(
    TAB_DIR,
    "20B9_patient_support_wide.csv"
  )
)

# ============================================================
# MAIN COMBINED FIGURE
#
# 10 panels:
# row 1: follicular niche + first fibroblast
# row 2: fibroblast + immune
# row 3: immune + vascular
#
# Fixed target order is retained.
# ============================================================

final <- wrap_plots(
  plot_list,
  ncol=4,
  byrow=TRUE
) +

  plot_annotation(
    title=
      "Cell-state-resolved transcriptional programs in AA and FFA",

    subtitle=
      paste0(
        "Patient-level L2 pseudobulk expression; ",
        "HC, FFA and AA shown using a common three-group TMM display normalization"
      ),

    caption=
      paste0(
        "Each point represents one eligible patient (>=20 cells in the indicated L2 state). ",
        "Brackets reproduce frozen 18C edgeR pseudobulk inference; ",
        "*** FDR<0.001, ** FDR<0.01, * FDR<0.05, ",
        "\u2020 FDR 0.05-0.10. ",
        "No statistical tests were recomputed from plotted values."
      ),

    theme=theme(
      plot.title=element_text(
        face="bold",
        size=16,
        hjust=0.5
      ),

      plot.subtitle=element_text(
        size=10,
        hjust=0.5
      ),

      plot.caption=element_text(
        size=8,
        hjust=1
      )
    )
  )

# ============================================================
# SAVE MAIN FIGURE
# ============================================================

MAIN_STEM <-
  "20B9_targeted_L2_patient_pseudobulk_boxplots"

ggsave(
  file.path(
    SVG_DIR,
    paste0(
      MAIN_STEM,
      ".svg"
    )
  ),
  final,
  width=14.2,
  height=10.0,
  units="in",
  device=svglite
)

ggsave(
  file.path(
    PDF_DIR,
    paste0(
      MAIN_STEM,
      ".pdf"
    )
  ),
  final,
  width=14.2,
  height=10.0,
  units="in"
)

ggsave(
  file.path(
    PNG_DIR,
    paste0(
      MAIN_STEM,
      ".png"
    )
  ),
  final,
  width=14.2,
  height=10.0,
  units="in",
  dpi=400
)

# ============================================================
# BIOLOGICAL-BLOCK FIGURES
#
# Also save separate figures for each biological block.
# These are useful for manuscript assembly without changing
# the master combined figure.
# ============================================================

block_order <- c(
  "Follicular niche",
  "Fibroblast remodeling",
  "Immune",
  "Vascular"
)

block_summary <- list()

for (bb in block_order) {

  idx <- which(
    targets$axis_group == bb
  )

  pp <- plot_list[
    idx
  ]

  if (length(pp) == 0) {
    next
  }

  ncol_block <- ifelse(
    length(pp) <= 2,
    length(pp),
    3
  )

  block_plot <- wrap_plots(
    pp,
    ncol=ncol_block
  ) +

    plot_annotation(
      title=bb,

      subtitle=
        "Patient-level L2 pseudobulk expression",

      caption=
        paste0(
          "Frozen 18C edgeR inference. ",
          "*** FDR<0.001; ** FDR<0.01; ",
          "* FDR<0.05; \u2020 FDR 0.05-0.10."
        ),

      theme=theme(
        plot.title=element_text(
          face="bold",
          size=14,
          hjust=0.5
        ),

        plot.subtitle=element_text(
          size=9,
          hjust=0.5
        ),

        plot.caption=element_text(
          size=7.5,
          hjust=1
        )
      )
    )

  stem_block <- paste0(
    "20B9_",
    str_replace_all(
      toupper(bb),
      "[^A-Z0-9]+",
      "_"
    ),
    "_L2_boxplots"
  )

  block_width <- ifelse(
    length(pp) == 1,
    4.3,
    ifelse(
      length(pp) == 2,
      7.4,
      10.5
    )
  )

  ggsave(
    file.path(
      SVG_DIR,
      paste0(
        stem_block,
        ".svg"
      )
    ),
    block_plot,
    width=block_width,
    height=4.8,
    units="in",
    device=svglite
  )

  ggsave(
    file.path(
      PDF_DIR,
      paste0(
        stem_block,
        ".pdf"
      )
    ),
    block_plot,
    width=block_width,
    height=4.8,
    units="in"
  )

  ggsave(
    file.path(
      PNG_DIR,
      paste0(
        stem_block,
        ".png"
      )
    ),
    block_plot,
    width=block_width,
    height=4.8,
    units="in",
    dpi=400
  )

  block_summary[[bb]] <- tibble(
    axis_group=bb,
    n_panels=length(pp),
    svg=file.path(
      SVG_DIR,
      paste0(
        stem_block,
        ".svg"
      )
    )
  )
}

block_summary_tbl <- bind_rows(
  block_summary
)

write_csv(
  block_summary_tbl,
  file.path(
    TAB_DIR,
    "20B9_biological_block_summary.csv"
  )
)

# ============================================================
# OUTPUT INVENTORY
# ============================================================

svg_files <- list.files(
  SVG_DIR,
  pattern="\\.svg$",
  full.names=TRUE
)

pdf_files <- list.files(
  PDF_DIR,
  pattern="\\.pdf$",
  full.names=TRUE
)

png_files <- list.files(
  PNG_DIR,
  pattern="\\.png$",
  full.names=TRUE
)

inventory_files <- c(
  svg_files,
  pdf_files,
  png_files
)

inventory <- tibble(
  file=inventory_files,
  bytes=file.info(
    inventory_files
  )$size,
  extension=tolower(
    tools::file_ext(
      inventory_files
    )
  )
)

write_csv(
  inventory,
  file.path(
    TAB_DIR,
    "20B9_output_inventory.csv"
  )
)

# ============================================================
# FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "ten_frozen_targets",
    "ten_plots_generated",
    "all_targets_have_expression",
    "all_display_patients_ge20_cells",
    "all_targets_min_three_patients_each_group",
    "four_biological_blocks",
    "main_SVG_exists",
    "main_PDF_exists",
    "main_PNG_exists",
    "five_SVG_total",
    "five_PDF_total",
    "five_PNG_total",
    "only_frozen_18C_statistics_used",
    "all_outputs_nonzero"
  ),

  pass=c(
    nrow(targets) == 10,

    length(plot_list) == 10,

    n_distinct(
      expression_tbl$panel
    ) == 10,

    all(
      expression_tbl$n_cells >=
        MIN_CELLS
    ),

    all(
      support_wide$all_groups_ge3
    ),

    n_distinct(
      targets$axis_group
    ) == 4,

    file.exists(
      file.path(
        SVG_DIR,
        paste0(
          MAIN_STEM,
          ".svg"
        )
      )
    ),

    file.exists(
      file.path(
        PDF_DIR,
        paste0(
          MAIN_STEM,
          ".pdf"
        )
      )
    ),

    file.exists(
      file.path(
        PNG_DIR,
        paste0(
          MAIN_STEM,
          ".png"
        )
      )
    ),

    length(svg_files) == 5,

    length(pdf_files) == 5,

    length(png_files) == 5,

    all(
      stats_tbl$contrast %in%
        contrast_order
    ),

    all(
      inventory$bytes > 0
    )
  )
)

write_csv(
  checks,
  file.path(
    TAB_DIR,
    "20B9_formal_checks.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B9 — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nTargets:",
  nrow(targets),
  "\n"
)

cat(
  "Biological blocks:",
  n_distinct(
    targets$axis_group
  ),
  "\n"
)

cat(
  "\nDisplay unit: patient x frozen L2 state\n"
)

cat(
  "Minimum cells per displayed patient-L2:",
  MIN_CELLS,
  "\n"
)

cat(
  "Disease order: HC -> FFA -> AA\n"
)

cat(
  "Common HC+FFA+AA TMM normalization for display: YES\n"
)

cat(
  "New differential-expression testing: NO\n"
)

cat(
  "Frozen 18C edgeR statistics used for brackets: YES\n"
)

cat(
  "FFA-vs-HC 18C/18D3-adjusted equivalence audited previously: YES\n"
)

cat(
  "\nFDR display threshold:",
  FDR_DISPLAY,
  "\n"
)

cat(
  "*** FDR<0.001\n"
)

cat(
  "**  FDR<0.01\n"
)

cat(
  "*   FDR<0.05\n"
)

cat(
  "\u2020   FDR 0.05-0.10\n"
)

cat(
  "\nFrozen target panel:\n"
)

print(
  targets %>%
    select(
      panel,
      axis_group,
      L2,
      gene
    ),
  n=Inf
)

cat(
  "\nPatient support:\n"
)

print(
  support_wide %>%
    select(
      panel,
      L2,
      gene,
      n_patients_HC,
      n_patients_FFA,
      n_patients_AA,
      all_groups_ge3,
      all_groups_ge4
    ),
  n=Inf
)

cat(
  "\nDisplayed significant comparisons:\n"
)

print(
  display_ledger %>%
    filter(
      display_significant
    ) %>%
    select(
      panel,
      final_L2_label,
      gene,
      contrast,
      logFC,
      FDR_BH,
      symbol
    ),
  n=Inf
)

cat(
  "\nOutputs:\n"
)

cat(
  "SVG:",
  length(svg_files),
  "\n"
)

cat(
  "PDF:",
  length(pdf_files),
  "\n"
)

cat(
  "PNG:",
  length(png_files),
  "\n"
)

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  checks,
  n=Inf
)

cat(
  "\nOriginal 18B pseudobulks modified: NO\n"
)

cat(
  "Frozen 18C DE modified: NO\n"
)

cat(
  "18D3 batch audit modified: NO\n"
)

cat(
  "Frozen L2 annotation modified: NO\n"
)

if (!all(
  checks$pass
)) {

  cat(
    "\nSTATUS=FAIL\n"
  )

  print(
    checks %>%
      filter(
        !pass
      )
  )

  stop(
    "20B9 formal checks failed."
  )
}

cat(
  "\nSTATUS=PASS\n"
)
