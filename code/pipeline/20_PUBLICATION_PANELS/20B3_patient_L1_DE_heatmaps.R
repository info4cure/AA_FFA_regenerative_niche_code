#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

cat("============================================================\n")
cat("20B3 — PATIENT-LEVEL L1 DIFFERENTIAL EXPRESSION HEATMAPS\n")
cat("PROTOTYPE: FIBROBLAST\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

COUNTS_FILE <- file.path(
  ROOT,
  "results/disease_biology/18_PSEUDOBULK_DE/18B_pseudobulk/objects",
  "18B_09_FIBROBLAST_L2_pseudobulk_counts.rds"
)

META_FILE <- file.path(
  ROOT,
  "results/disease_biology/18_PSEUDOBULK_DE/18B_pseudobulk/objects",
  "18B_pseudobulk_metadata_ALL.rds"
)

OUT <- file.path(
  ROOT,
  "results/publication_panels/20B3_patient_L1_DE_heatmaps",
  "Fibroblast"
)

FIG_SVG <- file.path(OUT, "figures/SVG")
FIG_PDF <- file.path(OUT, "figures/PDF")
FIG_PNG <- file.path(OUT, "figures/PNG")
DE_DIR  <- file.path(OUT, "DE")
SEL_DIR <- file.path(OUT, "selected_genes")
SRC_DIR <- file.path(OUT, "source_data")
TAB_DIR <- file.path(OUT, "tables")

for (d in c(
  FIG_SVG, FIG_PDF, FIG_PNG,
  DE_DIR, SEL_DIR, SRC_DIR, TAB_DIR
)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

stopifnot(
  file.exists(COUNTS_FILE),
  file.exists(META_FILE)
)

MODULE <- "09_FIBROBLAST_L2"
COMPARTMENT <- "Fibroblast"

DISEASE_COLS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

# ============================================================
# 2. LOAD L2 PSEUDOBULKS
# ============================================================

cat("\nSTEP 1/8 — loading frozen L2 pseudobulks\n")

mat <- readRDS(COUNTS_FILE)

meta <- readRDS(META_FILE) %>%
  filter(module == MODULE) %>%
  as.data.frame()

stopifnot(
  inherits(mat, "Matrix") || is.matrix(mat),
  nrow(mat) == 23412,
  setequal(colnames(mat), meta$pseudobulk_id)
)

meta <- meta[
  match(colnames(mat), meta$pseudobulk_id),
  ,
  drop=FALSE
]

stopifnot(
  identical(colnames(mat), meta$pseudobulk_id)
)

# ============================================================
# 3. COLLAPSE L2 -> PATIENT L1
# ============================================================

cat("STEP 2/8 — collapsing L2 pseudobulks to patient-level L1\n")

patient_info <- meta %>%
  group_by(patient_id) %>%
  summarise(
    disease = first(disease),
    batch = first(batch),
    n_cells = sum(n_cells),
    n_L2_present = n(),
    .groups="drop"
  )

stopifnot(
  nrow(patient_info) == 16,
  all(patient_info$n_cells >= 20)
)

disease_order <- c("HC","FFA","AA")

patient_info <- patient_info %>%
  mutate(
    disease_order = match(disease, disease_order)
  ) %>%
  arrange(
    disease_order,
    patient_id
  )

patient_levels <- patient_info$patient_id

patient_counts <- do.call(
  cbind,
  lapply(
    patient_levels,
    function(pid) {
      idx <- which(meta$patient_id == pid)
      Matrix::rowSums(
        mat[, idx, drop=FALSE]
      )
    }
  )
)

rownames(patient_counts) <- rownames(mat)
colnames(patient_counts) <- patient_levels

patient_counts <- as.matrix(patient_counts)

stopifnot(
  ncol(patient_counts) == 16,
  identical(
    colnames(patient_counts),
    patient_levels
  )
)

write.csv(
  patient_counts,
  file.path(
    SRC_DIR,
    "20B3_Fibroblast_patient_L1_raw_counts.csv"
  ),
  quote=FALSE
)

write_csv(
  patient_info,
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_patient_metadata.csv"
  )
)

# ============================================================
# 4. GLOBAL 16-PATIENT NORMALIZATION FOR DISPLAY
# ============================================================

cat("STEP 3/8 — TMM normalization for 16-patient display matrix\n")

y_all <- DGEList(
  counts=patient_counts
)

y_all <- calcNormFactors(
  y_all,
  method="TMM"
)

logCPM_all <- cpm(
  y_all,
  log=TRUE,
  prior.count=2
)

write.csv(
  logCPM_all,
  file.path(
    SRC_DIR,
    "20B3_Fibroblast_patient_logCPM_TMM.csv"
  ),
  quote=FALSE
)

# ============================================================
# 5. CONTRAST SPECIFICATION — FROZEN FROM 18C
# ============================================================

contrast_spec <- tribble(
  ~contrast, ~group1, ~group0, ~adjust_batch,
  "AA_vs_HC",  "AA",  "HC", FALSE,
  "FFA_vs_HC", "FFA", "HC", TRUE
)

# ============================================================
# 6. L1 edgeR QL
# ============================================================

cat("STEP 4/8 — L1 edgeR quasi-likelihood DE\n")

run_L1_DE <- function(
  contrast_name,
  group1,
  group0,
  adjust_batch
) {

  md <- patient_info %>%
    filter(
      disease %in% c(group0, group1),
      n_cells >= 20
    ) %>%
    mutate(
      disease_model=factor(
        disease,
        levels=c(group0,group1)
      )
    )

  n0 <- sum(md$disease_model == group0)
  n1 <- sum(md$disease_model == group1)

  stopifnot(
    n0 >= 3,
    n1 >= 3
  )

  counts <- patient_counts[
    ,
    md$patient_id,
    drop=FALSE
  ]

  if (adjust_batch) {

    md$batch_model <- factor(md$batch)

    design <- model.matrix(
      ~ batch_model + disease_model,
      data=md
    )

  } else {

    design <- model.matrix(
      ~ disease_model,
      data=md
    )
  }

  stopifnot(
    qr(design)$rank == ncol(design)
  )

  coef_name <- paste0(
    "disease_model",
    group1
  )

  coef_idx <- which(
    colnames(design) == coef_name
  )

  stopifnot(length(coef_idx) == 1)

  y <- DGEList(
    counts=counts
  )

  keep <- filterByExpr(
    y,
    design=design
  )

  stopifnot(sum(keep) > 0)

  y <- y[
    keep,
    ,
    keep.lib.sizes=FALSE
  ]

  y <- calcNormFactors(
    y,
    method="TMM"
  )

  y <- estimateDisp(
    y,
    design=design,
    robust=TRUE
  )

  fit <- glmQLFit(
    y,
    design=design,
    robust=TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef=coef_idx
  )

  tt <- topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table

  res <- tt %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    rename(FDR_BH=FDR) %>%
    mutate(
      contrast=contrast_name,
      numerator=group1,
      denominator=group0,
      compartment=COMPARTMENT,
      n_numerator=n1,
      n_denominator=n0,
      batch_adjusted=adjust_batch,
      n_genes_tested=sum(keep),
      .before=1
    )

  list(
    result=res,
    design=design,
    keep=keep
  )
}

DE <- list()

for (i in seq_len(nrow(contrast_spec))) {

  cc <- contrast_spec[i,]

  cat(
    "  ",
    cc$contrast,
    " | batch adjusted=",
    cc$adjust_batch,
    "\n",
    sep=""
  )

  z <- run_L1_DE(
    cc$contrast,
    cc$group1,
    cc$group0,
    cc$adjust_batch
  )

  DE[[cc$contrast]] <- z$result

  write_csv(
    z$result,
    file.path(
      DE_DIR,
      paste0(
        "20B3_Fibroblast_",
        cc$contrast,
        "_full_DE.csv"
      )
    )
  )
}

# ============================================================
# 7. SELECT TOP50 UP + TOP50 DOWN
# ============================================================

cat("STEP 5/8 — selecting top-ranked UP/DOWN genes\n")

select_top <- function(res) {

  up <- res %>%
    filter(logFC > 0) %>%
    arrange(PValue, desc(abs(logFC))) %>%
    slice_head(n=50) %>%
    mutate(
      direction="UP",
      direction_order=1,
      rank_direction=row_number()
    )

  down <- res %>%
    filter(logFC < 0) %>%
    arrange(PValue, desc(abs(logFC))) %>%
    slice_head(n=50) %>%
    mutate(
      direction="DOWN",
      direction_order=2,
      rank_direction=row_number()
    )

  bind_rows(up,down)
}

# ============================================================
# 8. HEATMAP
# ============================================================

make_heatmap <- function(
  contrast_name,
  res
) {

  selected <- select_top(res)

  stopifnot(
    nrow(selected) == 100,
    n_distinct(selected$gene) == 100
  )

  write_csv(
    selected,
    file.path(
      SEL_DIR,
      paste0(
        "20B3_Fibroblast_",
        contrast_name,
        "_top50_UP_DOWN.csv"
      )
    )
  )

  genes <- selected$gene

  stopifnot(
    all(genes %in% rownames(logCPM_all))
  )

  m <- logCPM_all[
    genes,
    patient_levels,
    drop=FALSE
  ]

  # Row-wise z-score.
  z <- t(
    scale(
      t(m),
      center=TRUE,
      scale=TRUE
    )
  )

  # Genes with zero variance would generate NA.
  z[!is.finite(z)] <- 0

  # Clip for graphical robustness.
  z <- pmax(
    pmin(z,2.5),
    -2.5
  )

  write.csv(
    z,
    file.path(
      SRC_DIR,
      paste0(
        "20B3_Fibroblast_",
        contrast_name,
        "_heatmap_zscore_matrix.csv"
      )
    ),
    quote=FALSE
  )

  row_key <- selected %>%
    transmute(
      gene,
      direction,
      rank_direction,
      row_order=row_number()
    )

  plot_df <- as.data.frame(z) %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      -gene,
      names_to="patient_id",
      values_to="z"
    ) %>%
    left_join(
      row_key,
      by="gene"
    ) %>%
    left_join(
      patient_info %>%
        select(
          patient_id,
          disease
        ),
      by="patient_id"
    ) %>%
    mutate(
      patient_id=factor(
        patient_id,
        levels=patient_levels
      ),
      # ggplot draws the first factor level at the bottom.
      # Reverse the selected-gene order so that the first
      # selected block (UP) is displayed at the TOP and
      # the second block (DOWN) at the BOTTOM.
      gene=factor(
        gene,
        levels=rev(genes)
      )
    )

  # Number of FDR-significant genes among displayed genes.
  n_sig_up <- selected %>%
    filter(
      direction=="UP",
      FDR_BH < 0.05
    ) %>%
    nrow()

  n_sig_down <- selected %>%
    filter(
      direction=="DOWN",
      FDR_BH < 0.05
    ) %>%
    nrow()

  title_txt <- gsub(
    "_",
    " ",
    contrast_name
  )

  subtitle_txt <- paste0(
    "Top-ranked 50 UP + 50 DOWN genes; ",
    "FDR<0.05 among displayed: ",
    n_sig_up,
    " UP / ",
    n_sig_down,
    " DOWN"
  )

  p <- ggplot(
    plot_df,
    aes(
      x=patient_id,
      y=gene,
      fill=z
    )
  ) +
    geom_tile() +

    # Group separators HC | FFA | AA.
    geom_vline(
      xintercept=c(5.5,10.5),
      linewidth=0.55,
      colour="white"
    ) +

    # UP / DOWN separator.
    geom_hline(
      yintercept=50.5,
      linewidth=0.7,
      colour="black"
    ) +

    scale_fill_gradient2(
      low="#2166AC",
      mid="#F7F7F7",
      high="#B2182B",
      midpoint=0,
      limits=c(-2.5,2.5),
      name="Row z-score"
    ) +

    labs(
      title=paste0(
        "Fibroblast — ",
        title_txt
      ),
      subtitle=subtitle_txt,
      x=NULL,
      y=NULL
    ) +

    theme_minimal(
      base_size=8
    ) +

    theme(
      panel.grid=element_blank(),

      axis.text.x=element_text(
        angle=45,
        hjust=1,
        vjust=1,
        size=7
      ),

      axis.text.y=element_text(
        size=5.5,
        colour="black"
      ),

      plot.title=element_text(
        face="bold",
        size=11
      ),

      plot.subtitle=element_text(
        size=7.5
      ),

      legend.position="right",

      plot.margin=margin(
        8,8,8,8
      )
    )

  # Disease annotation strip as a separate ggplot.
  ann <- patient_info %>%
    mutate(
      patient_id=factor(
        patient_id,
        levels=patient_levels
      ),
      disease=factor(
        disease,
        levels=c("HC","FFA","AA")
      ),
      y=1
    )

  p_ann <- ggplot(
    ann,
    aes(
      x=patient_id,
      y=y,
      fill=disease
    )
  ) +
    geom_tile() +
    scale_fill_manual(
      values=DISEASE_COLS,
      drop=FALSE
    ) +
    scale_x_discrete(
      expand=c(0,0)
    ) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=margin(0,0,0,0)
    )

  if (!requireNamespace(
    "patchwork",
    quietly=TRUE
  )) {
    stop("Package 'patchwork' required.")
  }

  final <- p_ann / p +
    patchwork::plot_layout(
      heights=c(0.035,0.965)
    )

  stem <- paste0(
    "20B3_Fibroblast_",
    contrast_name,
    "_top50_UP_DOWN_heatmap"
  )

  ggsave(
    file.path(
      FIG_SVG,
      paste0(stem,".svg")
    ),
    final,
    width=7.4,
    height=10.0,
    units="in",
    device=svglite::svglite
  )

  ggsave(
    file.path(
      FIG_PDF,
      paste0(stem,".pdf")
    ),
    final,
    width=7.4,
    height=10.0,
    units="in"
  )

  ggsave(
    file.path(
      FIG_PNG,
      paste0(stem,".png")
    ),
    final,
    width=7.4,
    height=10.0,
    units="in",
    dpi=400
  )

  tibble(
    contrast=contrast_name,
    n_displayed=nrow(selected),
    n_up=sum(selected$direction=="UP"),
    n_down=sum(selected$direction=="DOWN"),
    n_FDR005_up=n_sig_up,
    n_FDR005_down=n_sig_down
  )
}

cat("STEP 6/8 — generating patient heatmaps\n")

heatmap_summary <- bind_rows(
  lapply(
    names(DE),
    function(cc) {
      make_heatmap(
        cc,
        DE[[cc]]
      )
    }
  )
)

write_csv(
  heatmap_summary,
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_heatmap_summary.csv"
  )
)

# ============================================================
# FINAL AUDIT
# ============================================================

cat("STEP 7/8 —
final audit\n")

# ============================================================
# 9. FORMAL CHECKS
# ============================================================

svg_files <- list.files(
  FIG_SVG,
  pattern="\\.svg$",
  full.names=TRUE
)

pdf_files <- list.files(
  FIG_PDF,
  pattern="\\.pdf$",
  full.names=TRUE
)

png_files <- list.files(
  FIG_PNG,
  pattern="\\.png$",
  full.names=TRUE
)

de_files <- list.files(
  DE_DIR,
  pattern="_full_DE\\.csv$",
  full.names=TRUE
)

selected_files <- list.files(
  SEL_DIR,
  pattern="_top50_UP_DOWN\\.csv$",
  full.names=TRUE
)

matrix_files <- list.files(
  SRC_DIR,
  pattern="_heatmap_zscore_matrix\\.csv$",
  full.names=TRUE
)

checks <- tibble(
  check = c(
    "sixteen_patient_L1_pseudobulks",
    "HC_n5",
    "FFA_n5",
    "AA_n6",
    "all_patients_ge20_cells",
    "two_DE_contrasts",
    "FFA_HC_batch_adjusted",
    "AA_HC_not_batch_adjusted",
    "two_top50_gene_tables",
    "100_genes_each_heatmap",
    "50_UP_each_heatmap",
    "50_DOWN_each_heatmap",
    "two_zscore_matrices",
    "two_SVG_exports",
    "two_PDF_exports",
    "two_PNG_exports",
    "patient_order_HC_FFA_AA"
  ),
  pass = c(
    ncol(patient_counts) == 16,
    sum(patient_info$disease == "HC") == 5,
    sum(patient_info$disease == "FFA") == 5,
    sum(patient_info$disease == "AA") == 6,
    all(patient_info$n_cells >= 20),
    length(DE) == 2,
    all(DE$FFA_vs_HC$batch_adjusted),
    all(!DE$AA_vs_HC$batch_adjusted),
    length(selected_files) == 2,
    all(heatmap_summary$n_displayed == 100),
    all(heatmap_summary$n_up == 50),
    all(heatmap_summary$n_down == 50),
    length(matrix_files) == 2,
    length(svg_files) == 2,
    length(pdf_files) == 2,
    length(png_files) == 2,
    identical(
      as.character(patient_info$disease),
      c(
        rep("HC",5),
        rep("FFA",5),
        rep("AA",6)
      )
    )
  )
)

write_csv(
  checks,
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_formal_checks.csv"
  )
)

print(checks)

if (!all(checks$pass)) {
  stop(
    "20B3 formal checks failed: ",
    paste(
      checks$check[!checks$pass],
      collapse=", "
    )
  )
}

# ============================================================
# 10. ANALYSIS / DISPLAY SPECIFICATION
# ============================================================

specification <- tibble(
  item = c(
    "aggregation_level",
    "minimum_cells_per_patient_L1",
    "FFA_vs_HC_model",
    "AA_vs_HC_model",
    "gene_filter",
    "normalization",
    "dispersion",
    "DE_framework",
    "multiple_testing",
    "heatmap_gene_selection",
    "genes_per_direction",
    "heatmap_expression",
    "heatmap_scaling",
    "heatmap_clipping",
    "column_clustering",
    "column_order",
    "third_disease_gene_selection_role"
  ),
  value = c(
    "patient x L1 compartment; L2 pseudobulks summed within patient",
    "20",
    "~ batch_model + disease_model",
    "~ disease_model; batch not identifiable because of confounding",
    "edgeR filterByExpr within each L1 contrast",
    "TMM",
    "edgeR estimateDisp robust=TRUE",
    "edgeR quasi-likelihood; glmQLFit robust=TRUE; glmQLFTest",
    "Benjamini-Hochberg within each L1 contrast",
    "top-ranked genes by PValue separately within positive and negative logFC",
    "50 UP + 50 DOWN",
    "TMM-normalized logCPM across all 16 patient-level L1 pseudobulks",
    "row z-score across patients",
    "[-2.5, 2.5]",
    "NO",
    "HC then FFA then AA; patient ID within disease",
    "display only; does not participate in gene selection for the focal contrast"
  )
)

write_csv(
  specification,
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_analysis_specification.csv"
  )
)

# ============================================================
# 11. HASH OUTPUTS
# ============================================================

cat("STEP 8/8 — writing audit and hashes\n")

all_output_files <- c(
  svg_files,
  pdf_files,
  png_files,
  de_files,
  selected_files,
  matrix_files,
  file.path(
    SRC_DIR,
    "20B3_Fibroblast_patient_L1_raw_counts.csv"
  ),
  file.path(
    SRC_DIR,
    "20B3_Fibroblast_patient_logCPM_TMM.csv"
  ),
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_patient_metadata.csv"
  ),
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_heatmap_summary.csv"
  ),
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_formal_checks.csv"
  ),
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_analysis_specification.csv"
  )
)

all_output_files <- unique(
  all_output_files[
    file.exists(all_output_files)
  ]
)

hashes <- tools::md5sum(
  all_output_files
)

hash_tbl <- tibble(
  file=names(hashes),
  md5=unname(hashes)
)

write_csv(
  hash_tbl,
  file.path(
    TAB_DIR,
    "20B3_Fibroblast_output_hashes.csv"
  )
)

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("20B3 — FINAL REPORT\n")
cat("============================================================\n")

cat("\nPatient-level L1 pseudobulks:", ncol(patient_counts), "\n")

cat(
  "HC / FFA / AA:",
  sum(patient_info$disease=="HC"), "/",
  sum(patient_info$disease=="FFA"), "/",
  sum(patient_info$disease=="AA"), "\n"
)

cat("\nDE MODELS\n")
cat("FFA vs HC: batch + disease\n")
cat("AA vs HC : disease only; batch confounded\n")

cat("\nHEATMAPS\n")
print(heatmap_summary)

cat("\nEXPORTS\n")
cat("SVG:", length(svg_files), "\n")
cat("PDF:", length(pdf_files), "\n")
cat("PNG:", length(png_files), "\n")

cat("\nGene selection: top-ranked by PValue within logFC direction\n")
cat("UP genes per contrast: 50\n")
cat("DOWN genes per contrast: 50\n")
cat("Columns clustered: NO\n")
cat("Column order: HC -> FFA -> AA\n")
cat("Heatmap values: patient-level L1 logCPM row z-scores\n")
cat("Third disease influences gene selection: NO\n")

cat("\nDE recalculated at L1 level: YES\n")
cat("L2 DE results modified: NO\n")
cat("Frozen L2 annotations modified: NO\n")
cat("Original 18B pseudobulks modified: NO\n")

cat("\nSTATUS=PASS\n")
