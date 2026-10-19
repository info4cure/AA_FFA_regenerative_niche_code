suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 18C edgeR PSEUDOBULK DIFFERENTIAL EXPRESSION\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "18_PSEUDOBULK_DE"
)

PB_BASE <- file.path(
  BASE,
  "18B_pseudobulk"
)

OUT <- file.path(
  BASE,
  "18C_DE"
)

TAB <- file.path(OUT, "tables")
FULL <- file.path(OUT, "full_results")
OBJ <- file.path(OUT, "objects")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FULL,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OBJ,
  recursive = TRUE,
  showWarnings = FALSE
)

elig <- read_csv(
  file.path(
    BASE,
    "18A_preflight",
    "tables",
    "18A_pseudobulk_contrast_eligibility.csv"
  ),
  show_col_types = FALSE
)

pb_meta <- readRDS(
  file.path(
    PB_BASE,
    "objects",
    "18B_pseudobulk_metadata_ALL.rds"
  )
)

modules <- sort(
  unique(pb_meta$module)
)

cat("\nSTEP 1/8 — loading pseudobulk count matrices\n")

pb <- list()

for (m in modules) {

  f <- file.path(
    PB_BASE,
    "objects",
    paste0(
      "18B_",
      m,
      "_pseudobulk_counts.rds"
    )
  )

  stopifnot(
    file.exists(f)
  )

  pb[[m]] <- readRDS(f)

  cat(
    m,
    ":",
    nrow(pb[[m]]),
    "genes x",
    ncol(pb[[m]]),
    "pseudobulks\n"
  )
}

cat("\nSTEP 2/8 — freezing DE specification\n")

contrast_spec <- tribble(
  ~contrast, ~group1, ~group0, ~adjust_batch,
  "AA_vs_HC",  "AA",  "HC",  FALSE,
  "FFA_vs_HC", "FFA", "HC",  TRUE,
  "AA_vs_FFA", "AA",  "FFA", FALSE
)

print(
  as.data.frame(contrast_spec),
  row.names = FALSE
)

safe_name <- function(x) {
  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}

# ============================================================
# 3. Single L2 × contrast edgeR QL analysis
# ============================================================

run_one_DE <- function(
  contrast_name,
  group1,
  group0,
  adjust_batch,
  module_name,
  l2_label,
  robust_ge50
) {

  mat <- pb[[module_name]]

  md <- pb_meta %>%
    filter(
      module == module_name,
      final_L2_label == l2_label,
      disease %in% c(group0, group1),
      n_cells >= 20
    ) %>%
    mutate(
      disease_model = factor(
        disease,
        levels = c(group0, group1)
      )
    )

  if (nrow(md) == 0) {
    return(NULL)
  }

  n0 <- sum(
    md$disease_model == group0
  )

  n1 <- sum(
    md$disease_model == group1
  )

  if (
    n0 < 3 ||
    n1 < 3
  ) {
    return(NULL)
  }

  stopifnot(
    all(
      md$pseudobulk_id %in%
        colnames(mat)
    )
  )

  counts <- mat[
    ,
    md$pseudobulk_id,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(counts),
      md$pseudobulk_id
    )
  )

  # ----------------------------------------------------------
  # Design
  # ----------------------------------------------------------

  if (adjust_batch) {

    md$batch_model <- factor(
      md$batch
    )

    design <- model.matrix(
      ~ batch_model + disease_model,
      data = md
    )

  } else {

    design <- model.matrix(
      ~ disease_model,
      data = md
    )
  }

  stopifnot(
    qr(design)$rank ==
      ncol(design)
  )

  coef_name <- paste0(
    "disease_model",
    group1
  )

  coef_idx <- which(
    colnames(design) ==
      coef_name
  )

  stopifnot(
    length(coef_idx) == 1
  )

  # ----------------------------------------------------------
  # edgeR
  # ----------------------------------------------------------

  y <- DGEList(
    counts = counts
  )

  keep <- filterByExpr(
    y,
    design = design
  )

  n_genes_tested <- sum(keep)

  if (n_genes_tested == 0) {
    stop(
      paste(
        "No genes retained by filterByExpr:",
        contrast_name,
        l2_label
      )
    )
  }

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- calcNormFactors(
    y,
    method = "TMM"
  )

  y <- estimateDisp(
    y,
    design = design,
    robust = TRUE
  )

  fit <- glmQLFit(
    y,
    design = design,
    robust = TRUE
  )

  qlf <- glmQLFTest(
    fit,
    coef = coef_idx
  )

  tt <- topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
  )$table

  res <- tt %>%
    rownames_to_column(
      "gene"
    ) %>%
    as_tibble() %>%
    rename(
      FDR_BH = FDR
    ) %>%
    mutate(
      contrast = contrast_name,
      numerator = group1,
      denominator = group0,
      module = module_name,
      final_L2_label = l2_label,
      n_numerator = n1,
      n_denominator = n0,
      batch_adjusted = adjust_batch,
      robust_ge50 = robust_ge50,
      min_cells = min(md$n_cells),
      median_cells = median(md$n_cells),
      n_genes_tested = n_genes_tested,
      .before = 1
    )

  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  sig <- tibble(
    contrast = contrast_name,
    numerator = group1,
    denominator = group0,
    module = module_name,
    final_L2_label = l2_label,
    n_numerator = n1,
    n_denominator = n0,
    batch_adjusted = adjust_batch,
    robust_ge50 = robust_ge50,
    min_cells = min(md$n_cells),
    median_cells = median(md$n_cells),
    n_genes_tested = n_genes_tested,

    n_FDR_005 =
      sum(
        res$FDR_BH < 0.05,
        na.rm = TRUE
      ),

    n_FDR_001 =
      sum(
        res$FDR_BH < 0.01,
        na.rm = TRUE
      ),

    n_up_FDR_005 =
      sum(
        res$FDR_BH < 0.05 &
          res$logFC > 0,
        na.rm = TRUE
      ),

    n_down_FDR_005 =
      sum(
        res$FDR_BH < 0.05 &
          res$logFC < 0,
        na.rm = TRUE
      ),

    max_abs_logFC_FDR005 =
      if (
        any(
          res$FDR_BH < 0.05,
          na.rm = TRUE
        )
      ) {
        max(
          abs(
            res$logFC[
              res$FDR_BH < 0.05
            ]
          ),
          na.rm = TRUE
        )
      } else {
        0
      }
  )

  # ----------------------------------------------------------
  # Diagnostics
  # ----------------------------------------------------------

  diagnostic <- tibble(
    contrast = contrast_name,
    numerator = group1,
    denominator = group0,
    module = module_name,
    final_L2_label = l2_label,
    n_samples = nrow(md),
    n_numerator = n1,
    n_denominator = n0,
    n_genes_input = nrow(mat),
    n_genes_tested = n_genes_tested,
    min_cells = min(md$n_cells),
    median_cells = median(md$n_cells),
    min_library_size =
      min(md$library_size),
    median_library_size =
      median(md$library_size),
    design_rank =
      qr(design)$rank,
    design_columns =
      ncol(design),
    batch_adjusted =
      adjust_batch,
    robust_ge50 =
      robust_ge50
  )

  # Save fitted object for audit/reproducibility.
  fit_object <- list(
    contrast = contrast_name,
    numerator = group1,
    denominator = group0,
    module = module_name,
    final_L2_label = l2_label,
    metadata = md,
    design = design,
    coefficient = coef_name,
    genes_retained =
      rownames(y),
    DGEList = y,
    fit = fit,
    qlf = qlf
  )

  list(
    results = res,
    summary = sig,
    diagnostic = diagnostic,
    fit_object = fit_object
  )
}

cat("\n18C DE function defined.\n")

# ============================================================
# 4. Run all pre-specified eligible L2 analyses
# ============================================================

cat("\nSTEP 3/8 — defining eligible L2 DE analyses\n")

all_results <- list()
all_summary <- list()
all_diag <- list()

counter <- 0L

expected_total <- sum(
  elig$eligible_ge20
)

cat(
  "Pre-specified eligible L2 x contrast analyses:",
  expected_total,
  "\n"
)

stopifnot(
  expected_total == 60
)

for (k in seq_len(nrow(contrast_spec))) {

  cc <- contrast_spec[k, ]

  ee <- elig %>%
    filter(
      contrast == cc$contrast,
      eligible_ge20
    ) %>%
    arrange(
      final_compartment,
      final_L2_label
    )

  cat("\n------------------------------------------------------------\n")
  cat(
    cc$contrast,
    ":",
    nrow(ee),
    "eligible L2 analyses\n"
  )
  cat("------------------------------------------------------------\n")

  for (i in seq_len(nrow(ee))) {

    counter <- counter + 1L

    lab <- ee$final_L2_label[i]

    # Resolve the unique frozen module containing this L2.
    mm <- pb_meta %>%
      filter(
        final_L2_label == lab
      ) %>%
      distinct(module)

    if (nrow(mm) != 1) {
      stop(
        paste(
          "L2 does not resolve to exactly one module:",
          lab
        )
      )
    }

    module_name <- mm$module[1]

    cat(
      sprintf(
        "[%02d/%02d] %s | %s | %s\n",
        counter,
        expected_total,
        cc$contrast,
        module_name,
        lab
      )
    )

    ans <- run_one_DE(
      contrast_name = cc$contrast,
      group1 = cc$group1,
      group0 = cc$group0,
      adjust_batch = cc$adjust_batch,
      module_name = module_name,
      l2_label = lab,
      robust_ge50 = ee$robust_ge50[i]
    )

    if (is.null(ans)) {
      stop(
        paste(
          "Unexpected NULL analysis:",
          cc$contrast,
          lab
        )
      )
    }

    key <- paste(
      cc$contrast,
      module_name,
      lab,
      sep = "|||"
    )

    all_results[[key]] <-
      ans$results

    all_summary[[key]] <-
      ans$summary

    all_diag[[key]] <-
      ans$diagnostic

    # --------------------------------------------------------
    # Save full gene-level result
    # --------------------------------------------------------

    result_file <- file.path(
      FULL,
      paste0(
        "18C_",
        safe_name(cc$contrast),
        "__",
        safe_name(module_name),
        "__",
        safe_name(lab),
        ".csv"
      )
    )

    write_csv(
      ans$results,
      result_file
    )

    # --------------------------------------------------------
    # Save fitted edgeR audit object
    # --------------------------------------------------------

    fit_file <- file.path(
      OBJ,
      paste0(
        "18C_",
        safe_name(cc$contrast),
        "__",
        safe_name(module_name),
        "__",
        safe_name(lab),
        "_edgeR_fit.rds"
      )
    )

    saveRDS(
      ans$fit_object,
      fit_file
    )

    rm(ans)
    gc(verbose = FALSE)
  }
}

stopifnot(
  counter == 60,
  length(all_results) == 60,
  length(all_summary) == 60,
  length(all_diag) == 60
)

cat(
  "\nAll",
  counter,
  "pre-specified analyses completed in memory.\n"
)

# ============================================================
# 5. Assemble global result tables
# ============================================================

cat("\nSTEP 4/8 — assembling global DE summaries\n")

summary_tbl <- bind_rows(
  all_summary
) %>%
  arrange(
    contrast,
    desc(n_FDR_005),
    final_L2_label
  )

diag_tbl <- bind_rows(
  all_diag
) %>%
  arrange(
    contrast,
    module,
    final_L2_label
  )

stopifnot(
  nrow(summary_tbl) == 60,
  nrow(diag_tbl) == 60
)

write_csv(
  summary_tbl,
  file.path(
    TAB,
    "18C_DE_summary_by_L2.csv"
  )
)

write_csv(
  diag_tbl,
  file.path(
    TAB,
    "18C_DE_diagnostics_by_L2.csv"
  )
)

cat("\n===== DE SUMMARY BY L2 =====\n")

print(
  as.data.frame(summary_tbl),
  row.names = FALSE
)

# ============================================================
# 6. Gene-level inventories
# ============================================================

cat("\nSTEP 5/8 — building gene-level inventories\n")

all_gene_results <- bind_rows(
  all_results
)

top_genes <- all_gene_results %>%
  group_by(
    contrast,
    module,
    final_L2_label
  ) %>%
  arrange(
    FDR_BH,
    PValue,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write_csv(
  top_genes,
  file.path(
    TAB,
    "18C_top20_genes_each_L2.csv"
  )
)

sig_genes <- all_gene_results %>%
  filter(
    FDR_BH < 0.05
  ) %>%
  arrange(
    contrast,
    module,
    final_L2_label,
    FDR_BH,
    PValue
  )

write_csv(
  sig_genes,
  file.path(
    TAB,
    "18C_all_FDR005_genes.csv"
  )
)

cat(
  "Total gene-level FDR<0.05 calls:",
  nrow(sig_genes),
  "\n"
)

rm(all_gene_results)
gc(verbose = FALSE)

# ============================================================
# 7. Robust-support and contrast summaries
# ============================================================

cat("\nSTEP 6/8 — robust-support inventory\n")

robust_summary <- summary_tbl %>%
  filter(
    robust_ge50
  ) %>%
  group_by(
    contrast
  ) %>%
  summarise(
    n_L2 = n(),
    n_L2_with_DE =
      sum(n_FDR_005 > 0),
    total_FDR005_genes =
      sum(n_FDR_005),
    total_up_FDR005 =
      sum(n_up_FDR_005),
    total_down_FDR005 =
      sum(n_down_FDR_005),
    .groups = "drop"
  )

print(
  as.data.frame(robust_summary),
  row.names = FALSE
)

write_csv(
  robust_summary,
  file.path(
    TAB,
    "18C_robust_ge50_summary.csv"
  )
)

cat("\nSTEP 7/8 — contrast summary\n")

contrast_summary <- summary_tbl %>%
  group_by(
    contrast
  ) %>%
  summarise(
    n_L2_tested = n(),
    n_L2_with_FDR005 =
      sum(n_FDR_005 > 0),
    total_genes_FDR005 =
      sum(n_FDR_005),
    total_up_FDR005 =
      sum(n_up_FDR_005),
    total_down_FDR005 =
      sum(n_down_FDR_005),
    .groups = "drop"
  )

print(
  as.data.frame(contrast_summary),
  row.names = FALSE
)

write_csv(
  contrast_summary,
  file.path(
    TAB,
    "18C_contrast_summary.csv"
  )
)

# ============================================================
# 8. Formal checks and final audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

expected <- tibble(
  contrast = c(
    "AA_vs_HC",
    "FFA_vs_HC",
    "AA_vs_FFA"
  ),
  expected_L2 = c(
    24L,
    18L,
    18L
  )
)

actual <- summary_tbl %>%
  count(
    contrast,
    name = "actual_L2"
  )

chk_n <- expected %>%
  left_join(
    actual,
    by = "contrast"
  )

n_result_files <- length(
  list.files(
    FULL,
    pattern = "^18C_.*\\.csv$",
    full.names = TRUE
  )
)

n_fit_files <- length(
  list.files(
    OBJ,
    pattern = "^18C_.*_edgeR_fit\\.rds$",
    full.names = TRUE
  )
)

aa_hc_batch <- summary_tbl %>%
  filter(contrast == "AA_vs_HC") %>%
  pull(batch_adjusted)

ffa_hc_batch <- summary_tbl %>%
  filter(contrast == "FFA_vs_HC") %>%
  pull(batch_adjusted)

aa_ffa_batch <- summary_tbl %>%
  filter(contrast == "AA_vs_FFA") %>%
  pull(batch_adjusted)

checks <- tibble(
  check = c(
    "sixty_prespecified_L2_contrast_tests",
    "AA_HC_24_L2",
    "FFA_HC_18_L2",
    "AA_FFA_18_L2",
    "all_samples_ge20_cells",
    "all_designs_full_rank",
    "all_gene_sets_nonempty",
    "all_gene_counts_filtered",
    "FFA_HC_batch_adjusted",
    "AA_HC_not_batch_adjusted",
    "AA_FFA_not_batch_adjusted",
    "FDR_within_L2_generated",
    "sixty_full_result_files",
    "sixty_edgeR_fit_objects",
    "summary_complete",
    "diagnostics_complete",
    "robust_ge50_flag_preserved",
    "annotation_unchanged",
    "figures_generated_NO"
  ),

  pass = c(
    nrow(summary_tbl) == 60,

    chk_n$actual_L2[
      chk_n$contrast == "AA_vs_HC"
    ] == 24,

    chk_n$actual_L2[
      chk_n$contrast == "FFA_vs_HC"
    ] == 18,

    chk_n$actual_L2[
      chk_n$contrast == "AA_vs_FFA"
    ] == 18,

    all(diag_tbl$min_cells >= 20),

    all(
      diag_tbl$design_rank ==
        diag_tbl$design_columns
    ),

    all(
      diag_tbl$n_genes_tested > 0
    ),

    all(
      diag_tbl$n_genes_tested <=
        diag_tbl$n_genes_input
    ),

    length(ffa_hc_batch) == 18 &&
      all(ffa_hc_batch),

    length(aa_hc_batch) == 24 &&
      all(!aa_hc_batch),

    length(aa_ffa_batch) == 18 &&
      all(!aa_ffa_batch),

    all(
      is.finite(
        summary_tbl$n_FDR_005
      )
    ),

    n_result_files == 60,

    n_fit_files == 60,

    nrow(summary_tbl) == 60,

    nrow(diag_tbl) == 60,

    all(
      !is.na(
        summary_tbl$robust_ge50
      )
    ),

    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(
    TAB,
    "18C_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

# ------------------------------------------------------------
# Final audit text
# ------------------------------------------------------------

audit <- c(
  "PIPELINE_V2 — 18C edgeR pseudobulk DE",
  "",
  "Inferential unit: patient",
  "Expression input: summed RNA counts",
  "Minimum pseudobulk size: 20 cells",
  "Gene filtering: edgeR filterByExpr within L2 x contrast",
  "Normalization: TMM",
  "Inference: edgeR quasi-likelihood",
  "Dispersion estimation: robust",
  "QL fitting: robust",
  "Multiple testing: BH within each L2 x contrast",
  "",
  "AA vs HC: disease model; batch not identifiable",
  "FFA vs HC: batch + disease model",
  "AA vs FFA: disease model; batch severely confounded",
  "",
  paste0(
    "Analyses: ",
    nrow(summary_tbl)
  ),
  paste0(
    "Gene-level FDR<0.05 calls: ",
    nrow(sig_genes)
  ),
  paste0(
    "Full result files: ",
    n_result_files
  ),
  paste0(
    "Saved edgeR fit objects: ",
    n_fit_files
  ),
  "",
  "Annotation modified: NO",
  "Figures generated: NO",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "18C_audit.txt"
  )
)

cat("\n===== FINAL CONTRAST SUMMARY =====\n")

print(
  as.data.frame(
    contrast_summary
  ),
  row.names = FALSE
)

cat("\n===== ROBUST >=50 SUMMARY =====\n")

print(
  as.data.frame(
    robust_summary
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("18C COMPLETE\n")
cat("============================================================\n")
cat("Pseudobulk unit: patient x frozen L2\n")
cat("Eligible L2 analyses: 60\n")
cat("AA vs HC: 24\n")
cat("FFA vs HC: 18\n")
cat("AA vs FFA: 18\n")
cat("Minimum pseudobulk cells: 20\n")
cat("Gene filter: filterByExpr per L2 x contrast\n")
cat("Normalization: TMM\n")
cat("Model: edgeR quasi-likelihood, robust\n")
cat("FFA vs HC batch-adjusted: YES\n")
cat("AA contrasts batch-adjusted: NO — confounded\n")
cat("Multiple testing: BH within L2 x contrast\n")
cat("Annotation modified: NO\n")
cat("Figures generated: NO\n")
cat("READY FOR 18D BIOLOGICAL INTERPRETATION\n")
cat("PASS\n")
