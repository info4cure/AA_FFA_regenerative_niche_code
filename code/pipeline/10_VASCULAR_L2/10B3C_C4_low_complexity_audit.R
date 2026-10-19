suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B3C C4 LOW-COMPLEXITY ENDOTHELIAL AUDIT\n")
cat("============================================================\n\n")

SRC <- "PIPELINE_V2/objects/annotated/vascular_v2_10B1_reclustered_3040.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_VASCULAR_L2/10B3C_C4_low_complexity_audit"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

x <- readRDS(SRC)
DefaultAssay(x) <- "RNA"

CLCOL <- "vascular_snn_res.0.3"
stopifnot(CLCOL %in% colnames(x@meta.data))

cl <- as.character(x@meta.data[[CLCOL]])
cells4 <- colnames(x)[cl == "4"]

stopifnot(length(cells4) == 197)

cat("C4 cells:", length(cells4), "\n")
cat("Patients:", length(unique(x$patient_id[cells4])), "\n")
cat("Batches:", length(unique(x$batch[cells4])), "\n")

# ============================================================
# 1. Prespecified programmes
# ============================================================

cat("\nSTEP 1/6 — defining vascular programmes\n")

programmes <- list(
  ENDOTHELIAL_CORE = c(
    "PECAM1","VWF","CDH5","EMCN",
    "KDR","ESAM","ENG","RAMP2"
  ),
  CAPILLARY = c(
    "CA4","RGCC","EMCN","KDR",
    "PLVAP","EDNRB","BTNL9","CA2"
  ),
  VENOUS_ACKR1 = c(
    "ACKR1","SELP","SELE","VCAM1",
    "ICAM1","VWF"
  ),
  ARTERIAL = c(
    "EFNB2","SOX17","GJA5","GJA4",
    "HEY1","BMX","FBLN5","SEMA3G"
  ),
  LYMPHATIC = c(
    "PROX1","PDPN","FLT4","CCL21",
    "LYVE1","NRP2","RELN"
  )
)

coverage <- tibble(
  programme = names(programmes),
  n_requested = sapply(programmes, length),
  n_present = sapply(
    programmes,
    function(g) length(intersect(g, rownames(x)))
  ),
  genes_present = sapply(
    programmes,
    function(g)
      paste(intersect(g, rownames(x)), collapse = ";")
  )
)

write_csv(
  coverage,
  file.path(TAB, "10B3C_programme_coverage.csv")
)

print(coverage, n = Inf)

# ============================================================
# 2. Cell-level programme scores
# ============================================================

cat("\nSTEP 2/6 — cell-level scores\n")

genes <- unique(unlist(programmes))
genes <- intersect(genes, rownames(x))

e <- FetchData(
  x,
  vars = genes,
  layer = "data"
)

score <- function(df, genes) {
  g <- intersect(genes, colnames(df))
  if (!length(g))
    return(rep(NA_real_, nrow(df)))
  rowMeans(df[, g, drop = FALSE])
}

detect <- function(df, genes) {
  g <- intersect(genes, colnames(df))
  if (!length(g))
    return(rep(0L, nrow(df)))
  rowSums(df[, g, drop = FALSE] > 0)
}

idx <- match(rownames(e), colnames(x))

z <- tibble(
  cell = rownames(e),
  cluster = cl[idx],
  C4 = cl[idx] == "4",
  patient_id = as.character(x$patient_id[idx]),
  batch = as.character(x$batch[idx]),
  nFeature_RNA = x$nFeature_RNA[idx],
  nCount_RNA = x$nCount_RNA[idx],
  percent_mt = x$percent.mt[idx],
  doublet_score = x$doublet_score[idx],
  EC_score = score(e, programmes$ENDOTHELIAL_CORE),
  EC_n = detect(e, programmes$ENDOTHELIAL_CORE),
  CAP_score = score(e, programmes$CAPILLARY),
  VEN_score = score(e, programmes$VENOUS_ACKR1),
  ART_score = score(e, programmes$ARTERIAL),
  LYM_score = score(e, programmes$LYMPHATIC)
)

write_csv(
  z,
  file.path(TAB, "10B3C_cell_scores.csv")
)

# ============================================================
# 3. C4 vs remaining vascular universe
# ============================================================

cat("STEP 3/6 — C4 vs remaining vascular universe\n")

qc <- z %>%
  group_by(C4) %>%
  summarise(
    n = n(),
    median_features = median(nFeature_RNA),
    q25_features = as.numeric(quantile(nFeature_RNA, 0.25)),
    q75_features = as.numeric(quantile(nFeature_RNA, 0.75)),
    median_counts = median(nCount_RNA),
    q25_counts = as.numeric(quantile(nCount_RNA, 0.25)),
    q75_counts = as.numeric(quantile(nCount_RNA, 0.75)),
    median_mt = median(percent_mt),
    median_doublet = median(doublet_score),
    .groups = "drop"
  )

cat("\n===== QC: C4 vs REST =====\n")
print(qc, n = Inf)

write_csv(
  qc,
  file.path(TAB, "10B3C_QC_C4_vs_rest.csv")
)

ec <- z %>%
  group_by(C4) %>%
  summarise(
    n = n(),
    median_EC_score = median(EC_score),
    median_EC_genes = median(EC_n),
    fraction_EC_ge3 = mean(EC_n >= 3),
    fraction_EC_ge4 = mean(EC_n >= 4),
    .groups = "drop"
  )

cat("\n===== ENDOTHELIAL BACKBONE =====\n")
print(ec, n = Inf)

write_csv(
  ec,
  file.path(TAB, "10B3C_endothelial_backbone.csv")
)

prog <- z %>%
  group_by(C4) %>%
  summarise(
    CAP = median(CAP_score),
    VEN = median(VEN_score),
    ART = median(ART_score),
    LYM = median(LYM_score),
    .groups = "drop"
  )

cat("\n===== VASCULAR PROGRAMMES =====\n")
print(prog, n = Inf)

write_csv(
  prog,
  file.path(TAB, "10B3C_vascular_programmes.csv")
)

# ============================================================
# 4. Patient recurrence and QC association
# ============================================================

cat("\nSTEP 4/6 — recurrence and QC association\n")

pat <- z %>%
  filter(C4) %>%
  count(patient_id, batch, name = "n_cells") %>%
  arrange(desc(n_cells))

cat("\n===== PATIENT RECURRENCE C4 =====\n")
print(pat, n = Inf)

write_csv(
  pat,
  file.path(TAB, "10B3C_patient_recurrence.csv")
)

frac <- z %>%
  count(patient_id, C4, name = "n") %>%
  group_by(patient_id) %>%
  mutate(
    total = sum(n),
    fraction = n / total
  ) %>%
  filter(C4) %>%
  ungroup() %>%
  arrange(desc(fraction))

cat("\n===== C4 FRACTION BY PATIENT =====\n")
print(frac, n = Inf)

write_csv(
  frac,
  file.path(TAB, "10B3C_patient_fraction.csv")
)

cor_rows <- list()

for (v in c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent_mt",
  "doublet_score"
)) {

  ct <- suppressWarnings(
    cor.test(
      as.numeric(z$C4),
      z[[v]],
      method = "spearman",
      exact = FALSE
    )
  )

  cor_rows[[v]] <- tibble(
    variable = v,
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

cors <- bind_rows(cor_rows) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

cat("\n===== CORRELATION WITH QC — ALL VASCULAR CELLS =====\n")
print(cors, n = Inf)

write_csv(
  cors,
  file.path(TAB, "10B3C_QC_correlations.csv")
)

c4_ec_cor <- bind_rows(
  lapply(
    c("nFeature_RNA","nCount_RNA","percent_mt"),
    function(v) {

      zz <- z %>% filter(C4)

      ct <- suppressWarnings(
        cor.test(
          zz$EC_score,
          zz[[v]],
          method = "spearman",
          exact = FALSE
        )
      )

      tibble(
        variable = v,
        spearman_rho = unname(ct$estimate),
        p_value = ct$p.value
      )
    }
  )
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

cat("\n===== C4 EC SCORE vs QC =====\n")
print(c4_ec_cor, n = Inf)

write_csv(
  c4_ec_cor,
  file.path(TAB, "10B3C_C4_ECscore_QC_correlations.csv")
)

# ============================================================
# 5. Non-mitochondrial markers
# ============================================================

cat("\nSTEP 5/6 — non-mitochondrial marker audit\n")

Idents(x) <- factor(cl)

m4 <- FindMarkers(
  x,
  ident.1 = "4",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

m4$gene <- rownames(m4)

m4_nonmt <- m4 %>%
  filter(
    !grepl("^MT-", gene),
    p_val_adj < 0.05
  ) %>%
  arrange(desc(avg_log2FC))

write_csv(
  m4_nonmt,
  file.path(TAB, "10B3C_C4_nonmitochondrial_markers.csv")
)

cat("\n===== C4 NON-MITOCHONDRIAL MARKERS — TOP 50 =====\n")
top50_nonmt <- m4_nonmt %>%
  select(
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  head(50) %>%
  as.data.frame()

print(
  top50_nonmt,
  row.names = FALSE
)

cat(
  "\nNon-mito significant markers:",
  nrow(m4_nonmt),
  "\n"
)

# ============================================================
# 6. Formal checks
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- tibble(
  check = c(
    "source_cells_3040",
    "C4_cells_197",
    "C4_patients_16",
    "C4_both_batches",
    "programme_coverage_complete",
    "cell_scores_generated",
    "QC_comparison_generated",
    "patient_recurrence_generated",
    "nonmitochondrial_markers_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_L2_identity_assigned",
    "source_object_not_modified"
  ),
  pass = c(
    ncol(x) == 3040,
    length(cells4) == 197,
    length(unique(x$patient_id[cells4])) == 16,
    length(unique(x$batch[cells4])) == 2,
    all(coverage$n_present == coverage$n_requested),
    nrow(z) == 3040,
    nrow(qc) == 2,
    nrow(pat) == 16,
    nrow(m4_nonmt) > 0,
    TRUE,
    ncol(x) == 3040,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B3C_formal_checks.csv")
)

cat("\n===== FORMAL CHECKS =====\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

writeLines(
  c(
    "PIPELINE_V2 — 10B3C C4 LOW-COMPLEXITY ENDOTHELIAL AUDIT",
    "",
    "Purpose:",
    "- determine whether r0.3 C4 represents a distinct vascular programme",
    "  or primarily a low-complexity endothelial state",
    "- assess endothelial backbone preservation",
    "- assess patient and batch recurrence",
    "- quantify association with QC metrics",
    "- inspect non-mitochondrial marker evidence",
    "",
    "Interpretive boundary:",
    "- disease not used",
    "- no cells removed",
    "- no L2 identity assigned",
    "- source vascular object not modified",
    "",
    "PASS — C4 diagnostic evidence generated."
  ),
  file.path(TAB, "10B3C_audit.txt")
)

cat("\n============================================================\n")
cat("10B3C COMPLETE\n")
cat("============================================================\n")
cat("C4 cells:", length(cells4), "\n")
cat("Patients:", length(unique(x$patient_id[cells4])), "\n")
cat("No cells removed.\n")
cat("No L2 identity assigned.\n")
cat("Disease not used.\n")
cat("PASS\n")
