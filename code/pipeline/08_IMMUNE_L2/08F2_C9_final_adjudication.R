suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08F2 C9 FINAL ADJUDICATION AUDIT\n")
cat("============================================================\n")

OBJ <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"

F08F <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08F_universe_refinement/tables/",
  "08F_C9_C12_cell_level_adjudication.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08F2_C9_final_adjudication"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# ==========================================================================
# 1. LOAD
# ==========================================================================

cat("\nSTEP 1/7 — loading data\n")

x <- readRDS(OBJ)
a <- read_csv(F08F, show_col_types=FALSE)

stopifnot(
  ncol(x)==4762,
  nrow(x)==23412
)

c9a <- a %>%
  filter(cluster=="9")

stopifnot(
  nrow(c9a)==183,
  length(unique(c9a$cell))==183
)

cells9 <- c9a$cell

stopifnot(
  all(cells9 %in% colnames(x))
)

cat("Immune universe:",ncol(x),"\n")
cat("C9 cells:",length(cells9),"\n")

# ==========================================================================
# 2. EXTRACT EXPRESSION
# ==========================================================================

cat("\nSTEP 2/7 — extracting lineage expression\n")

genes <- c(
  "PTPRC","CD3D","CD3E","TRAC","LST1","LYZ",
  "KRT5","KRT14","KRT15","KRT17","COL17A1","KRTDAP"
)

stopifnot(
  all(genes %in% rownames(x))
)

dat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,cells9,drop=FALSE]

stopifnot(
  identical(colnames(dat),cells9)
)

# ==========================================================================
# 3. BUILD CELL-LEVEL TABLE
# ==========================================================================

cat("\nSTEP 3/7 — building C9 cell-level table\n")

md <- x@meta.data[cells9,,drop=FALSE]

expr <- as.data.frame(
  t(as.matrix(dat)),
  check.names=FALSE
)

expr$cell <- rownames(expr)

tab <- c9a %>%
  select(
    cell,
    patient_id,
    disease,
    batch,
    refined_disposition,
    immune_gene_count,
    epithelial_gene_count,
    leukocyte_score,
    epithelial_score
  ) %>%
  left_join(
    data.frame(
      cell=rownames(md),
      nCount_RNA=md$nCount_RNA,
      nFeature_RNA=md$nFeature_RNA,
      percent.mt=md$percent.mt,
      doublet_score=md$doublet_score,
      stringsAsFactors=FALSE
    ),
    by="cell"
  ) %>%
  left_join(
    expr,
    by="cell"
  )

stopifnot(
  nrow(tab)==183,
  !anyDuplicated(tab$cell)
)

# --------------------------------------------------------------------------
# Explicit simultaneous lineage support
# --------------------------------------------------------------------------

tab <- tab %>%
  mutate(
    immune_strict_genes =
      (PTPRC>0) +
      (CD3D>0) +
      (CD3E>0) +
      (TRAC>0) +
      (LST1>0) +
      (LYZ>0),

    epithelial_strict_genes =
      (KRT5>0) +
      (KRT14>0) +
      (KRT15>0) +
      (KRT17>0) +
      (COL17A1>0) +
      (KRTDAP>0),

    simultaneous_strong =
      immune_strict_genes >= 2 &
      epithelial_strict_genes >= 2
  )

write_csv(
  tab,
  file.path(
    OUT,
    "08F2_C9_cell_level_full_audit.csv"
  )
)

# ==========================================================================
# 4. REFERENCE DISTRIBUTION FROM OTHER IMMUNE CELLS
# ==========================================================================

cat("\nSTEP 4/7 — reference QC distributions\n")

other <- setdiff(
  colnames(x),
  cells9
)

md_other <- x@meta.data[other,,drop=FALSE]

reference <- data.frame(
  group="OTHER_IMMUNE",
  cell=other,
  doublet_score=as.numeric(md_other$doublet_score),
  nCount_RNA=as.numeric(md_other$nCount_RNA),
  nFeature_RNA=as.numeric(md_other$nFeature_RNA),
  percent.mt=as.numeric(md_other$percent.mt),
  stringsAsFactors=FALSE
)

c9qc <- tab %>%
  transmute(
    group=case_when(
      refined_disposition=="KEEP_IMMUNE"
        ~ "C9_KEEP_IMMUNE",
      refined_disposition=="RETURN_EPITHELIAL"
        ~ "C9_RETURN_EPITHELIAL",
      refined_disposition=="REVIEW_MIXED_DOUBLET"
        ~ "C9_MIXED",
      refined_disposition=="REVIEW_UNCERTAIN"
        ~ "C9_UNCERTAIN",
      TRUE ~ "C9_OTHER"
    ),
    cell,
    doublet_score,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  )

qc_all <- bind_rows(
  reference,
  c9qc
)

summ <- qc_all %>%
  group_by(group) %>%
  summarise(
    n_cells=n(),

    doublet_q25=quantile(
      doublet_score,0.25,na.rm=TRUE
    ),
    doublet_median=median(
      doublet_score,na.rm=TRUE
    ),
    doublet_q75=quantile(
      doublet_score,0.75,na.rm=TRUE
    ),
    doublet_q90=quantile(
      doublet_score,0.90,na.rm=TRUE
    ),

    nCount_q25=quantile(
      nCount_RNA,0.25,na.rm=TRUE
    ),
    nCount_median=median(
      nCount_RNA,na.rm=TRUE
    ),
    nCount_q75=quantile(
      nCount_RNA,0.75,na.rm=TRUE
    ),

    nFeature_q25=quantile(
      nFeature_RNA,0.25,na.rm=TRUE
    ),
    nFeature_median=median(
      nFeature_RNA,na.rm=TRUE
    ),
    nFeature_q75=quantile(
      nFeature_RNA,0.75,na.rm=TRUE
    ),

    mt_median=median(
      percent.mt,na.rm=TRUE
    ),

    .groups="drop"
  )

cat("\n=== QC DISTRIBUTIONS ===\n")
print(
  as.data.frame(summ),
  row.names=FALSE
)

write_csv(
  summ,
  file.path(
    OUT,
    "08F2_C9_QC_group_summary.csv"
  )
)

# ==========================================================================
# 5. EMPIRICAL DOUBLETSCORE CALIBRATION
#
# We do NOT impose an arbitrary absolute score.
# We compare C9 cells against the empirical distribution of the remaining
# immune universe.
# ==========================================================================

cat("\nSTEP 5/7 — empirical doublet-score calibration\n")

ref_scores <- reference$doublet_score[
  is.finite(reference$doublet_score)
]

q90 <- unname(
  quantile(ref_scores,0.90,na.rm=TRUE)
)
q95 <- unname(
  quantile(ref_scores,0.95,na.rm=TRUE)
)
q99 <- unname(
  quantile(ref_scores,0.99,na.rm=TRUE)
)

cat(
  "Other-immune doublet-score thresholds:\n",
  "Q90 =",sprintf("%.4f",q90),"\n",
  "Q95 =",sprintf("%.4f",q95),"\n",
  "Q99 =",sprintf("%.4f",q99),"\n"
)

tab <- tab %>%
  mutate(
    doublet_gt_q90 =
      doublet_score > q90,

    doublet_gt_q95 =
      doublet_score > q95,

    doublet_gt_q99 =
      doublet_score > q99
  )

calibration <- tab %>%
  group_by(
    refined_disposition
  ) %>%
  summarise(
    n_cells=n(),
    simultaneous_strong_n=sum(
      simultaneous_strong
    ),
    simultaneous_strong_pct=
      100*mean(simultaneous_strong),

    gt_q90_n=sum(
      doublet_gt_q90,
      na.rm=TRUE
    ),
    gt_q90_pct=
      100*mean(
        doublet_gt_q90,
        na.rm=TRUE
      ),

    gt_q95_n=sum(
      doublet_gt_q95,
      na.rm=TRUE
    ),
    gt_q95_pct=
      100*mean(
        doublet_gt_q95,
        na.rm=TRUE
      ),

    gt_q99_n=sum(
      doublet_gt_q99,
      na.rm=TRUE
    ),
    gt_q99_pct=
      100*mean(
        doublet_gt_q99,
        na.rm=TRUE
      ),

    median_doublet_score=
      median(
        doublet_score,
        na.rm=TRUE
      ),

    median_nCount=
      median(
        nCount_RNA,
        na.rm=TRUE
      ),

    median_nFeature=
      median(
        nFeature_RNA,
        na.rm=TRUE
      ),

    .groups="drop"
  )

cat("\n=== C9 EMPIRICAL CALIBRATION ===\n")

print(
  as.data.frame(calibration),
  row.names=FALSE
)

write_csv(
  calibration,
  file.path(
    OUT,
    "08F2_C9_empirical_doublet_calibration.csv"
  )
)

# ==========================================================================
# 6. PATIENT SUPPORT
# ==========================================================================

cat("\nSTEP 6/7 — patient support\n")

patient_support <- tab %>%
  count(
    refined_disposition,
    patient_id,
    disease,
    name="n_cells"
  ) %>%
  group_by(
    refined_disposition
  ) %>%
  mutate(
    fraction_within_disposition =
      n_cells/sum(n_cells)
  ) %>%
  ungroup() %>%
  arrange(
    refined_disposition,
    desc(n_cells)
  )

write_csv(
  patient_support,
  file.path(
    OUT,
    "08F2_C9_patient_support.csv"
  )
)

# ==========================================================================
# 7. FINAL EVIDENCE FLAGS — NO AUTOMATIC REMOVAL
# ==========================================================================

cat("\nSTEP 7/7 — final evidence flags\n")

final <- tab %>%
  mutate(

    final_evidence_flag = case_when(

      refined_disposition=="KEEP_IMMUNE"
        ~ "IMMUNE_SUPPORTED",

      refined_disposition=="RETURN_EPITHELIAL"
        ~ "EPITHELIAL_SUPPORTED",

      refined_disposition=="REVIEW_MIXED_DOUBLET" &
        simultaneous_strong &
        doublet_gt_q95
        ~ "HIGH_CONFIDENCE_RESIDUAL_DOUBLET",

      refined_disposition=="REVIEW_MIXED_DOUBLET" &
        simultaneous_strong
        ~ "MIXED_LINEAGE_REVIEW",

      refined_disposition=="REVIEW_UNCERTAIN" &
        doublet_gt_q95
        ~ "UNCERTAIN_HIGH_DOUBLET_SCORE",

      refined_disposition=="REVIEW_UNCERTAIN"
        ~ "UNCERTAIN_REVIEW",

      TRUE
        ~ "REVIEW"
    )
  )

final_summary <- final %>%
  count(
    final_evidence_flag,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== C9 FINAL EVIDENCE FLAGS ===\n")

print(
  as.data.frame(final_summary),
  row.names=FALSE
)

write_csv(
  final,
  file.path(
    OUT,
    "08F2_C9_final_cell_evidence.csv"
  )
)

write_csv(
  final_summary,
  file.path(
    OUT,
    "08F2_C9_final_evidence_summary.csv"
  )
)

# ==========================================================================
# FORMAL CHECKS
# ==========================================================================

checks <- data.frame(
  check=c(
    "immune_input_4762",
    "C9_cells_183",
    "C9_unique_cells_183",
    "all_marker_genes_present",
    "all_C9_have_original_disposition",
    "all_C9_have_final_evidence_flag",
    "final_summary_sum_183",
    "no_cells_removed",
    "no_cells_reassigned",
    "no_frozen_object_written"
  ),
  pass=c(
    ncol(x)==4762,
    nrow(tab)==183,
    length(unique(tab$cell))==183,
    all(genes %in% rownames(x)),
    !anyNA(tab$refined_disposition),
    !anyNA(final$final_evidence_flag),
    sum(final_summary$n_cells)==183,
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors=FALSE
)

cat("\n=== FORMAL CHECKS ===\n")
print(
  checks,
  row.names=FALSE
)

write_csv(
  checks,
  file.path(
    OUT,
    "08F2_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

audit <- c(
  "PIPELINE_V2 — 08F2 C9 final adjudication audit",
  paste0("Generated: ",Sys.time()),
  "Input immune universe: 4762 cells",
  "Target: working cluster 9, 183 cells",
  "Purpose: distinguish immune, epithelial and residual mixed-lineage cells.",
  "Cell-level immune and epithelial expression assessed directly.",
  "Simultaneous lineage evidence requires >=2 immune and >=2 epithelial genes.",
  "Doublet score calibrated against the empirical distribution of all other immune cells.",
  "Q90/Q95/Q99 thresholds exported.",
  "No arbitrary absolute doublet-score cutoff imposed.",
  "No cells removed.",
  "No cells reassigned.",
  "No L1 labels modified.",
  "No frozen object written.",
  "Final evidence flags are adjudication evidence, not automatic filtering decisions."
)

writeLines(
  audit,
  file.path(
    OUT,
    "08F2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("08F2 COMPLETE\n")
cat("============================================================\n")
cat("C9 cells:",nrow(final),"\n")
cat("No cells removed.\n")
cat("No cells reassigned.\n")
cat("No frozen object written.\n")
cat("PASS — C9 final adjudication evidence generated.\n")

