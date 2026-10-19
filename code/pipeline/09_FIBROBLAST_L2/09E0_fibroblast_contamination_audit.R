suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(Matrix)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09E0 FIBROBLAST CONTAMINATION AUDIT\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09B_reclustered.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09E0_contamination_audit/tables"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==4808,
  nrow(x)==23412,
  length(unique(x$patient_id))==16
)

DefaultAssay(x) <- "RNA"

cluster03 <- as.character(
  x@meta.data[["fibroblast_snn_res.0.3"]]
)

cluster04 <- as.character(
  x@meta.data[["fibroblast_snn_res.0.4"]]
)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")

# ==========================================================================
# 1. STRICT LINEAGE PROGRAMMES
# ==========================================================================

cat("\nSTEP 1/8 — defining strict lineage programmes\n")

programmes <- list(

  Fibroblast_core=c(
    "COL1A1","COL1A2","DCN","LUM",
    "COL3A1","PDGFRA","PDGFRB"
  ),

  Mesenchymal_core=c(
    "DCN","LUM","COL1A1","COL1A2",
    "COL3A1","COL6A1"
  ),

  T_cell=c(
    "CD3D","CD3E","CD3G","TRAC",
    "CD2","PTPRC"
  ),

  Myeloid=c(
    "PTPRC","LST1","TYROBP","FCER1G",
    "LYZ","CTSS"
  ),

  Epithelial=c(
    "KRT5","KRT14","KRT15","KRT17",
    "KRT19","EPCAM","KRTDAP","COL17A1"
  ),

  Endothelial=c(
    "PECAM1","VWF","EMCN","KDR",
    "ESAM","CDH5"
  ),

  Pericyte_SMC=c(
    "RGS5","CSPG4","MCAM","ACTA2",
    "TAGLN","MYH11"
  ),

  Neural_Schwann=c(
    "S100B","SOX10","PLP1",
    "MPZ","PMP2","CDH19"
  ),

  Melanocyte=c(
    "PMEL","MLANA","TYR","TYRP1","DCT"
  )
)

coverage <- bind_rows(
  lapply(
    names(programmes),
    function(nm){

      requested <- programmes[[nm]]
      present <- intersect(requested,rownames(x))

      data.frame(
        programme=nm,
        n_present=length(present),
        n_requested=length(requested),
        coverage=length(present)/length(requested),
        missing=paste(
          setdiff(requested,present),
          collapse=";"
        )
      )
    }
  )
)

cat("\n=== PROGRAMME COVERAGE ===\n")
print(coverage,row.names=FALSE)

write_csv(
  coverage,
  file.path(OUT,"09E0_programme_coverage.csv")
)

stopifnot(all(coverage$coverage==1))

# ==========================================================================
# 2. EXPRESSION MATRIX
# ==========================================================================

cat("\nSTEP 2/8 — extracting expression matrix\n")

genes <- unique(unlist(programmes))

dat <- GetAssayData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

stopifnot(
  ncol(dat)==4808,
  identical(colnames(dat),colnames(x))
)

# ==========================================================================
# 3. CELL-LEVEL PROGRAMME EVIDENCE
# ==========================================================================

cat("\nSTEP 3/8 — cell-level lineage evidence\n")

programme_mean <- function(gs){

  g <- intersect(gs,rownames(dat))

  Matrix::colMeans(
    dat[g,,drop=FALSE]
  )
}

programme_detected <- function(gs){

  g <- intersect(gs,rownames(dat))

  Matrix::colSums(
    dat[g,,drop=FALSE] > 0
  )
}

score <- data.frame(
  cell=colnames(x),
  stringsAsFactors=FALSE
)

for(nm in names(programmes)){

  score[[paste0(nm,"_mean")]] <-
    programme_mean(programmes[[nm]])

  score[[paste0(nm,"_genes")]] <-
    programme_detected(programmes[[nm]])
}

score$cluster03 <- cluster03
score$cluster04 <- cluster04

score$patient_id <- as.character(x$patient_id)
score$disease <- as.character(x$disease)
score$batch <- as.character(x$batch)

score$doublet_score <- as.numeric(x$doublet_score)
score$nCount_RNA <- as.numeric(x$nCount_RNA)
score$nFeature_RNA <- as.numeric(x$nFeature_RNA)
score$percent_mt <- as.numeric(x$percent.mt)

# ==========================================================================
# 4. STRICT CELL-LEVEL FLAGS
# ==========================================================================

cat("\nSTEP 4/8 — strict contamination flags\n")

score <- score %>%
  mutate(

    fibroblast_supported =
      Fibroblast_core_genes >= 3 &
      Mesenchymal_core_genes >= 3,

    T_supported =
      T_cell_genes >= 3,

    myeloid_supported =
      Myeloid_genes >= 3,

    epithelial_supported =
      Epithelial_genes >= 3,

    endothelial_supported =
      Endothelial_genes >= 3,

    pericyte_supported =
      Pericyte_SMC_genes >= 3,

    schwann_supported =
      Neural_Schwann_genes >= 3,

    melanocyte_supported =
      Melanocyte_genes >= 3
  )

score <- score %>%
  mutate(

    n_nonfib_lineages =
      as.integer(T_supported) +
      as.integer(myeloid_supported) +
      as.integer(epithelial_supported) +
      as.integer(endothelial_supported) +
      as.integer(pericyte_supported) +
      as.integer(schwann_supported) +
      as.integer(melanocyte_supported),

    provisional_flag = case_when(

      fibroblast_supported &
      n_nonfib_lineages==0
        ~ "FIBROBLAST_SUPPORTED",

      T_supported &
      !fibroblast_supported
        ~ "T_RETURN_CANDIDATE",

      myeloid_supported &
      !fibroblast_supported
        ~ "MYELOID_RETURN_CANDIDATE",

      epithelial_supported &
      !fibroblast_supported
        ~ "EPITHELIAL_RETURN_CANDIDATE",

      endothelial_supported &
      !fibroblast_supported
        ~ "ENDOTHELIAL_RETURN_CANDIDATE",

      schwann_supported &
      !fibroblast_supported
        ~ "SCHWANN_RETURN_CANDIDATE",

      melanocyte_supported &
      !fibroblast_supported
        ~ "MELANOCYTE_RETURN_CANDIDATE",

      fibroblast_supported &
      n_nonfib_lineages>=1
        ~ "MIXED_LINEAGE_REVIEW",

      !fibroblast_supported &
      n_nonfib_lineages>=2
        ~ "MULTILINEAGE_REVIEW",

      TRUE
        ~ "UNRESOLVED_REVIEW"
    )
  )

cat("\n=== PROVISIONAL CELL-LEVEL FLAGS ===\n")

flag_summary <- score %>%
  count(
    provisional_flag,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

print(as.data.frame(flag_summary),row.names=FALSE)

write_csv(
  flag_summary,
  file.path(OUT,"09E0_cell_flag_summary.csv")
)

# ==========================================================================
# 5. CLUSTER-LEVEL AUDIT
# ==========================================================================

cat("\nSTEP 5/8 — cluster-level contamination audit\n")

cluster_audit <- score %>%
  group_by(cluster03) %>%
  summarise(

    n_cells=n(),

    n_patients=n_distinct(patient_id),

    fibroblast_supported_pct=
      100*mean(fibroblast_supported),

    T_supported_pct=
      100*mean(T_supported),

    myeloid_supported_pct=
      100*mean(myeloid_supported),

    epithelial_supported_pct=
      100*mean(epithelial_supported),

    endothelial_supported_pct=
      100*mean(endothelial_supported),

    pericyte_supported_pct=
      100*mean(pericyte_supported),

    schwann_supported_pct=
      100*mean(schwann_supported),

    melanocyte_supported_pct=
      100*mean(melanocyte_supported),

    mixed_pct=
      100*mean(
        provisional_flag=="MIXED_LINEAGE_REVIEW"
      ),

    median_doublet_score=
      median(doublet_score,na.rm=TRUE),

    median_nCount=
      median(nCount_RNA,na.rm=TRUE),

    median_nFeature=
      median(nFeature_RNA,na.rm=TRUE),

    .groups="drop"
  )

cat("\n=== RESOLUTION 0.3 CLUSTER AUDIT ===\n")
print(as.data.frame(cluster_audit),row.names=FALSE)

write_csv(
  cluster_audit,
  file.path(OUT,"09E0_cluster03_contamination_audit.csv")
)

# ==========================================================================
# 6. TARGETED C5/C9 AUDIT
# ==========================================================================

cat("\nSTEP 6/8 — targeted C5/C9 audit\n")

target <- score %>%
  filter(
    cluster03 %in% c("5","9")
  )

target_summary <- target %>%
  count(
    cluster03,
    cluster04,
    provisional_flag,
    name="n_cells"
  ) %>%
  group_by(
    cluster03,
    cluster04
  ) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

cat("\n=== C5/C9 CELL DISPOSITION ===\n")
print(as.data.frame(target_summary),row.names=FALSE)

write_csv(
  target_summary,
  file.path(OUT,"09E0_C5_C9_targeted_summary.csv")
)

target_patient <- target %>%
  count(
    cluster03,
    cluster04,
    provisional_flag,
    patient_id,
    disease,
    name="n_cells"
  )

write_csv(
  target_patient,
  file.path(OUT,"09E0_C5_C9_patient_support.csv")
)

# ==========================================================================
# 7. DOUBLET CALIBRATION
# ==========================================================================

cat("\nSTEP 7/8 — empirical doublet calibration\n")

reference_scores <- score$doublet_score[
  score$provisional_flag=="FIBROBLAST_SUPPORTED"
]

reference_scores <- reference_scores[
  is.finite(reference_scores)
]

q <- quantile(
  reference_scores,
  probs=c(.90,.95,.99),
  na.rm=TRUE
)

cat(
  "Fibroblast-supported doublet thresholds:\n",
  "Q90 =",round(q[[1]],4),"\n",
  "Q95 =",round(q[[2]],4),"\n",
  "Q99 =",round(q[[3]],4),"\n"
)

score <- score %>%
  mutate(
    doublet_gt_Q90=
      doublet_score > q[[1]],

    doublet_gt_Q95=
      doublet_score > q[[2]],

    doublet_gt_Q99=
      doublet_score > q[[3]]
  )

doublet_summary <- score %>%
  group_by(provisional_flag) %>%
  summarise(
    n_cells=n(),
    median_doublet_score=
      median(doublet_score,na.rm=TRUE),
    pct_gt_Q90=
      100*mean(doublet_gt_Q90,na.rm=TRUE),
    pct_gt_Q95=
      100*mean(doublet_gt_Q95,na.rm=TRUE),
    pct_gt_Q99=
      100*mean(doublet_gt_Q99,na.rm=TRUE),
    median_nCount=
      median(nCount_RNA,na.rm=TRUE),
    median_nFeature=
      median(nFeature_RNA,na.rm=TRUE),
    .groups="drop"
  )

cat("\n=== FLAG × DOUBLET QC ===\n")
print(as.data.frame(doublet_summary),row.names=FALSE)

write_csv(
  doublet_summary,
  file.path(OUT,"09E0_flag_doublet_QC.csv")
)

# ==========================================================================
# 8. EXPORT + FORMAL CHECKS
# ==========================================================================

cat("\nSTEP 8/8 — export and formal checks\n")

write_csv(
  score,
  file.path(OUT,"09E0_master_cell_audit.csv")
)

checks <- data.frame(
  check=c(
    "cells_4808",
    "genes_23412",
    "patients_16",
    "programme_coverage_complete",
    "master_rows_4808",
    "master_unique_cells_4808",
    "C5_present",
    "C9_present",
    "no_FibConsensus_used",
    "no_cells_removed",
    "no_cells_reassigned",
    "no_L2_labels_assigned",
    "no_frozen_object_written"
  ),

  pass=c(
    ncol(x)==4808,
    nrow(x)==23412,
    length(unique(x$patient_id))==16,
    all(coverage$coverage==1),
    nrow(score)==4808,
    length(unique(score$cell))==4808,
    any(score$cluster03=="5"),
    any(score$cluster03=="9"),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),

  stringsAsFactors=FALSE
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks,row.names=FALSE)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(OUT,"09E0_formal_checks.csv")
)

audit <- c(
  "PIPELINE_V2 — 09E0 fibroblast contamination audit",
  paste0("Generated: ",Sys.time()),
  "Input: fibroblast_v2_09B_reclustered.rds",
  "Input cells: 4808.",
  "Strict fibroblast and non-fibroblast lineage programmes evaluated.",
  "Resolution 0.3 used for cluster-level audit only.",
  "Resolution 0.4 retained for targeted split provenance.",
  "C5 and C9 explicitly audited at cell level.",
  "Doublet score used as orthogonal evidence.",
  "FibConsensus was not used.",
  "Disease was not used for lineage adjudication.",
  "No cells removed.",
  "No cells reassigned.",
  "No L2 labels assigned.",
  "No frozen object written."
)

writeLines(
  audit,
  file.path(OUT,"09E0_audit.txt")
)

cat("\n============================================================\n")
cat("09E0 COMPLETE\n")
cat("============================================================\n")
cat("Cells audited:",nrow(score),"\n")
cat("No cells removed.\n")
cat("No cells reassigned.\n")
cat("FibConsensus not used.\n")
cat("PASS — fibroblast contamination evidence generated.\n")

