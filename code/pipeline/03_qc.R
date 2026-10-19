# =============================================================================
# PIPELINE_V2 — 03_qc.R
#
# PURPOSE
#   Reproduce the historical A01 QC exactly.
#
# INPUT
#   objects/checkpoints/02_preQC_objects.rds
#
# FILTERS
#   nFeature_RNA >= 200
#   percent.mt < 50
#
# IMPORTANT
#   percent.mt and percent.ribo are calculated BEFORE filtering,
#   reproducing historical A00.
#
# NO:
#   doublet removal
#   normalization
#   integration
#   clustering
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

options(stringsAsFactors = FALSE)

set.seed(20260928)

source("PIPELINE_V2/config/paths.R")

cat("\n")
cat("============================================================\n")
cat("PIPELINE_V2 — 03 QC\n")
cat("============================================================\n\n")

# =============================================================================
# PATHS
# =============================================================================

IN_RDS <- file.path(
  V2_CHECKPOINTS,
  "02_preQC_objects.rds"
)

OUT_RDS <- file.path(
  V2_CHECKPOINTS,
  "03_postQC_objects.rds"
)

OUT_TABLE <- file.path(
  V2_TABLES,
  "03_QC_summary.csv"
)

QC_DIR <- file.path(
  V2_QC,
  "cell_filtering"
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_AUDIT <- file.path(
  QC_DIR,
  "03_QC_audit.txt"
)

THRESHOLD_FILE <- file.path(
  V2_ROOT,
  "config",
  "qc_thresholds.csv"
)

# Historical summary — used ONLY for reproducibility validation
HISTORICAL_QC <- file.path(
  LEGACY_SCRNA,
  "01_QC",
  "A01_QC_summary.tsv"
)

# =============================================================================
# LOAD
# =============================================================================

if (!file.exists(IN_RDS)) {
  stop(
    "Input checkpoint not found: ",
    IN_RDS
  )
}

objs <- readRDS(IN_RDS)

stopifnot(
  length(objs) == 16
)

# =============================================================================
# LOAD THRESHOLDS
# =============================================================================

if (!file.exists(THRESHOLD_FILE)) {
  stop(
    "QC threshold file not found: ",
    THRESHOLD_FILE
  )
}

thr <- read.csv(
  THRESHOLD_FILE,
  stringsAsFactors = FALSE
)

get_threshold <- function(name) {
  
  x <- thr$value[
    thr$parameter == name
  ]
  
  if (length(x) != 1) {
    stop(
      "Threshold must occur exactly once: ",
      name
    )
  }
  
  as.numeric(x)
}

MIN_FEATURES <- get_threshold(
  "min_features"
)

MAX_PERCENT_MT <- get_threshold(
  "max_percent_mt"
)

cat("QC thresholds:\n")
cat(
  "  nFeature_RNA >= ",
  MIN_FEATURES,
  "\n",
  sep = ""
)
cat(
  "  percent.mt < ",
  MAX_PERCENT_MT,
  "\n\n",
  sep = ""
)

# =============================================================================
# QC
# =============================================================================

results <- vector(
  "list",
  length(objs)
)

qc_objects <- vector(
  "list",
  length(objs)
)

names(qc_objects) <- names(objs)

for (i in seq_along(objs)) {
  
  nm <- names(objs)[i]
  
  x <- objs[[i]]
  
  # ---------------------------------------------------------------------------
  # Reproduce metrics calculated in historical A00.
  # ---------------------------------------------------------------------------
  
  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "^MT-"
  )
  
  x[["percent.ribo"]] <- PercentageFeatureSet(
    x,
    pattern = "^RP[SL]"
  )
  
  before <- ncol(x)
  
  # ---------------------------------------------------------------------------
  # Reproduce historical A01 exactly.
  # ---------------------------------------------------------------------------
  
  low_feature <- (
    x$nFeature_RNA < MIN_FEATURES
  )
  
  high_mt <- (
    x$percent.mt >= MAX_PERCENT_MT
  )
  
  keep <- (
    !low_feature &
      !high_mt
  )
  
  # Preserve QC decisions as metadata.
  x$QC_low_feature <- low_feature
  x$QC_high_mt <- high_mt
  x$QC_pass <- keep
  
  # ---------------------------------------------------------------------------
  # Summary before subsetting
  # ---------------------------------------------------------------------------
  
  results[[i]] <- data.frame(
    
    patient_id =
      unique(x$patient_id),
    
    disease =
      unique(x$disease),
    
    batch =
      unique(x$batch),
    
    cells_preQC =
      before,
    
    removed_low_feature =
      sum(low_feature),
    
    removed_high_mt =
      sum(high_mt),
    
    removed_both =
      sum(
        low_feature &
          high_mt
      ),
    
    removed_union =
      sum(!keep),
    
    cells_postQC =
      sum(keep),
    
    retained_pct =
      100 * sum(keep) / before,
    
    median_genes_pre =
      median(
        x$nFeature_RNA
      ),
    
    median_genes_post =
      median(
        x$nFeature_RNA[keep]
      ),
    
    median_UMI_pre =
      median(
        x$nCount_RNA
      ),
    
    median_UMI_post =
      median(
        x$nCount_RNA[keep]
      ),
    
    median_mt_pre =
      median(
        x$percent.mt
      ),
    
    median_mt_post =
      median(
        x$percent.mt[keep]
      ),
    
    median_ribo_pre =
      median(
        x$percent.ribo
      ),
    
    median_ribo_post =
      median(
        x$percent.ribo[keep]
      ),
    
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Subset
  # ---------------------------------------------------------------------------
  
  y <- x[, keep]
  
  y$pipeline_stage <- "V2_postQC"
  
  qc_objects[[i]] <- y
  
  cat(
    sprintf(
      "%-7s %-3s %-12s %5d -> %5d (%5.1f%% retained)\n",
      nm,
      unique(x$disease),
      unique(x$batch),
      before,
      ncol(y),
      100 * ncol(y) / before
    )
  )
  
  rm(x, y)
  gc(verbose = FALSE)
}

tab <- do.call(
  rbind,
  results
)

rownames(tab) <- NULL

# =============================================================================
# GLOBAL CHECKS
# =============================================================================

PRE_EXPECTED <- 99222L
POST_HISTORICAL_EXPECTED <- 89668L

pre_total <- sum(
  tab$cells_preQC
)

post_total <- sum(
  tab$cells_postQC
)

removed_total <- sum(
  tab$removed_union
)

stopifnot(
  pre_total == PRE_EXPECTED
)

cat("\n")
cat("Pre-QC cells: ", pre_total, "\n", sep = "")
cat("Post-QC cells: ", post_total, "\n", sep = "")
cat(
  "Removed: ",
  removed_total,
  " (",
  sprintf(
    "%.2f",
    100 * removed_total / pre_total
  ),
  "%)\n",
  sep = ""
)

# =============================================================================
# HISTORICAL A01 COMPARISON
# =============================================================================

tab$historical_postQC <- NA_integer_
tab$postQC_delta <- NA_integer_
tab$historical_match <- NA

historical_available <- file.exists(
  HISTORICAL_QC
)

if (historical_available) {
  
  hist <- read.delim(
    HISTORICAL_QC,
    check.names = FALSE
  )
  
  required_hist <- c(
    "patient_id",
    "cells_postQC"
  )
  
  if (
    all(
      required_hist %in%
      colnames(hist)
    )
  ) {
    
    idx <- match(
      tab$patient_id,
      hist$patient_id
    )
    
    tab$historical_postQC <-
      hist$cells_postQC[idx]
    
    tab$postQC_delta <-
      tab$cells_postQC -
      tab$historical_postQC
    
    tab$historical_match <-
      tab$postQC_delta == 0
    
  } else {
    
    warning(
      "Historical A01 table does not contain expected columns."
    )
  }
}

# =============================================================================
# WRITE TABLE BEFORE HARD STOP
# =============================================================================

write.csv(
  tab,
  OUT_TABLE,
  row.names = FALSE,
  quote = TRUE
)

# =============================================================================
# PRINT COMPARISON
# =============================================================================

cat("\nPER PATIENT\n")

print(
  tab[
    ,
    c(
      "patient_id",
      "disease",
      "batch",
      "cells_preQC",
      "cells_postQC",
      "retained_pct",
      "median_genes_post",
      "median_mt_post"
    )
  ],
  row.names = FALSE
)

historical_test <- any(
  !is.na(
    tab$historical_match
  )
)

if (historical_test) {
  
  cat(
    "\nHISTORICAL A01 REPRODUCIBILITY\n"
  )
  
  print(
    tab[
      ,
      c(
        "patient_id",
        "cells_postQC",
        "historical_postQC",
        "postQC_delta",
        "historical_match"
      )
    ],
    row.names = FALSE
  )
  
  if (
    any(
      is.na(
        tab$historical_match
      )
    ) ||
    !all(
      tab$historical_match
    )
  ) {
    
    stop(
      "V2 QC does not exactly reproduce historical A01."
    )
  }
}

# =============================================================================
# HARD REPRODUCIBILITY TARGET
# =============================================================================

if (
  post_total !=
  POST_HISTORICAL_EXPECTED
) {
  
  stop(
    "Expected historical post-QC total = ",
    POST_HISTORICAL_EXPECTED,
    "; obtained = ",
    post_total
  )
}

# =============================================================================
# DISEASE / BATCH SUMMARIES
# =============================================================================

cat("\nPOST-QC BY DISEASE\n")

print(
  aggregate(
    cells_postQC ~ disease,
    data = tab,
    sum
  ),
  row.names = FALSE
)

cat("\nRETENTION BY BATCH\n")

print(
  aggregate(
    cbind(
      cells_preQC,
      cells_postQC
    ) ~ batch,
    data = tab,
    sum
  ),
  row.names = FALSE
)

# =============================================================================
# SAVE ONLY AFTER ALL REPRODUCIBILITY CHECKS PASS
# =============================================================================

saveRDS(
  qc_objects,
  OUT_RDS,
  compress = FALSE
)

# =============================================================================
# AUDIT REPORT
# =============================================================================

audit <- c(
  
  "PIPELINE_V2 — 03 QC audit",
  
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste0(
    "Pre-QC cells: ",
    pre_total
  ),
  
  paste0(
    "Post-QC cells: ",
    post_total
  ),
  
  paste0(
    "Removed cells: ",
    removed_total
  ),
  
  paste0(
    "Removed percent: ",
    sprintf(
      "%.4f",
      100 * removed_total / pre_total
    )
  ),
  
  "",
  
  paste0(
    "Threshold nFeature_RNA >= ",
    MIN_FEATURES
  ),
  
  paste0(
    "Threshold percent.mt < ",
    MAX_PERCENT_MT
  ),
  
  "",
  
  "percent.mt pattern: ^MT-",
  "percent.ribo pattern: ^RP[SL]",
  
  "",
  
  paste0(
    "Historical A01 table available: ",
    historical_available
  ),
  
  paste0(
    "Historical patient-level comparison performed: ",
    historical_test
  ),
  
  paste0(
    "Historical total target: ",
    POST_HISTORICAL_EXPECTED
  ),
  
  paste0(
    "Historical total reproduced: ",
    post_total ==
      POST_HISTORICAL_EXPECTED
  ),
  
  "",
  
  paste0(
    "Checkpoint: ",
    OUT_RDS
  ),
  
  paste0(
    "Summary: ",
    OUT_TABLE
  )
)

writeLines(
  audit,
  OUT_AUDIT
)

# =============================================================================
# SUCCESS
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("V2 QC COMPLETE — EXACT HISTORICAL REPRODUCTION\n")
cat("============================================================\n")

cat(
  "99,222 pre-QC -> ",
  format(
    post_total,
    big.mark = ","
  ),
  " post-QC cells\n",
  sep = ""
)

cat(
  "Historical A01 patient-level counts: EXACT MATCH\n"
)

cat(
  "\nCheckpoint:\n  ",
  OUT_RDS,
  "\n"
)

cat(
  "Summary:\n  ",
  OUT_TABLE,
  "\n"
)

cat(
  "Audit:\n  ",
  OUT_AUDIT,
  "\n"
)

cat(
  "\nNo doublet removal was performed.\n"
)