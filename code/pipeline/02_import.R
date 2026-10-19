# =============================================================================
# PIPELINE_V2 — 02_import.R
#
# PURPOSE
#   Import the 16 frozen Cell Ranger filtered matrices into Seurat.
#
# INPUT
#   PIPELINE_V2/objects/frozen/scrnaseq_provenance_manifest.csv
#
# OUTPUT
#   PIPELINE_V2/objects/checkpoints/02_preQC_objects.rds
#   PIPELINE_V2/results/tables/02_import_summary.csv
#   PIPELINE_V2/results/QC/import/02_import_audit.txt
#
# IMPORTANT
#   - NO QC filtering
#   - NO doublet removal
#   - NO normalization
#   - NO integration
#   - NO clustering
#   - NO annotation
#
# REPRODUCIBILITY
#   The imported objects are compared with the historical A00 import whenever
#   the historical manifest is available.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

set.seed(20260928)

source("PIPELINE_V2/config/paths.R")

cat("\n")
cat("============================================================\n")
cat("PIPELINE_V2 — 02 IMPORT\n")
cat("============================================================\n\n")

# =============================================================================
# PATHS
# =============================================================================

PROVENANCE <- file.path(
  V2_FROZEN,
  "scrnaseq_provenance_manifest.csv"
)

OUT_RDS <- file.path(
  V2_CHECKPOINTS,
  "02_preQC_objects.rds"
)

OUT_TABLE <- file.path(
  V2_TABLES,
  "02_import_summary.csv"
)

QC_DIR <- file.path(
  V2_QC,
  "import"
)

dir.create(
  V2_CHECKPOINTS,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  V2_TABLES,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_AUDIT <- file.path(
  QC_DIR,
  "02_import_audit.txt"
)

# Historical A00 manifest for reproducibility check
HISTORICAL_MANIFEST <- file.path(
  LEGACY_SCRNA,
  "00_inputs",
  "A00_16_sample_manifest.tsv"
)

# =============================================================================
# LOAD FROZEN PROVENANCE
# =============================================================================

if (!file.exists(PROVENANCE)) {
  stop(
    "Frozen provenance manifest not found:\n",
    PROVENANCE
  )
}

m <- read.csv(
  PROVENANCE,
  check.names = FALSE
)

required <- c(
  "source_sample_id",
  "sample_id",
  "patient_id",
  "disease",
  "batch",
  "source_dataset",
  "matrix_path",
  "provenance_status"
)

missing_cols <- setdiff(
  required,
  colnames(m)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing provenance columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (!all(m$provenance_status == "PASS")) {
  stop(
    "Not all provenance records are PASS."
  )
}

stopifnot(
  nrow(m) == 16,
  length(unique(m$sample_id)) == 16,
  length(unique(m$patient_id)) == 16
)

# =============================================================================
# IMPORT
# =============================================================================

objects <- vector(
  "list",
  nrow(m)
)

names(objects) <- m$patient_id

summary_rows <- vector(
  "list",
  nrow(m)
)

for (i in seq_len(nrow(m))) {
  
  z <- m[i, ]
  
  cat(
    sprintf(
      "[%02d/16] %-7s %-3s %-12s ",
      i,
      z$patient_id,
      z$disease,
      z$batch
    )
  )
  
  # ---------------------------------------------------------------------------
  # Matrix existence
  # ---------------------------------------------------------------------------
  
  if (!dir.exists(z$matrix_path)) {
    stop(
      "\nMatrix directory disappeared since provenance freeze:\n",
      z$matrix_path
    )
  }
  
  # ---------------------------------------------------------------------------
  # Read exactly as historical A00:
  # Read10X(filtered_feature_bc_matrix, gene.column=2, unique.features=TRUE)
  # ---------------------------------------------------------------------------
  
  mat <- Read10X(
    data.dir = z$matrix_path,
    gene.column = 2,
    unique.features = TRUE
  )
  
  # Cell Ranger can return a list if multiple feature types are present.
  if (is.list(mat)) {
    
    if (!"Gene Expression" %in% names(mat)) {
      stop(
        "\nRead10X returned multiple modalities but no 'Gene Expression' matrix for ",
        z$patient_id
      )
    }
    
    mat <- mat[["Gene Expression"]]
  }
  
  if (!inherits(mat, "sparseMatrix")) {
    warning(
      "Matrix for ",
      z$patient_id,
      " is not a sparseMatrix."
    )
  }
  
  n_genes_raw <- nrow(mat)
  n_cells_raw <- ncol(mat)
  
  # ---------------------------------------------------------------------------
  # Create Seurat object — deliberately NO QC thresholds.
  # ---------------------------------------------------------------------------
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = z$patient_id,
    min.cells = 0,
    min.features = 0
  )
  
  # ---------------------------------------------------------------------------
  # Stable globally unique cell IDs
  # ---------------------------------------------------------------------------
  
  obj <- RenameCells(
    object = obj,
    add.cell.id = z$patient_id
  )
  
  # ---------------------------------------------------------------------------
  # Canonical metadata
  # ---------------------------------------------------------------------------
  
  obj$source_sample_id <- z$source_sample_id
  obj$sample_id        <- z$sample_id
  obj$patient_id       <- z$patient_id
  obj$disease          <- z$disease
  obj$batch            <- z$batch
  obj$source_dataset   <- z$source_dataset
  
  # Preserve import stage explicitly
  obj$pipeline_stage <- "V2_preQC"
  
  # ---------------------------------------------------------------------------
  # Integrity checks
  # ---------------------------------------------------------------------------
  
  stopifnot(
    ncol(obj) == n_cells_raw
  )
  
  stopifnot(
    nrow(obj) == n_genes_raw
  )
  
  stopifnot(
    length(unique(obj$patient_id)) == 1,
    unique(obj$patient_id) == z$patient_id
  )
  
  stopifnot(
    length(unique(obj$disease)) == 1,
    unique(obj$disease) == z$disease
  )
  
  # Cell names must begin with canonical patient ID
  expected_prefix <- paste0(
    z$patient_id,
    "_"
  )
  
  if (!all(
    startsWith(
      colnames(obj),
      expected_prefix
    )
  )) {
    stop(
      "Cell renaming failed for ",
      z$patient_id
    )
  }
  
  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------
  
  summary_rows[[i]] <- data.frame(
    
    source_sample_id =
      z$source_sample_id,
    
    sample_id =
      z$sample_id,
    
    patient_id =
      z$patient_id,
    
    disease =
      z$disease,
    
    batch =
      z$batch,
    
    source_dataset =
      z$source_dataset,
    
    n_genes =
      nrow(obj),
    
    n_cells_preQC =
      ncol(obj),
    
    matrix_path =
      z$matrix_path,
    
    stringsAsFactors = FALSE
  )
  
  objects[[i]] <- obj
  
  cat(
    sprintf(
      "%d genes | %d cells\n",
      nrow(obj),
      ncol(obj)
    )
  )
  
  # Release matrix reference
  rm(mat)
  gc(verbose = FALSE)
}

summary_df <- do.call(
  rbind,
  summary_rows
)

# =============================================================================
# GLOBAL INTEGRITY
# =============================================================================

stopifnot(
  length(objects) == 16
)

stopifnot(
  length(unique(names(objects))) == 16
)

all_cells <- unlist(
  lapply(
    objects,
    colnames
  ),
  use.names = FALSE
)

if (anyDuplicated(all_cells)) {
  stop(
    "Duplicated cell IDs detected across imported objects."
  )
}

total_cells <- sum(
  summary_df$n_cells_preQC
)

cat("\n")
cat("Total pre-QC cells:", total_cells, "\n")

# =============================================================================
# HISTORICAL A00 REPRODUCIBILITY CHECK
# =============================================================================

summary_df$historical_n_cells <- NA_integer_
summary_df$cell_count_delta <- NA_integer_
summary_df$historical_match <- NA

historical_available <- file.exists(
  HISTORICAL_MANIFEST
)

if (historical_available) {
  
  cat(
    "\nHistorical A00 manifest found.\n"
  )
  
  hist <- read.delim(
    HISTORICAL_MANIFEST,
    check.names = FALSE
  )
  
  # Determine historical cell-count column robustly.
  candidate_count_cols <- c(
    "n_cells",
    "cells",
    "ncell",
    "n_cells_raw",
    "n_cells_preQC"
  )
  
  count_col <- intersect(
    candidate_count_cols,
    colnames(hist)
  )
  
  if (length(count_col) > 0) {
    
    count_col <- count_col[[1]]
    
    hist_patient_col <- if (
      "patient_id" %in% colnames(hist)
    ) {
      "patient_id"
    } else if (
      "sample_id" %in% colnames(hist)
    ) {
      "sample_id"
    } else {
      NA_character_
    }
    
    if (!is.na(hist_patient_col)) {
      
      historical_counts <- data.frame(
        patient_id =
          as.character(
            hist[[hist_patient_col]]
          ),
        
        historical_n_cells =
          as.integer(
            hist[[count_col]]
          ),
        
        stringsAsFactors = FALSE
      )
      
      idx <- match(
        summary_df$patient_id,
        historical_counts$patient_id
      )
      
      summary_df$historical_n_cells <-
        historical_counts$historical_n_cells[idx]
      
      summary_df$cell_count_delta <-
        summary_df$n_cells_preQC -
        summary_df$historical_n_cells
      
      summary_df$historical_match <-
        summary_df$cell_count_delta == 0
      
    } else {
      
      warning(
        "Historical manifest has no recognizable patient/sample ID column."
      )
    }
    
  } else {
    
    warning(
      paste0(
        "Historical manifest exists but no recognized cell-count column was found. ",
        "Columns: ",
        paste(
          colnames(hist),
          collapse = ", "
        )
      )
    )
  }
  
} else {
  
  cat(
    "\nHistorical A00 manifest not found; ",
    "import will be validated internally only.\n"
  )
}

# =============================================================================
# WRITE SUMMARY BEFORE HARD HISTORICAL CHECK
# =============================================================================

write.csv(
  summary_df,
  OUT_TABLE,
  row.names = FALSE,
  quote = TRUE
)

# =============================================================================
# HISTORICAL HARD CHECK, WHEN COUNTS ARE AVAILABLE
# =============================================================================

historical_test_performed <- any(
  !is.na(
    summary_df$historical_match
  )
)

historical_all_match <- NA

if (historical_test_performed) {
  
  historical_all_match <- all(
    summary_df$historical_match,
    na.rm = FALSE
  )
  
  cat("\nHistorical import comparison:\n")
  
  print(
    summary_df[
      ,
      c(
        "patient_id",
        "n_cells_preQC",
        "historical_n_cells",
        "cell_count_delta",
        "historical_match"
      )
    ],
    row.names = FALSE
  )
  
  if (!historical_all_match) {
    
    bad <- summary_df[
      is.na(summary_df$historical_match) |
        !summary_df$historical_match,
      ,
      drop = FALSE
    ]
    
    cat(
      "\nNON-MATCHING SAMPLES:\n"
    )
    
    print(
      bad,
      row.names = FALSE
    )
    
    stop(
      "V2 import does not exactly reproduce historical A00 cell counts."
    )
  }
}

# =============================================================================
# SAVE CHECKPOINT
# =============================================================================

saveRDS(
  objects,
  OUT_RDS,
  compress = FALSE
)

# =============================================================================
# AUDIT REPORT
# =============================================================================

audit <- c(
  
  "PIPELINE_V2 — 02 import audit",
  
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste0(
    "Samples imported: ",
    length(objects)
  ),
  
  paste0(
    "Individuals: ",
    length(
      unique(
        summary_df$patient_id
      )
    )
  ),
  
  paste0(
    "Total pre-QC cells: ",
    total_cells
  ),
  
  "",
  
  paste0(
    "HC samples: ",
    sum(
      summary_df$disease == "HC"
    )
  ),
  
  paste0(
    "AA samples: ",
    sum(
      summary_df$disease == "AA"
    )
  ),
  
  paste0(
    "FFA samples: ",
    sum(
      summary_df$disease == "FFA"
    )
  ),
  
  "",
  
  paste0(
    "Historical A00 manifest available: ",
    historical_available
  ),
  
  paste0(
    "Historical cell-count test performed: ",
    historical_test_performed
  ),
  
  paste0(
    "Historical counts all match: ",
    historical_all_match
  ),
  
  "",
  
  "Import method:",
  "Read10X(filtered_feature_bc_matrix, gene.column=2, unique.features=TRUE)",
  "CreateSeuratObject(min.cells=0, min.features=0)",
  "RenameCells(add.cell.id=patient_id)",
  
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
cat("V2 IMPORT COMPLETE\n")
cat("============================================================\n")

cat(
  "Objects:",
  length(objects),
  "\n"
)

cat(
  "Total pre-QC cells:",
  total_cells,
  "\n"
)

if (historical_test_performed) {
  
  cat(
    "Historical A00 counts: EXACT MATCH\n"
  )
}

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
  "\nNo QC filtering was performed.\n"
)