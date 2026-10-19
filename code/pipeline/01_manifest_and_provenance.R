# =============================================================================
# PIPELINE_V2 — 01_manifest_and_provenance.R
#
# PURPOSE
#   Resolve and validate the 16 scRNA-seq Cell Ranger inputs.
#
# DOES
#   - reads config/sample_manifest.csv
#   - constructs source paths
#   - distinguishes symlink path from resolved physical target
#   - validates filtered_feature_bc_matrix
#   - validates matrix/features/barcodes files
#   - records optional H5 + metrics_summary.csv
#   - records file sizes and modification times
#   - freezes a canonical provenance manifest
#
# DOES NOT
#   - read expression matrices
#   - create Seurat objects
#   - perform QC
#   - modify raw data
# =============================================================================

options(stringsAsFactors = FALSE)

source("PIPELINE_V2/config/paths.R")

cat("\n")
cat("============================================================\n")
cat("PIPELINE_V2 — 01 MANIFEST AND PROVENANCE\n")
cat("============================================================\n\n")

# =============================================================================
# PATHS
# =============================================================================

CONFIG_MANIFEST <- file.path(
  V2_ROOT,
  "config",
  "sample_manifest.csv"
)

OUTDIR <- file.path(
  V2_ROOT,
  "objects",
  "frozen"
)

QCDIR <- file.path(
  V2_ROOT,
  "results",
  "QC",
  "provenance"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QCDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_MANIFEST <- file.path(
  OUTDIR,
  "scrnaseq_provenance_manifest.csv"
)

OUT_TSV <- file.path(
  OUTDIR,
  "scrnaseq_provenance_manifest.tsv"
)

OUT_AUDIT <- file.path(
  QCDIR,
  "01_provenance_audit.txt"
)

# =============================================================================
# LOAD CONFIG MANIFEST
# =============================================================================

if (!file.exists(CONFIG_MANIFEST)) {
  stop(
    "Config manifest not found: ",
    CONFIG_MANIFEST
  )
}

m <- read.csv(
  CONFIG_MANIFEST,
  check.names = FALSE
)

required_cols <- c(
  "source_sample_id",
  "sample_id",
  "patient_id",
  "disease",
  "batch",
  "source_dataset"
)

missing_cols <- setdiff(
  required_cols,
  colnames(m)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required manifest columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

# =============================================================================
# BASIC COHORT VALIDATION
# =============================================================================

stopifnot(
  nrow(m) == 16
)

stopifnot(
  length(unique(m$sample_id)) == 16
)

stopifnot(
  length(unique(m$patient_id)) == 16
)

stopifnot(
  sum(m$disease == "HC") == 5
)

stopifnot(
  sum(m$disease == "AA") == 6
)

stopifnot(
  sum(m$disease == "FFA") == 5
)

allowed_disease <- c(
  "HC",
  "AA",
  "FFA"
)

if (!all(m$disease %in% allowed_disease)) {
  stop(
    "Unexpected disease label(s): ",
    paste(
      setdiff(
        unique(m$disease),
        allowed_disease
      ),
      collapse = ", "
    )
  )
}

# =============================================================================
# SOURCE ROOTS
# =============================================================================

HISTORICAL_ROOT <- file.path(
  RAW_DATA,
  "scrnaseq",
  "01_iScience_GSE314972",
  "cellranger_links"
)

CABIMER_ROOT <- file.path(
  RAW_DATA,
  "scrnaseq",
  "02_CABIMER_2026",
  "cellranger"
)

if (!dir.exists(HISTORICAL_ROOT)) {
  stop(
    "Historical Cell Ranger root not found: ",
    HISTORICAL_ROOT
  )
}

if (!dir.exists(CABIMER_ROOT)) {
  stop(
    "CABIMER Cell Ranger root not found: ",
    CABIMER_ROOT
  )
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

resolve_source_root <- function(
    batch,
    source_sample_id
) {
  
  root <- switch(
    batch,
    
    historical = HISTORICAL_ROOT,
    
    CABIMER2026 = CABIMER_ROOT,
    
    stop(
      "Unknown batch: ",
      batch
    )
  )
  
  file.path(
    root,
    source_sample_id
  )
}


first_existing <- function(paths) {
  
  hit <- paths[
    file.exists(paths)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[[1]]
}


first_existing_dir <- function(paths) {
  
  hit <- paths[
    dir.exists(paths)
  ]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[[1]]
}


file_size_safe <- function(path) {
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    return(NA_real_)
  }
  
  as.numeric(
    file.info(path)$size
  )
}


mtime_safe <- function(path) {
  
  if (
    is.na(path) ||
    !file.exists(path)
  ) {
    return(NA_character_)
  }
  
  format(
    file.info(path)$mtime,
    "%Y-%m-%d %H:%M:%S"
  )
}


resolve_symlink_target <- function(path) {
  
  link <- Sys.readlink(path)
  
  if (
    length(link) == 0 ||
    is.na(link) ||
    link == ""
  ) {
    
    # Not a symlink.
    return(
      normalizePath(
        path,
        mustWork = FALSE
      )
    )
  }
  
  # Absolute symlink
  if (
    grepl(
      "^/",
      link
    )
  ) {
    
    return(
      normalizePath(
        link,
        mustWork = FALSE
      )
    )
  }
  
  # Relative symlink
  normalizePath(
    file.path(
      dirname(path),
      link
    ),
    mustWork = FALSE
  )
}

# =============================================================================
# RESOLVE ALL 16 INPUTS
# =============================================================================

rows <- vector(
  "list",
  nrow(m)
)

for (i in seq_len(nrow(m))) {
  
  z <- m[i, ]
  
  source_root <- resolve_source_root(
    batch = z$batch,
    source_sample_id = z$source_sample_id
  )
  
  source_exists <- file.exists(
    source_root
  ) || dir.exists(
    source_root
  )
  
  source_is_symlink <- (
    Sys.readlink(source_root) != ""
  )
  
  resolved_target <- if (
    source_exists
  ) {
    
    resolve_symlink_target(
      source_root
    )
    
  } else {
    
    NA_character_
  }
  
  # ---------------------------------------------------------------------------
  # Canonical input for V2:
  # filtered_feature_bc_matrix directory, matching historical A00 pipeline.
  # ---------------------------------------------------------------------------
  
  filtered_dir <- first_existing_dir(
    c(
      file.path(
        source_root,
        "filtered_feature_bc_matrix"
      ),
      
      file.path(
        source_root,
        "outs",
        "filtered_feature_bc_matrix"
      )
    )
  )
  
  raw_dir <- first_existing_dir(
    c(
      file.path(
        source_root,
        "raw_feature_bc_matrix"
      ),
      
      file.path(
        source_root,
        "outs",
        "raw_feature_bc_matrix"
      )
    )
  )
  
  filtered_h5 <- first_existing(
    c(
      file.path(
        source_root,
        "filtered_feature_bc_matrix.h5"
      ),
      
      file.path(
        source_root,
        "outs",
        "filtered_feature_bc_matrix.h5"
      )
    )
  )
  
  raw_h5 <- first_existing(
    c(
      file.path(
        source_root,
        "raw_feature_bc_matrix.h5"
      ),
      
      file.path(
        source_root,
        "outs",
        "raw_feature_bc_matrix.h5"
      )
    )
  )
  
  metrics_summary <- first_existing(
    c(
      file.path(
        source_root,
        "metrics_summary.csv"
      ),
      
      file.path(
        source_root,
        "outs",
        "metrics_summary.csv"
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Validate MEX files.
  # Support gzipped and uncompressed Cell Ranger outputs.
  # ---------------------------------------------------------------------------
  
  matrix_file <- if (
    !is.na(filtered_dir)
  ) {
    
    first_existing(
      c(
        file.path(
          filtered_dir,
          "matrix.mtx.gz"
        ),
        
        file.path(
          filtered_dir,
          "matrix.mtx"
        )
      )
    )
    
  } else {
    
    NA_character_
  }
  
  features_file <- if (
    !is.na(filtered_dir)
  ) {
    
    first_existing(
      c(
        file.path(
          filtered_dir,
          "features.tsv.gz"
        ),
        
        file.path(
          filtered_dir,
          "features.tsv"
        ),
        
        # Older 10x nomenclature
        file.path(
          filtered_dir,
          "genes.tsv.gz"
        ),
        
        file.path(
          filtered_dir,
          "genes.tsv"
        )
      )
    )
    
  } else {
    
    NA_character_
  }
  
  barcodes_file <- if (
    !is.na(filtered_dir)
  ) {
    
    first_existing(
      c(
        file.path(
          filtered_dir,
          "barcodes.tsv.gz"
        ),
        
        file.path(
          filtered_dir,
          "barcodes.tsv"
        )
      )
    )
    
  } else {
    
    NA_character_
  }
  
  # ---------------------------------------------------------------------------
  # PASS / FAIL
  # ---------------------------------------------------------------------------
  
  required_ok <- all(
    c(
      source_exists,
      !is.na(filtered_dir),
      !is.na(matrix_file),
      !is.na(features_file),
      !is.na(barcodes_file)
    )
  )
  
  status <- if (
    required_ok
  ) {
    "PASS"
  } else {
    "FAIL"
  }
  
  rows[[i]] <- data.frame(
    
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
    
    source_root =
      source_root,
    
    source_exists =
      source_exists,
    
    source_is_symlink =
      source_is_symlink,
    
    resolved_target =
      resolved_target,
    
    input_stage =
      "cellranger_filtered",
    
    input_format =
      "10x_MEX",
    
    matrix_path =
      filtered_dir,
    
    matrix_file =
      matrix_file,
    
    features_file =
      features_file,
    
    barcodes_file =
      barcodes_file,
    
    filtered_h5 =
      filtered_h5,
    
    raw_matrix_path =
      raw_dir,
    
    raw_h5 =
      raw_h5,
    
    metrics_summary =
      metrics_summary,
    
    matrix_size_bytes =
      file_size_safe(
        matrix_file
      ),
    
    features_size_bytes =
      file_size_safe(
        features_file
      ),
    
    barcodes_size_bytes =
      file_size_safe(
        barcodes_file
      ),
    
    matrix_mtime =
      mtime_safe(
        matrix_file
      ),
    
    provenance_status =
      status,
    
    stringsAsFactors = FALSE
  )
}

prov <- do.call(
  rbind,
  rows
)

# =============================================================================
# FINAL VALIDATION
# =============================================================================

cat("\nResolved provenance:\n\n")

print(
  prov[
    ,
    c(
      "source_sample_id",
      "sample_id",
      "disease",
      "batch",
      "source_is_symlink",
      "provenance_status"
    )
  ],
  row.names = FALSE
)

cat("\n")

n_pass <- sum(
  prov$provenance_status == "PASS"
)

n_fail <- sum(
  prov$provenance_status == "FAIL"
)

cat(
  "PASS:",
  n_pass,
  "/",
  nrow(prov),
  "\n"
)

cat(
  "FAIL:",
  n_fail,
  "/",
  nrow(prov),
  "\n"
)

# =============================================================================
# WRITE FROZEN MANIFEST
# =============================================================================

write.csv(
  prov,
  OUT_MANIFEST,
  row.names = FALSE,
  quote = TRUE
)

write.table(
  prov,
  OUT_TSV,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =============================================================================
# HUMAN-READABLE AUDIT
# =============================================================================

audit_lines <- c(
  
  "PIPELINE_V2 — scRNA-seq provenance audit",
  
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  
  "",
  
  paste0(
    "Samples: ",
    nrow(prov)
  ),
  
  paste0(
    "Individuals: ",
    length(
      unique(
        prov$patient_id
      )
    )
  ),
  
  paste0(
    "HC: ",
    sum(
      prov$disease == "HC"
    )
  ),
  
  paste0(
    "AA: ",
    sum(
      prov$disease == "AA"
    )
  ),
  
  paste0(
    "FFA: ",
    sum(
      prov$disease == "FFA"
    )
  ),
  
  "",
  
  paste0(
    "PASS: ",
    n_pass
  ),
  
  paste0(
    "FAIL: ",
    n_fail
  ),
  
  "",
  
  "Canonical V2 expression input:",
  "Cell Ranger filtered_feature_bc_matrix (10x MEX).",
  
  "",
  
  "Source ID normalization:",
  
  paste(
    paste0(
      prov$source_sample_id,
      " -> ",
      prov$sample_id
    ),
    collapse = "\n"
  ),
  
  "",
  
  "Symlink targets:",
  
  paste(
    paste0(
      prov$sample_id,
      " | ",
      prov$source_root,
      " -> ",
      prov$resolved_target
    ),
    collapse = "\n"
  )
)

writeLines(
  audit_lines,
  OUT_AUDIT
)

# =============================================================================
# HARD STOP IF ANY SAMPLE FAILS
# =============================================================================

if (n_fail > 0) {
  
  cat(
    "\nFAILED samples:\n"
  )
  
  print(
    prov[
      prov$provenance_status == "FAIL",
      ,
      drop = FALSE
    ]
  )
  
  stop(
    "Provenance validation failed. ",
    n_fail,
    " sample(s) are incomplete."
  )
}

# =============================================================================
# SUCCESS
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("PROVENANCE VALIDATION PASSED\n")
cat("============================================================\n")

cat(
  "16/16 samples resolved.\n"
)

cat(
  "Frozen manifest:\n  ",
  OUT_MANIFEST,
  "\n"
)

cat(
  "Audit:\n  ",
  OUT_AUDIT,
  "\n"
)

cat(
  "\nNo expression matrix was loaded and no raw file was modified.\n"
)