# =============================================================================
# PIPELINE_V2 — 04_doublets.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
})

options(stringsAsFactors = FALSE)


SEED <- 20260924
set.seed(SEED)

source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 04 DOUBLETS\n")
cat("============================================================\n\n")

IN_RDS <- file.path(V2_CHECKPOINTS, "03_postQC_objects.rds")
OUT_RDS <- file.path(V2_CHECKPOINTS, "04_singlet_objects.rds")
OUT_TABLE <- file.path(V2_TABLES, "04_doublet_summary.csv")

QC_DIR <- file.path(V2_QC, "doublets")
dir.create(QC_DIR, recursive=TRUE, showWarnings=FALSE)

OUT_AUDIT <- file.path(QC_DIR, "04_doublet_audit.txt")

HIST_SUMMARY <- file.path(
  LEGACY_SCRNA,
  "01_QC",
  "A02_doublet_summary.tsv"
)

HIST_RDS <- file.path(
  LEGACY_SCRNA,
  "01_QC",
  "A02_16_singlet_objects.rds"
)

if (!file.exists(IN_RDS))
  stop("Missing input: ", IN_RDS)

objs <- readRDS(IN_RDS)

stopifnot(length(objs) == 16)

# =============================================================================
# SOFTWARE PROVENANCE
# =============================================================================

pkg_versions <- c(
  R = R.version.string,
  Seurat = as.character(packageVersion("Seurat")),
  SeuratObject = as.character(packageVersion("SeuratObject")),
  SingleCellExperiment = as.character(packageVersion("SingleCellExperiment")),
  scDblFinder = as.character(packageVersion("scDblFinder"))
)

cat("Software:\n")
print(pkg_versions)
cat("\n")

# =============================================================================
# RUN scDblFinder EXACTLY AS HISTORICAL A02
# =============================================================================

clean <- vector("list", length(objs))
names(clean) <- names(objs)

summ <- vector("list", length(objs))

for (i in seq_along(objs)) {
  
  nm <- names(objs)[i]
  
  cat(sprintf("[%02d/16] Processing %-7s ... ", i, nm))
  
  # Reset seed per sample for deterministic per-sample execution.
  # This is deliberately explicit in V2.
  
  x <- objs[[i]]
  
  x[["RNA"]] <- JoinLayers(x[["RNA"]])
  
  counts <- LayerData(
    x,
    assay = "RNA",
    layer = "counts"
  )
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts)
  )
  
  colData(sce)$sample <- nm
  
  sce <- scDblFinder(
    sce,
    samples = "sample",
    verbose = FALSE
  )
  
  cls <- as.character(
    sce$scDblFinder.class
  )
  
  score <- sce$scDblFinder.score
  
  if (length(cls) != ncol(x))
    stop("Classification length mismatch: ", nm)
  
  x$doublet_class <- cls
  x$doublet_score <- score
  
  keep <- cls == "singlet"
  
  y <- x[, keep]
  
  y$pipeline_stage <- "V2_singlet"
  
  clean[[i]] <- y
  
  summ[[i]] <- data.frame(
    patient_id = nm,
    disease = unique(x$disease),
    batch = unique(x$batch),
    input_cells = ncol(x),
    singlets = sum(keep),
    doublets = sum(!keep),
    doublet_pct = 100 * sum(!keep) / ncol(x),
    stringsAsFactors = FALSE
  )
  
  cat(
    sum(keep), " singlets; ",
    sum(!keep), " doublets; ",
    sprintf("%.2f%%\n", 100*sum(!keep)/ncol(x)),
    sep=""
  )
  
  rm(x, y, counts, sce)
  gc(verbose=FALSE)
}

tab <- do.call(rbind, summ)
rownames(tab) <- NULL

# =============================================================================
# HISTORICAL COUNT COMPARISON
# =============================================================================

tab$historical_singlets <- NA_integer_
tab$historical_doublets <- NA_integer_
tab$singlet_delta <- NA_integer_
tab$count_match <- NA

if (file.exists(HIST_SUMMARY)) {
  
  hist <- read.delim(
    HIST_SUMMARY,
    check.names=FALSE
  )
  
  idx <- match(
    tab$patient_id,
    hist$patient_id
  )
  
  tab$historical_singlets <- hist$singlets[idx]
  tab$historical_doublets <- hist$doublets[idx]
  
  tab$singlet_delta <-
    tab$singlets -
    tab$historical_singlets
  
  tab$count_match <-
    tab$singlets ==
    tab$historical_singlets
}

write.csv(
  tab,
  OUT_TABLE,
  row.names=FALSE
)

cat("\nHISTORICAL COUNT COMPARISON\n")
print(
  tab[, c(
    "patient_id",
    "input_cells",
    "singlets",
    "historical_singlets",
    "singlet_delta",
    "count_match"
  )],
  row.names=FALSE
)

# =============================================================================
# CELL-IDENTITY COMPARISON
# =============================================================================

identity_check_performed <- FALSE
identity_exact <- NA
identity_rows <- list()

if (file.exists(HIST_RDS)) {
  
  cat("\nComparing singlet cell identities with historical A02...\n")
  
  hist_objs <- readRDS(HIST_RDS)
  
  identity_check_performed <- TRUE
  
  for (nm in names(clean)) {
    
    new_cells <- sort(colnames(clean[[nm]]))
    old_cells <- sort(colnames(hist_objs[[nm]]))
    
    same <- identical(
      new_cells,
      old_cells
    )
    
    identity_rows[[nm]] <- data.frame(
      patient_id=nm,
      v2_singlets=length(new_cells),
      historical_singlets=length(old_cells),
      shared=length(intersect(new_cells, old_cells)),
      only_v2=length(setdiff(new_cells, old_cells)),
      only_historical=length(setdiff(old_cells, new_cells)),
      exact_identity_match=same,
      stringsAsFactors=FALSE
    )
  }
  
  identity_tab <- do.call(
    rbind,
    identity_rows
  )
  
  rownames(identity_tab) <- NULL
  
  write.csv(
    identity_tab,
    file.path(
      QC_DIR,
      "04_singlet_identity_comparison.csv"
    ),
    row.names=FALSE
  )
  
  print(identity_tab, row.names=FALSE)
  
  identity_exact <- all(
    identity_tab$exact_identity_match
  )
  
  rm(hist_objs)
  gc(verbose=FALSE)
}

# =============================================================================
# GLOBAL TARGETS
# =============================================================================

TOTAL_INPUT_EXPECTED <- 89668L
TOTAL_SINGLET_EXPECTED <- 82893L

total_input <- sum(tab$input_cells)
total_singlets <- sum(tab$singlets)
total_doublets <- sum(tab$doublets)

cat("\nTOTALS\n")
cat("Input:   ", total_input, "\n")
cat("Singlets:", total_singlets, "\n")
cat("Doublets:", total_doublets, "\n")
cat(
  "Doublet %:",
  sprintf("%.2f", 100*total_doublets/total_input),
  "\n"
)

cat("\nSINGLETS BY DISEASE\n")
print(
  aggregate(
    singlets ~ disease,
    data=tab,
    sum
  ),
  row.names=FALSE
)

# =============================================================================
# HARD REPRODUCIBILITY CHECKS
# =============================================================================

stopifnot(
  total_input == TOTAL_INPUT_EXPECTED
)

if (
  any(is.na(tab$count_match)) ||
  !all(tab$count_match)
) {
  stop(
    "V2 scDblFinder counts do not reproduce historical A02 exactly."
  )
}

if (total_singlets != TOTAL_SINGLET_EXPECTED) {
  stop(
    "Expected 82,893 singlets; obtained ",
    total_singlets
  )
}

# If historical cell identities are available, require exact identity.
if (
  identity_check_performed &&
  !identity_exact
) {
  stop(
    paste0(
      "Singlet counts may match, but cell identities differ from historical A02. ",
      "Inspect 04_singlet_identity_comparison.csv before continuing."
    )
  )
}

# =============================================================================
# SAVE ONLY AFTER VALIDATION
# =============================================================================

saveRDS(
  clean,
  OUT_RDS,
  compress=FALSE
)

# =============================================================================
# AUDIT
# =============================================================================

audit <- c(
  "PIPELINE_V2 — 04 doublet audit",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0(names(pkg_versions), ": ", pkg_versions),
  "",
  paste0("Seed: ", SEED),
  "Method: scDblFinder per biological sample",
  "Input assay: RNA counts",
  "samples argument: sample",
  "",
  paste0("Input cells: ", total_input),
  paste0("Singlets: ", total_singlets),
  paste0("Doublets: ", total_doublets),
  paste0(
    "Overall doublet percent: ",
    sprintf("%.4f",100*total_doublets/total_input)
  ),
  "",
  paste0(
    "Historical count match all samples: ",
    all(tab$count_match)
  ),
  paste0(
    "Historical identity check performed: ",
    identity_check_performed
  ),
  paste0(
    "Historical cell identities exact match: ",
    identity_exact
  ),
  "",
  paste0("Checkpoint: ", OUT_RDS),
  paste0("Summary: ", OUT_TABLE)
)

writeLines(
  audit,
  OUT_AUDIT
)

cat("\n============================================================\n")
cat("V2 DOUBLETS COMPLETE — EXACT HISTORICAL REPRODUCTION\n")
cat("============================================================\n")
cat(
  format(total_input,big.mark=","),
  " -> ",
  format(total_singlets,big.mark=","),
  " singlets\n",
  sep=""
)

if (identity_check_performed)
  cat("Historical singlet cell identities: EXACT MATCH\n")

cat("\nCheckpoint:\n  ",OUT_RDS,"\n")
cat("Summary:\n  ",OUT_TABLE,"\n")
cat("Audit:\n  ",OUT_AUDIT,"\n")
