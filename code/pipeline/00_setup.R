source("PIPELINE_V2/config/paths.R")

dirs <- c(
  V2_CHECKPOINTS,
  V2_FROZEN,
  V2_TABLES,
  V2_QC,
  V2_FIGURES,
  V2_LOGS
)

invisible(
  lapply(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

manifest <- read.csv(
  "PIPELINE_V2/config/sample_manifest.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(manifest) == 16,
  length(unique(manifest$patient_id)) == 16,
  sum(manifest$disease == "HC") == 5,
  sum(manifest$disease == "AA") == 6,
  sum(manifest$disease == "FFA") == 5
)

ann <- read.csv(
  "PIPELINE_V2/config/annotations_L1.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(ann) == 27,
  length(unique(ann$global_cluster)) == 27
)

contrasts <- read.csv(
  "PIPELINE_V2/config/contrasts.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(contrasts) == 3
)

cat("\n========================================\n")
cat("PIPELINE_V2 SETUP OK\n")
cat("========================================\n")
cat("Individuals: 16\n")
cat("HC:  5\n")
cat("AA:  6\n")
cat("FFA: 5\n")
cat("Global clusters:", nrow(ann), "\n")
cat(
  "Contrasts:",
  paste(contrasts$contrast, collapse = ", "),
  "\n"
)
