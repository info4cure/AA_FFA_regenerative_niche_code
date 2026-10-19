suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H4A RETURN FIBROBLAST C9 TO EPITHELIAL\n")
cat("============================================================\n\n")

FIB_IN <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09F_refined.rds"
EPI_FROZEN <- "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds"

OUTDIR <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09H4A_C9_return_epithelial"
)
TAB <- file.path(OUTDIR, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create("PIPELINE_V2/objects/annotated",
           recursive=TRUE, showWarnings=FALSE)

FIB_OUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4A_C9_removed.rds"
)

EPI_OUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "epithelial_v2_09H4A_C9_returned_unadjudicated.rds"
)

cat("STEP 1/8 — loading immutable source objects\n")

fib <- readRDS(FIB_IN)
epi <- readRDS(EPI_FROZEN)

cat("Fibroblast source cells:", ncol(fib), "\n")
cat("Frozen epithelial cells:", ncol(epi), "\n")

stopifnot(ncol(fib) == 4714)

RES03 <- "fibroblast_refined_snn_res.0.3"

if (!RES03 %in% colnames(fib@meta.data)) {
  stop("Missing fibroblast r0.3 metadata.")
}

cat("\nSTEP 2/8 — identifying r0.3 cluster 9\n")

c9 <- colnames(fib)[
  as.character(fib@meta.data[[RES03]]) == "9"
]

cat("C9 cells:", length(c9), "\n")

stopifnot(length(c9) == 29)
stopifnot(length(unique(c9)) == 29)

meta_c9 <- fib@meta.data[c9, , drop=FALSE] %>%
  rownames_to_column("cell")

write_csv(
  meta_c9,
  file.path(TAB, "09H4A_C9_source_metadata.csv")
)

patient_summary <- meta_c9 %>%
  count(patient_id, disease, batch, name="n_cells") %>%
  arrange(disease, patient_id)

write_csv(
  patient_summary,
  file.path(TAB, "09H4A_C9_patient_distribution.csv")
)

print(patient_summary)

cat("\nSTEP 3/8 — formal return ledger\n")

ledger <- meta_c9 %>%
  transmute(
    cell = cell,
    patient_id = patient_id,
    disease = disease,
    batch = batch,
    source_L1 = "Fibroblast",
    source_object = FIB_IN,
    source_stage = "09F_refined",
    source_resolution = "0.3",
    source_cluster = "9",
    adjudication_stage = "09H4A",
    final_disposition = "RETURN_EPITHELIAL",
    return_reason = paste0(
      "Recurrent multi-patient cluster lacking positive ",
      "FibConsensus-family support and showing coherent ",
      "epithelial marker programme"
    ),
    epithelial_L2_status = "UNADJUDICATED_AFTER_RETURN"
  )

write_csv(
  ledger,
  file.path(TAB, "09H4A_C9_return_ledger.csv")
)

cat("\nSTEP 4/8 — checking relationship to frozen epithelial universe\n")

already_epi <- intersect(c9, colnames(epi))

cat("C9 cells already present in frozen epithelial object:",
    length(already_epi), "\n")

if (length(already_epi) > 0) {
  stop(
    "C9 cells unexpectedly already exist in frozen epithelial object."
  )
}

cat("\nSTEP 5/8 — constructing refined fibroblast universe\n")

fib_new <- subset(
  fib,
  cells=setdiff(colnames(fib), c9)
)

cat("Fibroblast cells after return:", ncol(fib_new), "\n")

stopifnot(ncol(fib_new) == 4685)

fib_new$fibroblast_09H4A_status <- "RETAIN_FIBROBLAST_UNIVERSE"
fib_new$fibroblast_analysis_stage <- "09H4A_C9_removed"

# Do NOT treat inherited geometry as final.
# It remains present only as provenance until 09H4B rebuild.
fib_new@misc[["09H4A"]] <- list(
  source_object = FIB_IN,
  source_cells = 4714,
  returned_epithelial_cells = 29,
  retained_fibroblast_cells = 4685,
  returned_source_cluster = "r0.3_C9",
  geometry_status = "STALE_REBUILD_REQUIRED_09H4B"
)

saveRDS(fib_new, FIB_OUT)

cat("\nSTEP 6/8 — constructing epithelial return object\n")

# Extract the 29 cells directly from the fibroblast source.
# This preserves their corrected expression matrix and provenance.
epi_return <- subset(fib, cells=c9)

epi_return$L1_previous <- "Fibroblast"
epi_return$L1 <- "Epithelial"
epi_return$L1_status <- "RETURNED_AFTER_FIBROBLAST_REVIEW"
epi_return$L1_flag <- "09H4A_RETURN_EPITHELIAL"
epi_return$L1_frozen <- FALSE

epi_return$epithelial_return_stage <- "09H4A"
epi_return$epithelial_return_source <- "fibroblast_r0.3_C9"
epi_return$epithelial_L2_status <- "UNADJUDICATED_AFTER_RETURN"

# Remove fibroblast geometry/graphs from the returned-cell object.
epi_return@reductions <- list()
epi_return@graphs <- list()
epi_return@neighbors <- list()

saveRDS(epi_return, EPI_OUT)

cat("Returned epithelial object cells:", ncol(epi_return), "\n")

stopifnot(ncol(epi_return) == 29)

cat("\nSTEP 7/8 — quantitative audit\n")

summary_tbl <- tibble(
  metric = c(
    "source_fibroblast_cells",
    "returned_epithelial_cells",
    "refined_fibroblast_cells",
    "frozen_epithelial_source_cells",
    "patients_in_return",
    "diseases_in_return"
  ),
  value = c(
    4714,
    length(c9),
    ncol(fib_new),
    ncol(epi),
    length(unique(meta_c9$patient_id)),
    length(unique(meta_c9$disease))
  )
)

write_csv(
  summary_tbl,
  file.path(TAB, "09H4A_summary.csv")
)

print(summary_tbl)

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "source_fibroblast_4714",
    "C9_exactly_29_cells",
    "C9_ten_patients",
    "C9_three_diseases",
    "C9_not_already_in_frozen_epithelial",
    "refined_fibroblast_4685",
    "returned_object_29",
    "all_return_cells_preserved",
    "frozen_epithelial_not_modified",
    "source_fibroblast_not_overwritten",
    "returned_L2_unadjudicated",
    "fibroblast_geometry_marked_stale"
  ),
  pass = c(
    ncol(fib) == 4714,
    length(c9) == 29,
    length(unique(meta_c9$patient_id)) == 10,
    length(unique(meta_c9$disease)) == 3,
    length(already_epi) == 0,
    ncol(fib_new) == 4685,
    ncol(epi_return) == 29,
    setequal(colnames(epi_return), c9),
    file.exists(EPI_FROZEN),
    file.exists(FIB_IN),
    all(epi_return$epithelial_L2_status ==
        "UNADJUDICATED_AFTER_RETURN"),
    identical(
      fib_new@misc[["09H4A"]]$geometry_status,
      "STALE_REBUILD_REQUIRED_09H4B"
    )
  )
)

print(checks)

write_csv(
  checks,
  file.path(TAB, "09H4A_formal_checks.csv")
)

if (!all(checks$pass)) {
  stop("09H4A formal checks failed.")
}

audit <- c(
  "PIPELINE_V2 — 09H4A C9 RETURN TO EPITHELIAL",
  "============================================",
  "",
  "Source fibroblast universe: 4,714 cells",
  "Returned cluster: refined r0.3 cluster 9",
  "Returned cells: 29",
  paste0(
    "Patients represented: ",
    length(unique(meta_c9$patient_id))
  ),
  paste0(
    "Diseases represented: ",
    paste(sort(unique(meta_c9$disease)), collapse=", ")
  ),
  "Refined fibroblast universe after return: 4,685 cells",
  "",
  paste0("Fibroblast derived object: ", FIB_OUT),
  paste0("Returned epithelial object: ", EPI_OUT),
  "",
  "The frozen epithelial L2 object was NOT modified.",
  "Returned epithelial cells remain L2-unadjudicated.",
  "The source 09F fibroblast object was NOT overwritten.",
  "Inherited fibroblast geometry is explicitly marked stale.",
  "A complete fibroblast rebuild is required in 09H4B.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(TAB, "09H4A_audit.txt")
)

cat("\n============================================================\n")
cat("09H4A COMPLETE\n")
cat("============================================================\n")
cat("Source fibroblast cells: 4714\n")
cat("Returned to epithelial:", length(c9), "\n")
cat("Refined fibroblast cells:", ncol(fib_new), "\n")
cat("Patients represented in return:",
    length(unique(meta_c9$patient_id)), "\n")
cat("Frozen epithelial object NOT modified.\n")
cat("Returned cells remain epithelial-L2 unadjudicated.\n")
cat("Fibroblast geometry must be rebuilt in 09H4B.\n")
cat("PASS — C9 formally returned to epithelial.\n")
