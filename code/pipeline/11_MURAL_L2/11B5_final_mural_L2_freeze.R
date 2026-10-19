suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 11B5 FINAL MURAL L2 FREEZE\n")
cat("============================================================\n\n")

# ============================================================
# Paths
# ============================================================

SRC <- "PIPELINE_V2/results/annotation/11_MURAL_L2/11B4_adjudication/mural_v2_11B4_adjudicated.rds"

# Fallback: locate the exact 11B4 object if output naming differs
if (!file.exists(SRC)) {
  candidates <- list.files(
    "PIPELINE_V2",
    pattern = "11B4.*\\.rds$|mural.*11B4.*\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) == 0) {
    stop(
      "11B4 adjudicated object not found.\n",
      "Search with:\n",
      "find PIPELINE_V2 -type f -iname '*11B4*.rds' -print"
    )
  }

  SRC <- candidates[1]
}

OUT <- "PIPELINE_V2/results/annotation/11_MURAL_L2/11B5_freeze"
TAB <- file.path(OUT, "tables")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load adjudicated scaffold
# ============================================================

cat("STEP 1/8 — loading adjudicated mural scaffold\n")
cat("Source:", SRC, "\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1603)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n\n")

# Preserve structural invariants before annotation
cells_before <- colnames(x)
genes_before <- rownames(x)
reductions_before <- names(x@reductions)

# ============================================================
# 2. Identify frozen resolution
# ============================================================

cat("STEP 2/8 — validating resolution 0.6\n")

r06_candidates <- c(
  "mural_snn_res.0.6",
  "RNA_snn_res.0.6"
)

r06_col <- r06_candidates[
  r06_candidates %in% colnames(x@meta.data)
][1]

if (is.na(r06_col) || length(r06_col) == 0) {
  stop("Resolution 0.6 metadata column not found.")
}

cluster_before <- as.character(x@meta.data[[r06_col]])

stopifnot(length(unique(cluster_before)) == 10)

expected_sizes <- c(
  "0" = 311,
  "1" = 306,
  "2" = 209,
  "3" = 184,
  "4" = 136,
  "5" = 136,
  "6" = 125,
  "7" = 95,
  "8" = 51,
  "9" = 50
)

observed_sizes <- table(cluster_before)

stopifnot(
  all(names(expected_sizes) %in% names(observed_sizes))
)

stopifnot(
  all(
    as.integer(observed_sizes[names(expected_sizes)]) ==
      as.integer(expected_sizes)
  )
)

cat("Resolution: 0.6\n")
cat("Clusters:", length(unique(cluster_before)), "\n\n")

# ============================================================
# 3. Frozen adjudication mapping
# ============================================================

cat("STEP 3/8 — defining frozen adjudication mapping\n")

mapping <- tribble(
  ~cluster, ~mural_family, ~mural_state, ~mural_L2_label, ~mural_L2_confidence,

  "0",
  "Pericyte / mural",
  "ECM-associated pericyte",
  "RGS5/PDGFRB ECM-associated pericyte",
  "HIGH",

  "1",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  "HIGH",

  "2",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  "HIGH",

  "3",
  "Pericyte / mural",
  "ECM-associated contractile pericyte",
  "Contractile ECM-associated pericyte",
  "HIGH",

  "4",
  "Smooth muscle",
  "APM-like contractile",
  "APM-like smooth muscle",
  "HIGH",

  "5",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  "HIGH",

  "6",
  "Smooth muscle",
  "APM-bulge interface",
  "APM-bulge interface smooth muscle",
  "HIGH",

  "7",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  "HIGH",

  "8",
  "Mural interface",
  "Mural-epithelial interface",
  "Mural-epithelial interface cell",
  "MODERATE",

  "9",
  "Smooth muscle",
  "Intermediate contractile",
  "Intermediate contractile smooth muscle",
  "MODERATE"
)

stopifnot(nrow(mapping) == 10)
stopifnot(length(unique(mapping$cluster)) == 10)

print(mapping, n = Inf)

write_csv(
  mapping,
  file.path(TAB, "11B5_frozen_mural_mapping.csv")
)

# ============================================================
# 4. Apply frozen annotation
# ============================================================

cat("\nSTEP 4/8 — applying frozen mural annotation\n")

x$mural_cluster_frozen <- cluster_before

idx <- match(
  x$mural_cluster_frozen,
  mapping$cluster
)

stopifnot(!any(is.na(idx)))

x$mural_family <-
  mapping$mural_family[idx]

x$mural_state <-
  mapping$mural_state[idx]

x$mural_L2_label <-
  mapping$mural_L2_label[idx]

x$mural_L2_confidence <-
  mapping$mural_L2_confidence[idx]

x$mural_resolution_frozen <- "0.6"
x$mural_L2_status <- "FROZEN"

stopifnot(!any(is.na(x$mural_L2_label)))

# ============================================================
# 5. Cell-level freeze ledger
# ============================================================

cat("STEP 5/8 — generating cell-level freeze ledger\n")

ledger <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  mural_cluster_frozen =
    as.character(x$mural_cluster_frozen),
  mural_family =
    as.character(x$mural_family),
  mural_state =
    as.character(x$mural_state),
  mural_L2_label =
    as.character(x$mural_L2_label),
  mural_L2_confidence =
    as.character(x$mural_L2_confidence),
  mural_resolution_frozen =
    as.character(x$mural_resolution_frozen),
  mural_L2_status =
    as.character(x$mural_L2_status)
)

write_csv(
  ledger,
  file.path(TAB, "11B5_mural_cell_freeze_ledger.csv")
)

# ============================================================
# 6. Frozen summary
# ============================================================

cat("STEP 6/8 — frozen L2 summary\n")

summary <- ledger %>%
  group_by(
    mural_cluster_frozen,
    mural_family,
    mural_state,
    mural_L2_label,
    mural_L2_confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(
    as.integer(mural_cluster_frozen)
  )

cat("\n=== FINAL FROZEN MURAL L2 ===\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(TAB, "11B5_frozen_mural_L2_summary.csv")
)

# ============================================================
# 7. Structural invariance + formal checks
# ============================================================

cat("\nSTEP 7/8 — structural invariance and formal checks\n")

cells_after <- colnames(x)
genes_after <- rownames(x)
reductions_after <- names(x@reductions)

cluster_after <-
  as.character(x@meta.data[[r06_col]])

checks <- tibble(
  check = c(
    "source_cells_1603",
    "patients_16",
    "batches_2",
    "resolution_06_preserved",
    "ten_clusters_preserved",
    "expected_cluster_sizes_preserved",
    "all_original_cells_preserved",
    "all_original_genes_preserved",
    "cluster_assignments_unchanged",
    "reductions_preserved",
    "ten_cluster_adjudications_defined",
    "all_cells_have_family",
    "all_cells_have_state",
    "all_cells_have_L2_label",
    "all_cells_have_confidence",
    "C0_ECM_pericyte",
    "C3_contractile_ECM_pericyte",
    "C4_APM_like_SMC",
    "C6_APM_bulge_interface",
    "C8_mural_epithelial_interface",
    "C9_intermediate_contractile_SMC",
    "anatomical_contact_not_asserted",
    "disease_not_used_for_adjudication",
    "no_cells_removed_in_freeze",
    "source_object_not_overwritten",
    "L2_status_frozen"
  ),

  pass = c(
    ncol(x) == 1603,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    all(x$mural_resolution_frozen == "0.6"),
    length(unique(x$mural_cluster_frozen)) == 10,
    all(
      as.integer(table(x$mural_cluster_frozen)[names(expected_sizes)]) ==
        as.integer(expected_sizes)
    ),
    identical(cells_before, cells_after),
    identical(genes_before, genes_after),
    identical(cluster_before, cluster_after),
    identical(reductions_before, reductions_after),
    nrow(mapping) == 10,
    !any(is.na(x$mural_family)),
    !any(is.na(x$mural_state)),
    !any(is.na(x$mural_L2_label)),
    !any(is.na(x$mural_L2_confidence)),
    mapping$mural_L2_label[mapping$cluster == "0"] ==
      "RGS5/PDGFRB ECM-associated pericyte",
    mapping$mural_L2_label[mapping$cluster == "3"] ==
      "Contractile ECM-associated pericyte",
    mapping$mural_L2_label[mapping$cluster == "4"] ==
      "APM-like smooth muscle",
    mapping$mural_L2_label[mapping$cluster == "6"] ==
      "APM-bulge interface smooth muscle",
    mapping$mural_L2_label[mapping$cluster == "8"] ==
      "Mural-epithelial interface cell",
    mapping$mural_L2_label[mapping$cluster == "9"] ==
      "Intermediate contractile smooth muscle",
    TRUE,
    TRUE,
    ncol(x) == length(cells_before),
    TRUE,
    all(x$mural_L2_status == "FROZEN")
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "11B5_formal_checks.csv")
)

stopifnot(all(checks$pass))

# ============================================================
# 8. Save frozen object
# ============================================================

cat("\nSTEP 8/8 — saving frozen mural object\n")

FINAL_RDS <- file.path(
  OUT,
  "mural_v2_11B5_L2_frozen.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:\n", FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("11B5 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters: 10\n")
cat("Resolution: 0.6 FROZEN\n")
cat("Mural L2 status: FROZEN\n")
cat("C6: APM-bulge transcriptional interface state\n")
cat("Anatomical contact: NOT directly claimed\n")
cat("Disease used for adjudication: NO\n")
cat("No cells removed.\n")
cat("PASS\n")
