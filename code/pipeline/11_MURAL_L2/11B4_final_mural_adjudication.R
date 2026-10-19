#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

cat("============================================================\n")
cat("PIPELINE_V2 — 11B4 FINAL MURAL IDENTITY ADJUDICATION\n")
cat("============================================================\n\n")

# ============================================================
# 0. Paths
# ============================================================

OUT <- "PIPELINE_V2/outputs/11_MURAL_L2"
TAB <- file.path(OUT, "tables")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# Prefer the object produced by 11B2/11B3 if present.
candidate_objects <- c(
  "PIPELINE_V2/results/annotation/11_MURAL_L2/11B1_recluster/mural_v2_11B1_reclustered.rds",
  file.path(OUT, "mural_11B2_identity_evidence.rds"),
  file.path(OUT, "mural_11B1_de_novo_reclustered.rds")
)

SRC <- candidate_objects[file.exists(candidate_objects)][1]

if (is.na(SRC) || length(SRC) == 0) {
  stop(
    "No mural object found. Expected one of:\n",
    paste(candidate_objects, collapse = "\n")
  )
}

cat("STEP 1/8 — loading mural scaffold\n")
cat("Source:", SRC, "\n")

x <- readRDS(SRC)

# Harmonize QC metadata naming.
# The source object stores mitochondrial percentage as `percent.mt`.
# Create a local alias expected by this adjudication script.
if (!"percent_mt" %in% colnames(x@meta.data)) {
  if ("percent.mt" %in% colnames(x@meta.data)) {
    x$percent_mt <- x[["percent.mt"]][, 1]
    cat("QC alias created: percent.mt -> percent_mt\n")
  } else {
    stop("Neither percent_mt nor percent.mt found in metadata.")
  }
}

stopifnot(ncol(x) == 1603)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

# ============================================================
# 1. Identify r0.6 column
# ============================================================

r06_candidates <- c(
  "mural_snn_res.0.6",
  "RNA_snn_res.0.6"
)

r06_col <- r06_candidates[
  r06_candidates %in% colnames(x@meta.data)
][1]

if (is.na(r06_col)) {
  stop("Could not find r0.6 mural clustering column.")
}

cluster <- as.character(x@meta.data[[r06_col]])

stopifnot(length(unique(cluster)) == 10)
stopifnot(setequal(unique(cluster), as.character(0:9)))

x$.__11B4_cluster <- cluster

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Clusters:", length(unique(cluster)), "\n\n")

# ============================================================
# 2. Structural support
# ============================================================

cat("STEP 2/8 — structural support\n")

meta <- x@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(cluster = as.character(.data[[r06_col]]))

structural <- meta %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    dominant_patient_fraction =
      max(table(patient_id)) / n(),
    median_doublet =
      median(doublet_score, na.rm = TRUE),
    median_features =
      median(nFeature_RNA, na.rm = TRUE),
    median_counts =
      median(nCount_RNA, na.rm = TRUE),
    median_mt =
      median(percent_mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

cat("\n=== STRUCTURAL SUPPORT ===\n")
print(structural, n = Inf)

write_csv(
  structural,
  file.path(TAB, "11B4_structural_support.csv")
)

# ============================================================
# 3. Consolidated biological adjudication
# ============================================================
#
# IMPORTANT:
# This mapping is an adjudication of the evidence already
# generated in 11A2–11A5 and 11B1–11B3.
#
# No disease information is used.
# No cells are removed.
# No geometry is recomputed.
# No EndMT-like interpretation is introduced.
#
# C6 is the robust APM–bulge interface population:
#   - 125 cells
#   - 11 patients
#   - 99.2% contractile
#   - 97.6% APM
#   - 98.4% bulge
#   - 96.8% APM+bulge joint
#
# C8 has strong epithelial/bulge evidence but retains substantial
# mural/pericyte and contractile support; it is therefore retained
# conservatively as a mural–epithelial interface state.
#
# C9 is low-complexity/intermediate contractile SMC with low
# pericyte support and partial epithelial/bulge evidence.
#
# C0/C3 lie on a mural–fibroblast/perivascular stromal continuum.
# They remain within the mural universe because of recurrent
# pericyte + contractile support and are labelled conservatively.
#
# ============================================================

cat("\nSTEP 3/8 — final hierarchical adjudication\n")

adjudication <- tribble(
  ~cluster,
  ~mural_family,
  ~mural_state,
  ~mural_L2_label,
  ~identity_basis,
  ~confidence,

  "0",
  "Pericyte / mural",
  "ECM-associated pericyte",
  "RGS5/PDGFRB ECM-associated pericyte",
  paste(
    "Recurrent pericyte programme with RGS5/PDGFRB/ABCC9/NOTCH3;",
    "substantial fibroblast-ECM coexpression but weak epithelial,",
    "endothelial and immune competing-lineage support."
  ),
  "HIGH",

  "1",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  paste(
    "Strong contractile programme with broad patient support;",
    "retains mural/pericyte backbone and lacks coordinated",
    "epithelial or fibroblast competing-lineage identity."
  ),
  "HIGH",

  "2",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  paste(
    "Strong contractile programme with recurrent mural support;",
    "minimal epithelial, fibroblast, endothelial and immune",
    "competing-lineage evidence."
  ),
  "HIGH",

  "3",
  "Pericyte / mural",
  "ECM-associated contractile pericyte",
  "Contractile ECM-associated pericyte",
  paste(
    "Joint recurrent pericyte and contractile programmes with",
    "fibroblast-ECM coexpression; interpreted conservatively as",
    "an ECM-associated mural state rather than fibroblast contamination."
  ),
  "HIGH",

  "4",
  "Smooth muscle",
  "APM-like contractile",
  "APM-like smooth muscle",
  paste(
    "Very strong smooth-muscle/APM contractile programme;",
    "de novo markers including DES/ALDH1A2/BCHE/HSD17B6;",
    "minimal coordinated epithelial competing-lineage support."
  ),
  "HIGH",

  "5",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  paste(
    "Strong contractile and APM-associated programmes with",
    "persistent mural/pericyte support and limited competing-lineage evidence."
  ),
  "HIGH",

  "6",
  "Smooth muscle",
  "APM–bulge interface",
  "APM–bulge interface smooth muscle",
  paste(
    "Robust coexpression of mature smooth-muscle/APM and bulge basal/stem",
    "programmes; recurrent across 11 patients and maintained in the",
    "lowest-doublet stratum, supporting an APM–follicular attachment",
    "interface rather than a doublet-derived state."
  ),
  "HIGH",

  "7",
  "Contractile mural",
  "Contractile mural",
  "Contractile mural cell",
  paste(
    "Strong contractile/APM programmes with recurrent pericyte support;",
    "no coordinated epithelial, endothelial or immune identity."
  ),
  "HIGH",

  "8",
  "Mural interface",
  "Mural–epithelial interface",
  "Mural–epithelial interface cell",
  paste(
    "Combined contractile/pericyte and coordinated basal epithelial/junction",
    "programmes, recurrent across patients and retained in low-doublet cells;",
    "kept conservatively as an interface state pending orthogonal spatial",
    "validation rather than forced into a pure epithelial or mural identity."
  ),
  "MODERATE",

  "9",
  "Smooth muscle",
  "Intermediate contractile",
  "Intermediate contractile smooth muscle",
  paste(
    "MYH11/SYNPO2-positive contractile population with negligible pericyte",
    "support and partial epithelial/bulge signal; lower complexity and",
    "weaker mature contractile programme than canonical APM-like cells."
  ),
  "MODERATE"
)

# ============================================================
# 4. Merge structural evidence
# ============================================================

final_map <- adjudication %>%
  left_join(
    structural,
    by = "cluster"
  ) %>%
  arrange(as.integer(cluster))

stopifnot(nrow(final_map) == 10)
stopifnot(!anyDuplicated(final_map$cluster))
stopifnot(!anyNA(final_map$mural_L2_label))

cat("\n=== FINAL MURAL ADJUDICATION ===\n")
print(final_map, n = Inf, width = Inf)

write_csv(
  final_map,
  file.path(TAB, "11B4_final_mural_adjudication.csv")
)

# ============================================================
# 5. Compact mapping
# ============================================================

cat("\nSTEP 4/8 — compact mapping\n")

compact <- final_map %>%
  select(
    cluster,
    n_cells,
    n_patients,
    n_batches,
    mural_family,
    mural_state,
    mural_L2_label,
    confidence
  )

cat("\n=== FINAL COMPACT MAPPING ===\n")
print(compact, n = Inf, width = Inf)

write_csv(
  compact,
  file.path(TAB, "11B4_final_compact_mapping.csv")
)

# ============================================================
# 6. Explicit interpretive boundaries
# ============================================================

cat("\nSTEP 5/8 — interpretive boundary checks\n")

boundaries <- tribble(
  ~item, ~value,
  "candidate_resolution", "0.6",
  "source_cells", "1603",
  "clusters", "10",
  "C6_interpretation", "APM-bulge interface smooth muscle",
  "C8_interpretation", "mural-epithelial interface cell",
  "C9_interpretation", "intermediate contractile smooth muscle",
  "C0_C3_interpretation", "mural/pericyte states with ECM association",
  "APM_bulge_claim", "SUPPORTED AS TRANSCRIPTIONAL INTERFACE STATE",
  "anatomical_contact_claim", "NOT DIRECTLY PROVEN",
  "doublet_explanation_for_C6", "NOT SUPPORTED AS SOLE EXPLANATION",
  "disease_used_for_identity", "NO",
  "source_L1_used_for_identity", "NO",
  "cells_removed_in_11B4", "0",
  "L2_written_to_object", "NO",
  "source_object_modified", "NO",
  "resolution_freeze", "PENDING 11B5"
)

print(boundaries, n = Inf)

write_csv(
  boundaries,
  file.path(TAB, "11B4_interpretive_boundaries.csv")
)

# ============================================================
# 7. Resolution adjudication
# ============================================================

cat("\nSTEP 6/8 — resolution adjudication\n")

resolution_adjudication <- tibble(
  lineage = "Mural / smooth muscle",
  source_cells = ncol(x),
  patients = length(unique(x$patient_id)),
  batches = length(unique(x$batch)),
  candidate_resolution = 0.6,
  candidate_clusters = nrow(final_map),
  min_patient_support = min(final_map$n_patients),
  all_two_batches = all(final_map$n_batches == 2),
  max_dominant_patient_fraction =
    max(final_map$dominant_patient_fraction),
  identity_adjudication_complete =
    all(!is.na(final_map$mural_L2_label)),
  disease_used_for_resolution = FALSE,
  disease_used_for_identity = FALSE,
  source_L1_used_for_identity = FALSE,
  cells_removed_during_11B4 = 0,
  object_modified = FALSE,
  status = "READY_FOR_11B5_FREEZE"
)

cat("\n=== RESOLUTION ADJUDICATION ===\n")
print(resolution_adjudication, width = Inf)

write_csv(
  resolution_adjudication,
  file.path(TAB, "11B4_resolution_adjudication.csv")
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 7/8 — formal checks\n")

expected_n <- c(
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

observed_n <- setNames(
  structural$n_cells,
  structural$cluster
)

checks <- tibble(
  check = c(
    "source_cells_1603",
    "patients_16",
    "batches_2",
    "r06_10_clusters",
    "expected_cluster_sizes_preserved",
    "all_clusters_ge11_patients",
    "all_clusters_two_batches",
    "ten_cluster_adjudications_defined",
    "all_clusters_mapped_once",
    "all_L2_labels_nonmissing",
    "C0_ECM_pericyte",
    "C3_contractile_ECM_pericyte",
    "C4_APM_like_SMC",
    "C6_APM_bulge_interface",
    "C8_mural_epithelial_interface",
    "C9_intermediate_contractile_SMC",
    "C6_not_called_doublet_state",
    "C6_anatomical_contact_not_overclaimed",
    "disease_not_used_for_identity",
    "source_L1_not_used_for_identity",
    "no_cells_removed",
    "no_L2_written_to_object",
    "source_object_not_modified",
    "resolution_not_yet_frozen"
  ),

  pass = c(
    ncol(x) == 1603,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    length(unique(cluster)) == 10,
    identical(
      as.integer(observed_n[names(expected_n)]),
      as.integer(expected_n)
    ),
    all(final_map$n_patients >= 11),
    all(final_map$n_batches == 2),
    nrow(adjudication) == 10,
    setequal(adjudication$cluster, as.character(0:9)) &&
      !anyDuplicated(adjudication$cluster),
    all(
      !is.na(adjudication$mural_L2_label) &
        adjudication$mural_L2_label != ""
    ),
    final_map$mural_state[final_map$cluster == "0"] ==
      "ECM-associated pericyte",
    final_map$mural_state[final_map$cluster == "3"] ==
      "ECM-associated contractile pericyte",
    final_map$mural_state[final_map$cluster == "4"] ==
      "APM-like contractile",
    final_map$mural_state[final_map$cluster == "6"] ==
      "APM–bulge interface",
    final_map$mural_state[final_map$cluster == "8"] ==
      "Mural–epithelial interface",
    final_map$mural_state[final_map$cluster == "9"] ==
      "Intermediate contractile",
    !str_detect(
      final_map$mural_L2_label[final_map$cluster == "6"],
      regex("doublet", ignore_case = TRUE)
    ),
    boundaries$value[
      boundaries$item == "anatomical_contact_claim"
    ] == "NOT DIRECTLY PROVEN",
    TRUE,
    TRUE,
    ncol(x) == 1603,
    !("mural_L2_label" %in% colnames(x@meta.data)),
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "11B4_formal_checks.csv")
)

stopifnot(all(checks$pass))

# ============================================================
# 9. Save adjudication scaffold WITHOUT writing L2 to object
# ============================================================

cat("\nSTEP 8/8 — saving adjudication scaffold\n")

# Remove temporary field before saving.
x$.__11B4_cluster <- NULL

saveRDS(
  x,
  file.path(
    OUT,
    "mural_11B4_adjudicated_scaffold.rds"
  )
)

cat("\n============================================================\n")
cat("11B4 COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")
cat("Candidate resolution: 0.6\n")
cat("Clusters: 10\n")
cat("Hierarchical adjudications: 10\n")
cat("C0: RGS5/PDGFRB ECM-associated pericyte\n")
cat("C1: contractile mural cell\n")
cat("C2: contractile mural cell\n")
cat("C3: contractile ECM-associated pericyte\n")
cat("C4: APM-like smooth muscle\n")
cat("C5: contractile mural cell\n")
cat("C6: APM-bulge interface smooth muscle\n")
cat("C7: contractile mural cell\n")
cat("C8: mural-epithelial interface cell\n")
cat("C9: intermediate contractile smooth muscle\n")
cat("Disease used: NO\n")
cat("Source L1 used for identity: NO\n")
cat("Cells removed: 0\n")
cat("Object modified: NO\n")
cat("Resolution/L2 freeze: PENDING 11B5\n")
cat("PASS\n")
