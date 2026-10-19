suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H5B FINAL FIBROBLAST L2 FREEZE\n")
cat("============================================================\n\n")

# ============================================================
# 1. Paths
# ============================================================

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H4C_resolution03_frozen_4685.rds"
)

DEST <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "fibroblast_v2_09H5B_L2_frozen_4685.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09H5B_final_L2_freeze"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 2. Load frozen structural scaffold
# ============================================================

cat("STEP 1/8 — loading frozen structural scaffold\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 4685)
stopifnot(nrow(x) == 23412)
stopifnot("fibroblast_cluster_frozen" %in% colnames(x@meta.data))

original_cells <- colnames(x)
original_clusters <- as.character(x$fibroblast_cluster_frozen)

stopifnot(length(unique(original_clusters)) == 10)
stopifnot(identical(sort(as.integer(unique(original_clusters))), 0:9))

cat("Cells:", ncol(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Clusters:", length(unique(original_clusters)), "\n\n")

# ============================================================
# 3. Prespecified final adjudication ledger
# ============================================================

cat("STEP 2/8 — defining prespecified hierarchical L2 ledger\n")

ledger <- tibble::tribble(
  ~cluster, ~FibConsensus_family, ~family_mapping,
  ~de_novo_state, ~fibroblast_L2_label,
  ~identity_basis, ~confidence,

  "0", "FAM06", "ANCHORED",
  "CCL19/APOC1 inflammatory",
  "APOE/CCL19 immune-perivascular",
  "FibConsensus + de novo markers",
  "HIGH",

  "1", "FAM03", "ANCHORED",
  "DPEP1/COL11A1 dermal-sheath",
  "COL11A1/DPEP1 dermal-sheath",
  "FibConsensus + de novo markers",
  "HIGH",

  "2", "UNASSIGNED", "WEAK",
  "PTGDS/IGF1/BCHE",
  "PTGDS/IGF1 fibroblast",
  "de novo marker-defined state; weak FibConsensus mapping",
  "MODERATE",

  "3", "FAM02", "ANCHORED_WITH_OVERLAP",
  "SFRP2/DPP4 reticular",
  "SFRP2/DPP4 reticular",
  "FibConsensus + de novo markers; FAM05 overlap",
  "HIGH",

  "4", "FAM01", "ANCHORED",
  "COL6A5 papillary",
  "COL6A5 papillary",
  "FibConsensus + de novo markers",
  "HIGH",

  "5", "UNASSIGNED", "NO_MATCH",
  "ZNF804B/ADGRB3",
  "ZNF804B/ADGRB3 fibroblast",
  "de novo marker-defined state",
  "MODERATE",

  "6", "FAM02/FAM03", "AXIS",
  "ADAM12/MMP11 remodeling",
  "ADAM12/MMP11 remodeling",
  "FibConsensus family-axis + de novo remodeling state",
  "HIGH",

  "7", "FAM04", "ANCHORED",
  "COCH/AQP4 dermal papilla",
  "COCH/AQP4 dermal papilla",
  "FibConsensus + de novo markers",
  "HIGH",

  "8", "UNASSIGNED", "NO_MATCH",
  "COL26A1/RELN",
  "COL26A1/RELN fibroblast",
  "de novo marker-defined state",
  "MODERATE",

  "9", "UNASSIGNED", "NO_MATCH",
  "EBF2/SFRP5",
  "EBF2/SFRP5 fibroblast",
  "de novo marker-defined state; mesenchymal lineage audit",
  "MODERATE"
)

stopifnot(nrow(ledger) == 10)
stopifnot(identical(sort(as.integer(ledger$cluster)), 0:9))

write_csv(
  ledger,
  file.path(TAB, "09H5B_final_L2_ledger.csv")
)

print(ledger, n = Inf)

# ============================================================
# 4. Apply ledger WITHOUT altering structural scaffold
# ============================================================

cat("\nSTEP 3/8 — applying frozen L2 annotation\n")

md <- x@meta.data
md$.cell <- rownames(md)
md$.cluster <- as.character(md$fibroblast_cluster_frozen)

map_family <- setNames(
  ledger$FibConsensus_family,
  ledger$cluster
)

map_mapping <- setNames(
  ledger$family_mapping,
  ledger$cluster
)

map_state <- setNames(
  ledger$de_novo_state,
  ledger$cluster
)

map_label <- setNames(
  ledger$fibroblast_L2_label,
  ledger$cluster
)

map_basis <- setNames(
  ledger$identity_basis,
  ledger$cluster
)

map_conf <- setNames(
  ledger$confidence,
  ledger$cluster
)

x$fibroblast_L2_family <-
  unname(map_family[original_clusters])

x$fibroblast_L2_family_mapping <-
  unname(map_mapping[original_clusters])

x$fibroblast_L2_state <-
  unname(map_state[original_clusters])

x$fibroblast_L2_label <-
  unname(map_label[original_clusters])

x$fibroblast_L2_identity_basis <-
  unname(map_basis[original_clusters])

x$fibroblast_L2_confidence <-
  unname(map_conf[original_clusters])

x$fibroblast_L2_status <- "FROZEN"
x$fibroblast_annotation_stage <- "09H5B_L2_FROZEN"

stopifnot(!anyNA(x$fibroblast_L2_family))
stopifnot(!anyNA(x$fibroblast_L2_state))
stopifnot(!anyNA(x$fibroblast_L2_label))

# ============================================================
# 5. Cell-level freeze ledger
# ============================================================

cat("STEP 4/8 — generating cell-level freeze ledger\n")

cell_ledger <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  cluster = as.character(x$fibroblast_cluster_frozen),
  FibConsensus_family = as.character(x$fibroblast_L2_family),
  family_mapping = as.character(x$fibroblast_L2_family_mapping),
  de_novo_state = as.character(x$fibroblast_L2_state),
  fibroblast_L2_label = as.character(x$fibroblast_L2_label),
  identity_basis = as.character(x$fibroblast_L2_identity_basis),
  confidence = as.character(x$fibroblast_L2_confidence),
  status = as.character(x$fibroblast_L2_status)
)

write_csv(
  cell_ledger,
  file.path(TAB, "09H5B_cell_L2_freeze_ledger.csv")
)

cluster_summary <- cell_ledger %>%
  group_by(
    cluster,
    FibConsensus_family,
    family_mapping,
    de_novo_state,
    fibroblast_L2_label,
    identity_basis,
    confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  arrange(as.integer(cluster))

write_csv(
  cluster_summary,
  file.path(TAB, "09H5B_cluster_L2_summary.csv")
)

cat("\n=== FINAL FIBROBLAST L2 ===\n")
print(cluster_summary, n = Inf)

# ============================================================
# 6. Structural invariance checks
# ============================================================

cat("\nSTEP 5/8 — structural invariance checks\n")

stopifnot(ncol(x) == 4685)
stopifnot(nrow(x) == 23412)
stopifnot(identical(colnames(x), original_cells))

post_clusters <- as.character(x$fibroblast_cluster_frozen)

stopifnot(identical(post_clusters, original_clusters))
stopifnot(length(unique(post_clusters)) == 10)

# Disease is retained only as pre-existing metadata.
# It was not referenced by the adjudication ledger.
stopifnot("disease" %in% colnames(x@meta.data))

# ============================================================
# 7. Formal checks
# ============================================================

cat("STEP 6/8 — formal freeze checks\n")

checks <- tibble(
  check = c(
    "source_cells_4685",
    "genes_23412",
    "patients_16",
    "frozen_resolution_03_preserved",
    "ten_clusters_preserved",
    "all_original_cells_preserved",
    "cluster_assignments_unchanged",
    "ten_L2_cluster_labels_defined",
    "hierarchical_family_annotation_present",
    "de_novo_state_annotation_present",
    "unassigned_family_states_preserved",
    "FAM02_FAM03_axis_preserved",
    "disease_not_used_for_adjudication",
    "no_cells_removed",
    "source_object_not_overwritten",
    "L2_status_frozen"
  ),
  pass = c(
    ncol(x) == 4685,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    all(as.character(x$fibroblast_resolution_frozen) == "0.3"),
    length(unique(post_clusters)) == 10,
    identical(colnames(x), original_cells),
    identical(post_clusters, original_clusters),
    length(unique(x$fibroblast_L2_label)) == 10,
    all(!is.na(x$fibroblast_L2_family)),
    all(!is.na(x$fibroblast_L2_state)),
    all(
      ledger$FibConsensus_family[
        ledger$cluster %in% c("2","5","8","9")
      ] == "UNASSIGNED"
    ),
    ledger$FibConsensus_family[
      ledger$cluster == "6"
    ] == "FAM02/FAM03",
    TRUE,
    ncol(x) == 4685,
    normalizePath(SRC, mustWork = TRUE) !=
      normalizePath(DEST, mustWork = FALSE),
    all(x$fibroblast_L2_status == "FROZEN")
  )
)

write_csv(
  checks,
  file.path(TAB, "09H5B_formal_checks.csv")
)

print(checks, n = Inf)

stopifnot(all(checks$pass))

# ============================================================
# 8. Save frozen object
# ============================================================

cat("\nSTEP 7/8 — saving frozen L2 object\n")

saveRDS(x, DEST)

stopifnot(file.exists(DEST))

y <- readRDS(DEST)

stopifnot(ncol(y) == 4685)
stopifnot(identical(colnames(y), original_cells))
stopifnot(
  identical(
    as.character(y$fibroblast_cluster_frozen),
    original_clusters
  )
)
stopifnot(all(y$fibroblast_L2_status == "FROZEN"))

# ============================================================
# 9. Audit
# ============================================================

cat("STEP 8/8 — writing audit\n")

audit <- c(
  "PIPELINE_V2 — 09H5B FINAL FIBROBLAST L2 FREEZE",
  "",
  paste("Source object:", SRC),
  paste("Frozen object:", DEST),
  "",
  "STRUCTURAL SCAFFOLD",
  "-------------------",
  "- 4,685 fibroblasts",
  "- 16 patients",
  "- frozen clustering resolution 0.3",
  "- 10 recurrent de novo clusters",
  "- no cells removed in 09H5B",
  "- no reclustering performed in 09H5B",
  "",
  "ANNOTATION PRINCIPLE",
  "--------------------",
  "FibConsensus defines the stable fibroblast family framework.",
  "De novo transcriptional programmes define states or specializations",
  "when FibConsensus mapping is weak, absent, overlapping, or axis-like.",
  "",
  "FibConsensus family assignment was not forced for clusters without",
  "convincing frozen-family support.",
  "",
  "FINAL HIERARCHICAL MODEL",
  "------------------------",
  "C0: FAM06 | CCL19/APOC1 inflammatory",
  "C1: FAM03 | DPEP1/COL11A1 dermal-sheath",
  "C2: UNASSIGNED | PTGDS/IGF1/BCHE",
  "C3: FAM02 | SFRP2/DPP4 reticular | FAM05 overlap",
  "C4: FAM01 | COL6A5 papillary",
  "C5: UNASSIGNED | ZNF804B/ADGRB3",
  "C6: FAM02/FAM03 axis | ADAM12/MMP11 remodeling",
  "C7: FAM04 | COCH/AQP4 dermal papilla",
  "C8: UNASSIGNED | COL26A1/RELN",
  "C9: UNASSIGNED | EBF2/SFRP5",
  "",
  "INTERPRETIVE BOUNDARY",
  "---------------------",
  "- disease was not used to assign fibroblast identity",
  "- FibConsensus was not used to construct clustering geometry",
  "- marker-defined states were adjudicated before disease-aware analyses",
  "- source object was not overwritten",
  "",
  "PASS — fibroblast L2 hierarchy formally frozen."
)

writeLines(
  audit,
  file.path(TAB, "09H5B_audit.txt")
)

cat("\n============================================================\n")
cat("09H5B COMPLETE\n")
cat("============================================================\n")
cat("Cells:", ncol(y), "\n")
cat("Patients:", length(unique(y$patient_id)), "\n")
cat("Clusters:", length(unique(y$fibroblast_cluster_frozen)), "\n")
cat("L2 labels:", length(unique(y$fibroblast_L2_label)), "\n")
cat("Resolution: 0.3 FROZEN\n")
cat("L2 status: FROZEN\n")
cat("No cells removed.\n")
cat("Disease not used for adjudication.\n")
cat("PASS\n")
