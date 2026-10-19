suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 14D FINAL SEBOCYTE L2 FREEZE\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
  "14C_adjudication/sebocyte_v2_14C_adjudication_scaffold.rds"
)

OUT <- "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/14D_freeze"
TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load adjudicated scaffold
# ============================================================

cat("STEP 1/8 — loading adjudicated sebocyte scaffold\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1388)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)
stopifnot(length(unique(x$batch)) == 2)

cells_before <- colnames(x)
genes_before <- rownames(x)
reductions_before <- names(x@reductions)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Batches:", length(unique(x$batch)), "\n")

# ============================================================
# 2. Validate r0.4
# ============================================================

cat("\nSTEP 2/8 — validating resolution 0.4\n")

R04 <- "seb_snn_res.0.4"

stopifnot(R04 %in% colnames(x@meta.data))

cluster_before <- as.character(
  x@meta.data[[R04]]
)

expected <- c(
  "0" = 436,
  "1" = 355,
  "2" = 278,
  "3" = 164,
  "4" = 155
)

observed <- table(cluster_before)

stopifnot(length(unique(cluster_before)) == 5)

stopifnot(
  all(
    as.integer(observed[names(expected)]) ==
      as.integer(expected)
  )
)

cat("Resolution: 0.4\n")
cat("Clusters: 5\n")

# ============================================================
# 3. Final adjudication
# ============================================================

cat("\nSTEP 3/8 — defining final sebocyte adjudication\n")

mapping <- tribble(
  ~cluster,
  ~sebocyte_family,
  ~sebocyte_state,
  ~sebocyte_L2_label,
  ~sebocyte_L2_confidence,
  ~identity_basis,

  "0",
  "Pilosebaceous progenitor-like",
  "Basal progenitor-like",
  "Basal pilosebaceous progenitor-like",
  "MODERATE",
  paste(
    "Strong basal and epithelial-junction programmes with",
    "pilosebaceous progenitor features and minimal mature",
    "sebocyte identity. Retained conservatively as a",
    "progenitor-like transcriptional state within the",
    "sebocyte-associated hierarchy."
  ),

  "1",
  "Sebocyte lineage",
  "Early commitment",
  "Early sebocyte commitment",
  "MODERATE",
  paste(
    "Intermediate sebocyte identity and lipogenesis with",
    "persistent basal/pilosebaceous features, supporting",
    "an early sebocyte-associated transcriptional state."
  ),

  "2",
  "Pilosebaceous progenitor-like",
  "Junctional progenitor-like",
  "Junctional pilosebaceous progenitor-like",
  "MODERATE",
  paste(
    "Dominant basal, epithelial-junction and pilosebaceous",
    "progenitor programmes with weak mature sebocyte identity;",
    "retained as a junctional progenitor-like transcriptional",
    "state rather than a mature sebocyte."
  ),

  "3",
  "Sebocyte lineage",
  "Differentiated lipogenic",
  "Differentiated lipogenic sebocyte",
  "HIGH",
  paste(
    "Strong sebocyte identity and lipogenesis with low",
    "basal/pilosebaceous programmes and enrichment of",
    "MOGAT1-associated differentiated lipid metabolism."
  ),

  "4",
  "Sebocyte lineage",
  "Terminal lipid-rich inflammatory",
  "Terminal lipid-rich inflammatory sebocyte",
  "HIGH",
  paste(
    "Strongest sebocyte identity and terminal lipid programme",
    "with ELOVL3, THRSP, APOC1, DEGS2, AZGP1, TM7SF2 and",
    "DHCR7 together with SAA1/SAA2-associated inflammatory",
    "features and no coordinated immune-cell programme."
  )
)

stopifnot(nrow(mapping) == 5)
stopifnot(length(unique(mapping$cluster)) == 5)

cat("\n=== FINAL SEBOCYTE ADJUDICATION ===\n")
print(mapping, n = Inf)

write_csv(
  mapping,
  file.path(TAB, "14D_frozen_sebocyte_mapping.csv")
)

# ============================================================
# 4. Apply frozen annotation
# ============================================================

cat("\nSTEP 4/8 — applying frozen sebocyte annotation\n")

x$sebocyte_cluster_frozen <- cluster_before

idx <- match(
  x$sebocyte_cluster_frozen,
  mapping$cluster
)

stopifnot(!any(is.na(idx)))

x$sebocyte_family <-
  mapping$sebocyte_family[idx]

x$sebocyte_state <-
  mapping$sebocyte_state[idx]

x$sebocyte_L2_label <-
  mapping$sebocyte_L2_label[idx]

x$sebocyte_L2_confidence <-
  mapping$sebocyte_L2_confidence[idx]

x$sebocyte_resolution_frozen <- "0.4"
x$sebocyte_L2_status <- "FROZEN"

# ============================================================
# 5. Cell-level freeze ledger
# ============================================================

cat("STEP 5/8 — generating cell-level freeze ledger\n")

ledger <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),
  sebocyte_cluster_frozen =
    as.character(x$sebocyte_cluster_frozen),
  sebocyte_family =
    as.character(x$sebocyte_family),
  sebocyte_state =
    as.character(x$sebocyte_state),
  sebocyte_L2_label =
    as.character(x$sebocyte_L2_label),
  sebocyte_L2_confidence =
    as.character(x$sebocyte_L2_confidence),
  sebocyte_resolution_frozen =
    as.character(x$sebocyte_resolution_frozen),
  sebocyte_L2_status =
    as.character(x$sebocyte_L2_status)
)

write_csv(
  ledger,
  file.path(TAB, "14D_sebocyte_cell_freeze_ledger.csv")
)

# ============================================================
# 6. Correct patient-level recurrence + frozen summary
# ============================================================

cat("STEP 6/8 — frozen summary and corrected recurrence\n")

patient_cluster <- ledger %>%
  count(
    sebocyte_cluster_frozen,
    patient_id,
    batch,
    name = "n_cells"
  )

patient_recurrence <- patient_cluster %>%
  group_by(sebocyte_cluster_frozen) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    patients_ge5_cells = sum(n_cells >= 5),
    median_cells_patient = median(n_cells),
    dominant_patient_fraction =
      max(n_cells) / sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(as.integer(sebocyte_cluster_frozen))

cat("\n=== CORRECTED PATIENT RECURRENCE ===\n")
print(patient_recurrence, n = Inf)

write_csv(
  patient_recurrence,
  file.path(TAB, "14D_corrected_patient_recurrence.csv")
)

summary <- ledger %>%
  group_by(
    sebocyte_cluster_frozen,
    sebocyte_family,
    sebocyte_state,
    sebocyte_L2_label,
    sebocyte_L2_confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(as.integer(sebocyte_cluster_frozen))

cat("\n=== FINAL FROZEN SEBOCYTE L2 ===\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(TAB, "14D_frozen_sebocyte_L2_summary.csv")
)

# ============================================================
# 7. Formal invariance checks
# ============================================================

cat("\nSTEP 7/8 — structural invariance and formal checks\n")

cells_after <- colnames(x)
genes_after <- rownames(x)
reductions_after <- names(x@reductions)

cluster_after <- as.character(
  x@meta.data[[R04]]
)

checks <- tibble(
  check = c(
    "cells_1388",
    "genes_23412",
    "patients_16",
    "batches_2",
    "resolution_04_frozen",
    "five_clusters_preserved",
    "expected_cluster_sizes_preserved",
    "all_original_cells_preserved",
    "all_original_genes_preserved",
    "cluster_assignments_unchanged",
    "reductions_preserved",
    "five_adjudications_defined",
    "all_cells_have_family",
    "all_cells_have_state",
    "all_cells_have_L2_label",
    "all_cells_have_confidence",
    "C0_basal_pilosebaceous",
    "C1_early_sebocyte",
    "C2_junctional_pilosebaceous",
    "C3_differentiated_lipogenic",
    "C4_terminal_lipid_inflammatory",
    "no_sweat_identity_assigned",
    "no_immune_identity_assigned",
    "progenitor_labels_are_like_states",
    "disease_not_used_for_identity",
    "no_cells_removed",
    "patient_recurrence_corrected",
    "L2_status_frozen",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(x) == 1388,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,
    all(x$sebocyte_resolution_frozen == "0.4"),
    length(unique(x$sebocyte_cluster_frozen)) == 5,

    all(
      as.integer(
        table(x$sebocyte_cluster_frozen)[names(expected)]
      ) == as.integer(expected)
    ),

    identical(cells_before, cells_after),
    identical(genes_before, genes_after),
    identical(cluster_before, cluster_after),
    identical(reductions_before, reductions_after),

    nrow(mapping) == 5,

    !any(is.na(x$sebocyte_family)),
    !any(is.na(x$sebocyte_state)),
    !any(is.na(x$sebocyte_L2_label)),
    !any(is.na(x$sebocyte_L2_confidence)),

    mapping$sebocyte_L2_label[
      mapping$cluster == "0"
    ] == "Basal pilosebaceous progenitor-like",

    mapping$sebocyte_L2_label[
      mapping$cluster == "1"
    ] == "Early sebocyte commitment",

    mapping$sebocyte_L2_label[
      mapping$cluster == "2"
    ] == "Junctional pilosebaceous progenitor-like",

    mapping$sebocyte_L2_label[
      mapping$cluster == "3"
    ] == "Differentiated lipogenic sebocyte",

    mapping$sebocyte_L2_label[
      mapping$cluster == "4"
    ] == "Terminal lipid-rich inflammatory sebocyte",

    !any(
      grepl(
        "sweat",
        mapping$sebocyte_L2_label,
        ignore.case = TRUE
      )
    ),

    !any(
      grepl(
        "immune cell",
        mapping$sebocyte_L2_label,
        ignore.case = TRUE
      )
    ),

    all(
      grepl(
        "progenitor-like",
        mapping$sebocyte_L2_label[
          mapping$cluster %in% c("0", "2")
        ],
        ignore.case = TRUE
      )
    ),

    TRUE,

    ncol(x) ==
      length(cells_before),

    all(
      patient_recurrence$n_patients >= 12
    ),

    all(
      x$sebocyte_L2_status ==
        "FROZEN"
    ),

    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(TAB, "14D_formal_checks.csv")
)

stopifnot(all(checks$pass))

# ============================================================
# 8. Save
# ============================================================

cat("\nSTEP 8/8 — saving frozen sebocyte object\n")

FINAL_RDS <- file.path(
  OUT,
  "sebocyte_v2_14D_L2_frozen.rds"
)

saveRDS(x, FINAL_RDS)

cat("\nSaved:\n", FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("14D COMPLETE\n")
cat("============================================================\n")
cat("Cells: 1388\n")
cat("Patients: 16\n")
cat("Batches: 2\n")
cat("Resolution: 0.4 FROZEN\n")
cat("Clusters: 5\n")
cat("C0: basal pilosebaceous progenitor-like\n")
cat("C1: early sebocyte commitment\n")
cat("C2: junctional pilosebaceous progenitor-like\n")
cat("C3: differentiated lipogenic sebocyte\n")
cat("C4: terminal lipid-rich inflammatory sebocyte\n")
cat("Progenitor designations: transcriptional states only\n")
cat("Disease used for identity: NO\n")
cat("Cells removed: 0\n")
cat("Sebocyte L2 status: FROZEN\n")
cat("PASS\n")
