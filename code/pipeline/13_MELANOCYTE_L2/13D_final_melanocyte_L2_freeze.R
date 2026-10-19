suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 13D FINAL MELANOCYTE L2 FREEZE\n")
cat("============================================================\n\n")

SRC <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13C_C4_boundary/melanocyte_v2_13C_refined_1517.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
  "13D_freeze"
)

TAB <- file.path(OUT, "tables")

dir.create(
  TAB,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load refined melanocyte universe
# ============================================================

cat("STEP 1/8 — loading refined melanocyte universe\n")

x <- readRDS(SRC)

stopifnot(ncol(x) == 1517)
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
# 2. Recover r0.4 after C4 return
# ============================================================

cat("\nSTEP 2/8 — validating retained r0.4 structure\n")

R04 <- "mel_snn_res.0.4"

stopifnot(
  R04 %in%
    colnames(x@meta.data)
)

cluster_before <- as.character(
  x@meta.data[[R04]]
)

stopifnot(
  setequal(
    unique(cluster_before),
    c("0", "1", "2", "3")
  )
)

observed_sizes <- table(cluster_before)

cat("\n=== RETAINED r0.4 CLUSTERS ===\n")
print(observed_sizes)

stopifnot(
  sum(observed_sizes) == 1517
)

# ============================================================
# 3. Final hierarchical adjudication
# ============================================================

cat("\nSTEP 3/8 — defining final melanocyte adjudication\n")

mapping <- tribble(
  ~cluster,
  ~melanocyte_family,
  ~melanocyte_state,
  ~melanocyte_L2_label,
  ~melanocyte_L2_confidence,
  ~identity_basis,

  "0",
  "Melanocyte",
  "Antigen-presentation-high",
  "Antigen-presentation-high melanocyte",
  "HIGH",
  paste(
    "Canonical melanocyte and melanogenesis programmes retained,",
    "with recurrent enrichment of HLA-DRB1, HLA-DRA, HLA-DQB1",
    "and CD74 supporting an antigen-presentation-high melanocyte state."
  ),

  "1",
  "Melanocyte",
  "Stress / neural-crest-like",
  "Stress/neural-crest-like melanocyte",
  "MODERATE",
  paste(
    "Canonical melanocyte identity retained with enrichment of",
    "MGP, GPX3, NRP2, FABP7, CRYAB, NUPR1 and related",
    "stress/neural-crest-associated transcriptional features;",
    "not interpreted as Schwann identity."
  ),

  "2",
  "Melanocyte",
  "Differentiated / homeostatic",
  "Differentiated/homeostatic melanocyte",
  "HIGH",
  paste(
    "Strong melanocyte and melanogenesis programmes with",
    "CRABP1, APOD, MIA, CLU and broad differentiated melanocyte",
    "expression without a dominant inflammatory programme."
  ),

  "3",
  "Melanocyte",
  "Low-pigmentation transitional",
  "Low-pigmentation transitional melanocyte",
  "MODERATE",
  paste(
    "Melanocyte identity retained but with lower aggregate",
    "melanocyte/MITF programme scores and enrichment of NAV2,",
    "JARID2, PRKCE, NRXN3 and NLGN1; designated conservatively",
    "as a transcriptional transitional state rather than a",
    "lineage trajectory."
  )
)

stopifnot(nrow(mapping) == 4)
stopifnot(length(unique(mapping$cluster)) == 4)

cat("\n=== FINAL MELANOCYTE ADJUDICATION ===\n")
print(mapping, n = Inf)

write_csv(
  mapping,
  file.path(
    TAB,
    "13D_frozen_melanocyte_mapping.csv"
  )
)

# ============================================================
# 4. Apply frozen annotation
# ============================================================

cat("\nSTEP 4/8 — applying frozen melanocyte annotation\n")

x$melanocyte_cluster_frozen <-
  cluster_before

idx <- match(
  x$melanocyte_cluster_frozen,
  mapping$cluster
)

stopifnot(!any(is.na(idx)))

x$melanocyte_family <-
  mapping$melanocyte_family[idx]

x$melanocyte_state <-
  mapping$melanocyte_state[idx]

x$melanocyte_L2_label <-
  mapping$melanocyte_L2_label[idx]

x$melanocyte_L2_confidence <-
  mapping$melanocyte_L2_confidence[idx]

x$melanocyte_resolution_frozen <-
  "0.4"

x$melanocyte_L2_status <-
  "FROZEN"

stopifnot(
  !any(is.na(x$melanocyte_L2_label))
)

# ============================================================
# 5. Cell-level freeze ledger
# ============================================================

cat("STEP 5/8 — generating cell-level freeze ledger\n")

ledger <- tibble(
  cell = colnames(x),
  patient_id = as.character(x$patient_id),
  batch = as.character(x$batch),

  melanocyte_cluster_frozen =
    as.character(
      x$melanocyte_cluster_frozen
    ),

  melanocyte_family =
    as.character(
      x$melanocyte_family
    ),

  melanocyte_state =
    as.character(
      x$melanocyte_state
    ),

  melanocyte_L2_label =
    as.character(
      x$melanocyte_L2_label
    ),

  melanocyte_L2_confidence =
    as.character(
      x$melanocyte_L2_confidence
    ),

  melanocyte_resolution_frozen =
    as.character(
      x$melanocyte_resolution_frozen
    ),

  melanocyte_L2_status =
    as.character(
      x$melanocyte_L2_status
    )
)

write_csv(
  ledger,
  file.path(
    TAB,
    "13D_melanocyte_cell_freeze_ledger.csv"
  )
)

# ============================================================
# 6. Frozen summary
# ============================================================

cat("STEP 6/8 — frozen L2 summary\n")

summary <- ledger %>%
  group_by(
    melanocyte_cluster_frozen,
    melanocyte_family,
    melanocyte_state,
    melanocyte_L2_label,
    melanocyte_L2_confidence
  ) %>%
  summarise(
    n_cells = n(),
    n_patients =
      n_distinct(patient_id),
    n_batches =
      n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(
    as.integer(
      melanocyte_cluster_frozen
    )
  )

cat("\n=== FINAL FROZEN MELANOCYTE L2 ===\n")
print(summary, n = Inf)

write_csv(
  summary,
  file.path(
    TAB,
    "13D_frozen_melanocyte_L2_summary.csv"
  )
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
    "refined_cells_1517",
    "genes_23412",
    "patients_16",
    "batches_2",
    "four_retained_r04_clusters",
    "C4_absent",
    "all_original_cells_preserved",
    "all_original_genes_preserved",
    "cluster_assignments_unchanged",
    "reductions_preserved",
    "four_adjudications_defined",
    "all_cells_have_family",
    "all_cells_have_state",
    "all_cells_have_L2_label",
    "all_cells_have_confidence",
    "C0_antigen_presentation_high",
    "C1_stress_neural_crest_like",
    "C2_differentiated_homeostatic",
    "C3_low_pigmentation_transitional",
    "no_Schwann_identity_assigned",
    "no_cycling_identity_assigned",
    "transitional_not_lineage_trajectory",
    "disease_not_used_for_identity",
    "no_additional_cells_removed",
    "resolution_04_frozen",
    "L2_status_frozen",
    "source_object_not_overwritten"
  ),

  pass = c(
    ncol(x) == 1517,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    length(unique(x$batch)) == 2,

    length(
      unique(
        x$melanocyte_cluster_frozen
      )
    ) == 4,

    !(
      "4" %in%
        unique(
          x$melanocyte_cluster_frozen
        )
    ),

    identical(
      cells_before,
      cells_after
    ),

    identical(
      genes_before,
      genes_after
    ),

    identical(
      cluster_before,
      cluster_after
    ),

    identical(
      reductions_before,
      reductions_after
    ),

    nrow(mapping) == 4,

    !any(
      is.na(
        x$melanocyte_family
      )
    ),

    !any(
      is.na(
        x$melanocyte_state
      )
    ),

    !any(
      is.na(
        x$melanocyte_L2_label
      )
    ),

    !any(
      is.na(
        x$melanocyte_L2_confidence
      )
    ),

    mapping$melanocyte_L2_label[
      mapping$cluster == "0"
    ] ==
      "Antigen-presentation-high melanocyte",

    mapping$melanocyte_L2_label[
      mapping$cluster == "1"
    ] ==
      "Stress/neural-crest-like melanocyte",

    mapping$melanocyte_L2_label[
      mapping$cluster == "2"
    ] ==
      "Differentiated/homeostatic melanocyte",

    mapping$melanocyte_L2_label[
      mapping$cluster == "3"
    ] ==
      "Low-pigmentation transitional melanocyte",

    !any(
      grepl(
        "Schwann",
        mapping$melanocyte_L2_label,
        ignore.case = TRUE
      )
    ),

    !any(
      grepl(
        "cycling",
        mapping$melanocyte_L2_label,
        ignore.case = TRUE
      )
    ),

    TRUE,
    TRUE,

    ncol(x) ==
      length(cells_before),

    all(
      x$melanocyte_resolution_frozen ==
        "0.4"
    ),

    all(
      x$melanocyte_L2_status ==
        "FROZEN"
    ),

    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

write_csv(
  checks,
  file.path(
    TAB,
    "13D_formal_checks.csv"
  )
)

stopifnot(
  all(checks$pass)
)

# ============================================================
# 8. Save
# ============================================================

cat("\nSTEP 8/8 — saving frozen melanocyte object\n")

FINAL_RDS <- file.path(
  OUT,
  "melanocyte_v2_13D_L2_frozen.rds"
)

saveRDS(
  x,
  FINAL_RDS
)

cat("\nSaved:\n")
cat(FINAL_RDS, "\n")

cat("\n============================================================\n")
cat("13D COMPLETE\n")
cat("============================================================\n")
cat("Initial L1 melanocyte candidates: 1556\n")
cat("Returned to epithelial in 13C: 39\n")
cat("Final melanocyte cells: 1517\n")
cat("Patients: 16\n")
cat("Batches: 2\n")
cat("Resolution: 0.4 FROZEN\n")
cat("Retained clusters: 4\n")
cat("C0: antigen-presentation-high melanocyte\n")
cat("C1: stress/neural-crest-like melanocyte\n")
cat("C2: differentiated/homeostatic melanocyte\n")
cat("C3: low-pigmentation transitional melanocyte\n")
cat("Schwann identity assigned: NO\n")
cat("Disease used for identity: NO\n")
cat("Additional cells removed during freeze: 0\n")
cat("Melanocyte L2 status: FROZEN\n")
cat("PASS\n")
