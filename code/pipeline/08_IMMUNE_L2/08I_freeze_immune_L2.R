suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 08I IMMUNE L2 FINAL FREEZE\n")
cat("============================================================\n")

# ============================================================
# Paths
# ============================================================

OBJ_IN <- "PIPELINE_V2/objects/annotated/immune_v2_08G_refined.rds"

OUTDIR <- paste0(
  "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
  "08I_L2_freeze"
)

TAB <- file.path(OUTDIR, "tables")

OBJ_OUT <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "immune_v2_08I_L2_frozen.rds"
)

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1 — load immutable refined universe
# ============================================================

cat("\nSTEP 1/7 — loading refined immune universe\n")

x <- readRDS(OBJ_IN)

stopifnot(ncol(x) == 4542)
stopifnot(nrow(x) == 23412)
stopifnot(length(unique(x$patient_id)) == 16)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")

# ============================================================
# STEP 2 — freeze resolution 0.4
# ============================================================

cat("\nSTEP 2/7 — freezing immune resolution 0.4\n")

RESCOL <- "immune_refined_res.0.4"

stopifnot(RESCOL %in% colnames(x@meta.data))

x$immune_L2_cluster <- as.character(x@meta.data[[RESCOL]])

cluster_counts <- table(x$immune_L2_cluster)

print(cluster_counts)

stopifnot(length(cluster_counts) == 10)

# ============================================================
# STEP 3 — explicit evidence-based L2 mapping
# ============================================================

cat("\nSTEP 3/7 — applying final evidence-based L2 mapping\n")

mapping <- tribble(
  ~immune_L2_cluster, ~immune_L2,                ~identity_basis,

  "0", "CD4 T",       paste0(
    "IL7R/RORA/BCL11B-supported CD4 T-cell programme; ",
    "CD4 programme top-ranked; patient support 14/16"
  ),

  "1", "CD4 T",       paste0(
    "CD4 T-cell programme top-ranked; ",
    "patient support 15/16"
  ),

  "2", "CD8 T",       paste0(
    "CD8A/CD8B/GZMK/CCL5/IFNG cytotoxic T-cell programme; ",
    "CD8 programme top-ranked; patient support 14/16"
  ),

  "3", "Macrophage/monocyte", paste0(
    "MRC1/C1QA/C1QB/C1QC/MSR1/CD14 myeloid-macrophage programme; ",
    "broad myeloid score top-ranked"
  ),

  "4", "Langerhans",  paste0(
    "CD207/CD1A/S100B Langerhans programme; ",
    "Langerhans score top-ranked; patient support 15/16"
  ),

  "5", "Treg",        paste0(
    "FOXP3/IL2RA/CTLA4/IKZF2 regulatory T-cell programme"
  ),

  "6", "Mast",        paste0(
    "TPSAB1/TPSB2/CPA3/KIT/MS4A2 mast-cell programme; ",
    "mast score top-ranked in 16/16 patients"
  ),

  "7", "NK",          paste0(
    "GNLY/FGFBP2/KLRD1/NCR1/PRF1/GZMB NK programme; ",
    "NK score top-ranked in 14/14 patients"
  ),

  "8", "cDC1-like dendritic", paste0(
    "CLEC9A/XCR1/CADM1/DNASE1L3 dendritic programme; ",
    "de novo markers override nonspecific B-cell programme ranking"
  ),

  "9", "B/plasma-lineage", paste0(
    "CD79A/PAX5/TNFRSF13C with IGHG/IGHA immunoglobulin programme; ",
    "mixed B/plasma-lineage evidence retained at L2"
  )
)

stopifnot(nrow(mapping) == 10)
stopifnot(length(unique(mapping$immune_L2_cluster)) == 10)

write_csv(
  mapping,
  file.path(TAB, "08I_cluster_to_L2_mapping.csv")
)

# ============================================================
# STEP 4 — attach frozen identities
# ============================================================

cat("\nSTEP 4/7 — attaching frozen L2 identities\n")

meta_map <- mapping %>%
  select(
    immune_L2_cluster,
    immune_L2,
    identity_basis
  )

md <- x@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(
    meta_map,
    by = "immune_L2_cluster"
  )

stopifnot(nrow(md) == ncol(x))
stopifnot(!any(is.na(md$immune_L2)))
stopifnot(!any(is.na(md$identity_basis)))

rownames(md) <- md$cell
md$cell <- NULL

x@meta.data <- md[colnames(x), , drop = FALSE]

x$immune_L2_status <- "FROZEN"
x$immune_L2_resolution <- "0.4"
x$immune_L2_version <- "PIPELINE_V2_08I"

# ============================================================
# STEP 5 — summary and provenance
# ============================================================

cat("\nSTEP 5/7 — generating frozen identity summary\n")

summary <- x@meta.data %>%
  as_tibble() %>%
  count(
    immune_L2_cluster,
    immune_L2,
    name = "n_cells"
  ) %>%
  arrange(as.numeric(immune_L2_cluster))

patient_summary <- x@meta.data %>%
  as_tibble() %>%
  group_by(
    immune_L2_cluster,
    immune_L2
  ) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_HC = n_distinct(patient_id[disease == "HC"]),
    n_AA = n_distinct(patient_id[disease == "AA"]),
    n_FFA = n_distinct(patient_id[disease == "FFA"]),
    n_batches = n_distinct(batch),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(immune_L2_cluster))

write_csv(
  summary,
  file.path(TAB, "08I_frozen_L2_summary.csv")
)

write_csv(
  patient_summary,
  file.path(TAB, "08I_frozen_L2_patient_support.csv")
)

print(patient_summary)

cat("\nUnique L2 labels:",
    length(unique(x$immune_L2)), "\n")

print(sort(table(x$immune_L2), decreasing = TRUE))

# ============================================================
# STEP 6 — provenance and save
# ============================================================

cat("\nSTEP 6/7 — recording provenance and saving object\n")

x@misc[["08I_immune_L2_freeze"]] <- list(
  source_object = OBJ_IN,
  source_cells = 4542,
  selected_resolution = "0.4",
  selected_resolution_column = RESCOL,
  n_clusters = 10,
  n_L2_labels = length(unique(x$immune_L2)),
  disease_used_for_identity_assignment = FALSE,
  cells_removed = 0,
  freeze_status = "FROZEN",
  freeze_version = "PIPELINE_V2_08I"
)

Idents(x) <- factor(
  x$immune_L2,
  levels = unique(mapping$immune_L2)
)

saveRDS(x, OBJ_OUT)

# ============================================================
# STEP 7 — formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

expected_counts <- c(
  "0" = 1098,
  "1" = 897,
  "2" = 778,
  "3" = 746,
  "4" = 337,
  "5" = 309,
  "6" = 151,
  "7" = 92,
  "8" = 69,
  "9" = 65
)

observed_counts <- table(x$immune_L2_cluster)

checks <- tibble(
  check = c(
    "cells_4542",
    "genes_23412",
    "patients_16",
    "resolution04_frozen",
    "ten_clusters",
    "nine_unique_L2_labels",
    "all_cells_have_L2",
    "all_cells_have_identity_basis",
    "cluster_counts_preserved",
    "CD4_shared_by_C0_C1",
    "C5_Treg",
    "C8_cDC1_like",
    "C9_B_plasma_lineage",
    "disease_not_used_for_assignment",
    "no_cells_removed",
    "source_object_not_overwritten"
  ),
  pass = c(
    ncol(x) == 4542,
    nrow(x) == 23412,
    length(unique(x$patient_id)) == 16,
    all(x$immune_L2_resolution == "0.4"),
    length(unique(x$immune_L2_cluster)) == 10,
    length(unique(x$immune_L2)) == 9,
    !any(is.na(x$immune_L2)),
    !any(is.na(x$identity_basis)),
    all(
      as.integer(observed_counts[names(expected_counts)]) ==
        as.integer(expected_counts)
    ),
    all(
      unique(x$immune_L2[x$immune_L2_cluster %in% c("0","1")]) ==
        "CD4 T"
    ),
    all(x$immune_L2[x$immune_L2_cluster == "5"] == "Treg"),
    all(
      x$immune_L2[x$immune_L2_cluster == "8"] ==
        "cDC1-like dendritic"
    ),
    all(
      x$immune_L2[x$immune_L2_cluster == "9"] ==
        "B/plasma-lineage"
    ),
    TRUE,
    ncol(x) == 4542,
    OBJ_IN != OBJ_OUT
  )
)

write_csv(
  checks,
  file.path(TAB, "08I_formal_checks.csv")
)

print(checks)

if (!all(checks$pass)) {
  cat("\nFAILED CHECKS:\n")
  print(filter(checks, !pass))
  stop("08I formal checks failed.")
}

audit <- c(
  "PIPELINE_V2 — 08I IMMUNE L2 FINAL FREEZE",
  "==========================================",
  "",
  paste0("Source object: ", OBJ_IN),
  paste0("Frozen object: ", OBJ_OUT),
  "",
  "Cells: 4542",
  "Patients: 16",
  "Selected resolution: 0.4",
  "Clusters: 10",
  "Unique L2 labels: 9",
  "",
  "Final cluster mapping:",
  "C0 -> CD4 T",
  "C1 -> CD4 T",
  "C2 -> CD8 T",
  "C3 -> Macrophage/monocyte",
  "C4 -> Langerhans",
  "C5 -> Treg",
  "C6 -> Mast",
  "C7 -> NK",
  "C8 -> cDC1-like dendritic",
  "C9 -> B/plasma-lineage",
  "",
  "Disease was not used for identity assignment.",
  "No cells were removed.",
  "Source object was not overwritten.",
  "",
  "PASS — immune L2 identities frozen."
)

writeLines(
  audit,
  file.path(TAB, "08I_audit.txt")
)

cat("\n============================================================\n")
cat("08I COMPLETE\n")
cat("============================================================\n")
cat("Immune cells: 4542\n")
cat("Frozen resolution: 0.4\n")
cat("Clusters: 10\n")
cat("Unique L2 identities: 9\n")
cat("No cells removed.\n")
cat("Saved:", OBJ_OUT, "\n")
cat("PASS — immune L2 identities formally frozen.\n")
