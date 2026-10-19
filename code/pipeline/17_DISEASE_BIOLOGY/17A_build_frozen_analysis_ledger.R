suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 17A FROZEN DISEASE-ANALYSIS LEDGER\n")
cat("============================================================\n")

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17A_frozen_analysis_atlas/tables"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. GLOBAL SOURCE
# ============================================================

cat("\nSTEP 1/8 — global frozen source\n")

g <- readRDS(
  "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"
)

stopifnot(ncol(g) == 82893)
stopifnot(nrow(g) == 23412)

meta <- g@meta.data %>%
  rownames_to_column("cell")

needed <- c(
  "cell",
  "patient_id",
  "disease",
  "batch",
  "L1"
)

stopifnot(all(needed %in% colnames(meta)))

master <- meta %>%
  transmute(
    cell,
    patient_id = as.character(patient_id),
    disease = as.character(disease),
    batch = as.character(batch),
    L1_original = as.character(L1),

    final_compartment = as.character(L1),

    final_L2_family = NA_character_,
    final_L2_state = NA_character_,
    final_L2_label = NA_character_,
    final_L2_confidence = NA_character_,
    final_L2_status = NA_character_,

    annotation_module = NA_character_,

    disposition = "RETAIN",
    exclusion_reason = NA_character_,
    return_source = NA_character_,

    analysis_eligible_global = TRUE,
    analysis_eligible_L1 = TRUE,
    analysis_eligible_L2 = FALSE
  )

cat("Cells:", nrow(master), "\n")


# ============================================================
# helper
# ============================================================

apply_freeze <- function(
  master,
  ledger,
  family_col,
  state_col,
  label_col,
  confidence_col,
  module,
  final_compartment
) {

  stopifnot(!anyDuplicated(ledger$cell))

  idx <- match(
    ledger$cell,
    master$cell
  )

  stopifnot(!anyNA(idx))

  master$final_compartment[idx] <- final_compartment

  master$final_L2_family[idx] <-
    as.character(ledger[[family_col]])

  master$final_L2_state[idx] <-
    as.character(ledger[[state_col]])

  master$final_L2_label[idx] <-
    as.character(ledger[[label_col]])

  master$final_L2_confidence[idx] <-
    as.character(ledger[[confidence_col]])

  master$final_L2_status[idx] <- "FROZEN"
  master$annotation_module[idx] <- module
  master$analysis_eligible_L2[idx] <- TRUE

  master
}


# ============================================================
# 2. EPITHELIAL ORIGINAL FREEZE
# ============================================================

cat("\nSTEP 2/8 — epithelial frozen L2\n")

ep <- readRDS(
  "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds"
)

epm <- ep@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    family = as.character(L2_identity),
    state = as.character(L2_state),
    label = as.character(L2_identity),
    confidence = as.character(L2_confidence)
  )

stopifnot(nrow(epm) == 63893)

master <- apply_freeze(
  master,
  epm,
  "family",
  "state",
  "label",
  "confidence",
  "07_EPITHELIAL_L2",
  "Epithelial"
)


# ============================================================
# 3. IMMUNE FINAL FREEZE + DISPOSITIONS
# ============================================================

cat("\nSTEP 3/8 — immune accounting\n")

imm_disp <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
    "08G_refined_universe/tables/",
    "08G_master_adjudication_ledger.csv"
  ),
  show_col_types = FALSE
)

stopifnot(nrow(imm_disp) == 4762)
stopifnot(n_distinct(imm_disp$cell) == 4762)

# First process non-immune dispositions
for (i in seq_len(nrow(imm_disp))) {

  cell <- imm_disp$cell[i]
  action <- imm_disp$final_08G_disposition[i]

  j <- match(cell, master$cell)
  stopifnot(!is.na(j))

  if (action == "RETURN_EPITHELIAL") {

    master$final_compartment[j] <- "Epithelial"
    master$final_L2_status[j] <- "UNADJUDICATED_AFTER_RETURN"
    master$annotation_module[j] <- "08_IMMUNE_L2"
    master$disposition[j] <- "RETURN_EPITHELIAL"
    master$return_source[j] <- "Immune"
    master$analysis_eligible_L2[j] <- FALSE

  } else if (action == "RETURN_NEURAL_SCHWANN") {

    master$final_compartment[j] <- "Neural_Schwann"
    master$final_L2_status[j] <- "UNADJUDICATED_AFTER_RETURN"
    master$annotation_module[j] <- "08_IMMUNE_L2"
    master$disposition[j] <- "RETURN_NEURAL_SCHWANN"
    master$return_source[j] <- "Immune"
    master$analysis_eligible_L2[j] <- FALSE

  } else if (
    action %in%
    c(
      "EXCLUDE_MIXED_LINEAGE",
      "EXCLUDE_UNRESOLVED"
    )
  ) {

    master$final_compartment[j] <- "EXCLUDED"
    master$final_L2_status[j] <- "EXCLUDED"
    master$annotation_module[j] <- "08_IMMUNE_L2"
    master$disposition[j] <- action
    master$exclusion_reason[j] <- action
    master$analysis_eligible_global[j] <- FALSE
    master$analysis_eligible_L1[j] <- FALSE
    master$analysis_eligible_L2[j] <- FALSE
  }
}

imm <- readRDS(
  "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds"
)

cat("Frozen immune cells:", ncol(imm), "\n")

cat("Immune metadata candidates:\n")
print(
  grep(
    "family|state|L2|identity|confidence|cluster",
    colnames(imm@meta.data),
    value = TRUE,
    ignore.case = TRUE
  )
)


# ============================================================
# STOP HERE DELIBERATELY
# ============================================================

cat("\n============================================================\n")
cat("17A PHASE 1 COMPLETE\n")
cat("Master initialized: 82,893 cells\n")
cat("Epithelial original L2 applied: 63,893\n")
cat("Immune dispositions applied\n")
cat("No disease comparisons performed\n")
cat("No source objects modified\n")
cat("PASS\n")
cat("============================================================\n")

saveRDS(
  master,
  file.path(
    OUT,
    "17A_master_ledger_phase1.rds"
  )
)

write_csv(
  master,
  file.path(
    OUT,
    "17A_master_ledger_phase1.csv"
  )
)
