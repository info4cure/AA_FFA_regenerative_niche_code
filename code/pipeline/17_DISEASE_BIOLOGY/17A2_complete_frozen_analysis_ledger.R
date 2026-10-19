suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 17A2 COMPLETE FROZEN ANALYSIS LEDGER\n")
cat("============================================================\n")

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17A_frozen_analysis_atlas/tables"
)

master <- readRDS(
  file.path(OUT, "17A_master_ledger_phase1.rds")
)

stopifnot(nrow(master) == 82893)
stopifnot(n_distinct(master$cell) == 82893)

apply_l2 <- function(
  x, cells, compartment, label, module,
  family = NULL, state = NULL, confidence = NULL
) {
  stopifnot(!anyDuplicated(cells))
  idx <- match(cells, x$cell)
  stopifnot(!anyNA(idx))

  x$final_compartment[idx] <- compartment
  x$final_L2_label[idx] <- as.character(label)
  x$final_L2_status[idx] <- "FROZEN"
  x$annotation_module[idx] <- module
  x$analysis_eligible_L2[idx] <- TRUE

  if (!is.null(family))
    x$final_L2_family[idx] <- as.character(family)

  if (!is.null(state))
    x$final_L2_state[idx] <- as.character(state)

  if (!is.null(confidence))
    x$final_L2_confidence[idx] <- as.character(confidence)

  x
}

return_cells <- function(
  x, cells, compartment, source, action, module
) {
  cells <- unique(cells)
  idx <- match(cells, x$cell)
  stopifnot(!anyNA(idx))

  x$final_compartment[idx] <- compartment
  x$final_L2_family[idx] <- NA_character_
  x$final_L2_state[idx] <- NA_character_
  x$final_L2_label[idx] <- NA_character_
  x$final_L2_confidence[idx] <- NA_character_
  x$final_L2_status[idx] <- "UNADJUDICATED_AFTER_RETURN"

  x$annotation_module[idx] <- module
  x$disposition[idx] <- action
  x$return_source[idx] <- source

  x$analysis_eligible_global[idx] <- TRUE
  x$analysis_eligible_L1[idx] <- TRUE
  x$analysis_eligible_L2[idx] <- FALSE

  x
}

exclude_cells <- function(
  x, cells, module, reason
) {
  cells <- unique(cells)
  idx <- match(cells, x$cell)
  stopifnot(!anyNA(idx))

  x$final_compartment[idx] <- "EXCLUDED"
  x$final_L2_family[idx] <- NA_character_
  x$final_L2_state[idx] <- NA_character_
  x$final_L2_label[idx] <- NA_character_
  x$final_L2_confidence[idx] <- NA_character_
  x$final_L2_status[idx] <- "EXCLUDED"

  x$annotation_module[idx] <- module
  x$disposition[idx] <- "EXCLUDE"
  x$exclusion_reason[idx] <- reason

  x$analysis_eligible_global[idx] <- FALSE
  x$analysis_eligible_L1[idx] <- FALSE
  x$analysis_eligible_L2[idx] <- FALSE

  x
}

cat("\nSTEP 1/9 — immune frozen L2\n")

imm <- readRDS(
  "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds"
)

stopifnot(ncol(imm) == 4542)

im <- imm@meta.data %>%
  rownames_to_column("cell")

stopifnot(
  all(
    c("cell", "immune_L2", "immune_L2_status") %in%
      colnames(im)
  )
)

master <- apply_l2(
  master,
  im$cell,
  "Immune",
  im$immune_L2,
  "08_IMMUNE_L2"
)

cat("Immune L2:", nrow(im), "\n")


cat("\nSTEP 2/9 — fibroblast frozen L2\n")

fib <- readRDS(
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "fibroblast_v2_09H5B_L2_frozen_4685.rds"
  )
)

stopifnot(ncol(fib) == 4685)

fm <- fib@meta.data %>%
  rownames_to_column("cell")

label_col <- intersect(
  c("fibroblast_L2_label", "fibroblast_L2", "L2_label"),
  colnames(fm)
)

stopifnot(length(label_col) == 1)

family_col <- intersect(
  c("FibConsensus_family", "fibroblast_family"),
  colnames(fm)
)

state_col <- intersect(
  c("de_novo_state", "fibroblast_state"),
  colnames(fm)
)

confidence_col <- intersect(
  c("confidence", "fibroblast_L2_confidence"),
  colnames(fm)
)

master <- apply_l2(
  master,
  fm$cell,
  "Fibroblast",
  fm[[label_col]],
  "09_FIBROBLAST_L2",
  family = if (length(family_col) == 1) fm[[family_col]] else NULL,
  state = if (length(state_col) == 1) fm[[state_col]] else NULL,
  confidence =
    if (length(confidence_col) == 1)
      fm[[confidence_col]]
    else NULL
)

cat("Fibroblast L2:", nrow(fm), "\n")

cat("\nSTEP 3/9 — vascular/lymphatic frozen L2\n")

vas <- readRDS(
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "vascular_v2_10B6_L2_frozen_3023.rds"
  )
)

stopifnot(ncol(vas) == 3023)

vm <- vas@meta.data %>%
  rownames_to_column("cell")

master <- apply_l2(
  master,
  vm$cell,
  "Endothelial",
  vm$vascular_L2_label,
  "10_VASCULAR_L2",
  family = vm$vascular_family,
  state = vm$vascular_state,
  confidence = vm$vascular_L2_confidence
)

cat("Endothelial L2:", nrow(vm), "\n")


cat("\nSTEP 4/9 — mural frozen L2\n")

mur <- readRDS(
  paste0(
    "PIPELINE_V2/results/annotation/11_MURAL_L2/",
    "11B5_freeze/mural_v2_11B5_L2_frozen.rds"
  )
)

stopifnot(ncol(mur) == 1603)

mm <- mur@meta.data %>%
  rownames_to_column("cell")

master <- apply_l2(
  master,
  mm$cell,
  "Mural",
  mm$mural_L2_label,
  "11_MURAL_L2",
  family = mm$mural_family,
  state = mm$mural_state,
  confidence = mm$mural_L2_confidence
)

cat("Mural L2:", nrow(mm), "\n")


cat("\nSTEP 5/9 — sweat-gland frozen L2\n")

sw <- readRDS(
  paste0(
    "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
    "12D_freeze/sweat_gland_v2_12D_L2_frozen.rds"
  )
)

stopifnot(ncol(sw) == 1843)

sm <- sw@meta.data %>%
  rownames_to_column("cell")

master <- apply_l2(
  master,
  sm$cell,
  "Sweat_gland",
  sm$sweat_L2_label,
  "12_SWEAT_GLAND_L2",
  family = sm$sweat_family,
  state = sm$sweat_state,
  confidence = sm$sweat_L2_confidence
)

cat("Sweat L2:", nrow(sm), "\n")


cat("\nSTEP 6/9 — melanocyte frozen L2\n")

mel <- readRDS(
  paste0(
    "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
    "13D_freeze/melanocyte_v2_13D_L2_frozen.rds"
  )
)

stopifnot(ncol(mel) == 1517)

mem <- mel@meta.data %>%
  rownames_to_column("cell")

master <- apply_l2(
  master,
  mem$cell,
  "Melanocyte",
  mem$melanocyte_L2_label,
  "13_MELANOCYTE_L2",
  family = mem$melanocyte_family,
  state = mem$melanocyte_state,
  confidence = mem$melanocyte_L2_confidence
)

cat("Melanocyte L2:", nrow(mem), "\n")


cat("\nSTEP 7/9 — sebocyte frozen L2\n")

seb <- readRDS(
  paste0(
    "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
    "14D_freeze/sebocyte_v2_14D_L2_frozen.rds"
  )
)

stopifnot(ncol(seb) == 1388)

sem <- seb@meta.data %>%
  rownames_to_column("cell")

master <- apply_l2(
  master,
  sem$cell,
  "Sebocyte",
  sem$sebocyte_L2_label,
  "14_SEBOCYTE_L2",
  family = sem$sebocyte_family,
  state = sem$sebocyte_state,
  confidence = sem$sebocyte_L2_confidence
)

cat("Sebocyte L2:", nrow(sem), "\n")

cat("\nSTEP 8/9 — final returns and exclusions\n")

# IMMUNE
idisp <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/08_IMMUNE_L2/",
    "08G_refined_universe/tables/",
    "08G_master_adjudication_ledger.csv"
  ),
  show_col_types = FALSE
)

master <- return_cells(
  master,
  idisp$cell[idisp$final_08G_disposition == "RETURN_EPITHELIAL"],
  "Epithelial",
  "Immune",
  "RETURN_EPITHELIAL",
  "08_IMMUNE_L2"
)

master <- return_cells(
  master,
  idisp$cell[idisp$final_08G_disposition == "RETURN_NEURAL_SCHWANN"],
  "Neural_Schwann",
  "Immune",
  "RETURN_NEURAL_SCHWANN",
  "08_IMMUNE_L2"
)

master <- exclude_cells(
  master,
  idisp$cell[idisp$final_08G_disposition == "EXCLUDE_MIXED_LINEAGE"],
  "08_IMMUNE_L2",
  "IMMUNE_MIXED_LINEAGE"
)

master <- exclude_cells(
  master,
  idisp$cell[idisp$final_08G_disposition == "EXCLUDE_UNRESOLVED"],
  "08_IMMUNE_L2",
  "IMMUNE_UNRESOLVED"
)

# FIBROBLAST EARLY
fdisp <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
    "09E1_final_adjudication/tables/",
    "09E1_master_adjudication.csv"
  ),
  show_col_types = FALSE
)

stopifnot(
  all(
    c("cell", "proposed_disposition") %in%
      colnames(fdisp)
  )
)

stopifnot(
  sum(fdisp$proposed_disposition == "EXCLUDE_DOUBLET") == 85,
  sum(fdisp$proposed_disposition == "RETURN_EPITHELIAL") == 8,
  sum(fdisp$proposed_disposition == "RETURN_T") == 1
)

master <- exclude_cells(
  master,
  fdisp$cell[fdisp$proposed_disposition == "EXCLUDE_DOUBLET"],
  "09_FIBROBLAST_L2",
  "FIBROBLAST_MIXED_HIGH_DOUBLET"
)

master <- return_cells(
  master,
  fdisp$cell[fdisp$proposed_disposition == "RETURN_EPITHELIAL"],
  "Epithelial",
  "Fibroblast",
  "RETURN_EPITHELIAL",
  "09_FIBROBLAST_L2"
)

master <- return_cells(
  master,
  fdisp$cell[fdisp$proposed_disposition == "RETURN_T"],
  "Immune",
  "Fibroblast",
  "RETURN_T",
  "09_FIBROBLAST_L2"
)

# FIBROBLAST LATE
flate <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
    "09H4A_C9_return_epithelial/tables/",
    "09H4A_C9_return_ledger.csv"
  ),
  show_col_types = FALSE
)

stopifnot(n_distinct(flate$cell) == 29)

master <- return_cells(
  master,
  flate$cell,
  "Epithelial",
  "Fibroblast",
  "RETURN_EPITHELIAL",
  "09_FIBROBLAST_L2"
)

# VASCULAR
vex <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/10_VASCULAR_L2/",
    "10B4A_C8_mixed_EC_T/tables/",
    "10B4A_C8_cell_adjudication_ledger.csv"
  ),
  show_col_types = FALSE
)

stopifnot(n_distinct(vex$cell) == 17)

master <- exclude_cells(
  master,
  vex$cell,
  "10_VASCULAR_L2",
  "MIXED_EC_T_PROBABLE_RESIDUAL_DOUBLET"
)

# MELANOCYTE -> EPITHELIAL
mret <- read_csv(
  paste0(
    "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
    "13C_C4_boundary/tables/",
    "13C_return_to_epithelial_ledger.csv"
  ),
  show_col_types = FALSE
)

stopifnot(n_distinct(mret$cell) == 39)

master <- return_cells(
  master,
  mret$cell,
  "Epithelial",
  "Melanocyte",
  "RETURN_EPITHELIAL",
  "13_MELANOCYTE_L2"
)

cat("Returns/exclusions applied.\n")

cat("\nSTEP 9/9 — conservation and analysis eligibility\n")

eligible <- master %>%
  filter(analysis_eligible_global)

excluded <- master %>%
  filter(!analysis_eligible_global)

cat("\n===== GLOBAL ACCOUNTING =====\n")
cat("Original:", nrow(master), "\n")
cat("Eligible global:", nrow(eligible), "\n")
cat("Excluded:", nrow(excluded), "\n")
cat("Conservation:", nrow(eligible) + nrow(excluded), "\n")

cat("\n===== FINAL L1 COUNTS =====\n")

final_L1 <- master %>%
  count(
    final_compartment,
    analysis_eligible_L1,
    sort = TRUE
  )

print(as.data.frame(final_L1), row.names = FALSE)

cat("\n===== DISPOSITIONS =====\n")

disp <- master %>%
  count(
    disposition,
    return_source,
    exclusion_reason,
    sort = TRUE
  )

print(as.data.frame(disp), row.names = FALSE)

cat("\n===== L2 ACCOUNTING =====\n")

l2_summary <- master %>%
  filter(analysis_eligible_global) %>%
  count(
    final_compartment,
    final_L2_status,
    analysis_eligible_L2,
    sort = TRUE
  )

print(as.data.frame(l2_summary), row.names = FALSE)

cat("\n===== DISEASE × PATIENT =====\n")

patient_summary <- master %>%
  filter(analysis_eligible_global) %>%
  count(
    disease,
    patient_id,
    name = "n_cells"
  ) %>%
  arrange(
    factor(disease, levels = c("HC", "FFA", "AA")),
    patient_id
  )

print(as.data.frame(patient_summary), row.names = FALSE)

# ============================================================
# Formal checks
# ============================================================

n_L2 <- function(compartment) {
  sum(
    master$final_compartment == compartment &
    master$analysis_eligible_L2,
    na.rm = TRUE
  )
}

n_return <- function(source, action) {
  sum(
    master$return_source == source &
    master$disposition == action,
    na.rm = TRUE
  )
}

checks <- tibble(
  check = c(
    "master_cells_82893",
    "unique_cells_82893",
    "excluded_213",
    "eligible_global_82680",
    "conservation_82893",

    "immune_frozen_4542",
    "fibroblast_frozen_4685",
    "endothelial_frozen_3023",
    "mural_frozen_1603",
    "sweat_frozen_1843",
    "melanocyte_frozen_1517",
    "sebocyte_frozen_1388",
    "epithelial_original_L2_63893",

    "immune_to_epithelial_36",
    "immune_to_neural_73",
    "fibroblast_to_epithelial_37",
    "fibroblast_to_T_1",
    "melanocyte_to_epithelial_39",

    "all_excluded_not_L1_eligible",
    "all_excluded_not_L2_eligible",
    "all_returns_L1_eligible",
    "all_returns_not_L2_eligible",

    "patients_16",
    "diseases_AA_FFA_HC",
    "no_duplicate_cells"
  ),

  pass = c(
    nrow(master) == 82893,
    n_distinct(master$cell) == 82893,
    nrow(excluded) == 213,
    nrow(eligible) == 82680,
    nrow(eligible) + nrow(excluded) == 82893,

    n_L2("Immune") == 4542,
    n_L2("Fibroblast") == 4685,
    n_L2("Endothelial") == 3023,
    n_L2("Mural") == 1603,
    n_L2("Sweat_gland") == 1843,
    n_L2("Melanocyte") == 1517,
    n_L2("Sebocyte") == 1388,
    n_L2("Epithelial") == 63893,

    n_return("Immune", "RETURN_EPITHELIAL") == 36,
    n_return("Immune", "RETURN_NEURAL_SCHWANN") == 73,
    n_return("Fibroblast", "RETURN_EPITHELIAL") == 37,
    n_return("Fibroblast", "RETURN_T") == 1,
    n_return("Melanocyte", "RETURN_EPITHELIAL") == 39,

    all(!master$analysis_eligible_L1[
      !master$analysis_eligible_global
    ]),

    all(!master$analysis_eligible_L2[
      !master$analysis_eligible_global
    ]),

    all(master$analysis_eligible_L1[
      master$disposition %in%
        c(
          "RETURN_EPITHELIAL",
          "RETURN_NEURAL_SCHWANN",
          "RETURN_T"
        )
    ]),

    all(!master$analysis_eligible_L2[
      master$disposition %in%
        c(
          "RETURN_EPITHELIAL",
          "RETURN_NEURAL_SCHWANN",
          "RETURN_T"
        )
    ]),

    n_distinct(master$patient_id) == 16,
    setequal(unique(master$disease), c("AA", "FFA", "HC")),
    !anyDuplicated(master$cell)
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(as.data.frame(checks), row.names = FALSE)

if (!all(checks$pass)) {
  stop("17A2 formal checks failed")
}

# ============================================================
# Outputs
# ============================================================

write_csv(
  master,
  file.path(
    OUT,
    "17A_master_cell_ledger_FINAL.csv"
  )
)

saveRDS(
  master,
  file.path(
    OUT,
    "17A_master_cell_ledger_FINAL.rds"
  )
)

write_csv(
  final_L1,
  file.path(
    OUT,
    "17A_final_L1_counts.csv"
  )
)

write_csv(
  l2_summary,
  file.path(
    OUT,
    "17A_L2_eligibility_summary.csv"
  )
)

write_csv(
  patient_summary,
  file.path(
    OUT,
    "17A_patient_cell_counts.csv"
  )
)

write_csv(
  checks,
  file.path(
    OUT,
    "17A_formal_checks.csv"
  )
)

cat("\n============================================================\n")
cat("17A COMPLETE\n")
cat("============================================================\n")
cat("Original cells: 82,893\n")
cat("Formal exclusions: 213\n")
cat("Final global analysis universe: 82,680\n")
cat("Patient is the downstream inferential unit\n")
cat("Returned cells retained for L1 analyses\n")
cat("Returned cells excluded from L2-specific inference\n")
cat("Frozen L2 identities unchanged\n")
cat("Disease used to redefine identities: NO\n")
cat("READY FOR 17B COMPOSITION / DIFFERENTIAL ABUNDANCE\n")
cat("PASS\n")
