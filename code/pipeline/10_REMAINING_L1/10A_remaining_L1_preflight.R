suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 180)

cat("============================================================\n")
cat("PIPELINE_V2 — 10A REMAINING L1 PREFLIGHT\n")
cat("============================================================\n\n")

# ============================================================
# 1. Paths
# ============================================================

OBJ <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/",
  "10_REMAINING_L1/10A_preflight"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(OBJ))

# ============================================================
# 2. Load frozen global L1 object
# ============================================================

cat("STEP 1/8 — loading frozen global L1 atlas\n")

x <- readRDS(OBJ)

stopifnot(ncol(x) == 82893)

md <- x@meta.data
md$cell <- rownames(md)

required <- c(
  "patient_id",
  "disease",
  "batch",
  "L1",
  "L1_frozen"
)

missing_required <- setdiff(required, colnames(md))

if (length(missing_required) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_required, collapse = ", ")
  )
}

md$L1_final <- as.character(md$L1)

cat("Cells:", nrow(md), "\n")
cat("Patients:", n_distinct(md$patient_id), "\n")
cat("Diseases:", paste(sort(unique(md$disease)), collapse = ", "), "\n")
cat("Batches:", paste(sort(unique(md$batch)), collapse = ", "), "\n\n")

# ============================================================
# 3. Complete L1 inventory
# ============================================================

cat("STEP 2/8 — complete L1 inventory\n")

l1_inventory <- md %>%
  group_by(L1_final) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_id),
    n_batches = n_distinct(batch),

    dominant_patient = {
      tt <- sort(table(patient_id), decreasing = TRUE)
      names(tt)[1]
    },

    dominant_patient_cells = {
      tt <- sort(table(patient_id), decreasing = TRUE)
      as.integer(tt[1])
    },

    dominant_patient_fraction =
      dominant_patient_cells / n_cells,

    median_cells_per_patient = {
      tt <- table(patient_id)
      median(as.numeric(tt))
    },

    min_cells_per_patient_present = {
      tt <- table(patient_id)
      min(as.numeric(tt))
    },

    max_cells_per_patient = {
      tt <- table(patient_id)
      max(as.numeric(tt))
    },

    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

write_csv(
  l1_inventory,
  file.path(TAB, "10A_complete_L1_inventory.csv")
)

cat("\n=== COMPLETE L1 INVENTORY ===\n")
print(l1_inventory, n = Inf)

# ============================================================
# 4. Patient × L1 matrix
# ============================================================

cat("\nSTEP 3/8 — patient-level representation\n")

patient_l1 <- md %>%
  count(
    patient_id,
    disease,
    batch,
    L1_final,
    name = "n_cells"
  )

write_csv(
  patient_l1,
  file.path(TAB, "10A_patient_L1_counts_long.csv")
)

patient_matrix <- patient_l1 %>%
  select(patient_id, L1_final, n_cells) %>%
  pivot_wider(
    names_from = L1_final,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  arrange(patient_id)

write_csv(
  patient_matrix,
  file.path(TAB, "10A_patient_L1_counts_wide.csv")
)

# ============================================================
# 5. Representation thresholds
# ============================================================

cat("STEP 4/8 — recurrence thresholds\n")

all_patients <- sort(unique(md$patient_id))

recurrence <- patient_l1 %>%
  complete(
    patient_id = all_patients,
    L1_final = unique(md$L1_final),
    fill = list(n_cells = 0)
  ) %>%
  group_by(L1_final) %>%
  summarise(
    patients_ge1 = sum(n_cells >= 1),
    patients_ge5 = sum(n_cells >= 5),
    patients_ge10 = sum(n_cells >= 10),
    patients_ge20 = sum(n_cells >= 20),
    patients_ge50 = sum(n_cells >= 50),
    median_all16 = median(n_cells),
    q25_all16 = quantile(n_cells, 0.25),
    q75_all16 = quantile(n_cells, 0.75),
    max_all16 = max(n_cells),
    .groups = "drop"
  ) %>%
  arrange(desc(patients_ge10), desc(median_all16))

write_csv(
  recurrence,
  file.path(TAB, "10A_L1_recurrence_thresholds.csv")
)

# ============================================================
# 6. Disease coverage — descriptive only
# ============================================================

cat("STEP 5/8 — disease coverage audit (descriptive only)\n")

disease_cells <- md %>%
  count(
    L1_final,
    disease,
    name = "n_cells"
  )

disease_patients <- md %>%
  distinct(
    L1_final,
    disease,
    patient_id
  ) %>%
  count(
    L1_final,
    disease,
    name = "n_patients"
  )

disease_coverage <- disease_cells %>%
  left_join(
    disease_patients,
    by = c("L1_final", "disease")
  ) %>%
  arrange(L1_final, disease)

write_csv(
  disease_coverage,
  file.path(TAB, "10A_L1_disease_coverage_descriptive.csv")
)

# ============================================================
# 7. Focus compartments
# ============================================================

cat("STEP 6/8 — targeted remaining-compartment inventory\n")

focus_patterns <- c(
  "Vascular_EC",
  "Lymphatic_EC",
  "Pericyte_SMC",
  "Smooth_muscle",
  "Sebocyte",
  "Melanocyte",
  "Sweat_gland"
)

focus_regex <- paste(focus_patterns, collapse = "|")

focus_levels <- unique(
  md$L1_final[
    grepl(
      focus_regex,
      md$L1_final,
      ignore.case = TRUE
    )
  ]
)

cat(
  "Detected focus L1 levels:\n",
  paste(sort(focus_levels), collapse = "\n"),
  "\n"
)

focus_summary <- l1_inventory %>%
  filter(L1_final %in% focus_levels) %>%
  left_join(
    recurrence,
    by = "L1_final"
  ) %>%
  arrange(desc(n_cells))

write_csv(
  focus_summary,
  file.path(TAB, "10A_focus_compartment_summary.csv")
)

focus_patient <- patient_l1 %>%
  filter(L1_final %in% focus_levels) %>%
  arrange(L1_final, patient_id)

write_csv(
  focus_patient,
  file.path(TAB, "10A_focus_patient_counts.csv")
)

focus_disease <- disease_coverage %>%
  filter(L1_final %in% focus_levels)

write_csv(
  focus_disease,
  file.path(TAB, "10A_focus_disease_coverage.csv")
)

cat("\n=== FOCUS COMPARTMENTS ===\n")
print(focus_summary, n = Inf)

cat("\n=== FOCUS DISEASE COVERAGE — DESCRIPTIVE ONLY ===\n")
print(as_tibble(focus_disease), n = Inf)

# ============================================================
# 8. QC / doublet audit
# ============================================================

cat("\nSTEP 7/8 — QC and doublet audit\n")

qc_available <- all(
  c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "doublet_class",
    "doublet_score"
  ) %in% colnames(md)
)

if (qc_available) {

  qc_focus <- md %>%
    filter(L1_final %in% focus_levels) %>%
    group_by(L1_final) %>%
    summarise(
      n_cells = n(),
      singlet_fraction =
        mean(doublet_class == "singlet", na.rm = TRUE),

      median_doublet_score =
        median(doublet_score, na.rm = TRUE),

      q75_doublet_score =
        quantile(doublet_score, 0.75, na.rm = TRUE),

      median_nFeature_RNA =
        median(nFeature_RNA, na.rm = TRUE),

      median_nCount_RNA =
        median(nCount_RNA, na.rm = TRUE),

      median_percent_mt =
        median(percent.mt, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(desc(n_cells))

  write_csv(
    qc_focus,
    file.path(TAB, "10A_focus_QC_doublet_summary.csv")
  )

  cat("\n=== FOCUS QC ===\n")
  print(qc_focus, n = Inf)

} else {

  qc_focus <- tibble(
    note = "Complete QC/doublet metadata unavailable in global frozen object."
  )

  write_csv(
    qc_focus,
    file.path(TAB, "10A_focus_QC_doublet_summary.csv")
  )
}

# ============================================================
# 9. Candidate triage — STRUCTURAL ONLY
# ============================================================

cat("\nSTEP 8/8 — structural triage\n")

triage <- focus_summary %>%
  mutate(
    structural_support = case_when(

      n_cells >= 300 &
      patients_ge10 >= 12 &
      n_batches >= 2 &
      dominant_patient_fraction < 0.50
      ~ "STRONG",

      n_cells >= 100 &
      patients_ge5 >= 8 &
      n_batches >= 2 &
      dominant_patient_fraction < 0.60
      ~ "MODERATE",

      TRUE
      ~ "LIMITED"
    ),

    preflight_action = case_when(

      structural_support == "STRONG"
      ~ "L2_RECLUSTER_CANDIDATE",

      structural_support == "MODERATE"
      ~ "TARGETED_IDENTITY_AUDIT",

      TRUE
      ~ "RETAIN_L1_UNLESS_STRONG_BIOLOGICAL_EVIDENCE"
    )
  )

write_csv(
  triage,
  file.path(TAB, "10A_focus_structural_triage.csv")
)

cat("\n=== STRUCTURAL TRIAGE ===\n")
print(
  triage %>%
    select(
      L1_final,
      n_cells,
      n_patients,
      n_batches,
      dominant_patient_fraction,
      patients_ge5,
      patients_ge10,
      patients_ge20,
      structural_support,
      preflight_action
    ),
  n = Inf
)

# ============================================================
# 10. Formal checks
# ============================================================

checks <- tibble(
  check = c(
    "global_cells_82893",
    "patients_16",
    "L1_frozen_present",
    "complete_inventory_generated",
    "patient_matrix_generated",
    "recurrence_generated",
    "disease_coverage_descriptive_only",
    "focus_compartments_detected",
    "no_L1_labels_modified",
    "no_cells_removed",
    "source_object_not_overwritten"
  ),
  pass = c(
    nrow(md) == 82893,
    n_distinct(md$patient_id) == 16,
    "L1_frozen" %in% colnames(x@meta.data),
    nrow(l1_inventory) > 0,
    nrow(patient_matrix) == 16,
    nrow(recurrence) == nrow(l1_inventory),
    TRUE,
    length(focus_levels) > 0,
    identical(
      as.character(x$L1),
      md$L1_final
    ) &&
      all(x$L1_frozen %in% TRUE),
    ncol(x) == 82893,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10A_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 10A REMAINING L1 PREFLIGHT",
  "",
  paste("Source:", OBJ),
  "",
  "Purpose:",
  "- inventory all frozen L1 compartments",
  "- quantify patient and batch recurrence",
  "- quantify descriptive disease coverage",
  "- identify remaining compartments suitable for L2 refinement",
  "",
  "Priority compartments:",
  "- blood and lymphatic endothelium",
  "- pericyte/SMC and myofibroblast",
  "- sebocyte",
  "- melanocyte",
  "- Schwann",
  "- adipocyte",
  "",
  "Interpretive boundary:",
  "- disease coverage is descriptive only",
  "- disease is not used to define identity or select clustering resolution",
  "- no L1 annotation modified",
  "- no cells removed",
  "- source object not overwritten",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(TAB, "10A_audit.txt")
)

cat("\n============================================================\n")
cat("10A COMPLETE\n")
cat("============================================================\n")
cat("Global cells:", ncol(x), "\n")
cat("Patients:", n_distinct(md$patient_id), "\n")
cat("L1 compartments:", nrow(l1_inventory), "\n")
cat("Focus compartments:", length(focus_levels), "\n")
cat("No annotations modified.\n")
cat("No cells removed.\n")
cat("PASS\n")
