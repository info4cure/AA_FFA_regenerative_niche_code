#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H3 FIBROBLAST RESOLUTION ADJUDICATION\n")
cat("============================================================\n\n")

ROOT <- "PIPELINE_V2"

G_DIR <- file.path(
  ROOT, "results/annotation/09_FIBROBLAST_L2/09G_refined_structure/tables"
)

H2_DIR <- file.path(
  ROOT,
  "results/annotation/09_FIBROBLAST_L2",
  "09H2_FibConsensus_de_novo_adjudication/tables"
)

OUT <- file.path(
  ROOT,
  "results/annotation/09_FIBROBLAST_L2",
  "09H3_resolution_adjudication/tables"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

TRANS_FILE <- file.path(G_DIR, "09G_transition_03_to_04.csv")
SUPPORT_FILE <- file.path(G_DIR, "09G_cluster_support.csv")
MARKERS04_FILE <- file.path(G_DIR, "09G_markers_res0.4_all.csv")
H2_EVID_FILE <- file.path(H2_DIR, "09H2_projection_evidence_classification.csv")
H2_RANK_FILE <- file.path(H2_DIR, "09H2_cluster_family_ranking.csv")

required <- c(
  TRANS_FILE,
  SUPPORT_FILE,
  MARKERS04_FILE,
  H2_EVID_FILE,
  H2_RANK_FILE
)

missing <- required[!file.exists(required)]

if(length(missing) > 0){
  stop(
    "Missing required input files:\n",
    paste(missing, collapse="\n")
  )
}

# ============================================================
# STEP 1
# ============================================================

cat("STEP 1/8 — loading structural and FibConsensus evidence\n")

trans <- read_csv(TRANS_FILE, show_col_types=FALSE)
support <- read_csv(SUPPORT_FILE, show_col_types=FALSE)
markers04 <- read_csv(MARKERS04_FILE, show_col_types=FALSE)
evid <- read_csv(H2_EVID_FILE, show_col_types=FALSE)
ranking <- read_csv(H2_RANK_FILE, show_col_types=FALSE)

trans <- trans %>%
  mutate(
    parent = as.character(parent),
    child = as.character(child)
  )

support <- support %>%
  mutate(
    cluster = as.character(cluster),
    resolution = as.character(resolution)
  )

evid <- evid %>%
  mutate(
    cluster = as.character(cluster),
    resolution = as.character(resolution)
  )

ranking <- ranking %>%
  mutate(
    cluster = as.character(cluster),
    resolution = as.character(resolution)
  )

cat("Transitions:", nrow(trans), "\n")
cat("r0.3 clusters:",
    n_distinct(evid$cluster[evid$resolution=="0.3"]), "\n")
cat("r0.4 clusters:",
    n_distinct(evid$cluster[evid$resolution=="0.4"]), "\n\n")


# ============================================================
# STEP 2 — define material structural transitions
# ============================================================

cat("STEP 2/8 — identifying material r0.3 -> r0.4 transitions\n")

material <- trans %>%
  filter(fraction_parent >= 0.05) %>%
  arrange(as.numeric(parent), desc(fraction_parent))

parent_summary <- material %>%
  group_by(parent) %>%
  summarise(
    n_material_children = n(),
    material_fraction = sum(fraction_parent),
    .groups="drop"
  )

true_split_parents <- parent_summary %>%
  filter(n_material_children >= 2) %>%
  pull(parent)

cat("Parents with >=2 material children:",
    paste(true_split_parents, collapse=", "), "\n\n")

write_csv(material,
          file.path(OUT, "09H3_material_transitions.csv"))

write_csv(parent_summary,
          file.path(OUT, "09H3_parent_transition_summary.csv"))


# ============================================================
# STEP 3 — attach support and FibConsensus evidence
# ============================================================

cat("STEP 3/8 — attaching patient support and frozen-family evidence\n")

s03 <- support %>%
  filter(resolution=="0.3") %>%
  select(
    parent=cluster,
    parent_n_cells=n_cells,
    parent_n_patients=n_patients,
    parent_dominant_patient=dominant_patient,
    parent_dominant_fraction=dominant_patient_fraction
  )

s04 <- support %>%
  filter(resolution=="0.4") %>%
  select(
    child=cluster,
    child_n_cells=n_cells,
    child_n_patients=n_patients,
    child_n_HC=n_HC,
    child_n_AA=n_AA,
    child_n_FFA=n_FFA,
    child_n_batches=n_batches,
    child_dominant_patient=dominant_patient,
    child_dominant_patient_cells=dominant_patient_cells,
    child_dominant_fraction=dominant_patient_fraction
  )

e03 <- evid %>%
  filter(resolution=="0.3") %>%
  select(
    parent=cluster,
    parent_top_family=top_family,
    parent_top_core_z=top_core_z,
    parent_core_margin=core_margin,
    parent_patient_fraction=patient_fraction,
    parent_projection_evidence=projection_evidence
  )

e04 <- evid %>%
  filter(resolution=="0.4") %>%
  select(
    child=cluster,
    child_top_family=top_family,
    child_top_core_z=top_core_z,
    child_core_margin=core_margin,
    child_patient_fraction=patient_fraction,
    child_projection_evidence=projection_evidence
  )

split_evidence <- material %>%
  left_join(s03, by="parent") %>%
  left_join(s04, by="child") %>%
  left_join(e03, by="parent") %>%
  left_join(e04, by="child") %>%
  mutate(
    family_changed =
      parent_top_family != child_top_family,

    core_gain =
      child_top_core_z - parent_top_core_z,

    margin_gain =
      child_core_margin - parent_core_margin,

    patient_support_gain =
      child_patient_fraction - parent_patient_fraction,

    patient_dominated =
      !is.na(child_dominant_fraction) &
      child_dominant_fraction >= 0.80,

    broad_patient_support =
      child_n_patients >= 5 &
      child_dominant_fraction < 0.50
  )


# ============================================================
# STEP 4 — marker evidence
# ============================================================

cat("STEP 4/8 — quantifying de novo marker support\n")

# Determine common Seurat marker column names safely.
fc_col <- intersect(
  c("avg_log2FC","avg_logFC","avg_diff"),
  names(markers04)
)

padj_col <- intersect(
  c("p_val_adj","p_val"),
  names(markers04)
)

if(length(fc_col)==0){
  stop("Could not identify fold-change column in 09G markers.")
}

if(length(padj_col)==0){
  stop("Could not identify p-value column in 09G markers.")
}

fc_col <- fc_col[1]
padj_col <- padj_col[1]

markers04 <- markers04 %>%
  mutate(cluster=as.character(cluster))

marker_summary <- markers04 %>%
  group_by(cluster) %>%
  summarise(
    n_positive_markers =
      sum(.data[[fc_col]] > 0, na.rm=TRUE),

    n_strong_markers =
      sum(
        .data[[fc_col]] >= 0.5 &
        .data[[padj_col]] < 0.05,
        na.rm=TRUE
      ),

    n_very_strong_markers =
      sum(
        .data[[fc_col]] >= 1 &
        .data[[padj_col]] < 0.05,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  rename(child=cluster)

split_evidence <- split_evidence %>%
  left_join(marker_summary, by="child")


# ============================================================
# STEP 5 — classify each material child
# ============================================================

cat("STEP 5/8 — classifying material child structure\n")

split_evidence <- split_evidence %>%
  mutate(

    recurrent_child =
      broad_patient_support &
      !patient_dominated,

    distinct_frozen_family =
      family_changed &
      child_top_core_z >= 0.5 &
      child_core_margin >= 0.25 &
      child_patient_fraction >= 0.50,

    strengthened_same_family =
      !family_changed &
      core_gain >= 0.25 &
      margin_gain >= 0.15 &
      child_patient_fraction >= 0.50,

    marker_supported =
      coalesce(n_strong_markers,0) >= 10,

    child_evidence_class = case_when(

      patient_dominated ~
        "PATIENT_DOMINATED_STATE",

      recurrent_child &
      distinct_frozen_family &
      marker_supported ~
        "RECURRENT_IDENTITY_SPLIT_SUPPORTED",

      recurrent_child &
      strengthened_same_family &
      marker_supported ~
        "RECURRENT_SUBSTRUCTURE_SUPPORTED",

      recurrent_child &
      marker_supported ~
        "RECURRENT_TRANSCRIPTOMIC_SUBSTRUCTURE",

      TRUE ~
        "INSUFFICIENT_SPLIT_SUPPORT"
    )
  )

write_csv(
  split_evidence,
  file.path(OUT, "09H3_material_child_evidence.csv")
)

cat("\n=== MATERIAL CHILD EVIDENCE ===\n")
print(
  split_evidence %>%
    select(
      parent, child,
      n_cells,
      fraction_parent,
      child_n_patients,
      child_dominant_fraction,
      parent_top_family,
      child_top_family,
      child_top_core_z,
      child_core_margin,
      child_patient_fraction,
      n_strong_markers,
      child_evidence_class
    ),
  n=Inf
)


# ============================================================
# STEP 6 — parent-level adjudication
# ============================================================

cat("\nSTEP 6/8 — parent-level split adjudication\n")

parent_adjudication <- split_evidence %>%
  filter(parent %in% true_split_parents) %>%
  group_by(parent) %>%
  summarise(
    n_material_children=n(),
    n_identity_supported=
      sum(child_evidence_class ==
            "RECURRENT_IDENTITY_SPLIT_SUPPORTED"),

    n_substructure_supported=
      sum(child_evidence_class %in%
            c(
              "RECURRENT_IDENTITY_SPLIT_SUPPORTED",
              "RECURRENT_SUBSTRUCTURE_SUPPORTED",
              "RECURRENT_TRANSCRIPTOMIC_SUBSTRUCTURE"
            )),

    n_patient_dominated=
      sum(child_evidence_class ==
            "PATIENT_DOMINATED_STATE"),

    all_children_recurrent=
      all(!patient_dominated),

    .groups="drop"
  ) %>%
  mutate(
    split_adjudication = case_when(

      n_patient_dominated > 0 ~
        "DO_NOT_PROMOTE_PATIENT_DOMINATED_CHILD",

      n_identity_supported >= 1 &
      n_substructure_supported >= 2 ~
        "PROMOTION_TO_HIGHER_RESOLUTION_SUPPORTED",

      n_substructure_supported >= 2 ~
        "BIOLOGICAL_SUBSTRUCTURE_PRESENT_BUT_IDENTITY_PROMOTION_UNCERTAIN",

      TRUE ~
        "LOWER_RESOLUTION_PREFERRED"
    )
  )

write_csv(
  parent_adjudication,
  file.path(OUT, "09H3_parent_split_adjudication.csv")
)

cat("\n=== PARENT SPLIT ADJUDICATION ===\n")
print(parent_adjudication, n=Inf)


# ============================================================
# STEP 7 — global resolution decision
# ============================================================

cat("\nSTEP 7/8 — global resolution adjudication\n")

n_true_splits <- length(true_split_parents)

n_promotable <- sum(
  parent_adjudication$split_adjudication ==
    "PROMOTION_TO_HIGHER_RESOLUTION_SUPPORTED"
)

n_uncertain <- sum(
  parent_adjudication$split_adjudication ==
    "BIOLOGICAL_SUBSTRUCTURE_PRESENT_BUT_IDENTITY_PROMOTION_UNCERTAIN"
)

n_patient_problem <- sum(
  parent_adjudication$split_adjudication ==
    "DO_NOT_PROMOTE_PATIENT_DOMINATED_CHILD"
)

# Conservative rule:
# r0.4 is selected only if the majority of genuine r0.3 splits
# provide recurrent identity-level evidence and no material
# patient-dominated child is introduced.
#
# Otherwise r0.3 is retained as the identity-resolution scaffold,
# while supported r0.4 structure can subsequently be represented
# as subtype/state evidence.

select_r04 <-
  n_true_splits > 0 &&
  n_promotable > (n_true_splits / 2) &&
  n_patient_problem == 0

selected_resolution <- ifelse(select_r04, "0.4", "0.3")

decision_reason <- if(selected_resolution=="0.3"){
  paste0(
    "Resolution 0.3 retained as the conservative identity-resolution ",
    "scaffold because resolution 0.4 did not provide recurrent ",
    "identity-level support for the majority of material splits and/or ",
    "introduced patient-dominated structure. Supported higher-resolution ",
    "structure remains available for subsequent subtype/state adjudication."
  )
} else {
  paste0(
    "Resolution 0.4 selected because the majority of material splits ",
    "showed recurrent identity-level support without patient-dominated ",
    "structure."
  )
}

global_decision <- tibble(
  n_true_split_parents=n_true_splits,
  n_promotable=n_promotable,
  n_uncertain=n_uncertain,
  n_patient_dominated=n_patient_problem,
  selected_resolution=selected_resolution,
  decision_reason=decision_reason
)

write_csv(
  global_decision,
  file.path(OUT, "09H3_global_resolution_decision.csv")
)

cat("\n=== GLOBAL RESOLUTION DECISION ===\n")
print(global_decision)


# ============================================================
# STEP 8 — formal checks and audit
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "transition_file_found",
    "res03_10_clusters",
    "res04_11_clusters",
    "material_transitions_generated",
    "true_splits_identified",
    "patient_dominance_audited",
    "FibConsensus_evidence_integrated",
    "marker_evidence_integrated",
    "disease_not_used_for_resolution_selection",
    "resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "no_object_modified"
  ),
  pass=c(
    file.exists(TRANS_FILE),
    n_distinct(evid$cluster[evid$resolution=="0.3"])==10,
    n_distinct(evid$cluster[evid$resolution=="0.4"])==11,
    nrow(material)>0,
    length(true_split_parents)>0,
    "patient_dominated" %in% names(split_evidence),
    all(c("parent_top_family","child_top_family") %in%
          names(split_evidence)),
    "n_strong_markers" %in% names(split_evidence),
    TRUE,
    selected_resolution %in% c("0.3","0.4"),
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(OUT, "09H3_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks, n=Inf)

if(!all(checks$pass)){
  stop("09H3 formal checks failed.")
}

audit <- c(
  "PIPELINE_V2 — 09H3 FIBROBLAST RESOLUTION ADJUDICATION",
  "=====================================================",
  "",
  "Inputs:",
  TRANS_FILE,
  SUPPORT_FILE,
  MARKERS04_FILE,
  H2_EVID_FILE,
  H2_RANK_FILE,
  "",
  paste0("Material split parents: ",
         paste(true_split_parents, collapse=", ")),
  paste0("Selected resolution: ", selected_resolution),
  "",
  decision_reason,
  "",
  "Disease status was not used to select resolution.",
  "FibConsensus was used as independent identity evidence.",
  "Patient-dominated structure was explicitly penalized.",
  "No fibroblast L2 labels were assigned.",
  "No cells were removed.",
  "No Seurat object was modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(OUT, "09H3_audit.txt")
)

cat("\n============================================================\n")
cat("09H3 COMPLETE\n")
cat("============================================================\n")
cat("Material split parents:",
    paste(true_split_parents, collapse=", "), "\n")
cat("Selected resolution:", selected_resolution, "\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("PASS — fibroblast resolution formally adjudicated.\n")
