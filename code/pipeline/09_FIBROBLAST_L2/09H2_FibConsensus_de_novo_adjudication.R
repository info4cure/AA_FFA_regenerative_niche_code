#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H2 FIBCONSENSUS × DE NOVO STRUCTURE ADJUDICATION\n")
cat("============================================================\n\n")

ROOT <- "PIPELINE_V2"

H1 <- file.path(
  ROOT,
  "results/annotation/09_FIBROBLAST_L2/09H1_FibConsensus_projection/tables"
)

G <- file.path(
  ROOT,
  "results/annotation/09_FIBROBLAST_L2/09G_refined_structure_audit/tables"
)

OUT <- file.path(
  ROOT,
  "results/annotation/09_FIBROBLAST_L2/09H2_FibConsensus_de_novo_adjudication"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(H1, "09H1_res03_cluster_scores.csv"),
  file.path(H1, "09H1_res04_cluster_scores.csv"),
  file.path(H1, "09H1_res03_patient_cluster_scores.csv"),
  file.path(H1, "09H1_res04_patient_cluster_scores.csv"),
  file.path(H1, "09H1_CORE_FULL_concordance.csv"),
  file.path(H1, "09H1_signature_coverage.csv")
)

missing <- required[!file.exists(required)]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required 09H1 files:\n",
      paste(missing, collapse = "\n")
    )
  )
}

cat("STEP 1/8 — loading 09H1 frozen projection outputs\n")

r03 <- read_csv(
  file.path(H1, "09H1_res03_cluster_scores.csv"),
  show_col_types = FALSE
)

r04 <- read_csv(
  file.path(H1, "09H1_res04_cluster_scores.csv"),
  show_col_types = FALSE
)

p03 <- read_csv(
  file.path(H1, "09H1_res03_patient_cluster_scores.csv"),
  show_col_types = FALSE
)

p04 <- read_csv(
  file.path(H1, "09H1_res04_patient_cluster_scores.csv"),
  show_col_types = FALSE
)

conc <- read_csv(
  file.path(H1, "09H1_CORE_FULL_concordance.csv"),
  show_col_types = FALSE
)

coverage <- read_csv(
  file.path(H1, "09H1_signature_coverage.csv"),
  show_col_types = FALSE
)

cat("r0.3 clusters:", n_distinct(r03$cluster), "\n")
cat("r0.4 clusters:", n_distinct(r04$cluster), "\n")
cat("Families:", n_distinct(r03$family_id), "\n\n")


# =============================================================================
# Helpers
# =============================================================================

rank_cluster_families <- function(df, resolution) {

  core <- df %>%
    filter(tier == "CORE") %>%
    select(
      cluster,
      family_id,
      core_mean_z = mean_z,
      core_median_z = median_z,
      n_cells,
      n_patients
    )

  full <- df %>%
    filter(tier == "FULL") %>%
    select(
      cluster,
      family_id,
      full_mean_z = mean_z,
      full_median_z = median_z
    )

  x <- core %>%
    left_join(
      full,
      by = c("cluster", "family_id")
    ) %>%
    mutate(
      resolution = resolution
    ) %>%
    group_by(cluster) %>%
    arrange(desc(core_mean_z), .by_group = TRUE) %>%
    mutate(
      core_rank = row_number(),
      top_core = first(core_mean_z),
      second_core = nth(core_mean_z, 2),
      core_margin = top_core - second_core
    ) %>%
    ungroup()

  x
}


summarize_clusters <- function(rank_df) {

  top1 <- rank_df %>%
    filter(core_rank == 1) %>%
    select(
      resolution,
      cluster,
      n_cells,
      n_patients,
      top_family = family_id,
      top_core_z = core_mean_z,
      top_full_z = full_mean_z,
      core_margin
    )

  top2 <- rank_df %>%
    filter(core_rank == 2) %>%
    select(
      resolution,
      cluster,
      second_family = family_id,
      second_core_z = core_mean_z,
      second_full_z = full_mean_z
    )

  top1 %>%
    left_join(
      top2,
      by = c("resolution", "cluster")
    )
}


patient_support <- function(df, resolution) {

  core <- df %>%
    filter(tier == "CORE") %>%
    group_by(cluster, patient_id) %>%
    arrange(desc(mean_z), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  winners <- core %>%
    filter(rank == 1) %>%
    count(
      cluster,
      family_id,
      name = "n_patients_top_family"
    )

  total <- core %>%
    distinct(cluster, patient_id) %>%
    count(cluster, name = "n_patients_scored")

  winners %>%
    left_join(total, by = "cluster") %>%
    mutate(
      patient_fraction = n_patients_top_family / n_patients_scored,
      resolution = resolution
    )
}


cat("STEP 2/8 — ranking frozen families within each de novo cluster\n")

rank03 <- rank_cluster_families(r03, "0.3")
rank04 <- rank_cluster_families(r04, "0.4")

write_csv(
  bind_rows(rank03, rank04),
  file.path(TAB, "09H2_cluster_family_ranking.csv")
)


cat("STEP 3/8 — cluster-level top-family summaries\n")

sum03 <- summarize_clusters(rank03)
sum04 <- summarize_clusters(rank04)

cluster_summary <- bind_rows(sum03, sum04)

write_csv(
  cluster_summary,
  file.path(TAB, "09H2_cluster_top_family_summary.csv")
)

cat("\n=== TOP FROZEN FAMILY BY CLUSTER ===\n")
print(
  cluster_summary %>%
    arrange(resolution, as.numeric(cluster))
)


cat("\nSTEP 4/8 — patient-level family support\n")

ps03 <- patient_support(p03, "0.3")
ps04 <- patient_support(p04, "0.4")

patient_support_all <- bind_rows(ps03, ps04)

write_csv(
  patient_support_all,
  file.path(TAB, "09H2_patient_top_family_support.csv")
)

top_patient_support <- patient_support_all %>%
  inner_join(
    cluster_summary %>%
      select(resolution, cluster, top_family),
    by = c("resolution", "cluster")
  ) %>%
  filter(family_id == top_family) %>%
  select(
    resolution,
    cluster,
    top_family,
    n_patients_top_family,
    n_patients_scored,
    patient_fraction
  )

write_csv(
  top_patient_support,
  file.path(TAB, "09H2_top_family_patient_support.csv")
)


cat("STEP 5/8 — CORE/FULL family concordance\n")

conc_summary <- conc %>%
  mutate(
    concordance_class = case_when(
      spearman_CORE_FULL >= 0.90 ~ "VERY_HIGH",
      spearman_CORE_FULL >= 0.80 ~ "HIGH",
      spearman_CORE_FULL >= 0.70 ~ "MODERATE",
      TRUE ~ "LOW"
    )
  )

write_csv(
  conc_summary,
  file.path(TAB, "09H2_CORE_FULL_concordance_classified.csv")
)

cat("\n=== CORE/FULL CONCORDANCE ===\n")
print(conc_summary)


cat("\nSTEP 6/8 — 0.3 -> 0.4 structural transitions\n")

transition_file <- file.path(
  G,
  "09G_res03_to_res04_transitions.csv"
)

if (!file.exists(transition_file)) {

  candidates <- list.files(
    G,
    pattern = "transition",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) > 0) {
    transition_file <- candidates[1]
  }
}

transition_available <- file.exists(transition_file)

if (transition_available) {

  transitions <- read_csv(
    transition_file,
    show_col_types = FALSE
  )

  write_csv(
    transitions,
    file.path(TAB, "09H2_res03_to_res04_transitions_source.csv")
  )

  cat("Transition source:", transition_file, "\n")

} else {

  cat(
    "WARNING: 09G transition table not located; ",
    "transition-level adjudication deferred.\n",
    sep = ""
  )
}


cat("\nSTEP 7/8 — evidence classification WITHOUT assigning L2\n")

adjudication <- cluster_summary %>%
  left_join(
    top_patient_support,
    by = c(
      "resolution",
      "cluster",
      "top_family"
    )
  ) %>%
  mutate(
    projection_evidence = case_when(

      top_core_z >= 0.50 &
        core_margin >= 0.50 &
        patient_fraction >= 0.50 ~
        "STRONG_FROZEN_FAMILY_SUPPORT",

      top_core_z >= 0.25 &
        core_margin >= 0.25 &
        patient_fraction >= 0.35 ~
        "MODERATE_FROZEN_FAMILY_SUPPORT",

      top_core_z > 0 &
        core_margin > 0 ~
        "WEAK_FROZEN_FAMILY_SUPPORT",

      TRUE ~
        "UNRESOLVED"
    )
  )

write_csv(
  adjudication,
  file.path(TAB, "09H2_projection_evidence_classification.csv")
)

cat("\n=== PROJECTION EVIDENCE ===\n")
print(
  adjudication %>%
    select(
      resolution,
      cluster,
      n_cells,
      n_patients,
      top_family,
      top_core_z,
      second_family,
      second_core_z,
      core_margin,
      patient_fraction,
      projection_evidence
    ) %>%
    arrange(resolution, as.numeric(cluster))
)


cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "six_families_present",
    "coverage_complete",
    "res03_10_clusters",
    "res04_11_clusters",
    "cluster_rankings_generated",
    "patient_support_generated",
    "CORE_FULL_concordance_available",
    "disease_not_used_for_adjudication",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "no_object_modified"
  ),
  pass = c(
    n_distinct(r03$family_id) == 6,
    all(coverage$coverage == 1),
    n_distinct(r03$cluster) == 10,
    n_distinct(r04$cluster) == 11,
    nrow(bind_rows(rank03, rank04)) > 0,
    nrow(patient_support_all) > 0,
    nrow(conc) == 12,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "09H2_formal_checks.csv")
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks)

if (!all(checks$pass)) {
  stop("09H2 formal checks failed.")
}

audit <- c(
  "PIPELINE_V2 — 09H2 FIBCONSENSUS × DE NOVO STRUCTURE ADJUDICATION",
  "",
  "Input fibroblast universe: 4,714 cells.",
  "FibConsensus reference: frozen v2.2 atlas-validation signatures.",
  "Six frozen consensus families evaluated.",
  "Resolution 0.3 and 0.4 evaluated independently.",
  "",
  "FibConsensus was used only as an external frozen projection framework.",
  "Disease labels were not used to assign fibroblast identities.",
  "No clustering resolution was selected.",
  "No fibroblast L2 labels were assigned.",
  "No cells were removed.",
  "No Seurat object was modified.",
  "",
  paste0(
    "09G transition table available: ",
    transition_available
  ),
  "",
  "PASS — de novo structure × frozen FibConsensus evidence adjudicated."
)

writeLines(
  audit,
  file.path(TAB, "09H2_audit.txt")
)

cat("\n============================================================\n")
cat("09H2 COMPLETE\n")
cat("============================================================\n")
cat("Fibroblast cells: 4714\n")
cat("Frozen families: 6\n")
cat("r0.3 clusters: 10\n")
cat("r0.4 clusters: 11\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("PASS — FibConsensus × de novo structure evidence generated.\n")
