suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

cat("============================================================\n")
cat("21G3B1 — CROSS-MODAL EVIDENCE MATRIX\n")
cat("============================================================\n")

BASE <- "PIPELINE_V2/results/disease_biology"

COMM_DIR <- file.path(
  BASE,
  "21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION"
)

OUT <- file.path(
  COMM_DIR,
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(OUT, "tables")
OBJ <- file.path(OUT, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

COMM_FILE <- file.path(
  COMM_DIR,
  "21G3_BIOLOGICAL_READOUT",
  "tables",
  "21G2C_primary_A_all.csv"
)

DECOMP_FILE <- file.path(
  COMM_DIR,
  "21G3_BIOLOGICAL_READOUT",
  "tables",
  "21G2C_primary_signal_decomposition_wide.csv"
)

BRIDGE_FILE <- file.path(
  TAB,
  "21G3A4B_frozen_communication_pathway_axis_bridge.csv"
)

DA_FILE <- file.path(
  BASE,
  "17B_composition_DA",
  "confirmatory",
  "17B3_integrated_DA_evidence.csv"
)

DE_FILE <- file.path(
  BASE,
  "18_PSEUDOBULK_DE",
  "18C_DE",
  "tables",
  "18C_DE_summary_by_L2.csv"
)

PATH_FILE <- file.path(
  BASE,
  "20_PATHWAY_ANALYSIS",
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2B_integrated_prioritization",
  "tables",
  "20F2B_integrated_state_prioritization.csv"
)

stopifnot(
  file.exists(COMM_FILE),
  file.exists(DECOMP_FILE),
  file.exists(BRIDGE_FILE),
  file.exists(DA_FILE),
  file.exists(DE_FILE),
  file.exists(PATH_FILE)
)

cat("\nSTEP 1/8 — loading frozen sources\n")

comm <- read_csv(COMM_FILE, show_col_types=FALSE)
decomp <- read_csv(DECOMP_FILE, show_col_types=FALSE)
bridge <- read_csv(BRIDGE_FILE, show_col_types=FALSE)
da <- read_csv(DA_FILE, show_col_types=FALSE)
de <- read_csv(DE_FILE, show_col_types=FALSE)
path <- read_csv(PATH_FILE, show_col_types=FALSE)

cat("COMM rows:", nrow(comm), "\n")
cat("DA rows:", nrow(da), "\n")
cat("DE rows:", nrow(de), "\n")
cat("PATHWAY rows:", nrow(path), "\n")
cat("BRIDGE rows:", nrow(bridge), "\n")

# ------------------------------------------------------------
# 2. Freeze analytical anchor
# ------------------------------------------------------------

cat("\nSTEP 2/8 — freezing communication anchor\n")

anchor <- comm %>%
  transmute(
    receiver,
    macrocompartment,
    contrast,

    communication_program =
      biological_program,

    communication_tier =
      analysis_tier,

    communication_mechanistic_class =
      mechanistic_class,

    communication_action =
      action,

    communication_axis =
      axis,

    communication_candidate_status =
      candidate_status,

    communication_mean_difference =
      mean_difference,

    communication_log2_mean_ratio =
      log2_mean_ratio,

    communication_p_exact =
      p_exact,

    communication_q_BH =
      q_BH,

    communication_evidence_class =
      evidence_class,

    communication_direction =
      direction
  )

# Full frozen communication analytical key.
# `communication_action` must remain part of the key because INCLUDE and
# REVIEW are distinct prespecified records in the frozen 21G2 universe.
# They must not be collapsed, averaged, filtered, or selected post hoc.
anchor_key_check <- anchor %>%
  count(
    receiver,
    contrast,
    communication_program,
    communication_tier,
    communication_mechanistic_class,
    communication_action,
    communication_axis,
    name="n"
  ) %>%
  filter(n > 1)

cat("Communication anchor rows:", nrow(anchor), "\n")
cat("Duplicate anchor keys:", nrow(anchor_key_check), "\n")

stopifnot(
  nrow(anchor_key_check) == 0
)

# ------------------------------------------------------------
# 3. Add frozen bridge
# ------------------------------------------------------------

cat("\nSTEP 3/8 — attaching frozen axis bridge\n")

master <- anchor %>%
  left_join(
    bridge %>%
      select(
        communication_axis,
        pathway_axis,
        bridge_class,
        bridge_interpretation
      ),
    by="communication_axis"
  )

stopifnot(
  nrow(master) == nrow(anchor),
  all(!is.na(master$bridge_class))
)

# ------------------------------------------------------------
# 4. Add DA evidence
# ------------------------------------------------------------

cat("\nSTEP 4/8 — attaching frozen DA evidence\n")

da_small <- da %>%
  transmute(
    receiver =
      final_L2_label,

    contrast,

    DA_beta_log_odds =
      beta_log_odds,

    DA_OR =
      OR,

    DA_FDR_compartment =
      FDR_BH_compartment,

    DA_FDR_global =
      FDR_BH_global,

    DA_direction =
      confirmatory_direction,

    DA_wilcox_direction_concordant =
      wilcox_confirmatory_direction_concordant
  )

stopifnot(
  !anyDuplicated(
    da_small[
      c("receiver","contrast")
    ]
  )
)

master <- master %>%
  left_join(
    da_small,
    by=c(
      "receiver",
      "contrast"
    )
  )

# ------------------------------------------------------------
# 5. Add DE evidence
# ------------------------------------------------------------

cat("\nSTEP 5/8 — attaching frozen DE evidence\n")

de_small <- de %>%
  transmute(
    receiver =
      final_L2_label,

    contrast,

    DE_batch_adjusted =
      batch_adjusted,

    DE_robust_ge50 =
      robust_ge50,

    DE_n_genes_tested =
      n_genes_tested,

    DE_n_FDR005 =
      n_FDR_005,

    DE_n_FDR001 =
      n_FDR_001,

    DE_n_up_FDR005 =
      n_up_FDR_005,

    DE_n_down_FDR005 =
      n_down_FDR_005,

    DE_max_abs_logFC_FDR005 =
      max_abs_logFC_FDR005
  )

stopifnot(
  !anyDuplicated(
    de_small[
      c("receiver","contrast")
    ]
  )
)

master <- master %>%
  left_join(
    de_small,
    by=c(
      "receiver",
      "contrast"
    )
  )

# ------------------------------------------------------------
# 6. Add pathway evidence without many-to-many expansion
# ------------------------------------------------------------

cat("\nSTEP 6/8 — attaching frozen pathway evidence\n")

# Convert the disease-specific columns of 20F2B into
# contrast-resolved records. AA-vs-FFA remains the direct
# triangulation evidence already frozen in 20F2B.

path_AA <- path %>%
  transmute(
    receiver =
      final_L2_label,

    pathway_axis =
      axis,

    contrast =
      "AA_vs_HC",

    pathway_status =
      AA_primary_status,

    pathway_direction =
      AA_directions,

    pathway_n_significant =
      AA_n_significant_pathways,

    pathway_median_NES =
      AA_median_database_NES,

    pathway_min_FDR =
      AA_min_FDR,

    pathway_architecture =
      architecture,

    pathway_evidence_tier =
      integrated_evidence_tier,

    pathway_publication_role =
      publication_role
  )

path_FFA <- path %>%
  transmute(
    receiver =
      final_L2_label,

    pathway_axis =
      axis,

    contrast =
      "FFA_vs_HC",

    pathway_status =
      FFA_primary_status,

    pathway_direction =
      FFA_directions,

    pathway_n_significant =
      FFA_n_significant_pathways,

    pathway_median_NES =
      FFA_median_database_NES,

    pathway_min_FDR =
      FFA_min_FDR,

    pathway_architecture =
      architecture,

    pathway_evidence_tier =
      integrated_evidence_tier,

    pathway_publication_role =
      publication_role
  )

path_direct <- path %>%
  transmute(
    receiver =
      final_L2_label,

    pathway_axis =
      axis,

    contrast =
      "AA_vs_FFA",

    pathway_status =
      direct_status,

    pathway_direction =
      direct_status,

    pathway_n_significant =
      direct_n_significant_pathways,

    pathway_median_NES =
      direct_median_NES,

    pathway_min_FDR =
      direct_min_FDR,

    pathway_architecture =
      architecture,

    pathway_evidence_tier =
      integrated_evidence_tier,

    pathway_publication_role =
      publication_role
  )

path_long <- bind_rows(
  path_AA,
  path_FFA,
  path_direct
)

path_key_duplicates <- path_long %>%
  count(
    receiver,
    pathway_axis,
    contrast,
    name="n"
  ) %>%
  filter(n > 1)

cat(
  "Duplicate pathway receiver × axis × contrast keys:",
  nrow(path_key_duplicates),
  "\n"
)

stopifnot(
  nrow(path_key_duplicates) == 0
)

pre_path_n <- nrow(master)

master <- master %>%
  left_join(
    path_long,
    by=c(
      "receiver",
      "pathway_axis",
      "contrast"
    )
  )

stopifnot(
  nrow(master) == pre_path_n
)

master <- master %>%
  mutate(
    pathway_bridge_status =
      case_when(
        bridge_class ==
          "NO_PRESPECIFIED_BRIDGE" ~
          "NO_PRESPECIFIED_BRIDGE",

        !is.na(pathway_status) ~
          "PATHWAY_RECORD_AVAILABLE",

        TRUE ~
          "NO_MATCHING_20F2B_STATE_AXIS_RECORD"
      )
  )

# ------------------------------------------------------------
# 7. Add decomposition only where it already exists
# ------------------------------------------------------------

cat("\nSTEP 7/8 — attaching frozen communication decomposition\n")

# Decomposition is a property of the complete frozen communication
# analytical record. Preserve mechanistic_class and action explicitly so
# decomposition from an INCLUDE record cannot propagate to a distinct
# REVIEW record sharing the same receiver/program/axis.

decomp_small <- decomp %>%
  transmute(
    receiver,
    macrocompartment,
    contrast,

    communication_program =
      biological_program,

    communication_tier =
      analysis_tier,

    communication_mechanistic_class =
      mechanistic_class,

    communication_action =
      action,

    communication_axis =
      axis,

    decomposition_pattern,

    density_mean_difference =
      mean_difference__pooled_D_axis,

    density_q_BH =
      q_BH__pooled_D_axis,

    active_intensity_mean_difference =
      mean_difference__mean_sender_active_intensity,

    active_intensity_q_BH =
      q_BH__mean_sender_active_intensity
  )

# The frozen decomposition source must be unique on the complete
# communication analytical key.
decomp_duplicates <- decomp_small %>%
  count(
    receiver,
    macrocompartment,
    contrast,
    communication_program,
    communication_tier,
    communication_mechanistic_class,
    communication_action,
    communication_axis,
    name="n"
  ) %>%
  filter(n > 1)

cat(
  "Duplicate decomposition full analytical keys:",
  nrow(decomp_duplicates),
  "\n"
)

stopifnot(
  nrow(decomp_duplicates) == 0
)

pre_decomp_n <- nrow(master)

master <- master %>%
  left_join(
    decomp_small,
    by=c(
      "receiver",
      "macrocompartment",
      "contrast",
      "communication_program",
      "communication_tier",
      "communication_mechanistic_class",
      "communication_action",
      "communication_axis"
    )
  )

stopifnot(
  nrow(master) == pre_decomp_n
)

# Frozen 21G2C contains exactly 16 decomposition records. The full-key
# join must therefore annotate exactly those 16 communication records
# and must not create decomposition evidence on REVIEW records.
stopifnot(
  sum(!is.na(master$decomposition_pattern)) ==
    nrow(decomp_small),

  sum(
    !is.na(master$decomposition_pattern) &
      master$communication_action == "REVIEW"
  ) == 0
)

# ------------------------------------------------------------
# 8. Coverage + freeze
# ------------------------------------------------------------

cat("\nSTEP 8/8 — auditing cross-modal coverage\n")

coverage <- master %>%
  summarise(
    n_anchor_rows =
      n(),

    n_receivers =
      n_distinct(receiver),

    n_contrasts =
      n_distinct(contrast),

    DA_available =
      sum(!is.na(DA_beta_log_odds)),

    DE_available =
      sum(!is.na(DE_n_genes_tested)),

    pathway_bridge_defined =
      sum(
        bridge_class !=
          "NO_PRESPECIFIED_BRIDGE"
      ),

    pathway_record_available =
      sum(
        pathway_bridge_status ==
          "PATHWAY_RECORD_AVAILABLE"
      ),

    no_prespecified_bridge =
      sum(
        pathway_bridge_status ==
          "NO_PRESPECIFIED_BRIDGE"
      ),

    no_matching_20F2B_record =
      sum(
        pathway_bridge_status ==
          "NO_MATCHING_20F2B_STATE_AXIS_RECORD"
      ),

    decomposition_available =
      sum(!is.na(decomposition_pattern))
  )

cat("\n===== CROSS-MODAL COVERAGE =====\n")
print(coverage, width=Inf)

status_table <- master %>%
  count(
    contrast,
    macrocompartment,
    bridge_class,
    pathway_bridge_status,
    name="n"
  ) %>%
  arrange(
    contrast,
    macrocompartment,
    bridge_class,
    pathway_bridge_status
  )

cat("\n===== PATHWAY BRIDGE STATUS =====\n")
print(status_table, n=Inf, width=Inf)

formal_checks <- tribble(
  ~check, ~pass,

  "communication_anchor_retained",
  nrow(master) == nrow(anchor),

  "twenty_eight_receivers_retained",
  n_distinct(master$receiver) == 28,

  "three_contrasts_retained",
  n_distinct(master$contrast) == 3,

  "nine_macrocompartments_retained",
  n_distinct(master$macrocompartment) == 9,

  "all_bridge_classes_present",
  all(!is.na(master$bridge_class)),

  "DA_join_no_row_expansion",
  TRUE,

  "DE_join_no_row_expansion",
  TRUE,

  "pathway_join_no_row_expansion",
  TRUE,

  "decomposition_join_no_row_expansion",
  TRUE,

  "no_crossmodal_score",
  TRUE,

  "no_concordance_classification",
  TRUE,

  "no_combined_p_value",
  TRUE,

  "no_combined_q_value",
  TRUE,

  "no_new_statistical_test",
  TRUE,

  "no_DE_refit",
  TRUE,

  "no_DA_refit",
  TRUE,

  "no_pathway_refit",
  TRUE,

  "no_communication_refit",
  TRUE
)

write_csv(
  master,
  file.path(
    TAB,
    "21G3B1_crossmodal_evidence_matrix.csv"
  )
)

saveRDS(
  master,
  file.path(
    OBJ,
    "21G3B1_crossmodal_evidence_matrix.rds"
  )
)

write_csv(
  coverage,
  file.path(
    TAB,
    "21G3B1_crossmodal_coverage.csv"
  )
)

write_csv(
  status_table,
  file.path(
    TAB,
    "21G3B1_pathway_bridge_status.csv"
  )
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G3B1_formal_checks.csv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")
print(formal_checks, n=Inf, width=Inf)

stopifnot(
  all(formal_checks$pass)
)

cat("\n============================================================\n")
cat("21G3B1 COMPLETE\n")
cat("============================================================\n")
cat("Cross-modal evidence matrix constructed\n")
cat("Communication retained as anchor modality\n")
cat("DA and DE attached by exact receiver × contrast\n")
cat("Pathway evidence attached only through frozen axis bridge\n")
cat("No-prespecified-bridge states retained explicitly\n")
cat("Decomposition evidence retained without expansion\n")
cat("No cross-modal concordance classification\n")
cat("No cross-modal score\n")
cat("No combined p/q values\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
