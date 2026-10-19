suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B2 PATIENT MECHANISTIC LR-AXIS PROJECTION\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

F1 <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

F2 <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES"
)

OUT <- file.path(
  F2,
  "21F2B2_PATIENT_AXIS_PROJECTION"
)

TAB <- file.path(OUT, "tables")
OBJ <- file.path(OUT, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

NETWORK_FILE <- file.path(
  F1,
  "objects",
  "21F1_patient_level_LR_network.rds"
)

STRUCT_FILE <- file.path(
  F1,
  "objects",
  "21F1B_structural_sender_receiver_LR_universe.rds"
)

EDGE_FILE <- file.path(
  F1,
  "objects",
  "21F1B_patient_sender_receiver_metrics.rds"
)

ADJ_FILE <- file.path(
  F2,
  "tables",
  "21F2B1B_pair_level_adjudication.csv"
)

stopifnot(
  file.exists(NETWORK_FILE),
  file.exists(STRUCT_FILE),
  file.exists(EDGE_FILE),
  file.exists(ADJ_FILE)
)

# ============================================================
# 1. Recover frozen inputs
# ============================================================

cat("\nSTEP 1/11 — recovering frozen inputs\n")

net <- readRDS(NETWORK_FILE)
struct <- readRDS(STRUCT_FILE)
edge <- readRDS(EDGE_FILE)

adj <- read_csv(
  ADJ_FILE,
  show_col_types=FALSE
)

cat("Patient LR rows       :", nrow(net), "\n")
cat("Structural LR rows    :", nrow(struct), "\n")
cat("Evaluable edge rows   :", nrow(edge), "\n")
cat("Adjudicated axis pairs:", nrow(adj), "\n")

stopifnot(
  nrow(net) == 2939108,
  nrow(struct) == 562941,
  nrow(edge) == 9883,
  nrow(adj) == 822
)

# ============================================================
# 2. Freeze projection universe
# ============================================================

cat("\nSTEP 2/11 — freezing projection universe\n")

pair_annotation <- adj %>%
  select(
    axis,
    ligand,
    receptor,
    action,
    adjudication_status,
    mechanistic_class,
    polarity,
    interpretive_scope,
    super_axis,
    directional_interpretation,
    selected_rule_id,
    reason
  ) %>%
  distinct()

stopifnot(
  !anyDuplicated(
    pair_annotation[
      c(
        "axis",
        "ligand",
        "receptor"
      )
    ]
  )
)

projection_inventory <- pair_annotation %>%
  count(
    action,
    adjudication_status,
    mechanistic_class,
    name="n_pairs"
  ) %>%
  arrange(
    action,
    mechanistic_class
  )

cat("\n===== PROJECTION INVENTORY =====\n")
print(
  projection_inventory,
  n=Inf,
  width=Inf
)

write_csv(
  projection_inventory,
  file.path(
    TAB,
    "21F2B2_projection_inventory.csv"
  )
)

# EXCLUDE pairs remain documented but are not projected as
# active mechanistic communication.
projectable <- pair_annotation %>%
  filter(
    action %in%
      c(
        "INCLUDE",
        "REVIEW"
      )
  )

excluded <- pair_annotation %>%
  filter(
    action == "EXCLUDE"
  )

write_csv(
  excluded,
  file.path(
    TAB,
    "21F2B2_excluded_pair_inventory.csv"
  )
)

# ============================================================
# 3. Construct axis-specific structural denominators
# ============================================================

cat("\nSTEP 3/11 — constructing axis-specific structural denominators\n")

struct_axis <- struct %>%
  inner_join(
    projectable,
    by=c(
      "ligand",
      "receptor"
    ),
    relationship="many-to-many"
  ) %>%
  distinct(
    sender,
    receiver,
    axis,
    mechanistic_class,
    action,
    ligand,
    receptor,
    .keep_all=TRUE
  )

M_axis <- struct_axis %>%
  count(
    sender,
    receiver,
    axis,
    mechanistic_class,
    action,
    name="M_axis"
  )

cat(
  "Structural axis LR rows:",
  nrow(struct_axis),
  "\n"
)

cat(
  "Structural sender-receiver-axis denominators:",
  nrow(M_axis),
  "\n"
)

stopifnot(
  all(M_axis$M_axis > 0)
)

write_csv(
  M_axis,
  file.path(
    TAB,
    "21F2B2_axis_structural_denominator.csv"
  )
)

saveRDS(
  struct_axis,
  file.path(
    OBJ,
    "21F2B2_axis_structural_LR_universe.rds"
  )
)

# ============================================================
# 4. Project active patient-level network
# ============================================================

cat("\nSTEP 4/11 — projecting active patient-level LR network\n")

active_axis <- net %>%
  inner_join(
    projectable,
    by=c(
      "ligand",
      "receptor"
    ),
    relationship="many-to-many"
  ) %>%
  distinct(
    patient_id,
    disease,
    sender,
    receiver,
    axis,
    mechanistic_class,
    action,
    ligand,
    receptor,
    .keep_all=TRUE
  )

cat(
  "Projected active LR rows:",
  nrow(active_axis),
  "\n"
)

stopifnot(
  all(
    active_axis$communication_score > 0
  )
)

# ============================================================
# 5. Active axis metrics N, S and intensity
# ============================================================

cat("\nSTEP 5/11 — calculating active axis metrics\n")

active_metrics <- active_axis %>%
  group_by(
    patient_id,
    disease,
    sender,
    receiver,
    axis,
    mechanistic_class,
    action
  ) %>%
  summarise(
    N_axis=n(),
    S_axis=sum(
      communication_score
    ),
    meanC_axis=mean(
      communication_score
    ),
    medianC_axis=median(
      communication_score
    ),
    minC_axis=min(
      communication_score
    ),
    maxC_axis=max(
      communication_score
    ),
    sender_n_cells=first(
      sender_n_cells
    ),
    receiver_n_cells=first(
      receiver_n_cells
    ),
    .groups="drop"
  )

# ============================================================
# 6. Construct complete evaluable patient-edge-axis universe
# ============================================================

cat("\nSTEP 6/11 — constructing evaluable patient-edge-axis universe\n")

# 21F1B edge rows define which patient × sender × receiver
# combinations are evaluable. Disease-independent M_axis then
# defines which mechanistic axes are structurally possible.

evaluable_axis <- edge %>%
  select(
    patient_id,
    disease,
    sender,
    sender_n_cells,
    receiver,
    receiver_n_cells
  ) %>%
  inner_join(
    M_axis,
    by=c(
      "sender",
      "receiver"
    ),
    relationship="many-to-many"
  ) %>%
  left_join(
    active_metrics,
    by=c(
      "patient_id",
      "disease",
      "sender",
      "receiver",
      "axis",
      "mechanistic_class",
      "action"
    ),
    suffix=c("", "_active")
  ) %>%
  mutate(
    N_axis=replace_na(
      N_axis,
      0L
    ),
    S_axis=replace_na(
      S_axis,
      0
    ),

    # Active-only intensity is undefined when N=0.
    meanC_axis=if_else(
      N_axis > 0,
      meanC_axis,
      NA_real_
    ),
    medianC_axis=if_else(
      N_axis > 0,
      medianC_axis,
      NA_real_
    ),
    minC_axis=if_else(
      N_axis > 0,
      minC_axis,
      NA_real_
    ),
    maxC_axis=if_else(
      N_axis > 0,
      maxC_axis,
      NA_real_
    ),

    D_axis=
      N_axis /
      M_axis,

    A_axis=
      S_axis /
      M_axis,

    axis_evaluable=TRUE,
    axis_active=
      N_axis > 0
  ) %>%
  select(
    patient_id,
    disease,
    sender,
    sender_n_cells,
    receiver,
    receiver_n_cells,
    axis,
    mechanistic_class,
    action,
    M_axis,
    N_axis,
    D_axis,
    meanC_axis,
    medianC_axis,
    minC_axis,
    maxC_axis,
    S_axis,
    A_axis,
    axis_evaluable,
    axis_active
  )

# Formal metric identities
stopifnot(
  all(
    abs(
      evaluable_axis$D_axis -
        evaluable_axis$N_axis /
        evaluable_axis$M_axis
    ) < 1e-12
  ),
  all(
    abs(
      evaluable_axis$A_axis -
        evaluable_axis$S_axis /
        evaluable_axis$M_axis
    ) < 1e-12
  )
)

# A = D * meanC only for active axes.
active_check <- evaluable_axis %>%
  filter(
    N_axis > 0
  )

stopifnot(
  all(
    abs(
      active_check$A_axis -
        active_check$D_axis *
        active_check$meanC_axis
    ) < 1e-10
  )
)

# ============================================================
# 7. Freeze metric dictionary
# ============================================================

cat("\nSTEP 7/11 — freezing mechanistic-axis metric dictionary\n")

metric_dictionary <- tribble(
  ~metric,
  ~symbol,
  ~formula,
  ~denominator,
  ~interpretation,

  "Structural axis opportunity",
  "M_axis",
  "Number of structurally possible adjudicated LR pairs in mechanistic class a for sender s and receiver r",
  "none",
  "Disease-independent mechanistic LR opportunity",

  "Active axis richness",
  "N_axis",
  "sum_k I(C_psrk > 0), k in mechanistic class a",
  "none",
  "Number of active LR channels within the mechanistic class",

  "Axis LR density",
  "D_axis",
  "N_axis / M_axis",
  "M_axis",
  "Fraction of structurally possible mechanistic LR channels active",

  "Mean active axis intensity",
  "meanC_axis",
  "(1 / N_axis) * sum_{k:C>0} C_psrk",
  "N_axis active LR pairs",
  "Mean communication score among active mechanistic LR channels",

  "Aggregate axis score",
  "S_axis",
  "sum_{k:C>0} C_psrk",
  "none",
  "Total active communication score within the mechanistic class",

  "Normalized axis exposure",
  "A_axis",
  "S_axis / M_axis = D_axis * meanC_axis when N_axis > 0",
  "M_axis",
  "Structural-opportunity-normalized mechanistic communication exposure"
)

write_csv(
  metric_dictionary,
  file.path(
    TAB,
    "21F2B2_metric_dictionary.csv"
  )
)

# ============================================================
# 8. Follicular receiver panel
# ============================================================

cat("\nSTEP 8/11 — freezing follicular receiver projection\n")

follicular_receivers <- c(
  "Bulge_Isthmus",
  "ORS",
  "ORS_basal"
)

follicular <- evaluable_axis %>%
  filter(
    receiver %in%
      follicular_receivers
  )

cat(
  "Follicular patient-sender-receiver-axis rows:",
  nrow(follicular),
  "\n"
)

write_csv(
  follicular,
  file.path(
    TAB,
    "21F2B2_follicular_sender_receiver_axis_metrics.csv"
  )
)

# ============================================================
# 9. Patient × receiver × axis incoming aggregation
# ============================================================

cat("\nSTEP 9/11 — aggregating incoming mechanistic exposure\n")

incoming_axis <- evaluable_axis %>%
  group_by(
    patient_id,
    disease,
    receiver,
    axis,
    mechanistic_class,
    action
  ) %>%
  summarise(
    n_evaluable_senders=n(),
    n_active_senders=sum(
      axis_active
    ),
    total_M_axis=sum(
      M_axis
    ),
    total_N_axis=sum(
      N_axis
    ),
    pooled_D_axis=
      total_N_axis /
      total_M_axis,
    mean_sender_D_axis=mean(
      D_axis
    ),
    median_sender_D_axis=median(
      D_axis
    ),
    mean_sender_active_intensity=
      if (
        any(N_axis > 0)
      )
        mean(
          meanC_axis[
            N_axis > 0
          ]
        )
      else
        NA_real_,
    median_sender_active_intensity=
      if (
        any(N_axis > 0)
      )
        median(
          meanC_axis[
            N_axis > 0
          ]
        )
      else
        NA_real_,
    total_S_axis=sum(
      S_axis
    ),
    mean_incoming_A_axis=mean(
      A_axis
    ),
    median_incoming_A_axis=median(
      A_axis
    ),
    total_normalized_axis_exposure=sum(
      A_axis
    ),
    .groups="drop"
  )

follicular_incoming <- incoming_axis %>%
  filter(
    receiver %in%
      follicular_receivers
  )

write_csv(
  incoming_axis,
  file.path(
    TAB,
    "21F2B2_patient_receiver_axis_metrics.csv"
  )
)

write_csv(
  follicular_incoming,
  file.path(
    TAB,
    "21F2B2_follicular_incoming_axis_metrics.csv"
  )
)

# ============================================================
# 10. Save frozen objects and descriptive audit
# ============================================================

cat("\nSTEP 10/11 — writing frozen projection objects\n")

saveRDS(
  evaluable_axis,
  file.path(
    OBJ,
    "21F2B2_patient_sender_receiver_axis_metrics.rds"
  )
)

saveRDS(
  incoming_axis,
  file.path(
    OBJ,
    "21F2B2_patient_receiver_axis_metrics.rds"
  )
)

write_csv(
  evaluable_axis,
  file.path(
    TAB,
    "21F2B2_patient_sender_receiver_axis_metrics.csv"
  )
)

summary_axis <- evaluable_axis %>%
  group_by(
    axis,
    mechanistic_class,
    action
  ) %>%
  summarise(
    n_patient_edges=n(),
    n_patients=n_distinct(
      patient_id
    ),
    n_senders=n_distinct(
      sender
    ),
    n_receivers=n_distinct(
      receiver
    ),
    median_M_axis=median(
      M_axis
    ),
    median_N_axis=median(
      N_axis
    ),
    median_D_axis=median(
      D_axis
    ),
    median_meanC_axis=median(
      meanC_axis,
      na.rm=TRUE
    ),
    median_A_axis=median(
      A_axis
    ),
    .groups="drop"
  )

write_csv(
  summary_axis,
  file.path(
    TAB,
    "21F2B2_axis_projection_summary.csv"
  )
)

cat("\n===== MECHANISTIC AXIS PROJECTION SUMMARY =====\n")
print(
  summary_axis,
  n=Inf,
  width=Inf
)

# ============================================================
# 11. Formal checks
# ============================================================

cat("\nSTEP 11/11 — formal checks\n")

checks <- tibble(
  check=c(
    "21F1_network_present",
    "21F1B_structural_universe_present",
    "21F1B_evaluable_edges_present",
    "21F2B1B_adjudication_present",
    "2939108_active_LR_rows_recovered",
    "562941_structural_LR_rows_recovered",
    "9883_evaluable_edges_recovered",
    "822_adjudicated_pairs_recovered",
    "excluded_pairs_not_projected",
    "review_pairs_retained_contextually",
    "axis_structural_denominator_generated",
    "axis_denominator_positive",
    "patient_axis_network_generated",
    "evaluable_inactive_axes_zero_N",
    "active_intensity_NA_when_inactive",
    "D_equals_N_over_M",
    "A_equals_S_over_M",
    "A_equals_D_times_meanC_when_active",
    "metric_dictionary_generated",
    "follicular_receivers_retained",
    "sender_identity_retained",
    "incoming_axis_aggregation_generated",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_threshold_refit",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(NETWORK_FILE),
    file.exists
(STRUCT_FILE),
    file.exists(EDGE_FILE),
    file.exists(ADJ_FILE),

    nrow(net) == 2939108,
    nrow(struct) == 562941,
    nrow(edge) == 9883,
    nrow(adj) == 822,

    !any(
      active_axis$action == "EXCLUDE"
    ),

    all(
      projectable$action %in%
        c(
          "INCLUDE",
          "REVIEW"
        )
    ) &&
      any(
        projectable$action == "REVIEW"
      ),

    nrow(M_axis) > 0,

    all(
      M_axis$M_axis > 0
    ),

    nrow(evaluable_axis) > 0,

    all(
      evaluable_axis$N_axis[
        !evaluable_axis$axis_active
      ] == 0
    ),

    all(
      is.na(
        evaluable_axis$meanC_axis[
          !evaluable_axis$axis_active
        ]
      )
    ),

    all(
      abs(
        evaluable_axis$D_axis -
          evaluable_axis$N_axis /
            evaluable_axis$M_axis
      ) < 1e-12
    ),

    all(
      abs(
        evaluable_axis$A_axis -
          evaluable_axis$S_axis /
            evaluable_axis$M_axis
      ) < 1e-12
    ),

    all(
      abs(
        active_check$A_axis -
          active_check$D_axis *
            active_check$meanC_axis
      ) < 1e-10
    ),

    file.exists(
      file.path(
        TAB,
        "21F2B2_metric_dictionary.csv"
      )
    ),

    all(
      c(
        "Bulge_Isthmus",
        "ORS",
        "ORS_basal"
      ) %in%
        unique(
          follicular$receiver
        )
    ),

    "sender" %in%
      names(follicular),

    nrow(incoming_axis) > 0,

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  checks,
  n=Inf,
  width=Inf
)

stopifnot(
  length(checks$check) ==
    length(checks$pass),
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21F2B2_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F2B2 patient mechanistic LR-axis projection",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  "INPUTS",
  "------",
  paste0(
    "Patient-level active LR rows: ",
    nrow(net)
  ),
  paste0(
    "Structural sender-receiver LR rows: ",
    nrow(struct)
  ),
  paste0(
    "Evaluable patient sender-receiver edges: ",
    nrow(edge)
  ),
  paste0(
    "Pair-level mechanistic adjudications: ",
    nrow(adj)
  ),
  "",
  "PROJECTION",
  "----------",
  paste0(
    "Projectable INCLUDE/REVIEW LR pairs: ",
    nrow(projectable)
  ),
  paste0(
    "Mechanistically excluded LR pairs: ",
    nrow(excluded)
  ),
  paste0(
    "Structural mechanistic LR rows: ",
    nrow(struct_axis)
  ),
  paste0(
    "Sender-receiver-axis denominators: ",
    nrow(M_axis)
  ),
  paste0(
    "Projected active LR rows: ",
    nrow(active_axis)
  ),
  paste0(
    "Evaluable patient sender-receiver-axis rows: ",
    nrow(evaluable_axis)
  ),
  paste0(
    "Follicular sender-receiver-axis rows: ",
    nrow(follicular)
  ),
  "",
  "FROZEN METRICS",
  "--------------",
  "M_axis = structurally possible LR pairs within mechanistic class a for sender s and receiver r.",
  "N_axis = active LR richness within mechanistic class a.",
  "D_axis = N_axis / M_axis.",
  "meanC_axis = mean communication score among active LR pairs.",
  "S_axis = sum of active LR communication scores.",
  "A_axis = S_axis / M_axis.",
  "For active axes, A_axis = D_axis * meanC_axis.",
  "",
  "ZERO / NON-EVALUABLE SEMANTICS",
  "------------------------------",
  paste(
    "Patient sender-receiver combinations inherit evaluability",
    "from the frozen 21F1B edge universe."
  ),
  paste(
    "Within an evaluable edge, a structurally possible mechanistic",
    "axis with no active LR pair receives N_axis=0, S_axis=0,",
    "D_axis=0 and A_axis=0."
  ),
  paste(
    "Active-only intensity metrics remain NA when N_axis=0."
  ),
  paste(
    "Structurally impossible sender-receiver-axis combinations",
    "are absent rather than encoded as biological zero."
  ),
  "",
  "INTERPRETATION BOUNDARIES",
  "-------------------------",
  "EXCLUDE pairs are not projected into mechanistic communication metrics.",
  "REVIEW pairs remain explicitly contextual and are not silently promoted to primary biology.",
  "No disease comparison was performed.",
  "No AA-vs-FFA test was performed.",
  "No expression threshold was refit.",
  "No communication score was refit.",
  "No DE model was refit.",
  "No pathway model was refit.",
  "No annotation was modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21F2B2_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B2 COMPLETE\n")
cat("============================================================\n")
cat("Patient-level mechanistic LR-axis projection frozen\n")
cat("Axis-specific structural denominators generated\n")
cat("Richness, density, intensity and normalized exposure retained\n")
cat("Evaluable inactive axes encoded as zero exposure\n")
cat("Non-evaluable/structurally impossible axes not encoded as zero\n")
cat("Follicular receiver panel retained with sender identity\n")
cat("INCLUDE and REVIEW biology kept distinguishable\n")
cat("EXCLUDE pairs omitted from quantitative projection\n")
cat("No disease comparison\n")
cat("No analytical refit\n")
cat("PASS\n")
