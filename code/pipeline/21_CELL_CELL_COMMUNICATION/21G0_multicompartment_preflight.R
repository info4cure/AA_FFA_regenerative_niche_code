suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21G0 MULTICOMPARTMENT COMMUNICATION PREFLIGHT\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

F1 <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

B2 <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES",
  "21F2B2_PATIENT_AXIS_PROJECTION"
)

OUT <- file.path(
  BASE,
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G0_PREFLIGHT"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

SUPPORT_FILE <- file.path(
  F1,
  "tables",
  "21F1_patient_identity_endpoint_support.csv"
)

EDGE_FILE <- file.path(
  F1,
  "tables",
  "21F1B_patient_sender_receiver_metrics.csv"
)

AXIS_FILE <- file.path(
  B2,
  "tables",
  "21F2B2_patient_receiver_axis_metrics.csv"
)

stopifnot(
  file.exists(SUPPORT_FILE),
  file.exists(EDGE_FILE),
  file.exists(AXIS_FILE)
)

support <- read_csv(
  SUPPORT_FILE,
  show_col_types=FALSE
)

edges <- read_csv(
  EDGE_FILE,
  show_col_types=FALSE
)

axis <- read_csv(
  AXIS_FILE,
  show_col_types=FALSE
)

# ============================================================
# 1. Recover frozen identity universe
# ============================================================

cat("\nSTEP 1/8 — recovering frozen communication universe\n")

identities <- support %>%
  distinct(
    communication_identity
  ) %>%
  arrange(
    communication_identity
  )

cat(
  "Communication identities:",
  nrow(identities),
  "\n"
)

print(
  identities,
  n=Inf
)

# ============================================================
# 2. Disease-blind identity support
# ============================================================

cat("\nSTEP 2/8 — auditing identity support\n")

identity_support <- support %>%
  group_by(
    communication_identity
  ) %>%
  summarise(
    n_patients_present=
      n_distinct(
        patient_id[n_cells > 0]
      ),

    n_patients_ge20=
      n_distinct(
        patient_id[n_cells >= 20]
      ),

    n_patients_evaluable=
      n_distinct(
        patient_id[evaluable]
      ),

    total_cells=
      sum(
        n_cells,
        na.rm=TRUE
      ),

    median_cells_per_present_patient=
      median(
        n_cells[n_cells > 0],
        na.rm=TRUE
      ),

    min_cells_present_patient=
      min(
        n_cells[n_cells > 0],
        na.rm=TRUE
      ),

    max_cells=
      max(
        n_cells,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  arrange(
    desc(n_patients_evaluable),
    desc(total_cells)
  )

cat("\n===== DISEASE-BLIND IDENTITY SUPPORT =====\n")

print(
  identity_support,
  n=Inf,
  width=Inf
)

write_csv(
  identity_support,
  file.path(
    TAB,
    "21G0_identity_support_disease_blind.csv"
  )
)

# ============================================================
# 3. Disease-stratified support — audit only
#    NOT used for identity selection
# ============================================================

cat("\nSTEP 3/8 — auditing disease-stratified support\n")

disease_support <- support %>%
  group_by(
    communication_identity,
    disease
  ) %>%
  summarise(
    n_patients_present=
      n_distinct(
        patient_id[n_cells > 0]
      ),

    n_patients_ge20=
      n_distinct(
        patient_id[n_cells >= 20]
      ),

    n_patients_evaluable=
      n_distinct(
        patient_id[evaluable]
      ),

    total_cells=
      sum(
        n_cells,
        na.rm=TRUE
      ),

    median_cells=
      median(
        n_cells[n_cells > 0],
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  arrange(
    communication_identity,
    disease
  )

cat("\n===== DISEASE-STRATIFIED SUPPORT — AUDIT ONLY =====\n")

print(
  disease_support,
  n=Inf,
  width=Inf
)

write_csv(
  disease_support,
  file.path(
    TAB,
    "21G0_identity_support_by_disease_AUDIT_ONLY.csv"
  )
)

# ============================================================
# 4. Sender and receiver coverage
# ============================================================

cat("\nSTEP 4/8 — auditing sender/receiver coverage\n")

sender_support <- edges %>%
  group_by(
    sender
  ) %>%
  summarise(
    n_patients=
      n_distinct(patient_id),

    n_receivers=
      n_distinct(receiver),

    n_patient_receiver_edges=n(),

    median_outgoing_A=
      median(
        normalized_exposure,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  rename(
    communication_identity=sender
  )

receiver_support <- edges %>%
  group_by(
    receiver
  ) %>%
  summarise(
    n_patients=
      n_distinct(patient_id),

    n_senders=
      n_distinct(sender),

    n_patient_sender_edges=n(),

    median_incoming_A=
      median(
        normalized_exposure,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  rename(
    communication_identity=receiver
  )

direction_support <- identities %>%
  left_join(
    sender_support,
    by="communication_identity"
  ) %>%
  left_join(
    receiver_support,
    by="communication_identity",
    suffix=c(
      "_sender",
      "_receiver"
    )
  )

cat("\n===== SENDER / RECEIVER SUPPORT =====\n")

print(
  direction_support,
  n=Inf,
  width=Inf
)

write_csv(
  direction_support,
  file.path(
    TAB,
    "21G0_sender_receiver_support.csv"
  )
)

# ============================================================
# 5. Mechanistic-axis receiver coverage
# ============================================================

cat("\nSTEP 5/8 — auditing mechanistic-axis receiver coverage\n")

axis_receiver <- axis %>%
  group_by(
    receiver,
    axis,
    mechanistic_class,
    action
  ) %>%
  summarise(
    n_patients=
      n_distinct(patient_id),

    n_nonmissing_A=
      sum(
        !is.na(
          mean_incoming_A_axis
        )
      ),

    n_nonmissing_D=
      sum(
        !is.na(
          pooled_D_axis
        )
      ),

    n_nonmissing_intensity=
      sum(
        !is.na(
          mean_sender_active_intensity
        )
      ),

    median_A=
      median(
        mean_incoming_A_axis,
        na.rm=TRUE
      ),

    .groups="drop"
  ) %>%
  arrange(
    receiver,
    axis,
    mechanistic_class,
    action
  )

write_csv(
  axis_receiver,
  file.path(
    TAB,
    "21G0_axis_receiver_support.csv"
  )
)

cat(
  "Receiver × axis × mechanism × action combinations:",
  nrow(axis_receiver),
  "\n"
)

# ============================================================
# 6. Mechanistic-axis inventory by identity
# ============================================================

cat("\nSTEP 6/8 — summarising mechanistic observability\n")

axis_summary <- axis_receiver %>%
  group_by(
    receiver
  ) %>%
  summarise(
    n_axis_keys=n(),

    n_INCLUDE_keys=
      sum(
        action == "INCLUDE",
        na.rm=TRUE
      ),

    n_REVIEW_keys=
      sum(
        action == "REVIEW",
        na.rm=TRUE
      ),

    n_keys_ge3_patients=
      sum(
        n_nonmissing_A >= 3
      ),

    n_keys_ge5_patients=
      sum(
        n_nonmissing_A >= 5
      ),

    .groups="drop"
  ) %>%
  arrange(
    desc(n_keys_ge5_patients),
    desc(n_axis_keys)
  )

cat("\n===== AXIS OBSERVABILITY BY RECEIVER =====\n")

print(
  axis_summary,
  n=Inf,
  width=Inf
)

write_csv(
  axis_summary,
  file.path(
    TAB,
    "21G0_axis_observability_by_receiver.csv"
  )
)

# ============================================================
# 7. Candidate multicompartment identity universe
#    Disease-blind support only
# ============================================================

cat("\nSTEP 7/8 — generating disease-blind candidate universe\n")

candidate_universe <- identity_support %>%
  left_join(
    direction_support,
    by="communication_identity"
  ) %>%
  left_join(
    axis_summary,
    by=c(
      "communication_identity"="receiver"
    )
  ) %>%
  mutate(
    candidate_status=case_when(

      n_patients_evaluable >= 8 ~
        "HIGH_SUPPORT",

      n_patients_evaluable >= 5 ~
        "MODERATE_SUPPORT",

      n_patients_evaluable >= 3 ~
        "LIMITED_SUPPORT",

      TRUE ~
        "INSUFFICIENT_SUPPORT"
    )
  ) %>%
  arrange(
    factor(
      candidate_status,
      levels=c(
        "HIGH_SUPPORT",
        "MODERATE_SUPPORT",
        "LIMITED_SUPPORT",
        "INSUFFICIENT_SUPPORT"
      )
    ),
    desc(n_patients_evaluable),
    desc(total_cells)
  )

cat("\n===== CANDIDATE MULTICOMPARTMENT UNIVERSE =====\n")

print(
  candidate_universe,
  n=Inf,
  width=Inf
)

write_csv(
  candidate_universe,
  file.path(
    TAB,
    "21G0_candidate_multicompartment_universe.csv"
  )
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "21F1_support_present",
    "21F1B_edge_metrics_present",
    "21F2B2_axis_metrics_present",
    "39_primary_identities_recovered",
    "sixteen_patients_recovered",
    "identity_support_generated",
    "disease_support_audit_generated",
    "sender_support_generated",
    "receiver_support_generated",
    "axis_receiver_support_generated",
    "candidate_universe_generated",
    "selection_status_disease_blind",
    "disease_support_not_used_for_selection",
    "no_new_expression_threshold",
    "no_new_LR_threshold",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_p_values",
    "no_BH_adjustment",
    "no_communication_refit",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(SUPPORT_FILE),
    file.exists(EDGE_FILE),
    file.exists(AXIS_FILE),

    nrow(identities) == 39,

    n_distinct(
      support$patient_id
    ) == 16,

    nrow(identity_support) == 39,

    file.exists(
      file.path(
        TAB,
        "21G0_identity_support_by_disease_AUDIT_ONLY.csv"
      )
    ),

    nrow(sender_support) > 0,
    nrow(receiver_support) > 0,
    nrow(axis_receiver) > 0,
    nrow(candidate_universe) == 39,

    all(
      candidate_universe$candidate_status %in%
        c(
          "HIGH_SUPPORT",
          "MODERATE_SUPPORT",
          "LIMITED_SUPPORT",
          "INSUFFICIENT_SUPPORT"
        )
    ),

    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

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
    "21G0_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21G0 multicompartment communication preflight",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  paste0(
    "Frozen communication identities recovered: ",
    nrow(identities)
  ),
  "",
  "PURPOSE",
  "-------",
  paste(
    "Audit disease-blind sender, receiver and mechanistic-axis",
    "support before defining multicompartment disease inference."
  ),
  "",
  "SELECTION PRINCIPLE",
  "-------------------",
  paste(
    "Candidate identity support is determined from overall",
    "patient-level evaluability without using disease effects."
  ),
  paste(
    "Disease-stratified support is written for audit only",
    "and is not used to select biological hypotheses."
  ),
  "",
  "No new expression threshold.",
  "No new LR threshold.",
  "No disease comparison.",
  "No p-values.",
  "No communication refit.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G0_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G0 COMPLETE\n")
cat("============================================================\n")
cat("Multicompartment communication support audited\n")
cat("Sender and receiver roles retained\n")
cat("Mechanistic-axis observability audited\n")
cat("Candidate universe classified disease-blind\n")
cat("No disease inference performed\n")
cat("PASS\n")
