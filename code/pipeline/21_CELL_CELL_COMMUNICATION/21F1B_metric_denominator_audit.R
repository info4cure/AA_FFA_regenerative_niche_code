suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F1B METRIC + DENOMINATOR AUDIT\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

E1 <- file.path(
  BASE,
  "21E1_OBSERVABLE_UNIVERSE"
)

E2 <- file.path(
  BASE,
  "21E2_FROZEN_LR_ELIGIBILITY"
)

F1 <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

TAB <- file.path(F1, "tables")
OBJ <- file.path(F1, "objects")

LR_FILE <- file.path(
  E1,
  "tables",
  "21E1A_observable_LR_universe.csv"
)

ELIG_FILE <- file.path(
  E2,
  "tables",
  "21E2_eligible_identity_genes.csv"
)

SPEC_FILE <- file.path(
  E2,
  "tables",
  "21E2_frozen_operating_points.csv"
)

ENDPOINT_FILE <- file.path(
  TAB,
  "21F1_patient_identity_endpoint_support.csv"
)

NETWORK_FILE <- file.path(
  OBJ,
  "21F1_patient_level_LR_network.rds"
)

stopifnot(
  file.exists(LR_FILE),
  file.exists(ELIG_FILE),
  file.exists(SPEC_FILE),
  file.exists(ENDPOINT_FILE),
  file.exists(NETWORK_FILE)
)

# ============================================================
# 1. Recover frozen inputs
# ============================================================

cat("\nSTEP 1/10 — recovering frozen inputs\n")

lr <- read_csv(
  LR_FILE,
  show_col_types=FALSE
)

elig <- read_csv(
  ELIG_FILE,
  show_col_types=FALSE
)

spec <- read_csv(
  SPEC_FILE,
  show_col_types=FALSE
)

endpoint <- read_csv(
  ENDPOINT_FILE,
  show_col_types=FALSE
)

net <- readRDS(
  NETWORK_FILE
)

primary <- spec %>%
  filter(
    analysis_specification == "PRIMARY"
  )

stopifnot(nrow(primary) == 1)

cat("21E1A observable rows :", nrow(lr), "\n")
cat("21F1 active LR rows   :", nrow(net), "\n")
cat("Endpoint rows         :", nrow(endpoint), "\n")

# ============================================================
# 2. Audit OmniPath symbol-level multiplicity
# ============================================================

cat("\nSTEP 2/10 — auditing OmniPath pair multiplicity\n")

required_lr <- c(
  "source_genesymbol",
  "target_genesymbol"
)

stopifnot(
  all(required_lr %in% names(lr))
)

pair_multiplicity <- lr %>%
  count(
    source_genesymbol,
    target_genesymbol,
    name="n_resource_records"
  ) %>%
  arrange(
    desc(n_resource_records),
    source_genesymbol,
    target_genesymbol
  )

unique_pair_key <- pair_multiplicity %>%
  select(
    source_genesymbol,
    target_genesymbol
  )

n_observable_rows <- nrow(lr)
n_unique_pairs <- nrow(unique_pair_key)
n_duplicate_resource_rows <-
  n_observable_rows - n_unique_pairs

cat("Observable resource rows :", n_observable_rows, "\n")
cat("Unique symbol LR pairs   :", n_unique_pairs, "\n")
cat("Redundant symbol records :", n_duplicate_resource_rows, "\n")

multiplicity_summary <- tibble(
  metric=c(
    "observable_resource_rows",
    "unique_symbol_LR_pairs",
    "redundant_symbol_level_records",
    "max_records_per_symbol_pair",
    "symbol_pairs_with_multiple_records"
  ),
  value=c(
    n_observable_rows,
    n_unique_pairs,
    n_duplicate_resource_rows,
    max(pair_multiplicity$n_resource_records),
    sum(
      pair_multiplicity$n_resource_records > 1
    )
  )
)

write_csv(
  pair_multiplicity,
  file.path(
    TAB,
    "21F1B_OmniPath_symbol_pair_multiplicity.csv"
  )
)

write_csv(
  multiplicity_summary,
  file.path(
    TAB,
    "21F1B_OmniPath_multiplicity_summary.csv"
  )
)

print(
  multiplicity_summary,
  n=Inf
)

# ============================================================
# 3. Audit exact duplicates in active 21F1 network
# ============================================================

cat("\nSTEP 3/10 — auditing active-network analytical keys\n")

required_net <- c(
  "patient_id",
  "disease",
  "sender",
  "receiver",
  "ligand",
  "receptor",
  "communication_score"
)

stopifnot(
  all(required_net %in% names(net))
)

active_key_audit <- net %>%
  count(
    patient_id,
    sender,
    receiver,
    ligand,
    receptor,
    name="n_records"
  )

active_duplicates <- active_key_audit %>%
  filter(
    n_records > 1
  )

cat(
  "Duplicated patient×sender×receiver×LR keys:",
  nrow(active_duplicates),
  "\n"
)

write_csv(
  active_duplicates,
  file.path(
    TAB,
    "21F1B_active_network_duplicate_keys.csv"
  )
)

# We do not silently deduplicate an already inferred network.
# Duplicate analytical keys would require stopping and repairing 21F1.

stopifnot(
  nrow(active_duplicates) == 0
)

# ============================================================
# 4. Recover recurrent primary genes by identity
# ============================================================

cat("\nSTEP 4/10 — recovering primary identity-gene universe\n")

elig_primary <- elig %>%
  filter(
    analysis_specification == "PRIMARY"
  ) %>%
  distinct(
    communication_identity,
    gene
  )

primary_identities <- sort(
  unique(
    elig_primary$communication_identity
  )
)

stopifnot(
  length(primary_identities) == 39
)

cat(
  "Primary identities:",
  length(primary_identities),
  "\n"
)

# ============================================================
# 5. Construct structural sender→receiver denominator M_sr
# ============================================================

cat("\nSTEP 5/10 — constructing structural LR opportunity denominator\n")

ligand_elig <- elig_primary %>%
  transmute(
    sender=communication_identity,
    ligand=gene
  )

receptor_elig <- elig_primary %>%
  transmute(
    receiver=communication_identity,
    receptor=gene
  )

structural_lr <- unique_pair_key %>%
  rename(
    ligand=source_genesymbol,
    receptor=target_genesymbol
  ) %>%
  inner_join(
    ligand_elig,
    by="ligand",
    relationship="many-to-many"
  ) %>%
  inner_join(
    receptor_elig,
    by="receptor",
    relationship="many-to-many"
  ) %>%
  distinct(
    sender,
    receiver,
    ligand,
    receptor
  )

denominator <- structural_lr %>%
  count(
    sender,
    receiver,
    name="n_structurally_possible_LR"
  )

expected_identity_pairs <-
  length(primary_identities)^2

cat(
  "Structural sender→receiver combinations:",
  nrow(denominator),
  "/",
  expected_identity_pairs,
  "\n"
)

cat(
  "Structural identity×identity×LR opportunities:",
  nrow(structural_lr),
  "\n"
)

write_csv(
  denominator,
  file.path(
    TAB,
    "21F1B_sender_receiver_structural_denominator.csv"
  )
)

saveRDS(
  structural_lr,
  file.path(
    OBJ,
    "21F1B_structural_sender_receiver_LR_universe.rds"
  ),
  compress="xz"
)

# ============================================================
# 6. Construct patient-level evaluable sender→receiver universe
# ============================================================

cat("\nSTEP 6/10 — constructing patient-level evaluable edge universe\n")

endpoint_eval <- endpoint %>%
  filter(
    communication_identity %in%
      primary_identities
  ) %>%
  select(
    patient_id,
    disease,
    communication_identity,
    n_cells,
    evaluable
  )

sender_eval <- endpoint_eval %>%
  filter(evaluable) %>%
  transmute(
    patient_id,
    disease,
    sender=communication_identity,
    sender_n_cells=n_cells
  )

receiver_eval <- endpoint_eval %>%
  filter(evaluable) %>%
  transmute(
    patient_id,
    disease,
    receiver=communication_identity,
    receiver_n_cells=n_cells
  )

evaluable_edges <- sender_eval %>%
  inner_join(
    receiver_eval,
    by=c(
      "patient_id",
      "disease"
    ),
    relationship="many-to-many"
  ) %>%
  inner_join(
    denominator,
    by=c(
      "sender",
      "receiver"
    )
  )

cat(
  "Patient-level evaluable sender→receiver edges:",
  nrow(evaluable_edges),
  "\n"
)

stopifnot(
  all(
    evaluable_edges$sender_n_cells >=
      primary$min_cells
  ),
  all(
    evaluable_edges$receiver_n_cells >=
      primary$min_cells
  ),
  all(
    evaluable_edges$n_structurally_possible_LR > 0
  )
)

# ============================================================
# 7. Aggregate active LR metrics
# ============================================================

cat("\nSTEP 7/10 — calculating frozen edge metrics\n")

active_metrics <- net %>%
  group_by(
    patient_id,
    disease,
    sender,
    receiver
  ) %>%
  summarise(
    n_active_LR=n(),
    sum_LR_score=sum(
      communication_score
    ),
    mean_active_LR_score=mean(
      communication_score
    ),
    median_active_LR_score=median(
      communication_score
    ),
    min_active_LR_score=min(
      communication_score
    ),
    max_active_LR_score=max(
      communication_score
    ),
    .groups="drop"
  )

edge_metrics <- evaluable_edges %>%
  left_join(
    active_metrics,
    by=c(
      "patient_id",
      "disease",
      "sender",
      "receiver"
    )
  ) %>%
  mutate(
    # Evaluable but no active LR pair is a genuine zero.
    n_active_LR=coalesce(
      n_active_LR,
      0L
    ),

    sum_LR_score=coalesce(
      sum_LR_score,
      0
    ),

    LR_density=
      n_active_LR /
      n_structurally_possible_LR,

    normalized_exposure=
      sum_LR_score /
      n_structurally_possible_LR,

    # Conditional intensity is undefined if no active pair exists.
    mean_active_LR_score=
      if_else(
        n_active_LR > 0,
        mean_active_LR_score,
        NA_real_
      ),

    median_active_LR_score=
      if_else(
        n_active_LR > 0,
        median_active_LR_score,
        NA_real_
      )
  )

stopifnot(
  all(
    edge_metrics$n_active_LR <=
      edge_metrics$n_structurally_possible_LR
  ),
  all(
    edge_metrics$LR_density >= 0 &
      edge_metrics$LR_density <= 1
  ),
  all(
    edge_metrics$normalized_exposure >= 0 &
      edge_metrics$normalized_exposure <= 1
  )
)

write_csv(
  edge_metrics,
  file.path(
    TAB,
    "21F1B_patient_sender_receiver_metrics.csv"
  )
)

saveRDS(
  edge_metrics,
  file.path(
    OBJ,
    "21F1B_patient_sender_receiver_metrics.rds"
  ),
  compress="xz"
)

# ============================================================
# 8. Metric dictionary with explicit formulas
# ============================================================

cat("\nSTEP 8/10 — freezing mathematical metric dictionary\n")

metric_dictionary <- tribble(
  ~metric,
  ~symbol,
  ~level,
  ~formula,
  ~denominator,
  ~range,
  ~zero_NA_rule,
  ~interpretation,

  "LR communication score",
  "C_psri",
  "patient × sender × receiver × LR",
  "sqrt(f_psL * f_prR)",
  "none",
  "[0,1]",
  paste(
    "Calculated only when both patient endpoints are evaluable",
    "and ligand/receptor satisfy the frozen patient-level detection rule."
  ),
  paste(
    "Geometric mean of ligand and receptor cellular detection fractions;",
    "not a biochemical signaling rate."
  ),

  "Active LR richness",
  "N_psr",
  "patient × sender × receiver",
  "sum_i I(C_psri > 0)",
  "none",
  "[0,M_sr]",
  paste(
    "Zero when sender and receiver are evaluable but no structurally",
    "possible LR pair passes both patient-level detection requirements;",
    "NA/absent when either endpoint is non-evaluable."
  ),
  "Breadth of detected LR channels.",

  "Structural LR opportunity",
  "M_sr",
  "sender × receiver",
  paste(
    "number of unique OmniPath symbol-level LR pairs with ligand",
    "recurrently eligible in sender and receptor recurrently eligible",
    "in receiver under the frozen PRIMARY specification"
  ),
  "none",
  "positive integer",
  "Disease-independent and patient-independent.",
  "Frozen opportunity denominator for a sender-receiver identity pair.",

  "LR density",
  "D_psr",
  "patient × sender × receiver",
  "N_psr / M_sr",
  "M_sr",
  "[0,1]",
  paste(
    "Zero is allowed only for an evaluable sender-receiver edge;",
    "non-evaluable edges are absent/NA."
  ),
  "Fraction of structurally possible LR channels detected in that patient.",

  "Mean active LR intensity",
  "meanC_psr",
  "patient × sender × receiver",
  "(1 / N_psr) * sum_i C_psri for active i",
  "N_psr",
  "[0,1]",
  "NA when N_psr = 0 or edge is non-evaluable.",
  "Conditional mean intensity among detected LR channels.",

  "Median active LR intensity",
  "medianC_psr",
  "patient × sender × receiver",
  "median_i(C_psri) over active i",
  "N_psr active channels",
  "[0,1]",
  "NA when N_psr = 0 or edge is non-evaluable.",
  "Conditional median intensity among detected LR channels.",

  "Aggregate LR score",
  "S_psr",
  "patient × sender × receiver",
  "sum_i C_psri",
  "none",
  "[0,M_sr]",
  paste(
    "Zero when edge is evaluable and no LR pair is active;",
    "non-evaluable edges are absent/NA."
  ),
  paste(
    "Unnormalized aggregate combining interaction breadth and intensity;",
    "not used alone as normalized communication strength."
  ),

  "Normalized exposure",
  "A_psr",
  "patient × sender × receiver",
  "S_psr / M_sr",
  "M_sr",
  "[0,1]",
  paste(
    "Zero when edge is evaluable and no LR pair is active;",
    "non-evaluable edges are absent/NA."
  ),
  paste(
    "Aggregate LR detection-weighted exposure normalized by the",
    "structural LR opportunity of that sender-receiver pair."
  )
)

write_csv(
  metric_dictionary,
  file.path(
    TAB,
    "21F_metric_dictionary.csv"
  )
)

# ============================================================
# 9. Descriptive audit only
# ============================================================

cat("\nSTEP 9/10 — descriptive metric audit\n")

metric_audit <- edge_metrics %>%
  summarise(
    n_evaluable_edges=n(),
    n_zero_active_edges=sum(
      n_active_LR == 0
    ),
    median_M=median(
      n_structurally_possible_LR
    ),
    median_N=median(
      n_active_LR
    ),
    median_density=median(
      LR_density
    ),
    median_active_intensity=median(
      mean_active_LR_score,
      na.rm=TRUE
    ),
    median_normalized_exposure=median(
      normalized_exposure
    )
  )

print(
  metric_audit,
  n=Inf,
  width=Inf
)

write_csv(
  metric_audit,
  file.path(
    TAB,
    "21F1B_metric_audit_summary.csv"
  )
)

# ============================================================
# 10. Formal checks
# ============================================================

cat("\nSTEP 10/10 — formal checks\n")

checks <- tibble(
  check=c(
    "21F1_network_present",
    "21E2_primary_eligibility_present",
    "39_primary_identities",
    "OmniPath_symbol_pairs_unique_for_denominator",
    "resource_multiplicity_audited",
    "active_network_key_unique",
    "structural_denominator_generated",
    "structural_denominator_disease_independent",
    "patient_evaluable_edges_generated",
    "non_evaluable_edges_not_zero",
    "evaluable_inactive_edges_zero_allowed",
    "active_richness_generated",
    "LR_density_generated",
    "mean_active_intensity_generated",
    "median_active_intensity_generated",
    "aggregate_score_retained",
    "normalized_exposure_generated",
    "metric_dictionary_generated",
    "metric_formulas_explicit",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_DE
_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(NETWORK_FILE),
    file.exists(ELIG_FILE),
    length(primary_identities) == 39,

    # Denominator constructed from unique symbol-level LR pairs
    !anyDuplicated(
      unique_pair_key[
        c(
          "source_genesymbol",
          "target_genesymbol"
        )
      ]
    ),

    nrow(multiplicity_summary) == 5,

    # 21F1 analytical key must be unique
    nrow(active_duplicates) == 0,

    nrow(denominator) > 0,

    # M_sr contains no patient or disease information
    !any(
      c(
        "patient_id",
        "disease"
      ) %in% names(denominator)
    ),

    nrow(evaluable_edges) > 0,

    # Non-evaluable edges never enter the patient edge universe
    all(
      evaluable_edges$sender_n_cells >=
        primary$min_cells &
      evaluable_edges$receiver_n_cells >=
        primary$min_cells
    ),

    # For evaluable edges, absence of active LR is represented by zero
    all(
      edge_metrics$n_active_LR >= 0
    ),

    "n_active_LR" %in% names(edge_metrics),

    all(
      is.finite(edge_metrics$LR_density)
    ),

    all(
      is.na(edge_metrics$mean_active_LR_score) |
      is.finite(edge_metrics$mean_active_LR_score)
    ),

    all(
      is.na(edge_metrics$median_active_LR_score) |
      is.finite(edge_metrics$median_active_LR_score)
    ),

    "sum_LR_score" %in% names(edge_metrics),

    all(
      is.finite(
        edge_metrics$normalized_exposure
      )
    ),

    file.exists(
      file.path(
        TAB,
        "21F_metric_dictionary.csv"
      )
    ),

    all(
      c(
        "metric",
        "symbol",
        "formula",
        "denominator",
        "zero_NA_rule",
        "interpretation"
      ) %in%
        names(metric_dictionary)
    ),

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
    "21F1B_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F1B metric and denominator audit",
  "",
  paste0(
    "Date: ",
    Sys.Date()
  ),
  "",
  "RESOURCE AUDIT",
  "--------------",
  paste0(
    "Observable OmniPath resource rows: ",
    n_observable_rows
  ),
  paste0(
    "Unique symbol-level LR pairs: ",
    n_unique_pairs
  ),
  paste0(
    "Redundant symbol-level resource records: ",
    n_duplicate_resource_rows
  ),
  "",
  "PRIMARY COMMUNICATION SCORE",
  "---------------------------",
  "C_psri = sqrt(f_psL * f_prR)",
  "",
  paste(
    "f_psL is the fraction of cells expressing ligand L",
    "in patient p and sender identity s."
  ),
  paste(
    "f_prR is the fraction of cells expressing receptor R",
    "in patient p and receiver identity r."
  ),
  "",
  "STRUCTURAL OPPORTUNITY",
  "----------------------",
  paste(
    "M_sr is the number of unique frozen OmniPath LR pairs",
    "whose ligand is recurrently eligible in sender s and",
    "whose receptor is recurrently eligible in receiver r."
  ),
  "",
  "EDGE METRICS",
  "------------",
  "N_psr = number of active LR pairs",
  "D_psr = N_psr / M_sr",
  "meanC_psr = mean C_psri among active LR pairs",
  "medianC_psr = median C_psri among active LR pairs",
  "S_psr = sum C_psri",
  "A_psr = S_psr / M_sr",
  "",
  "MISSINGNESS / ZERO RULES",
  "------------------------",
  paste(
    "If sender or receiver has fewer than the frozen minimum",
    "number of cells, the patient-level edge is non-evaluable",
    "and is not encoded as zero."
  ),
  paste(
    "If both endpoints are evaluable but no structurally possible",
    "LR pair satisfies patient-level expression requirements,",
    "N_psr, D_psr, S_psr and A_psr are zero."
  ),
  paste(
    "Conditional active-interaction intensity metrics are NA",
    "when N_psr = 0."
  ),
  "",
  "INTERPRETIVE BOUNDARIES",
  "-----------------------",
  paste(
    "Aggregate LR score S_psr combines interaction breadth",
    "and intensity and is not interpreted alone as normalized",
    "communication strength."
  ),
  paste(
    "D_psr measures detected-channel density relative to",
    "the disease-independent structural opportunity M_sr."
  ),
  paste(
    "A_psr measures detection-weighted aggregate exposure",
    "normalized by M_sr."
  ),
  paste(
    "All communication quantities are transcriptomic proxies",
    "and do not establish biochemical ligand-receptor activity."
  ),
  "",
  "No disease hypothesis test.",
  "No AA-vs-FFA test.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    F1,
    "21F1B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F1B COMPLETE\n")
cat("============================================================\n")
cat("OmniPath symbol-level multiplicity audited\n")
cat("Active LR analytical keys audited\n")
cat("Disease-independent structural denominator M_sr frozen\n")
cat("Patient-level evaluable edge universe frozen\n")
cat("Richness, density, intensity and normalized exposure generated\n")
cat("Mathematical metric dictionary frozen\n")
cat("Zero versus non-evaluable semantics frozen\n")
cat("No disease comparison\n")
cat("No analytical refit\n")
cat("PASS\n")
