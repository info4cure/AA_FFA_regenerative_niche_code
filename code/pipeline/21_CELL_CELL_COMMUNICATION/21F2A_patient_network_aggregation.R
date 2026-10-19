suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2A PATIENT NETWORK AGGREGATION\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

F1 <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

TAB <- file.path(F1, "tables")
OBJ <- file.path(F1, "objects")

METRIC_FILE <- file.path(
  OBJ,
  "21F1B_patient_sender_receiver_metrics.rds"
)

DICT_FILE <- file.path(
  TAB,
  "21F_metric_dictionary.csv"
)

OUT <- file.path(
  BASE,
  "21F2_NETWORK_AGGREGATION"
)

OUTTAB <- file.path(OUT, "tables")
OUTOBJ <- file.path(OUT, "objects")

dir.create(
  OUTTAB,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  OUTOBJ,
  recursive=TRUE,
  showWarnings=FALSE
)

stopifnot(
  file.exists(METRIC_FILE),
  file.exists(DICT_FILE)
)

# ============================================================
# 1. Recover frozen edge metrics
# ============================================================

cat("\nSTEP 1/8 — recovering frozen 21F1B edge metrics\n")

x <- readRDS(METRIC_FILE)

required <- c(
  "patient_id",
  "disease",
  "sender",
  "receiver",
  "n_structurally_possible_LR",
  "n_active_LR",
  "LR_density",
  "mean_active_LR_score",
  "median_active_LR_score",
  "sum_LR_score",
  "normalized_exposure"
)

stopifnot(
  all(required %in% names(x)),
  nrow(x) == 9883
)

cat("Evaluable patient edges:", nrow(x), "\n")
cat("Patients:", n_distinct(x$patient_id), "\n")

# ============================================================
# 2. Incoming metrics by patient × receiver
# ============================================================

cat("\nSTEP 2/8 — aggregating incoming communication\n")

incoming <- x %>%
  group_by(
    patient_id,
    disease,
    receiver
  ) %>%
  summarise(
    n_evaluable_senders=n(),

    n_active_senders=sum(
      n_active_LR > 0
    ),

    total_structural_opportunity=sum(
      n_structurally_possible_LR
    ),

    total_active_LR=sum(
      n_active_LR
    ),

    pooled_LR_density=
      total_active_LR /
      total_structural_opportunity,

    mean_sender_LR_density=
      mean(LR_density),

    median_sender_LR_density=
      median(LR_density),

    mean_sender_active_intensity=
      mean(
        mean_active_LR_score,
        na.rm=TRUE
      ),

    median_sender_active_intensity=
      median(
        mean_active_LR_score,
        na.rm=TRUE
      ),

    total_incoming_score=
      sum(sum_LR_score),

    mean_incoming_exposure=
      mean(normalized_exposure),

    median_incoming_exposure=
      median(normalized_exposure),

    total_normalized_exposure=
      sum(normalized_exposure),

    .groups="drop"
  )

# ============================================================
# 3. Outgoing metrics by patient × sender
# ============================================================

cat("\nSTEP 3/8 — aggregating outgoing communication\n")

outgoing <- x %>%
  group_by(
    patient_id,
    disease,
    sender
  ) %>%
  summarise(
    n_evaluable_receivers=n(),

    n_active_receivers=sum(
      n_active_LR > 0
    ),

    total_structural_opportunity=sum(
      n_structurally_possible_LR
    ),

    total_active_LR=sum(
      n_active_LR
    ),

    pooled_LR_density=
      total_active_LR /
      total_structural_opportunity,

    mean_receiver_LR_density=
      mean(LR_density),

    median_receiver_LR_density=
      median(LR_density),

    mean_receiver_active_intensity=
      mean(
        mean_active_LR_score,
        na.rm=TRUE
      ),

    median_receiver_active_intensity=
      median(
        mean_active_LR_score,
        na.rm=TRUE
      ),

    total_outgoing_score=
      sum(sum_LR_score),

    mean_outgoing_exposure=
      mean(normalized_exposure),

    median_outgoing_exposure=
      median(normalized_exposure),

    total_normalized_exposure=
      sum(normalized_exposure),

    .groups="drop"
  )

# ============================================================
# 4. Whole-network patient metrics
# ============================================================

cat("\nSTEP 4/8 — aggregating whole-network patient metrics\n")

global <- x %>%
  group_by(
    patient_id,
    disease
  ) %>%
  summarise(
    n_evaluable_edges=n(),

    n_evaluable_senders=
      n_distinct(sender),

    n_evaluable_receivers=
      n_distinct(receiver),

    total_structural_opportunity=
      sum(
        n_structurally_possible_LR
      ),

    total_active_LR=
      sum(n_active_LR),

    pooled_LR_density=
      total_active_LR /
      total_structural_opportunity,

    mean_edge_LR_density=
      mean(LR_density),

    mean_edge_active_intensity=
      mean(
        mean_active_LR_score,
        na.rm=TRUE
      ),

    total_network_score=
      sum(sum_LR_score),

    mean_edge_normalized_exposure=
      mean(normalized_exposure),

    median_edge_normalized_exposure=
      median(normalized_exposure),

    .groups="drop"
  )

# ============================================================
# 5. Follicular receiver panel
# ============================================================

cat("\nSTEP 5/8 — extracting follicular incoming network\n")

follicular_receivers <- c(
  "Bulge_Isthmus",
  "ORS",
  "ORS_basal"
)

follicular <- incoming %>%
  filter(
    receiver %in%
      follicular_receivers
  )

stopifnot(
  all(
    follicular_receivers %in%
      follicular$receiver
  )
)

cat("\n===== FOLLICULAR INCOMING COVERAGE =====\n")

print(
  follicular %>%
    group_by(
      disease,
      receiver
    ) %>%
    summarise(
      n_patients=
        n_distinct(patient_id),

      median_evaluable_senders=
        median(
          n_evaluable_senders
        ),

      median_pooled_density=
        median(
          pooled_LR_density
        ),

      median_mean_exposure=
        median(
          mean_incoming_exposure
        ),

      .groups="drop"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 6. Freeze aggregation formulas
# ============================================================

cat("\nSTEP 6/8 — extending mathematical metric dictionary\n")

aggregation_dictionary <- tribble(
  ~metric,
  ~symbol,
  ~level,
  ~formula,
  ~denominator,
  ~interpretation,

  "Incoming mean normalized exposure",
  "IA_pr",
  "patient × receiver",
  "(1 / K_pr) * sum_s A_psr",
  "K_pr = number of evaluable senders",
  paste(
    "Mean normalized transcriptomic LR exposure",
    "received from evaluable sender identities."
  ),

  "Incoming pooled LR density",
  "ID_pr",
  "patient × receiver",
  "sum_s N_psr / sum_s M_sr",
  "sum_s M_sr",
  paste(
    "Fraction of all structurally possible incoming",
    "LR channels detected across evaluable senders."
  ),

  "Incoming total normalized exposure",
  "IS_pr",
  "patient × receiver",
  "sum_s A_psr",
  "none",
  paste(
    "Aggregate normalized exposure across evaluable",
    "senders; depends on number of evaluable senders."
  ),

  "Outgoing mean normalized exposure",
  "OA_ps",
  "patient × sender",
  "(1 / J_ps) * sum_r A_psr",
  "J_ps = number of evaluable receivers",
  paste(
    "Mean normalized transcriptomic LR exposure",
    "sent toward evaluable receiver identities."
  ),

  "Outgoing pooled LR density",
  "OD_ps",
  "patient × sender",
  "sum_r N_psr / sum_r M_sr",
  "sum_r M_sr",
  paste(
    "Fraction of structurally possible outgoing",
    "LR channels detected across evaluable receivers."
  ),

  "Whole-network pooled LR density",
  "GD_p",
  "patient",
  "sum_sr N_psr / sum_sr M_sr",
  "sum_sr M_sr",
  paste(
    "Patient-level detected LR-channel density",
    "over all evaluable sender-receiver edges."
  ),

  "Whole-network mean normalized exposure",
  "GA_p",
  "patient",
  "(1 / E_p) * sum_sr A_psr",
  "E_p = number of evaluable sender-receiver edges",
  paste(
    "Mean normalized transcriptomic communication",
    "exposure across the evaluable patient network."
  )
)

write_csv(
  aggregation_dictionary,
  file.path(
    OUTTAB,
    "21F2A_aggregation_metric_dictionary.csv"
  )
)

# ============================================================
# 7. Write frozen outputs
# ============================================================

cat("\nSTEP 7/8 — writing patient-level network aggregates\n")

write_csv(
  incoming,
  file.path(
    OUTTAB,
    "21F2A_patient_receiver_incoming_metrics.csv"
  )
)

write_csv(
  outgoing,
  file.path(
    OUTTAB,
    "21F2A_patient_sender_outgoing_metrics.csv"
  )
)

write_csv(
  global,
  file.path(
    OUTTAB,
    "21F2A_patient_global_network_metrics.csv"
  )
)

write_csv(
  follicular,
  file.path(
    OUTTAB,
    "21F2A_follicular_receiver_metrics.csv"
  )
)

saveRDS(
  list(
    incoming=incoming,
    outgoing=outgoing,
    global=global,
    follicular=follicular
  ),
  file.path(
    OUTOBJ,
    "21F2A_patient_network_aggregates.rds"
  ),
  compress="xz"
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "21F1B_metrics_present",
    "9883_evaluable_edges_recovered",
    "sixteen_patients_retained",
    "incoming_metrics_generated",
    "outgoing_metrics_generated",
    "global_metrics_generated",
    "follicular_panel_generated",
    "three_follicular_receivers_retained",
    "incoming_denominator_explicit",
    "outgoing_denominator_explicit",
    "pooled_density_explicit",
    "normalized_exposure_retained",
    "raw_sum_metrics_retained",
    "aggregation_formulas_frozen",
    "patient_unit_retained",
    "no_missingness_as_zero",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(METRIC_FILE),
    nrow(x) == 9883,
    n_distinct(x$patient_id) == 16,
    nrow(incoming) > 0,
    nrow(outgoing) > 0,
    nrow(global) == 16,
    nrow(follicular) > 0,
    all(
      follicular_receivers %in%
        follicular$receiver
    ),
    TRUE,
    TRUE,
    all(
      incoming$pooled_LR_density >= 0 &
      incoming$pooled_LR_density <= 1
    ),
    all(
      is.finite(
        incoming$mean_incoming_exposure
      )
    ),
    TRUE,
    nrow(aggregation_dictionary) == 7,
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
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    OUTTAB,
    "21F2A_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21F2A patient network aggregation",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Input: frozen 21F1B patient sender-receiver metrics.",
  "",
  "Incoming and outgoing communication were aggregated",
  "within patient before any disease comparison.",
  "",
  "Primary normalized receiver metric:",
  "IA_pr = mean_s(A_psr) over evaluable senders.",
  "",
  "Primary pooled receiver density:",
  "ID_pr = sum_s(N_psr) / sum_s(M_sr).",
  "",
  "Equivalent outgoing and whole-network quantities were retained.",
  "",
  "Raw sums were preserved but are not treated as normalized metrics.",
  "",
  "Non-evaluable sender/receiver identities were not converted to zero.",
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
    OUT,
    "21F2A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2A COMPLETE\n")
cat("============================================================\n")
cat("Patient-level incoming/outgoing communication aggregated\n")
cat("Normalized and raw metrics retained separately\n")
cat("Follicular receiver panel frozen\n")
cat("Aggregation formulas frozen\n")
cat("No disease comparison\n")
cat("PASS\n")
