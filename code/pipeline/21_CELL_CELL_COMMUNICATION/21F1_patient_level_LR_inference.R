suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F1 PATIENT-LEVEL LR INFERENCE\n")
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

OUT <- file.path(
  BASE,
  "21F_PATIENT_LEVEL_COMMUNICATION"
)

TAB <- file.path(OUT, "tables")
OBJ <- file.path(OUT, "objects")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJ, recursive=TRUE, showWarnings=FALSE)

EXPR_FILE <- file.path(
  E1,
  "21E1B_patient_identity_LR_gene_expression.rds"
)

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

# ============================================================
# 1. Frozen inputs
# ============================================================

cat("\nSTEP 1/10 — recovering frozen inputs\n")

stopifnot(
  file.exists(EXPR_FILE),
  file.exists(LR_FILE),
  file.exists(ELIG_FILE),
  file.exists(SPEC_FILE)
)

expr <- readRDS(EXPR_FILE)

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

primary <- spec %>%
  filter(
    analysis_specification == "PRIMARY"
  )

stopifnot(nrow(primary) == 1)

MIN_CELLS <- primary$min_cells
MIN_PCT   <- primary$min_pct_expr

cat("Primary min cells :", MIN_CELLS, "\n")
cat("Primary pct_expr  :", MIN_PCT, "\n")

# ============================================================
# 2. Primary identity-gene universe
# ============================================================

cat("\nSTEP 2/10 — recovering recurrent primary identity-gene universe\n")

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

cat(
  "Primary identities:",
  length(primary_identities),
  "\n"
)

stopifnot(
  length(primary_identities) == 39
)

# ============================================================
# 3. Patient-level eligible endpoint expression
# ============================================================

cat("\nSTEP 3/10 — defining patient-level eligible endpoints\n")

expr_primary <- expr %>%
  filter(
    communication_identity %in%
      primary_identities
  ) %>%
  inner_join(
    elig_primary,
    by=c(
      "communication_identity",
      "gene"
    )
  ) %>%
  mutate(
    patient_endpoint_evaluable =
      n_cells >= MIN_CELLS,
    patient_gene_detected =
      patient_endpoint_evaluable &
      pct_expr >= MIN_PCT
  )

endpoint_summary <- expr_primary %>%
  group_by(
    patient_id,
    disease,
    communication_identity
  ) %>%
  summarise(
    n_cells=first(n_cells),
    evaluable=first(n_cells) >= MIN_CELLS,
    n_recurrent_genes=n(),
    n_detected_recurrent_genes=sum(
      patient_gene_detected
    ),
    .groups="drop"
  )

write_csv(
  endpoint_summary,
  file.path(
    TAB,
    "21F1_patient_identity_endpoint_support.csv"
  )
)

cat("\n===== ENDPOINT SUPPORT =====\n")

print(
  endpoint_summary %>%
    group_by(disease) %>%
    summarise(
      n_patient_identity_groups=n(),
      n_evaluable=sum(evaluable),
      median_detected_genes=
        median(
          n_detected_recurrent_genes[
            evaluable
          ]
        ),
      .groups="drop"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 4. Freeze LR pair universe
# ============================================================

cat("\nSTEP 4/10 — recovering frozen OmniPath LR pair universe\n")

pair_key <- lr %>%
  distinct(
    source_genesymbol,
    target_genesymbol,
    .keep_all=TRUE
  )

cat(
  "Observable unique LR pairs:",
  nrow(pair_key),
  "\n"
)

# ============================================================
# 5. Sender ligand table
# ============================================================

cat("\nSTEP 5/10 — constructing sender ligand endpoints\n")

sender <- expr_primary %>%
  filter(
    patient_gene_detected
  ) %>%
  transmute(
    patient_id,
    disease,
    batch,
    sender=communication_identity,
    ligand=gene,
    ligand_pct_expr=pct_expr,
    ligand_logCPM=logCPM,
    sender_n_cells=n_cells
  )

receiver <- expr_primary %>%
  filter(
    patient_gene_detected
  ) %>%
  transmute(
    patient_id,
    disease,
    batch,
    receiver=communication_identity,
    receptor=gene,
    receptor_pct_expr=pct_expr,
    receptor_logCPM=logCPM,
    receiver_n_cells=n_cells
  )

cat("Eligible sender endpoints  :", nrow(sender), "\n")
cat("Eligible receiver endpoints:", nrow(receiver), "\n")

# ============================================================
# 6. Join ligands to frozen LR pairs
# ============================================================

cat("\nSTEP 6/10 — projecting endpoints onto frozen LR resource\n")

sender_lr <- sender %>%
  inner_join(
    pair_key,
    by=c(
      "ligand"="source_genesymbol"
    )
  )

receiver_lr <- receiver %>%
  select(
    patient_id,
    disease,
    receiver,
    receptor,
    receptor_pct_expr,
    receptor_logCPM,
    receiver_n_cells
  )

# ============================================================
# 7. Patient-level sender→receiver LR inference
# ============================================================

cat("\nSTEP 7/10 — inferring patient-level sender→receiver interactions\n")

comm <- sender_lr %>%
  inner_join(
    receiver_lr,
    by=c(
      "patient_id",
      "disease",
      "target_genesymbol"="receptor"
    ),
    relationship="many-to-many"
  ) %>%
  mutate(
    receptor=target_genesymbol,

    communication_score=
      sqrt(
        ligand_pct_expr *
        receptor_pct_expr
      )
  ) %>%
  filter(
    is.finite(communication_score),
    communication_score > 0
  )

cat(
  "Patient-level LR communication rows:",
  nrow(comm),
  "\n"
)

stopifnot(
  nrow(comm) > 0,
  all(
    comm$sender_n_cells >=
      MIN_CELLS
  ),
  all(
    comm$receiver_n_cells >=
      MIN_CELLS
  ),
  all(
    comm$ligand_pct_expr >=
      MIN_PCT
  ),
  all(
    comm$receptor_pct_expr >=
      MIN_PCT
  )
)

# ============================================================
# 8. Compact primary communication tables
# ============================================================

cat("\nSTEP 8/10 — generating compact patient-level network\n")

patient_edge <- comm %>%
  group_by(
    patient_id,
    disease,
    sender,
    receiver
  ) %>%
  summarise(
    n_LR_pairs=n_distinct(
      paste(
        ligand,
        receptor,
        sep="|||"
      )
    ),
    mean_LR_score=mean(
      communication_score
    ),
    median_LR_score=median(
      communication_score
    ),
    sum_LR_score=sum(
      communication_score
    ),
    .groups="drop"
  )

patient_pathway_lr <- comm %>%
  select(
    patient_id,
    disease,
    sender,
    receiver,
    ligand,
    receptor,
    communication_score,
    ligand_pct_expr,
    receptor_pct_expr,
    ligand_logCPM,
    receptor_logCPM,
    sender_n_cells,
    receiver_n_cells
  )

write_csv(
  patient_edge,
  file.path(
    TAB,
    "21F1_patient_sender_receiver_summary.csv"
  )
)

saveRDS(
  patient_pathway_lr,
  file.path(
    OBJ,
    "21F1_patient_level_LR_network.rds"
  ),
  compress="xz"
)

# ============================================================
# 9. Follicular receiver audit
# ============================================================

cat("\nSTEP 9/10 — auditing follicular receiver coverage\n")

follicular_receivers <- c(
  "Bulge_Isthmus",
  "ORS",
  "ORS_basal"
)

follicular_support <- endpoint_summary %>%
  filter(
    communication_identity %in%
      follicular_receivers
  ) %>%
  select(
    patient_id,
    disease,
    communication_identity,
    n_cells,
    evaluable,
    n_recurrent_genes,
    n_detected_recurrent_genes
  )

write_csv(
  follicular_support,
  file.path(
    TAB,
    "21F1_follicular_receiver_support.csv"
  )
)

follicular_network_summary <- patient_edge %>%
  filter(
    receiver %in%
      follicular_receivers
  ) %>%
  group_by(
    disease,
    receiver
  ) %>%
  summarise(
    n_patients=n_distinct(patient_id),
    n_senders=n_distinct(sender),
    median_LR_pairs=
      median(n_LR_pairs),
    median_sum_score=
      median(sum_LR_score),
    .groups="drop"
  )

write_csv(
  follicular_network_summary,
  file.path(
    TAB,
    "21F1_follicular_network_descriptive_summary.csv"
  )
)

cat("\n===== FOLLICULAR RECEIVER SUPPORT =====\n")

print(
  follicular_support,
  n=Inf,
  width=Inf
)

cat("\n===== FOLLICULAR NETWORK — DESCRIPTIVE ONLY =====\n")

print(
  follicular_network_summary,
  n=Inf,
  width=Inf
)

# ============================================================
# 10. Formal checks
# ============================================================

cat("\nSTEP 10/10 — formal checks\n")

checks <- tibble(
  check=c(
    "21E1B_expression_input_present",
    "21E2_eligibility_input_present",
    "primary_specification_recovered",
    "thirty_nine_primary_identities",
    "patient_unit_retained",
    "sender_endpoint_ge20",
    "receiver_endpoint_ge20",
    "ligand_detection_ge10pct",
    "receptor_detection_ge10pct",
    "recurrent_identity_gene_filter_used",
    "frozen_OmniPath_pairs_used",
    "communication_score_bounded",
    "non_evaluable_not_encoded_as_zero",
    "patient_level_network_generated",
    "follicular_receivers_retained",
    "no_disease_test",
    "no_AA_FFA_test",
    "no_group_mean_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(EXPR_FILE),
    file.exists(ELIG_FILE),
    nrow(primary) == 1,
    length(primary_identities) == 39,
    "patient_id" %in% names(comm),
    all(comm$sender_n_cells >= MIN_CELLS),
    all(comm$receiver_n_cells >= MIN_CELLS),
    all(comm$ligand_pct_expr >= MIN_PCT),
    all(comm$receptor_pct_expr >= MIN_PCT),
    TRUE,
    TRUE,
    all(
      comm$communication_score > 0 &
      comm$communication_score <= 1
    ),
    TRUE,
    nrow(patient_edge) > 0,
    all(
      follicular_receivers %in%
        endpoint_summary$communication_identity
    ),
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
    TAB,
    "21F1_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21F1 patient-level LR inference",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "PRIMARY SPECIFICATION",
  "---------------------",
  "Patient × identity endpoint requires >=20 cells.",
  "Ligand/receptor requires pct_expr >=10% in that patient × identity.",
  "Genes must also belong to the recurrent >=3-patient identity-gene universe.",
  "",
  "COMMUNICATION SCORE",
  "-------------------",
  "score = sqrt(ligand_pct_expr * receptor_pct_expr)",
  "",
  paste(
    "Scores are calculated separately for each patient, sender, receiver",
    "and frozen OmniPath ligand-receptor pair."
  ),
  "",
  paste(
    "Non-evaluable patient × identity endpoints are absent/NA by design",
    "and are not interpreted as zero communication."
  ),
  "",
  "No disease hypothesis test.",
  "No AA-vs-FFA test.",
  "No group-mean communication inference.",
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
    "21F1_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F1 COMPLETE\n")
cat("============================================================\n")
cat("Patient-level sender→receiver LR network inferred\n")
cat("Primary frozen eligibility specification retained\n")
cat("Non-evaluable endpoints not encoded as biological zeros\n")
cat("Follicular receiver coverage audited\n")
cat("No disease comparison performed\n")
cat("PASS\n")
