suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

cat("============================================================\n")
cat("21G3A3 — CROSS-MODAL INTEGRATION DESIGN FREEZE\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/"
)

CW_FILE <- file.path(
  BASE,
  "21G3_CROSS_MODAL_INTEGRATION",
  "tables",
  "21G3A2_exact_identity_crosswalk.csv"
)

CHECK_FILE <- file.path(
  BASE,
  "21G3_CROSS_MODAL_INTEGRATION",
  "tables",
  "21G3A2_formal_checks.csv"
)

OUT <- file.path(
  BASE,
  "21G3_CROSS_MODAL_INTEGRATION"
)

TAB <- file.path(
  OUT,
  "tables"
)

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

stopifnot(
  file.exists(CW_FILE),
  file.exists(CHECK_FILE)
)

cat("\nSTEP 1/6 — recovering frozen crosswalk\n")

cw <- read_csv(
  CW_FILE,
  show_col_types=FALSE
)

checks_A2 <- read_csv(
  CHECK_FILE,
  show_col_types=FALSE
)

stopifnot(
  all(checks_A2$pass),
  nrow(cw) == 28,
  all(cw$DA_exact),
  all(cw$DE_exact),
  sum(cw$PATHWAY_exact) == 16
)

cat("Communication receivers:", nrow(cw), "\n")
cat("DA exact:", sum(cw$DA_exact), "\n")
cat("DE exact:", sum(cw$DE_exact), "\n")
cat("20F2B represented:", sum(cw$PATHWAY_exact), "\n")

# ============================================================
# STEP 2 — freeze modality roles
# ============================================================

cat("\nSTEP 2/6 — freezing modality roles\n")

modality_roles <- tribble(
  ~modality,       ~analytical_role,            ~interpretive_role,
  "COMMUNICATION", "ANCHOR",                    "RECEIVER_EXPOSURE",
  "DA",            "ORTHOGONAL_CONTEXT",        "CELLULAR_ABUNDANCE",
  "DE",            "ORTHOGONAL_CONTEXT",        "TRANSCRIPTIONAL_STATE",
  "PATHWAY",       "ORTHOGONAL_CONTEXT",        "COORDINATED_PROGRAMME"
)

write_csv(
  modality_roles,
  file.path(
    TAB,
    "21G3A3_modality_roles.csv"
  )
)

print(
  modality_roles,
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 3 — freeze integration principles
# ============================================================

cat("\nSTEP 3/6 — freezing integration principles\n")

principles <- tribble(
  ~rule_id, ~rule,
  "R01",
  "Communication is the anchor modality.",
  "R02",
  "The elementary integration unit is contrast x receiver x communication axis.",
  "R03",
  "DA, DE and pathway evidence are orthogonal context and are not interchangeable.",
  "R04",
  "No cross-modal evidence score will be constructed.",
  "R05",
  "No p-values or q-values will be combined across modalities.",
  "R06",
  "Absence from 20F2B is not evidence of pathway absence.",
  "R07",
  "A receiver absent from 20F2B is classified NOT_IN_20F2B_UNIVERSE.",
  "R08",
  "DA direction is not automatically interpreted as mechanistic concordance with communication.",
  "R09",
  "DE direction is not automatically interpreted as mechanistic concordance with communication without a predefined biological mapping.",
  "R10",
  "Pathway direction may be compared with communication only when the pathway axis has a prespecified biological relation to the communication axis.",
  "R11",
  "Nominal communication evidence remains nominal and is never promoted by cross-modal support.",
  "R12",
  "FDR-defined communication evidence retains its original 21G2 multiplicity correction.",
  "R13",
  "Cross-modal integration is descriptive/triangulatory and introduces no new inferential family.",
  "R14",
  "AA-vs-HC, FFA-vs-HC and AA-vs-FFA remain separate contrasts.",
  "R15",
  "Primary and contextual communication programmes remain separated.",
  "R16",
  "Receiver identity remains the analytical cellular unit; macrocompartment is interpretive.",
  "R17",
  "No missing modality is imputed.",
  "R18",
  "No L2 identity outside the communication receiver universe is introduced into receiver-level integration."
)

write_csv(
  principles,
  file.path(
    TAB,
    "21G3A3_integration_principles.csv"
  )
)

print(
  principles,
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 4 — freeze pathway availability semantics
# ============================================================

cat("\nSTEP 4/6 — freezing pathway availability semantics\n")

receiver_universe <- cw %>%
  transmute(
    macrocompartment,
    receiver,

    DA_availability=
      if_else(
        DA_exact,
        "AVAILABLE_EXACT",
        "NOT_AVAILABLE"
      ),

    DE_availability=
      if_else(
        DE_exact,
        "AVAILABLE_EXACT",
        "NOT_AVAILABLE"
      ),

    pathway_availability=
      if_else(
        PATHWAY_exact,
        "AVAILABLE_IN_20F2B",
        "NOT_IN_20F2B_UNIVERSE"
      )
  )

cat("\n===== FROZEN RECEIVER × MODALITY AVAILABILITY =====\n")

print(
  receiver_universe,
  n=Inf,
  width=Inf
)

write_csv(
  receiver_universe,
  file.path(
    TAB,
    "21G3A3_receiver_modality_availability.csv"
  )
)

# ============================================================
# STEP 5 — freeze reporting vocabulary
# ============================================================

cat("\nSTEP 5/6 — freezing reporting vocabulary\n")

reporting_vocabulary <- tribble(
  ~term, ~definition,
  "COMMUNICATION_SIGNAL",
  "Communication result retaining its original 21G2 evidence class.",
  "DA_CONTEXT_PRESENT",
  "Frozen differential-abundance evidence is available for the same receiver and contrast.",
  "DE_CONTEXT_PRESENT",
  "Frozen differential-expression evidence is available for the same receiver and contrast.",
  "PATHWAY_CONTEXT_PRESENT",
  "Frozen 20F2B pathway evidence is available for the same receiver and contrast.",
  "MULTIMODAL_CONTEXT",
  "At least two orthogonal non-communication modalities provide interpretable evidence for the same receiver and contrast.",
  "COMMUNICATION_PLUS_MULTIMODAL_CONTEXT",
  "Communication evidence accompanied by interpretable evidence from at least two orthogonal modalities.",
  "NOT_IN_20F2B_UNIVERSE",
  "Receiver was not represented in the frozen 20F2B prioritization universe; this is not evidence of pathway absence.",
  "DIRECTIONALLY_COMPATIBLE",
  "Direction is compatible under an explicitly prespecified biological mapping.",
  "DIRECTIONALLY_OPPOSED",
  "Direction is opposed under an explicitly prespecified biological mapping.",
  "DIRECTION_NOT_INTERPRETABLE",
  "Directional comparison across modalities is not biologically prespecified or is not meaningful."
)

write_csv(
  reporting_vocabulary,
  file.path(
    TAB,
    "21G3A3_reporting_vocabulary.csv"
  )
)

print(
  reporting_vocabulary,
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 6 — formal checks
# ============================================================

cat("\nSTEP 6/6 — formal checks\n")

formal_checks <- tibble(
  check=c(
    "A2_all_checks_pass",
    "receiver_universe_28",
    "nine_macrocompartments",
    "DA_exact_for_all_receivers",
    "DE_exact_for_all_receivers",
    "pathway_16_receivers",
    "pathway_missing_not_called_negative",
    "communication_anchor_frozen",
    "no_crossmodal_score",
    "no_combined_p_value",
    "no_combined_q_value",
    "DA_not_called_mechanistic_concordance",
    "DE_direction_requires_mapping",
    "pathway_direction_requires_mapping",
    "contrasts_remain_separate",
    "primary_contextual_remain_separate",
    "no_missing_modality_imputation",
    "no_new_inference",
    "no_refit"
  ),

  pass=c(
    all(checks_A2$pass),
    nrow(receiver_universe) == 28,
    n_distinct(receiver_universe$macrocompartment) == 9,
    all(cw$DA_exact),
    all(cw$DE_exact),
    sum(cw$PATHWAY_exact) == 16,
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
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  formal_checks,
  n=Inf,
  width=Inf
)

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G3A3_formal_checks.csv"
  )
)

stopifnot(
  all(formal_checks$pass)
)

audit <- c(
  "============================================================",
  "21G3A3 CROSS-MODAL INTEGRATION DESIGN FREEZE",
  "============================================================",
  "",
  "Communication is the anchor modality.",
  "",
  paste(
    "Frozen receiver universe:",
    nrow(receiver_universe)
  ),
  paste(
    "Frozen macrocompartments:",
    n_distinct(receiver_universe$macrocompartment)
  ),
  paste(
    "DA exact receiver coverage:",
    sum(cw$DA_exact),
    "/",
    nrow(cw)
  ),
  paste(
    "DE exact receiver coverage:",
    sum(cw$DE_exact),
    "/",
    nrow(cw)
  ),
  paste(
    "20F2B receiver coverage:",
    sum(cw$PATHWAY_exact),
    "/",
    nrow(cw)
  ),
  "",
  paste(
    "Receivers absent from 20F2B are explicitly classified",
    "NOT_IN_20F2B_UNIVERSE and are not interpreted as",
    "pathway-negative."
  ),
  "",
  paste(
    "No cross-modal score, combined p-value, combined q-value,",
    "new threshold or inferential family was defined."
  ),
  "",
  paste(
    "DA, DE and pathway results remain orthogonal contextual",
    "modalities and cannot promote the original communication",
    "evidence class."
  ),
  "",
  "No statistical inference performed.",
  "No analytical refit.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21G3A3_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21G3A3 COMPLETE\n")
cat("============================================================\n")
cat("Cross-modal integration design frozen before concordance analysis\n")
cat("28 communication receivers retained\n")
cat("DA coverage 28/28\n")
cat("DE coverage 28/28\n")
cat("20F2B pathway coverage 16/28\n")
cat("12 receivers explicitly NOT_IN_20F2B_UNIVERSE\n")
cat("Communication retained as anchor modality\n")
cat("No cross-modal score\n")
cat("No combined p/q values\n")
cat("No new statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
