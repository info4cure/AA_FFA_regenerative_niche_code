suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("21G3A2 — CROSS-MODAL IDENTITY CROSSWALK\n")
cat("============================================================\n")

BASE <- "PIPELINE_V2/results/disease_biology"

COMM_FILE <- file.path(
  BASE,
  "21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION",
  "21G3_BIOLOGICAL_READOUT",
  "tables",
  "21G2C_primary_A_all.csv"
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

OUT <- file.path(
  BASE,
  "21_CELL_CELL_COMMUNICATION",
  "21G_MULTICOMPARTMENT_COMMUNICATION",
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
  file.exists(COMM_FILE),
  file.exists(DA_FILE),
  file.exists(DE_FILE),
  file.exists(PATH_FILE)
)

cat("\nSTEP 1/7 — loading frozen sources\n")

comm <- read_csv(
  COMM_FILE,
  show_col_types=FALSE
)

da <- read_csv(
  DA_FILE,
  show_col_types=FALSE
)

de <- read_csv(
  DE_FILE,
  show_col_types=FALSE
)

path <- read_csv(
  PATH_FILE,
  show_col_types=FALSE
)

cat("Communication rows:", nrow(comm), "\n")
cat("DA rows:", nrow(da), "\n")
cat("DE rows:", nrow(de), "\n")
cat("Pathway rows:", nrow(path), "\n")

# ============================================================
# Helper: discover plausible identity columns
# ============================================================

identity_candidates <- c(
  "final_L2_label",
  "L2",
  "L2_label",
  "L2_identity",
  "identity",
  "cell_identity",
  "cell_type",
  "celltype",
  "receiver"
)

discover_identity_column <- function(x, source_name) {

  hits <- intersect(
    identity_candidates,
    names(x)
  )

  tibble(
    source=source_name,
    candidate_column=hits
  )
}

identity_column_audit <- bind_rows(
  discover_identity_column(comm, "COMM"),
  discover_identity_column(da,   "DA"),
  discover_identity_column(de,   "DE"),
  discover_identity_column(path, "PATHWAY")
)

cat("\n===== IDENTITY COLUMN AUDIT =====\n")

print(
  identity_column_audit,
  n=Inf,
  width=Inf
)

write_csv(
  identity_column_audit,
  file.path(
    TAB,
    "21G3A2_identity_column_audit.csv"
  )
)

# ============================================================
# Resolve only when unambiguous
# ============================================================

resolve_identity_column <- function(x, source_name) {

  hits <- intersect(
    identity_candidates,
    names(x)
  )

  if (length(hits) != 1) {

    stop(
      paste0(
        source_name,
        ": expected exactly one identity column; found ",
        length(hits),
        " [",
        paste(hits, collapse=", "),
        "]"
      )
    )
  }

  hits
}

comm_col <- resolve_identity_column(
  comm,
  "COMM"
)

da_col <- resolve_identity_column(
  da,
  "DA"
)

de_col <- resolve_identity_column(
  de,
  "DE"
)

path_col <- resolve_identity_column(
  path,
  "PATHWAY"
)

cat("\nResolved identity columns:\n")
cat("COMM   :", comm_col, "\n")
cat("DA     :", da_col, "\n")
cat("DE     :", de_col, "\n")
cat("PATHWAY:", path_col, "\n")

# ============================================================
# STEP 2 — frozen communication receiver universe
# ============================================================

cat("\nSTEP 2/7 — freezing communication receiver universe\n")

comm_ids <- comm %>%
  distinct(
    macrocompartment,
    receiver
  ) %>%
  arrange(
    macrocompartment,
    receiver
  )

stopifnot(
  nrow(comm_ids) == 28,
  n_distinct(comm_ids$macrocompartment) == 9,
  n_distinct(comm_ids$receiver) == 28
)

cat("Communication receivers:", nrow(comm_ids), "\n")
cat(
  "Macrocompartments:",
  n_distinct(comm_ids$macrocompartment),
  "\n"
)

# ============================================================
# STEP 3 — frozen identity universes
# ============================================================

cat("\nSTEP 3/7 — extracting frozen identity universes\n")

extract_ids <- function(x, col, source_name) {

  tibble(
    source=source_name,
    identity=as.character(x[[col]])
  ) %>%
    filter(
      !is.na(identity),
      identity != ""
    ) %>%
    distinct() %>%
    arrange(identity)
}

da_ids <- extract_ids(
  da,
  da_col,
  "DA"
)

de_ids <- extract_ids(
  de,
  de_col,
  "DE"
)

path_ids <- extract_ids(
  path,
  path_col,
  "PATHWAY"
)

cat("DA identities:", nrow(da_ids), "\n")
cat("DE identities:", nrow(de_ids), "\n")
cat("Pathway identities:", nrow(path_ids), "\n")

cat("\n===== DA IDENTITIES =====\n")
print(da_ids, n=Inf)

cat("\n===== DE IDENTITIES =====\n")
print(de_ids, n=Inf)

cat("\n===== PATHWAY IDENTITIES =====\n")
print(path_ids, n=Inf)

# ============================================================
# STEP 4 — exact matching only
# ============================================================

cat("\nSTEP 4/7 — exact identity matching\n")

crosswalk <- comm_ids %>%
  mutate(
    DA_exact=
      receiver %in%
        da_ids$identity,

    DE_exact=
      receiver %in%
        de_ids$identity,

    PATHWAY_exact=
      receiver %in%
        path_ids$identity,

    n_exact_modalities=
      as.integer(DA_exact) +
      as.integer(DE_exact) +
      as.integer(PATHWAY_exact),

    exact_crossmodal_status=
      case_when(
        n_exact_modalities == 3 ~
          "EXACT_ALL_THREE",

        n_exact_modalities == 2 ~
          "EXACT_TWO",

        n_exact_modalities == 1 ~
          "EXACT_ONE",

        TRUE ~
          "NO_EXACT_MATCH"
      )
  )

cat("\n===== EXACT CROSS-MODAL CROSSWALK =====\n")

print(
  crosswalk,
  n=Inf,
  width=Inf
)

cat("\n===== EXACT MATCH SUMMARY =====\n")

print(
  crosswalk %>%
    count(
      exact_crossmodal_status,
      name="n_receivers"
    ) %>%
    arrange(
      desc(n_exact_modalities <- case_when(
        exact_crossmodal_status == "EXACT_ALL_THREE" ~ 3L,
        exact_crossmodal_status == "EXACT_TWO" ~ 2L,
        exact_crossmodal_status == "EXACT_ONE" ~ 1L,
        TRUE ~ 0L
      ))
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 5 — unmatched identities for manual adjudication
# ============================================================

cat("\nSTEP 5/7 — isolating identities requiring adjudication\n")

unmatched <- crosswalk %>%
  filter(
    n_exact_modalities < 3
  )

cat("\n===== RECEIVERS REQUIRING CROSSWALK REVIEW =====\n")

print(
  unmatched,
  n=Inf,
  width=Inf
)

# Candidate lexical matches are diagnostic only.
# They MUST NOT be accepted automatically.

normalise_label <- function(x) {

  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "") %>%
    str_trim()
}

all_external_ids <- bind_rows(
  da_ids,
  de_ids,
  path_ids
) %>%
  mutate(
    normalized_identity=
      normalise_label(identity)
  )

lexical_candidates <- unmatched %>%
  mutate(
    normalized_receiver=
      normalise_label(receiver)
  ) %>%
  select(
    macrocompartment,
    receiver,
    normalized_receiver
  ) %>%
  inner_join(
    all_external_ids,
    by=c(
      "normalized_receiver" =
        "normalized_identity"
    )
  ) %>%
  distinct(
    macrocompartment,
    receiver,
    source,
    identity
  ) %>%
  arrange(
    macrocompartment,
    receiver,
    source,
    identity
  )

cat("\n===== NORMALIZED-STRING CANDIDATES =====\n")
cat("DIAGNOSTIC ONLY — NO AUTOMATIC MAPPING\n")

print(
  lexical_candidates,
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 6 — reverse audit
# ============================================================

cat("\nSTEP 6/7 — reverse identity coverage audit\n")

comm_receiver_set <-
  comm_ids$receiver

reverse_audit <- bind_rows(
  da_ids,
  de_ids,
  path_ids
) %>%
  mutate(
    exact_COMM_match=
      identity %in%
        comm_receiver_set
  ) %>%
  arrange(
    source,
    identity
  )

cat("\n===== EXTERNAL IDENTITIES NOT EXACTLY REPRESENTED IN COMM =====\n")

print(
  reverse_audit %>%
    filter(
      !exact_COMM_match
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# STEP 7 — freeze audit tables
# ============================================================

cat("\nSTEP 7/7 — writing crosswalk audit\n")

write_csv(
  comm_ids,
  file.path(
    TAB,
    "21G3A2_frozen_communication_receiver_universe.csv"
  )
)

write_csv(
  crosswalk,
  file.path(
    TAB,
    "21G3A2_exact_identity_crosswalk.csv"
  )
)

write_csv(
  unmatched,
  file.path(
    TAB,
    "21G3A2_receivers_requiring_adjudication.csv"
  )
)

write_csv(
  lexical_candidates,
  file.path(
    TAB,
    "21G3A2_normalized_string_candidates_DIAGNOSTIC_ONLY.csv"
  )
)

write_csv(
  reverse_audit,
  file.path(
    TAB,
    "21G3A2_reverse_identity_coverage.csv"
  )
)

formal_checks <- tibble(
  check=c(
    "communication_source_present",
    "DA_source_present",
    "DE_source_present",
    "pathway_source_present",
    "communication_28_receivers",
    "communication_9_macrocompartments",
    "communication_receiver_unique",
    "no_effect_recalculated",
    "no_p_value_recalculated",
    "no_BH_recalculated",
    "no_DE_refit",
    "no_DA_refit",
    "no_pathway_refit",
    "no_communication_refit",
    "no_automatic_fuzzy_mapping"
  ),

  pass=c(
    file.exists(COMM_FILE),
    file.exists(DA_FILE),
    file.exists(DE_FILE),
    file.exists(PATH_FILE),
    nrow(comm_ids) == 28,
    n_distinct(comm_ids$macrocompartment) == 9,
    n_distinct(comm_ids$receiver) == 28,
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

write_csv(
  formal_checks,
  file.path(
    TAB,
    "21G3A2_formal_checks.csv"
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  formal_checks,
  n=Inf,
  width=Inf
)

stopifnot(
  all(formal_checks$pass)
)

cat("\n============================================================\n")
cat("21G3A2 COMPLETE\n")
cat("============================================================\n")
cat("Communication receiver universe frozen at 28 identities\n")
cat("Nine communication macrocompartments retained\n")
cat("Exact cross-modal identity matching performed\n")
cat("Non-exact identities retained for manual adjudication\n")
cat("Normalized-string candidates are diagnostic only\n")
cat("No biological concordance evaluated\n")
cat("No threshold applied\n")
cat("No statistical inference performed\n")
cat("No analytical refit\n")
cat("PASS\n")
