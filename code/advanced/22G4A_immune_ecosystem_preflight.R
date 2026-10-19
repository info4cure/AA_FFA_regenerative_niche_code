
# ============================================================
# 22G4A — IMMUNE-CENTERED ECOSYSTEM PREFLIGHT
#
# Purpose:
# Define the frozen immune L2 universe and determine which
# immune states have sufficient patient-level support for
# downstream coupling to:
#
#   1. Bulge trajectory / functional state
#   2. Follicular epithelium
#   3. Superficial epithelium
#   4. Mesenchymal / vascular niche
#
# No inferential modelling is performed here.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4A — IMMUNE-CENTERED ECOSYSTEM PREFLIGHT\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

L2_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4A_PREFLIGHT"
)

TAB <- file.path(
  OUT,
  "tables"
)

AUD <- file.path(
  OUT,
  "audit"
)

dir.create(
  TAB,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  AUD,
  recursive=TRUE,
  showWarnings=FALSE
)

stopifnot(
  file.exists(
    L2_FILE
  )
)

# ============================================================
# 2. LOAD FROZEN L2 COMPOSITION
# ============================================================

x <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

IMMUNE_L1 <- c(
  "T_NK",
  "Myeloid_DC"
)

immune <- x |>
  filter(
    L1 %in% IMMUNE_L1
  )

cat("\n===== IMMUNE L1 =====\n")

print(
  immune |>
    count(
      L1,
      name="rows"
    )
)

# ============================================================
# 3. GLOBAL SUPPORT
# ============================================================

global_support <- immune |>
  group_by(
    L1,
    L2
  ) |>
  summarise(
    n_patients=
      n_distinct(
        patient_id
      ),

    total_cells=
      sum(
        n,
        na.rm=TRUE
      ),

    median_n=
      median(
        n,
        na.rm=TRUE
      ),

    min_n=
      min(
        n,
        na.rm=TRUE
      ),

    zero_patients=
      sum(
        n==0,
        na.rm=TRUE
      ),

    .groups="drop"
  )

# ============================================================
# 4. DISEASE-SPECIFIC SUPPORT
# ============================================================

disease_support <- immune |>
  mutate(
    positive=
      n > 0
  ) |>
  group_by(
    L1,
    L2,
    disease
  ) |>
  summarise(
    n_positive=
      sum(
        positive,
        na.rm=TRUE
      ),

    .groups="drop"
  ) |>
  pivot_wider(
    names_from=disease,
    values_from=n_positive,
    values_fill=0,
    names_prefix="positive_"
  )

support <- global_support |>
  left_join(
    disease_support,
    by=c(
      "L1",
      "L2"
    )
  ) |>
  mutate(
    feature=paste(
      L1,
      L2,
      sep=" :: "
    ),

    eligible_immune=
      total_cells >= 40 &
      (n_patients-zero_patients) >= 12 &
      positive_HC >= 3 &
      positive_FFA >= 3 &
      positive_AA >= 3,

    exclusion_reason=case_when(
      eligible_immune
        ~ "INCLUDE",

      total_cells < 40
        ~ "<40 total cells",

      (n_patients-zero_patients) < 12
        ~ "<12 positive patients",

      positive_HC < 3
        ~ "<3 HC positive",

      positive_FFA < 3
        ~ "<3 FFA positive",

      positive_AA < 3
        ~ "<3 AA positive",

      TRUE
        ~ "OTHER"
    )
  ) |>
  arrange(
    L1,
    L2
  )

# ============================================================
# 5. DISPLAY
# ============================================================

cat("\n===== FROZEN IMMUNE L2 =====\n")

print(
  support,
  n=Inf
)

cat("\n===== ELIGIBLE IMMUNE STATES =====\n")

eligible <- support |>
  filter(
    eligible_immune
  )

print(
  eligible,
  n=Inf
)

cat(
  "\nRaw immune states:",
  nrow(support),
  "\n"
)

cat(
  "Eligible immune states:",
  nrow(eligible),
  "\n"
)

# ============================================================
# 6. CELL COUNTS BY IMMUNE COMPARTMENT
# ============================================================

L1_counts <- immune |>
  group_by(
    L1
  ) |>
  summarise(
    total_cells=
      sum(
        n,
        na.rm=TRUE
      ),

    n_patients=
      n_distinct(
        patient_id
      ),

    .groups="drop"
  )

cat("\n===== IMMUNE L1 SUPPORT =====\n")

print(
  L1_counts,
  n=Inf
)

# ============================================================
# 7. PATIENT × IMMUNE L2 SOURCE
# ============================================================

patient_source <- immune |>
  dplyr::select(
    patient_id,
    disease,
    L1,
    L2,
    n,
    denominator_L1,
    eligible,
    frequency,
    percent
  )

# ============================================================
# 8. SAVE
# ============================================================

write.csv(
  support,
  file.path(
    TAB,
    "22G4A_immune_L2_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  eligible,
  file.path(
    TAB,
    "22G4A_eligible_immune_states.csv"
  ),
  row.names=FALSE
)

write.csv(
  L1_counts,
  file.path(
    TAB,
    "22G4A_immune_L1_support.csv"
  ),
  row.names=FALSE
)

write.csv(
  patient_source,
  file.path(
    TAB,
    "22G4A_patient_immune_L2_source.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 9. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "16 patients represented",
    "T_NK represented",
    "Myeloid_DC represented",
    "Immune L2 universe non-empty",
    "At least one eligible immune state",
    "HC represented",
    "FFA represented",
    "AA represented"
  ),

  pass=c(
    n_distinct(
      immune$patient_id
    )==16,

    "T_NK" %in%
      immune$L1,

    "Myeloid_DC" %in%
      immune$L1,

    nrow(support)>0,

    nrow(eligible)>0,

    "HC" %in%
      immune$disease,

    "FFA" %in%
      immune$disease,

    "AA" %in%
      immune$disease
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  checks,
  n=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G4A_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 10. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4A — IMMUNE-CENTERED ECOSYSTEM PREFLIGHT",
  "============================================================",
  "",
  paste0(
    "Patients: ",
    n_distinct(
      immune$patient_id
    )
  ),
  paste0(
    "Raw immune L2 states: ",
    nrow(support)
  ),
  paste0(
    "Eligible immune L2 states: ",
    nrow(eligible)
  ),
  "",
  "Eligibility:",
  ">=40 total cells",
  ">=12 positive patients",
  ">=3 positive HC",
  ">=3 positive FFA",
  ">=3 positive AA",
  "",
  "Planned downstream ecosystem layers:",
  "1. Bulge trajectory/function",
  "2. Follicular epithelium",
  "3. Superficial epithelium",
  "4. Mesenchymal/vascular niche",
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(checks$pass),
    "/",
    nrow(checks),
    " PASS"
  ),
  "",
  if(all(checks$pass)) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUD,
    "22G4A_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

if(!all(checks$pass)) {
  stop(
    "22G4A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

