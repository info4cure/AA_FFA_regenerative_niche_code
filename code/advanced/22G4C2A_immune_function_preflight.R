
# ============================================================
# 22G4C2A — AA vs HC IMMUNE FUNCTIONAL-STATE PREFLIGHT
#
# Goal:
# Audit frozen patient × L2 Hallmark scores for the nine
# AA/HC-eligible immune states.
#
# NO differential testing.
# NO pathway selection based on significance.
# NO re-estimation of L2 identities.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4C2A — AA vs HC IMMUNE FUNCTIONAL-STATE PREFLIGHT\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE",
  "22G4C2_IMMUNE_FUNCTION"
)

TAB <- file.path(OUT,"tables")
AUD <- file.path(OUT,"audit")

# ============================================================
# 1. FIND FROZEN HALLMARK SOURCE
# ============================================================

hallmark_hits <- list.files(
  ROOT,
  pattern="20C_patient_L2_Hallmark_GSVA_scores\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

cat("\n===== HALLMARK SOURCE SEARCH =====\n")
cat("Hits:",length(hallmark_hits),"\n")

print(
  hallmark_hits
)

stopifnot(
  length(hallmark_hits)==1
)

HALLMARK_FILE <- hallmark_hits[1]

IMMUNE_FILE <- file.path(
  dirname(OUT),
  "tables",
  "22G4C_AA_HC_eligible_immune_states.csv"
)

stopifnot(
  file.exists(HALLMARK_FILE),
  file.exists(IMMUNE_FILE)
)

# ============================================================
# 2. LOAD
# ============================================================

h <- read.csv(
  HALLMARK_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  as_tibble()

immune_def <- read.csv(
  IMMUNE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

cat("\n===== RAW HALLMARK TABLE =====\n")
cat("Rows:",nrow(h),"\n")
cat("Columns:",ncol(h),"\n")

cat("\nCOLUMN NAMES:\n")
print(
  names(h)
)

cat("\nHEAD:\n")
print(
  head(h,10),
  width=Inf
)

# ============================================================
# 3. IDENTIFY REQUIRED METADATA COLUMNS
# ============================================================

required_meta <- c(
  "patient_id",
  "disease",
  "L1",
  "L2"
)

missing_meta <- setdiff(
  required_meta,
  names(h)
)

cat("\n===== REQUIRED METADATA =====\n")

if(length(missing_meta)==0) {
  cat("All required metadata columns present.\n")
} else {
  cat(
    "Missing:",
    paste(
      missing_meta,
      collapse=", "
    ),
    "\n"
  )
}

stopifnot(
  length(missing_meta)==0
)

# ============================================================
# 4. AA + HC IMMUNE SUBSET
# ============================================================

immune_scores <- h |>
  filter(
    disease %in%
      c(
        "HC",
        "AA"
      )
  ) |>
  inner_join(
    immune_def |>
      dplyr::select(
        L1,
        L2,
        feature
      ),
    by=c(
      "L1",
      "L2"
    )
  )

cat("\n===== MATCHED AA+HC IMMUNE SCORE ROWS =====\n")

cat(
  "Rows:",
  nrow(immune_scores),
  "\n"
)

cat(
  "Patients:",
  n_distinct(
    immune_scores$patient_id
  ),
  "\n"
)

cat(
  "Immune states:",
  n_distinct(
    immune_scores$feature
  ),
  "\n"
)

# ============================================================
# 5. DETECT HALLMARK SCORE COLUMNS
# ============================================================

meta_candidates <- c(
  "patient_id",
  "disease",
  "L1",
  "L2",
  "feature",
  "n",
  "n_cells",
  "eligible"
)

numeric_cols <- names(
  immune_scores
)[
  vapply(
    immune_scores,
    is.numeric,
    logical(1)
  )
]

score_cols <- setdiff(
  numeric_cols,
  meta_candidates
)

# Prefer explicit HALLMARK columns if present.

explicit_hallmark <- grep(
  "^HALLMARK_",
  names(immune_scores),
  value=TRUE
)

if(length(explicit_hallmark)>0) {
  score_cols <- explicit_hallmark
}

cat("\n===== DETECTED SCORE COLUMNS =====\n")
cat(
  "n =",
  length(score_cols),
  "\n"
)

print(
  score_cols
)

stopifnot(
  length(score_cols)>0
)

# ============================================================
# 6. IMMUNE-STATE COVERAGE
# ============================================================

coverage <- immune_scores |>
  group_by(
    L1,
    L2,
    feature
  ) |>
  summarise(
    n_rows=n(),

    n_patients=
      n_distinct(
        patient_id
      ),

    n_HC=
      n_distinct(
        patient_id[
          disease=="HC"
        ]
      ),

    n_AA=
      n_distinct(
        patient_id[
          disease=="AA"
        ]
      ),

    .groups="drop"
  )

cat("\n===== IMMUNE STATE COVERAGE =====\n")

print(
  coverage,
  n=Inf
)

# ============================================================
# 7. SCORE COMPLETENESS
# ============================================================

score_completeness <- lapply(
  score_cols,
  function(v) {

    tibble(
      pathway=v,

      n_nonmissing=
        sum(
          is.finite(
            immune_scores[[v]]
          )
        ),

      n_missing=
        sum(
          !is.finite(
            immune_scores[[v]]
          )
        ),

      n_unique=
        n_distinct(
          immune_scores[[v]][
            is.finite(
              immune_scores[[v]]
            )
          ]
        )
    )
  }
) |>
  bind_rows()

cat("\n===== SCORE COMPLETENESS =====\n")

print(
  score_completeness,
  n=Inf
)

# ============================================================
# 8. PRE-SPECIFIED IMMUNE HALLMARK PANEL
# ============================================================

candidate_panel <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_APOPTOSIS"
)

candidate_presence <- tibble(
  pathway=
    candidate_panel,

  available=
    candidate_panel %in%
      score_cols
)

cat("\n===== PRE-SPECIFIED IMMUNE HALLMARK PANEL =====\n")

print(
  candidate_presence,
  n=Inf
)

# ============================================================
# 9. CYTOTOXICITY SOURCE AUDIT
# ============================================================

# Hallmark does not contain a dedicated cytotoxicity gene set.
# We therefore only record the intended frozen custom signature
# here. It will NOT be calculated until its expression source
# has been explicitly audited.

cytotoxicity_genes <- c(
  "NKG7",
  "GNLY",
  "PRF1",
  "GZMB",
  "GZMH",
  "GZMK",
  "CTSW",
  "CCL5"
)

cytotoxicity_manifest <- tibble(
  programme="CYTOTOXICITY",
  gene=cytotoxicity_genes
)

cat("\n===== CUSTOM CYTOTOXICITY SIGNATURE =====\n")

print(
  cytotoxicity_manifest,
  n=Inf
)

# ============================================================
# 10. SAVE
# ============================================================

write.csv(
  coverage,
  file.path(
    TAB,
    "22G4C2A_immune_state_Hallmark_coverage.csv"
  ),
  row.names=FALSE
)

write.csv(
  score_completeness,
  file.path(
    TAB,
    "22G4C2A_Hallmark_score_completeness.csv"
  ),
  row.names=FALSE
)

write.csv(
  candidate_presence,
  file.path(
    TAB,
    "22G4C2A_prespecified_Hallmark_panel.csv"
  ),
  row.names=FALSE
)

write.csv(
  cytotoxicity_manifest,
  file.path(
    TAB,
    "22G4C2A_cytotoxicity_signature_manifest.csv"
  ),
  row.names=FALSE
)

write.csv(
  immune_scores,
  file.path(
    TAB,
    "22G4C2A_AA_HC_immune_Hallmark_source.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 11. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "Unique Hallmark source",
    "11 AA+HC patients represented",
    "9 immune states represented",
    "Hallmark score columns detected",
    "IFNG Hallmark available",
    "IFNA Hallmark available",
    "TNFA/NFKB Hallmark available",
    "IL2/STAT5 Hallmark available",
    "Inflammatory response Hallmark available"
  ),

  pass=c(
    length(hallmark_hits)==1,

    n_distinct(
      immune_scores$patient_id
    )==11,

    n_distinct(
      immune_scores$feature
    )==9,

    length(
      score_cols
    )>0,

    "HALLMARK_INTERFERON_GAMMA_RESPONSE" %in%
      score_cols,

    "HALLMARK_INTERFERON_ALPHA_RESPONSE" %in%
      score_cols,

    "HALLMARK_TNFA_SIGNALING_VIA_NFKB" %in%
      score_cols,

    "HALLMARK_IL2_STAT5_SIGNALING" %in%
      score_cols,

    "HALLMARK_INFLAMMATORY_RESPONSE" %in%
      score_cols
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
    "22G4C2A_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 12. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C2A — AA vs HC IMMUNE FUNCTIONAL-STATE PREFLIGHT",
  "============================================================",
  "",
  paste0(
    "Hallmark source: ",
    HALLMARK_FILE
  ),
  paste0(
    "AA+HC patients represented: ",
    n_distinct(
      immune_scores$patient_id
    )
  ),
  paste0(
    "Eligible immune states represented: ",
    n_distinct(
      immune_scores$feature
    )
  ),
  paste0(
    "Detected Hallmark scores: ",
    length(
      score_cols
    )
  ),
  "",
  "Pre-specified functional hypotheses:",
  "- IFN-gamma response",
  "- IFN-alpha response",
  "- TNFA/NFKB signaling",
  "- IL2/STAT5 signaling",
  "- inflammatory response",
  "- complement",
  "- apoptosis/stress",
  "",
  "Custom cytotoxicity signature defined but NOT calculated.",
  "Expression source must be audited before custom scoring.",
  "",
  "No differential testing performed.",
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
    "22G4C2A_audit.txt"
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
    "22G4C2A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C2A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

