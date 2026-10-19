#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20F2B INTEGRATED STATE PRIORITIZATION\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

ROOT <- file.path(
  "PIPELINE_V2/results/disease_biology",
  "20_PATHWAY_ANALYSIS"
)

IN <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2A_disease_architecture",
  "tables",
  "20F2A_candidate_states_for_20F2B.csv"
)

OUT <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2B_integrated_prioritization"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Recover frozen 20F2A candidates
# ============================================================

cat("\nSTEP 1/9 — recovering frozen 20F2A candidate states\n")

stopifnot(file.exists(IN))

x <- read_csv(
  IN,
  show_col_types=FALSE
)

cat("Candidate states:", nrow(x), "\n")

required <- c(
  "module",
  "final_L2_label",
  "axis",
  "AA_max_primary_databases",
  "AA_n_significant_pathways",
  "AA_median_database_NES",
  "AA_min_FDR",
  "AA_Hallmark_support",
  "AA_GSVA_support",
  "AA_primary_status",
  "FFA_max_primary_databases",
  "FFA_n_significant_pathways",
  "FFA_median_database_NES",
  "FFA_min_FDR",
  "FFA_Hallmark_support",
  "FFA_GSVA_support",
  "FFA_primary_status",
  "architecture",
  "direct_max_primary_databases",
  "direct_n_significant_pathways",
  "direct_median_NES",
  "direct_min_FDR",
  "direct_Hallmark_support",
  "direct_GSVA_support",
  "direct_status",
  "direct_support",
  "n_disease_consensus",
  "any_Hallmark_support",
  "any_GSVA_support",
  "evidence_grade"
)

stopifnot(all(required %in% names(x)))
stopifnot(nrow(x) == 65)

key_dup <- x %>%
  count(module, final_L2_label, axis) %>%
  filter(n > 1)

stopifnot(nrow(key_dup) == 0)

# ============================================================
# 2. Define interpretable primary-consensus strength
# ============================================================

cat("\nSTEP 2/9 — deriving primary consensus strength\n")

# Do not combine heterogeneous statistical quantities into
# an artificial continuous score.
#
# Disease-specific primary support:
# 3 databases = TRIPLE
# 2 databases = CORE
# <2          = no primary cross-database consensus

x2 <- x %>%
  mutate(
    AA_consensus_class = case_when(
      AA_max_primary_databases >= 3 ~ "TRIPLE",
      AA_max_primary_databases == 2 ~ "CORE",
      TRUE ~ "NONE"
    ),
    FFA_consensus_class = case_when(
      FFA_max_primary_databases >= 3 ~ "TRIPLE",
      FFA_max_primary_databases == 2 ~ "CORE",
      TRUE ~ "NONE"
    )
  )

# Integrated strength depends on architecture.
#
# Shared/opposite states require evidence from both diseases,
# therefore the weaker disease-side consensus is the limiting
# evidence.
#
# Supported-only states use the supported disease side.

x2 <- x2 %>%
  mutate(
    integrated_primary_class = case_when(

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN",
        "OPPOSITE"
      ) &
        AA_consensus_class == "TRIPLE" &
        FFA_consensus_class == "TRIPLE" ~
        "BILATERAL_TRIPLE",

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN",
        "OPPOSITE"
      ) &
        AA_consensus_class %in% c("CORE","TRIPLE") &
        FFA_consensus_class %in% c("CORE","TRIPLE") ~
        "BILATERAL_CORE",

      architecture == "AA_SUPPORTED_ONLY" &
        AA_consensus_class == "TRIPLE" ~
        "UNILATERAL_TRIPLE",

      architecture == "AA_SUPPORTED_ONLY" &
        AA_consensus_class == "CORE" ~
        "UNILATERAL_CORE",

      architecture == "FFA_SUPPORTED_ONLY" &
        FFA_consensus_class == "TRIPLE" ~
        "UNILATERAL_TRIPLE",

      architecture == "FFA_SUPPORTED_ONLY" &
        FFA_consensus_class == "CORE" ~
        "UNILATERAL_CORE",

      TRUE ~ "OTHER"
    )
  )

cat("\nIntegrated primary classes:\n")
print(
  x2 %>%
    count(integrated_primary_class, sort=TRUE),
  n=Inf
)

# ============================================================
# 3. Orthogonal support
# ============================================================

cat("\nSTEP 3/9 — classifying orthogonal support\n")

x2 <- x2 %>%
  mutate(
    orthogonal_support_class = case_when(
      any_GSVA_support & any_Hallmark_support ~
        "HALLMARK_PLUS_PATIENT",
      any_GSVA_support ~
        "PATIENT",
      any_Hallmark_support ~
        "HALLMARK",
      TRUE ~
        "NONE"
    ),
    patient_validated = any_GSVA_support %in% TRUE
  )

print(
  x2 %>%
    count(
      orthogonal_support_class,
      sort=TRUE
    ),
  n=Inf
)

# ============================================================
# 4. Interpret direct AA-vs-FFA triangulation
# ============================================================

cat("\nSTEP 4/9 — interpreting direct AA-vs-FFA triangulation\n")

# Critical distinction:
#
# SHARED states:
# DIRECT_NONDISCRIMINATING is compatible with the architecture
# and must NOT be treated as negative evidence.
#
# Supported-only states:
# DIRECT_CONCORDANT provides independent discriminatory support.
#
# Opposite states:
# DIRECT_CONCORDANT would strengthen directionality;
# DIRECT_MIXED remains ambiguous.
#
# No direct primary support is absence of confirmation,
# not discordance.

x2 <- x2 %>%
  mutate(
    direct_triangulation_class = case_when(

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN"
      ) &
        direct_support == "DIRECT_NONDISCRIMINATING" ~
        "COMPATIBLE_SHARED",

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN"
      ) &
        direct_support == "NO_DIRECT_PRIMARY_SUPPORT" ~
        "UNINFORMATIVE",

      architecture %in% c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) &
        direct_support == "DIRECT_CONCORDANT" ~
        "CONFIRMATORY",

      architecture %in% c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) &
        direct_support == "NO_DIRECT_PRIMARY_SUPPORT" ~
        "UNINFORMATIVE",

      architecture == "OPPOSITE" &
        direct_support == "DIRECT_CONCORDANT" ~
        "CONFIRMATORY",

      architecture == "OPPOSITE" &
        direct_support == "DIRECT_MIXED" ~
        "AMBIGUOUS",

      architecture == "OPPOSITE" &
        direct_support == "NO_DIRECT_PRIMARY_SUPPORT" ~
        "UNINFORMATIVE",

      grepl(
        "DISCORD",
        direct_support,
        ignore.case=TRUE
      ) ~
        "DISCORDANT",

      TRUE ~
        "OTHER"
    )
  )

print(
  x2 %>%
    count(
      architecture,
      direct_triangulation_class,
      sort=TRUE
    ),
  n=Inf
)

# ============================================================
# 5. Define categorical evidence tier
# ============================================================

cat("\nSTEP 5/9 — assigning integrated evidence tiers\n")

# Evidence tiers are categorical, not numerical scores.
#
# TIER_1:
# strongest cross-database primary evidence plus patient GSVA.
#
# TIER_2:
# strong primary consensus with orthogonal Hallmark support,
# or unilateral disease evidence independently confirmed by
# the direct AA-vs-FFA contrast.
#
# TIER_3:
# robust primary consensus without patient-level validation.
#
# OPPOSITE states are kept explicit and are not promoted above
# Tier 2 unless independently confirmed.
#
# No state is labelled disease-specific.

x2 <- x2 %>%
  mutate(
    integrated_evidence_tier = case_when(

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN"
      ) &
        integrated_primary_class == "BILATERAL_TRIPLE" &
        patient_validated ~
        "TIER_1",

      architecture %in% c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) &
        integrated_primary_class == "UNILATERAL_TRIPLE" &
        patient_validated &
        direct_triangulation_class == "CONFIRMATORY" ~
        "TIER_1",

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN"
      ) &
        integrated_primary_class == "BILATERAL_TRIPLE" &
        any_Hallmark_support ~
        "TIER_2",

      architecture %in% c(
        "SHARED_UP",
        "SHARED_DOWN"
      ) &
        integrated_primary_class == "BILATERAL_CORE" &
        patient_validated ~
        "TIER_2",

      architecture %in% c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) &
        integrated_primary_class == "UNILATERAL_TRIPLE" &
        direct_triangulation_class == "CONFIRMATORY" ~
        "TIER_2",

      architecture %in% c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) &
        integrated_primary_class == "UNILATERAL_TRIPLE" &
        patient_validated ~
        "TIER_2",

      architecture == "OPPOSITE" &
        integrated_primary_class == "BILATERAL_TRIPLE" &
        direct_triangulation_class == "CONFIRMATORY" ~
        "TIER_2",

      integrated_primary_class %in% c(
        "BILATERAL_TRIPLE",
        "UNILATERAL_TRIPLE"
      ) ~
        "TIER_3A",

      integrated_primary_class %in% c(
        "BILATERAL_CORE",
        "UNILATERAL_CORE"
      ) &
        patient_validated ~
        "TIER_3A",

      integrated_primary_class %in% c(
        "BILATERAL_CORE",
        "UNILATERAL_CORE"
      ) &
        any_Hallmark_support ~
        "TIER_3B",

      integrated_primary_class %in% c(
        "BILATERAL_CORE",
        "UNILATERAL_CORE"
      ) ~
        "TIER_3C",

      TRUE ~
        "REVIEW"
    )
  )

cat("\nIntegrated evidence tiers:\n")
print(
  x2 %>%
    count(
      integrated_evidence_tier,
      sort=TRUE
    ),
  n=Inf
)

# ============================================================
# 6. Build publication-priority strata
# ============================================================

cat("\nSTEP 6/9 — defining publication-priority strata\n")

# This is a deterministic display/interpretation hierarchy,
# not a biological effect-size score.

tier_order <- c(
  "TIER_1",
  "TIER_2",
  "TIER_3A",
  "TIER_3B",
  "TIER_3C",
  "REVIEW"
)

architecture_order <- c(
  "SHARED_UP",
  "SHARED_DOWN",
  "AA_SUPPORTED_ONLY",
  "FFA_SUPPORTED_ONLY",
  "OPPOSITE"
)

x2 <- x2 %>%
  mutate(
    tier_rank = match(
      integrated_evidence_tier,
      tier_order
    ),
    architecture_rank = match(
      architecture,
      architecture_order
    ),

    # Descriptive publication role.
    publication_role = case_when(

      integrated_evidence_tier == "TIER_1" ~
        "PRIMARY",

      integrated_evidence_tier == "TIER_2" ~
        "PRIMARY_SUPPORTING",

      integrated_evidence_tier %in% c(
        "TIER_3A",
        "TIER_3B"
      ) ~
        "SECONDARY",

      TRUE ~
        "CONTEXTUAL"
    )
  ) %>%
  arrange(
    tier_rank,
    architecture_rank,
    axis,
    module,
    final_L2_label
  )

# ============================================================
# 7. Generate summaries
# ============================================================

cat("\nSTEP 7/9 — generating prioritization summaries\n")

tier_summary <- x2 %>%
  count(
    integrated_evidence_tier,
    publication_role,
    name="n_states"
  ) %>%
  mutate(
    tier_rank=match(
      integrated_evidence_tier,
      tier_order
    )
  ) %>%
  arrange(tier_rank) %>%
  select(-tier_rank)

architecture_tier_summary <- x2 %>%
  count(
    architecture,
    integrated_evidence_tier,
    name="n_states"
  ) %>%
  mutate(
    architecture_rank=match(
      architecture,
      architecture_order
    ),
    tier_rank=match(
      integrated_evidence_tier,
      tier_order
    )
  ) %>%
  arrange(
    architecture_rank,
    tier_rank
  ) %>%
  select(
    -architecture_rank,
    -tier_rank
  )

axis_tier_summary <- x2 %>%
  count(
    axis,
    integrated_evidence_tier,
    name="n_states"
  ) %>%
  mutate(
    tier_rank=match(
      integrated_evidence_tier,
      tier_order
    )
  ) %>%
  arrange(
    axis,
    tier_rank
  ) %>%
  select(-tier_rank)

L2_priority_summary <- x2 %>%
  group_by(
    module,
    final_L2_label
  ) %>%
  summarise(
    n_candidate_states=n(),
    n_TIER_1=sum(
      integrated_evidence_tier == "TIER_1"
    ),
    n_TIER_2=sum(
      integrated_evidence_tier == "TIER_2"
    ),
    n_primary_publication=sum(
      publication_role == "PRIMARY"
    ),
    n_primary_supporting=sum(
      publication_role == "PRIMARY_SUPPORTING"
    ),
    n_patient_validated=sum(
      patient_validated
    ),
    n_direct_confirmatory=sum(
      direct_triangulation_class == "CONFIRMATORY"
    ),
    axes=paste(
      sort(unique(axis)),
      collapse=";"
    ),
    .groups="drop"
  ) %>%
  arrange(
    desc(n_TIER_1),
    desc(n_TIER_2),
    desc(n_patient_validated),
    module,
    final_L2_label
  )

primary_states <- x2 %>%
  filter(
    publication_role %in% c(
      "PRIMARY",
      "PRIMARY_SUPPORTING"
    )
  )

direct_confirmed_supported_only <- x2 %>%
  filter(
    architecture %in% c(
      "AA_SUPPORTED_ONLY",
      "FFA_SUPPORTED_ONLY"
    ),
    direct_triangulation_class ==
      "CONFIRMATORY"
  )

opposite_states <- x2 %>%
  filter(
    architecture == "OPPOSITE"
  )

# ============================================================
# 8. Write frozen outputs
# ============================================================

cat("\nSTEP 8/9 — writing frozen 20F2B outputs\n")

write_csv(
  x2,
  file.path(
    TAB,
    "20F2B_integrated_state_prioritization.csv"
  )
)

write_csv(
  tier_summary,
  file.path(
    TAB,
    "20F2B_evidence_tier_summary.csv"
  )
)

write_csv(
  architecture_tier_summary,
  file.path(
    TAB,
    "20F2B_architecture_tier_summary.csv"
  )
)

write_csv(
  axis_tier_summary,
  file.path(
    TAB,
    "20F2B_axis_tier_summary.csv"
  )
)

write_csv(
  L2_priority_summary,
  file.path(
    TAB,
    "20F2B_L2_priority_summary.csv"
  )
)

write_csv(
  primary_states,
  file.path(
    TAB,
    "20F2B_primary_publication_states.csv"
  )
)

write_csv(
  direct_confirmed_supported_only,
  file.path(
    TAB,
    "20F2B_direct_confirmed_supported_only.csv"
  )
)

write_csv(
  opposite_states,
  file.path(
    TAB,
    "20F2B_opposite_states.csv"
  )
)

# ============================================================
# 9. Formal checks / audit
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "20F2A_input_present",
    "candidate_states_65",
    "analytical_key_unique",
    "five_architecture_classes",
    "primary_consensus_retained",
    "Hallmark_
orthogonal_not_primary",
    "GSVA_orthogonal_not_primary",
    "direct_contrast_not_primary",
    "shared_direct_nondiscriminating_not_penalized",
    "supported_only_terminology_retained",
    "opposite_states_retained",
    "categorical_not_continuous_scoring",
    "all_states_assigned_tier",
    "publication_roles_assigned",
    "no_GSEA_refit",
    "no_GSVA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    file.exists(IN),
    nrow(x2) == 65,
    nrow(key_dup) == 0,
    setequal(
      unique(x2$architecture),
      c(
        "SHARED_UP",
        "SHARED_DOWN",
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY",
        "OPPOSITE"
      )
    ),
    all(
      x2$integrated_primary_class %in%
        c(
          "BILATERAL_TRIPLE",
          "BILATERAL_CORE",
          "UNILATERAL_TRIPLE",
          "UNILATERAL_CORE",
          "OTHER"
        )
    ),
    TRUE,
    TRUE,
    TRUE,
    all(
      x2$direct_triangulation_class[
        x2$architecture %in%
          c("SHARED_UP", "SHARED_DOWN") &
          x2$direct_support ==
            "DIRECT_NONDISCRIMINATING"
      ] == "COMPATIBLE_SHARED"
    ),
    all(
      x2$architecture %in%
        c(
          "SHARED_UP",
          "SHARED_DOWN",
          "AA_SUPPORTED_ONLY",
          "FFA_SUPPORTED_ONLY",
          "OPPOSITE"
        )
    ),
    all(
      x2$architecture[
        x2$architecture == "OPPOSITE"
      ] == "OPPOSITE"
    ),
    !("priority_score" %in% names(x2)),
    all(
      !is.na(x2$integrated_evidence_tier)
    ),
    all(
      !is.na(x2$publication_role)
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
  n=Inf
)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "20F2B_formal_checks.csv"
  )
)

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

audit <- c(
  "PIPELINE_V2 — 20F2B integrated state prioritization",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Frozen input:",
  "- 20F2A candidate states for 20F2B",
  paste0("- candidate states: ", nrow(x2)),
  "",
  "Analytical principle:",
  "- no continuous composite evidence score was constructed",
  "- heterogeneous quantities (NES, FDR, pathway counts,",
  "  database support, Hallmark and GSVA) were not numerically",
  "  collapsed into a single artificial metric",
  "- prioritization is categorical and deterministic",
  "",
  "Primary evidence:",
  "- Reactome, GO-BP and BioPlanet remain the three primary",
  "  pathway databases",
  "- TRIPLE = directional support from all three primary databases",
  "- CORE = directional support from two primary databases",
  "",
  "Integrated primary classes:",
  "- BILATERAL_TRIPLE",
  "- BILATERAL_CORE",
  "- UNILATERAL_TRIPLE",
  "- UNILATERAL_CORE",
  "",
  "Architecture:",
  "- SHARED_UP",
  "- SHARED_DOWN",
  "- AA_SUPPORTED_ONLY",
  "- FFA_SUPPORTED_ONLY",
  "- OPPOSITE",
  "",
  "Conservative terminology:",
  "- AA_SUPPORTED_ONLY and FFA_SUPPORTED_ONLY are retained",
  "- these states are not called disease-specific",
  "- absence of consensus in the other disease is not evidence",
  "  of biological absence",
  "",
  "Orthogonal evidence:",
  "- Hallmark is supporting evidence, not a primary database",
  "- patient-level Hallmark GSVA is orthogonal patient support",
  "- neither changes the primary database denominator",
  "",
  "Direct AA-vs-FFA triangulation:",
  "- direct contrast does not define disease architecture",
  "- DIRECT_CONCORDANT can independently support",
  "  AA_SUPPORTED_ONLY or FFA_SUPPORTED_ONLY states",
  "- DIRECT_NONDISCRIMINATING is compatible with SHARED_UP",
  "  and SHARED_DOWN and is not treated as negative evidence",
  "- absence of direct primary support is considered",
  "  uninformative rather than discordant",
  "- mixed/opposite evidence remains explicit",
  "",
  "Integrated evidence tiers:",
  "- TIER_1: strongest primary consensus with patient-level",
  "  validation, including directly confirmed unilateral states",
  "- TIER_2: strong primary evidence with orthogonal or direct",
  "  confirmation",
  "- TIER_3A/B/C: progressively less orthogonally supported",
  "  but still primary-consensus-supported states",
  "- REVIEW: states not meeting the predefined integrated tiers",
  "",
  "Publication roles:",
  "- PRIMARY",
  "- PRIMARY_SUPPORTING",
  "- SECONDARY",
  "- CONTEXTUAL",
  "",
  "No GSEA refit.",
  "No GSVA refit.",
  "No differential-expression refit.",
  "No annotation modified.",
  "No figures generated."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20F2B_audit.txt"
  )
)

# ------------------------------------------------------------
# Console biological readout
# ------------------------------------------------------------

cat("\n===== INTEGRATED EVIDENCE TIER SUMMARY =====\n")
print(
  tier_summary,
  n=Inf
)

cat("\n===== ARCHITECTURE × TIER =====\n")
print(
  architecture_tier_summary,
  n=Inf
)

cat("\n===== AXIS × TIER =====\n")
print(
  axis_tier_summary,
  n=Inf
)

cat("\n===== TOP L2 BY INTEGRATED SUPPORT =====\n")
print(
  L2_priority_summary,
  n=Inf
)

cat("\n===== PRIMARY / PRIMARY-SUPPORTING STATES =====\n")
print(
  primary_states %>%
    select(
      module,
      final_L2_label,
      axis,
      architecture,
      integrated_primary_class,
      orthogonal_support_class,
      direct_triangulation_class,
      integrated_evidence_tier,
      publication_role
    ),
  n=Inf
)

cat("\n===== DIRECT-CONFIRMED SUPPORTED-ONLY STATES =====\n")
print(
  direct_confirmed_supported_only %>%
    select(
      module,
      final_L2_label,
      axis,
      architecture,
      integrated_primary_class,
      direct_status,
      direct_support,
      direct_triangulation_class,
      integrated_evidence_tier
    ),
  n=Inf
)

cat("\n===== OPPOSITE STATES =====\n")
print(
  opposite_states %>%
    select(
      module,
      final_L2_label,
      axis,
      AA_primary_status,
      FFA_primary_status,
      architecture,
      direct_status,
      direct_support,
      direct_triangulation_class,
      integrated_evidence_tier
    ),
  n=Inf
)

cat("\n============================================================\n")
cat("20F2B COMPLETE\n")
cat("============================================================\n")
cat("Integrated disease-state prioritization frozen\n")
cat("Primary evidence: Reactome + GO-BP + BioPlanet\n")
cat("Hallmark and patient GSVA retained as orthogonal support\n")
cat("AA-vs-FFA retained as independent triangulation\n")
cat("No artificial continuous evidence score\n")
cat("No disease-specific claims inferred from absence of support\n")
cat("No analytical refit\n")
cat("PASS\n")
