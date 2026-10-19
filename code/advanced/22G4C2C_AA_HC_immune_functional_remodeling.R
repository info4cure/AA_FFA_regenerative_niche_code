
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4C2C — AA vs HC IMMUNE FUNCTIONAL REMODELING\n")
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
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

SCORE_FILE <- file.path(
  TAB,
  "22G4C2B_targeted_immune_functional_scores.csv"
)

stopifnot(
  file.exists(SCORE_FILE)
)

x <- read.csv(
  SCORE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

# ============================================================
# 1. DEFINITIONS
# ============================================================

PRIMARY_POPULATIONS <- c(
  "CD4 T",
  "CD8 T",
  "NK",
  "Macrophage/monocyte"
)

PRIMARY_PROGRAMMES <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_APOPTOSIS",
  "CYTOTOXICITY"
)

# ============================================================
# 2. PRIMARY AA vs HC TESTS
# ============================================================

cat("\n===== PRIMARY AA vs HC TESTS =====\n")

primary_list <- list()

k <- 1L

for(pop in PRIMARY_POPULATIONS) {

  programmes <- x |>
    filter(
      final_L2_label==pop
    ) |>
    distinct(
      programme
    ) |>
    pull(
      programme
    )

  for(prog in programmes) {

    d <- x |>
      filter(
        final_L2_label==pop,
        programme==prog
      )

    n_HC <- sum(
      d$disease=="HC"
    )

    n_AA <- sum(
      d$disease=="AA"
    )

    med_HC <- median(
      d$score[
        d$disease=="HC"
      ],
      na.rm=TRUE
    )

    med_AA <- median(
      d$score[
        d$disease=="AA"
      ],
      na.rm=TRUE
    )

    mean_HC <- mean(
      d$score[
        d$disease=="HC"
      ],
      na.rm=TRUE
    )

    mean_AA <- mean(
      d$score[
        d$disease=="AA"
      ],
      na.rm=TRUE
    )

    wt <- suppressWarnings(
      wilcox.test(
        score ~ disease,
        data=d,
        exact=FALSE
      )
    )

    primary_list[[k]] <- tibble(
      final_L2_label=pop,
      programme=prog,
      n_HC=n_HC,
      n_AA=n_AA,
      median_HC=med_HC,
      median_AA=med_AA,
      delta_median_AA_HC=
        med_AA-med_HC,
      mean_HC=mean_HC,
      mean_AA=mean_AA,
      delta_mean_AA_HC=
        mean_AA-mean_HC,
      p=wt$p.value
    )

    k <- k+1L
  }
}

primary <- bind_rows(
  primary_list
) |>
  group_by(
    final_L2_label
  ) |>
  mutate(
    FDR_within_population=
      p.adjust(
        p,
        method="BH"
      )
  ) |>
  ungroup() |>
  mutate(
    FDR_global=
      p.adjust(
        p,
        method="BH"
      ),

    direction=case_when(
      delta_median_AA_HC > 0
        ~ "AA_HIGHER",
      delta_median_AA_HC < 0
        ~ "AA_LOWER",
      TRUE
        ~ "NO_CHANGE"
    )
  )

# ============================================================
# 3. CELL-SUPPORT SENSITIVITY
# ============================================================

cat("\n===== CELL-SUPPORT SENSITIVITY =====\n")

thresholds <- c(
  ALL=0,
  GE5=5,
  GE10=10,
  GE20=20
)

sens_list <- list()

k <- 1L

for(nm in names(thresholds)) {

  threshold <- thresholds[[nm]]

  for(pop in PRIMARY_POPULATIONS) {

    programmes <- x |>
      filter(
        final_L2_label==pop
      ) |>
      distinct(
        programme
      ) |>
      pull(
        programme
      )

    for(prog in programmes) {

      d <- x |>
        filter(
          final_L2_label==pop,
          programme==prog,
          n_cells >= threshold
        )

      n_HC <- sum(
        d$disease=="HC"
      )

      n_AA <- sum(
        d$disease=="AA"
      )

      if(
        n_HC < 2 ||
        n_AA < 2
      ) {
        next
      }

      wt <- suppressWarnings(
        wilcox.test(
          score ~ disease,
          data=d,
          exact=FALSE
        )
      )

      med_HC <- median(
        d$score[
          d$disease=="HC"
        ],
        na.rm=TRUE
      )

      med_AA <- median(
        d$score[
          d$disease=="AA"
        ],
        na.rm=TRUE
      )

      sens_list[[k]] <- tibble(
        sensitivity=nm,
        min_cells=threshold,
        final_L2_label=pop,
        programme=prog,
        n_HC=n_HC,
        n_AA=n_AA,
        median_HC=med_HC,
        median_AA=med_AA,
        delta_median_AA_HC=
          med_AA-med_HC,
        p=wt$p.value
      )

      k <- k+1L
    }
  }
}

sensitivity <- bind_rows(
  sens_list
) |>
  group_by(
    sensitivity,
    final_L2_label
  ) |>
  mutate(
    FDR=
      p.adjust(
        p,
        method="BH"
      )
  ) |>
  ungroup()

# ============================================================
# 4. CABIMER HC ANCHOR
# ============================================================

cat("\n===== CABIMER HC ANCHOR =====\n")

anchor <- x |>
  filter(
    patient_id=="CON276"
  ) |>
  dplyr::select(
    final_L2_label,
    programme,
    anchor_score=score,
    anchor_n_cells=n_cells
  )

anchor_comparison <- x |>
  filter(
    disease=="AA"
  ) |>
  group_by(
    final_L2_label,
    programme
  ) |>
  summarise(
    n_AA=n(),
    median_AA=
      median(
        score,
        na.rm=TRUE
      ),
    mean_AA=
      mean(
        score,
        na.rm=TRUE
      ),
    .groups="drop"
  ) |>
  left_join(
    anchor,
    by=c(
      "final_L2_label",
      "programme"
    )
  ) |>
  mutate(
    delta_AA_vs_CABIMER_HC=
      median_AA-anchor_score,

    anchor_direction=case_when(
      delta_AA_vs_CABIMER_HC > 0
        ~ "AA_HIGHER",
      delta_AA_vs_CABIMER_HC < 0
        ~ "AA_LOWER",
      TRUE
        ~ "NO_CHANGE"
    )
  )

# ============================================================
# 5. HISTORICAL-HC SENSITIVITY
# ============================================================

historical_comparison <- x |>
  filter(
    disease=="HC",
    batch=="historical"
  ) |>
  group_by(
    final_L2_label,
    programme
  ) |>
  summarise(
    n_historical_HC=n(),
    median_historical_HC=
      median(
        score,
        na.rm=TRUE
      ),
    .groups="drop"
  ) |>
  left_join(
    x |>
      filter(
        disease=="AA"
      ) |>
      group_by(
        final_L2_label,
        programme
      ) |>
      summarise(
        n_AA=n(),
        median_AA=
          median(
            score,
            na.rm=TRUE
          ),
        .groups="drop"
      ),
    by=c(
      "final_L2_label",
      "programme"
    )
  ) |>
  mutate(
    delta_AA_vs_historical_HC=
      median_AA-
      median_historical_HC,

    historical_direction=case_when(
      delta_AA_vs_historical_HC > 0
        ~ "AA_HIGHER",
      delta_AA_vs_historical_HC < 0
        ~ "AA_LOWER",
      TRUE
        ~ "NO_CHANGE"
    )
  )

# ============================================================
# 6. INTEGRATE ROBUSTNESS
# ============================================================

integrated <- primary |>
  left_join(
    anchor_comparison |>
      dplyr::select(
        final_L2_label,
        programme,
        anchor_score,
        anchor_n_cells,
        delta_AA_vs_CABIMER_HC,
        anchor_direction
      ),
    by=c(
      "final_L2_label",
      "programme"
    )
  ) |>
  left_join(
    historical_comparison |>
      dplyr::select(
        final_L2_label,
        programme,
        delta_AA_vs_historical_HC,
        historical_direction
      ),
    by=c(
      "final_L2_label",
      "programme"
    )
  ) |>
  mutate(
    direction_concordant=
      direction==
        anchor_direction &
      direction==
        historical_direction,

    functional_class=case_when(

      FDR_within_population < 0.05 &
        direction_concordant
        ~ "ROBUST_SUPPORTED",

      FDR_within_population < 0.10 &
        direction_concordant
        ~ "SUPPORTED",

      p < 0.05 &
        direction_concordant
        ~ "NOMINAL_CONCORDANT",

      p < 0.10 &
        direction_concordant
        ~ "TREND_CONCORDANT",

      direction_concordant &
        abs(delta_median_AA_HC) > 0
        ~ "DIRECTIONALLY_CONCORDANT",

      TRUE
        ~ "UNSUPPORTED"
    )
  )

# ============================================================
# 7. PRIMARY BIOLOGICAL PROGRAMS
# ============================================================

biological_priority <- integrated |>
  filter(
    programme %in%c(
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "CYTOTOXICITY"
    )
  ) |>
  arrange(
    factor(
      functional_class,
      levels=c(
        "ROBUST_SUPPORTED",
        "SUPPORTED",
        "NOMINAL_CONCORDANT",
        "TREND_CONCORDANT",
        "DIRECTIONALLY_CONCORDANT",
        "UNSUPPORTED"
      )
    ),
    FDR_within_population,
    p
  )

# ============================================================
# 8. REPORT
# ============================================================

cat("\n===== PRIMARY AA vs HC RESULTS =====\n")

print(
  integrated |>
    arrange(
      FDR_within_population,
      p
    ),
  n=Inf,
  width=Inf
)

cat("\n===== BIOLOGICALLY PRIORITIZED PROGRAMS =====\n")

print(
  biological_priority |>
    dplyr::select(
      final_L2_label,
      programme,
      n_HC,
      n_AA,
      median_HC,
      median_AA,
      delta_median_AA_HC,
      p,
      FDR_within_population,
      delta_AA_vs_CABIMER_HC,
      delta_AA_vs_historical_HC,
      direction_concordant,
      functional_class
    ),
  n=Inf,
  width=Inf
)

cat("\n===== CELL-SUPPORT SENSITIVITY — CD8/NK =====\n")

print(
  sensitivity |>
    filter(
      final_L2_label %in%c(
        "CD8 T",
        "NK"
      )
    ) |>
    arrange(
      final_L2_label,
      programme,
      min_cells
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 9. SAVE
# ============================================================

write.csv(
  primary,
  file.path(
    TAB,
    "22G4C2C_AA_HC_primary_functional_tests.csv"
  ),
  row.names=FALSE
)

write.csv(
  sensitivity,
  file.path(
    TAB,
    "22G4C2C_cell_support_sensitivity.csv"
  ),
  row.names=FALSE
)

write.csv(
  anchor_comparison,
  file.path(
    TAB,
    "22G4C2C_CABIMER_HC_anchor.csv"
  ),
  row.names=FALSE
)

write.csv(
  historical_comparison,
  file.path(
    TAB,
    "22G4C2C_historical_HC_sensitivity.csv"
  ),
  row.names=FALSE
)

write.csv(
  integrated,
  file.path(
    TAB,
    "22G4C2C_integrated_functional_remodeling.csv"
  ),
  row.names=FALSE
)

write.csv(
  biological_priority,
  file.path(
    TAB,
    "22G4C2C_biologically_prioritized_programs.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 10. FREEZE
# ============================================================

freeze <- list(
  module="22G4C2C",

  purpose=
    "AA-vs-HC targeted immune functional remodeling",

  statistical_unit=
    "patient pseudobulk",

  populations=
    PRIMARY_POPULATIONS,

  primary_results=
    primary,

  support_sensitivity=
    sensitivity,

  CABIMER_HC_anchor=
    anchor_comparison,

  historical_HC_sensitivity=
    historical_comparison,

  integrated_results=
    integrated,

  biological_priority=
    biological_priority,

  batch_constraint=paste(
    "Disease and batch are strongly confounded.",
    "Disease+batch regression is therefore not interpreted",
    "as an identifiable causal disease effect."
  )
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4C2C_AA_HC_immune_functional_remodeling_freeze.rds"
  )
)

# ============================================================
# 11. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "Four immune populations",
    "CD4 tested",
    "CD8 tested",
    "NK tested",
    "Macrophage tested",
    "Cytotoxicity tested in CD8",
    "Cytotoxicity tested in NK",
    "Cell-support sensitivity generated",
    "CABIMER HC anchor generated",
    "Historical HC sensitivity generated",
    "Integrated table generated",
    "Freeze exists"
  ),

  pass=c(
    n_distinct(
      integrated$final_L2_label
    )==4,

    "CD4 T" %in%
      integrated$final_L2_label,

    "CD8 T" %in%
      integrated$final_L2_label,

    "NK" %in%
      integrated$final_L2_label,

    "Macrophage/monocyte" %in%
      integrated$final_L2_label,

    any(
      integrated$final_L2_label=="CD8 T" &
      integrated$programme=="CYTOTOXICITY"
    ),

    any(
      integrated$final_L2_label=="NK" &
      integrated$programme=="CYTOTOXICITY"
    ),

    nrow(sensitivity)>0,

    nrow(anchor_comparison)>0,

    nrow(historical_comparison)>0,

    nrow(integrated)>0,

    file.exists(
      file.path(
        OBJ,
        "22G4C2C_AA_HC_immune_functional_remodeling_freeze.rds"
      )
    )
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
    "22G4C2C_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 12. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C2C — AA vs HC IMMUNE FUNCTIONAL REMODELING",
  "============================================================",
  "",
  "Primary question:",
  paste(
    "Do frozen immune populations exhibit functional-state",
    "remodeling in AA despite limited compositional change?"
  ),
  "",
  "
Disease universe: AA + HC.",
  "",
  "Populations tested:",
  "- CD4 T",
  "- CD8 T",
  "- NK",
  "- Macrophage/monocyte",
  "",
  "Pre-specified functional programmes:",
  "- IFN-gamma response",
  "- IFN-alpha response",
  "- TNFA/NFKB signaling",
  "- IL2/STAT5 signaling",
  "- inflammatory response",
  "- complement",
  "- apoptosis",
  "- cytotoxicity in CD8 T and NK",
  "",
  "Primary comparison:",
  "Patient-level AA versus HC functional scores.",
  "",
  "Multiple testing:",
  paste(
    "BH FDR is calculated within each immune population;",
    "a global FDR is also retained for audit."
  ),
  "",
  "Cell-support sensitivity:",
  paste(
    "Analyses were repeated using all pseudobulks and",
    "minimum cell thresholds of >=5, >=10 and >=20 cells."
  ),
  "",
  "Batch sensitivity:",
  paste(
    "AA is strongly confounded with CABIMER2026 batch;",
    "therefore disease + batch regression is not treated",
    "as an identifiable disease-effect model."
  ),
  paste(
    "CON276 is the only HC CABIMER2026 patient and is used",
    "as a directional same-batch anchor, not as an",
    "independent inferential comparison."
  ),
  paste(
    "Historical HC samples are retained as a complementary",
    "directional sensitivity analysis."
  ),
  "",
  "Interpretation hierarchy:",
  "ROBUST_SUPPORTED: within-population FDR <0.05 and directional concordance.",
  "SUPPORTED: within-population FDR <0.10 and directional concordance.",
  "NOMINAL_CONCORDANT: nominal P <0.05 and directional concordance.",
  "TREND_CONCORDANT: nominal P <0.10 and directional concordance.",
  paste(
    "DIRECTIONALLY_CONCORDANT: same direction versus pooled HC,",
    "CABIMER HC anchor and historical HC without nominal support."
  ),
  "UNSUPPORTED: criteria above not fulfilled.",
  "",
  paste(
    "NK is retained as a biologically prioritized Tier-2",
    "population despite low patient-level cell support."
  ),
  paste(
    "NK and CD8 findings require explicit evaluation of",
    "cell-support sensitivity before biological interpretation."
  ),
  "",
  paste(
    "Functional-state differences are interpreted separately",
    "from immune-cell abundance differences."
  ),
  "",
  paste0(
    "Functional tests: ",
    nrow(primary)
  ),
  paste0(
    "Integrated results: ",
    nrow(integrated)
  ),
  paste0(
    "Biologically prioritized rows: ",
    nrow(biological_priority)
  ),
  "",
  paste0(
    "ROBUST_SUPPORTED: ",
    sum(
      integrated$functional_class==
        "ROBUST_SUPPORTED",
      na.rm=TRUE
    )
  ),
  paste0(
    "SUPPORTED: ",
    sum(
      integrated$functional_class==
        "SUPPORTED",
      na.rm=TRUE
    )
  ),
  paste0(
    "NOMINAL_CONCORDANT: ",
    sum(
      integrated$functional_class==
        "NOMINAL_CONCORDANT",
      na.rm=TRUE
    )
  ),
  paste0(
    "TREND_CONCORDANT: ",
    sum(
      integrated$functional_class==
        "TREND_CONCORDANT",
      na.rm=TRUE
    )
  ),
  paste0(
    "DIRECTIONALLY_CONCORDANT: ",
    sum(
      integrated$functional_class==
        "DIRECTIONALLY_CONCORDANT",
      na.rm=TRUE
    )
  ),
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
    "22G4C2C_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

# ============================================================
# 13. FINAL STATUS
# ============================================================

if(!all(checks$pass)) {
  stop(
    "22G4C2C formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C2C COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

