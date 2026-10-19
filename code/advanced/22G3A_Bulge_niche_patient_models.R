
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
})

cat("\n============================================================\n")
cat("22G3A — BULGE-CENTERED MESENCHYMAL/VASCULAR NICHE COUPLING\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

CFILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F4_BULGE_MONOCLE3_TRAJECTORIES",
  "22F4C_TRUNK_PSEUDOTIME",
  "objects",
  "22F4C_trunk_definition_freeze.rds"
)

SCORE_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F2_FUNCTIONAL_STATE_DENSITY",
  "22F2A_SCORE_FREEZE",
  "tables",
  "22F2A_Bulge_patient_functional_scores.csv"
)

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
  "22G3_MESENCHYMAL_VASCULAR_NICHE",
  "22G3A_PATIENT_MODELLING"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

for(d in c(TAB,OBJ,AUD)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

stopifnot(
  file.exists(CFILE),
  file.exists(SCORE_FILE),
  file.exists(L2_FILE)
)

# ============================================================
# 2. LOAD FROZEN BULGE PHENOTYPE
# ============================================================

C <- readRDS(CFILE)

position <- as_tibble(
  C$patient_position
) |>
  dplyr::select(
    patient_id,
    disease,
    n_cells,
    median_position,
    mean_position,
    pct_T1,
    pct_T2,
    pct_T3,
    pct_stem_high
  )

scores <- read.csv(
  SCORE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  dplyr::select(
    patient_id,
    disease,
    n_bulge,
    stem_niche_mean,
    MHCII_mean,
    cycling_mean,
    CD200_z_mean
  )

bulge <- position |>
  left_join(
    scores,
    by=c(
      "patient_id",
      "disease"
    )
  ) |>
  mutate(
    disease=factor(
      disease,
      levels=c(
        "HC",
        "FFA",
        "AA"
      )
    )
  )

stopifnot(
  nrow(bulge)==16,
  !anyDuplicated(bulge$patient_id)
)

# ============================================================
# 3. NICHE UNIVERSE
# ============================================================

NICHE_L1 <- c(
  "Fibroblast",
  "Pericyte_SMC",
  "Smooth_muscle",
  "Vascular_EC",
  "Lymphatic_EC"
)

raw <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  filter(
    L1 %in% NICHE_L1,
    eligible
  )

cat("\n===== RAW NICHE UNIVERSE =====\n")
cat("L1:",paste(sort(unique(raw$L1)),collapse=", "),"\n")
cat("L2 identities:",n_distinct(paste(raw$L1,raw$L2)),"\n")

# ============================================================
# 4. PRE-SPECIFIED SUPPORT FILTER
# ============================================================

support_by_disease <- raw |>
  mutate(
    positive=n>0
  ) |>
  group_by(
    L1,
    L2,
    disease
  ) |>
  summarise(
    n_positive=sum(positive),
    .groups="drop"
  ) |>
  pivot_wider(
    names_from=disease,
    values_from=n_positive,
    values_fill=0,
    names_prefix="positive_"
  )

support <- raw |>
  group_by(
    L1,
    L2
  ) |>
  summarise(
    n_patients=n_distinct(patient_id),
    positive_patients=sum(n>0),
    total_cells=sum(n),
    median_cells=median(n),
    zero_patients=sum(n==0),
    .groups="drop"
  ) |>
  left_join(
    support_by_disease,
    by=c(
      "L1",
      "L2"
    )
  ) |>
  mutate(
    eligible_niche=
      total_cells >= 40 &
      positive_patients >= 12 &
      positive_HC >= 3 &
      positive_FFA >= 3 &
      positive_AA >= 3,

    exclusion_reason=case_when(
      eligible_niche ~ "INCLUDE",
      total_cells < 40 ~ "<40 total cells",
      positive_patients < 12 ~ "<12 positive patients",
      positive_HC < 3 ~ "<3 HC positive",
      positive_FFA < 3 ~ "<3 FFA positive",
      positive_AA < 3 ~ "<3 AA positive",
      TRUE ~ "OTHER"
    ),

    feature=paste(
      L1,
      L2,
      sep=" :: "
    )
  ) |>
  arrange(
    L1,
    L2
  )

write.csv(
  support,
  file.path(
    TAB,
    "22G3A_niche_support_and_eligibility.csv"
  ),
  row.names=FALSE
)

cat("\n===== SUPPORT / ELIGIBILITY =====\n")

print(
  support,
  n=Inf
)

eligible <- support |>
  filter(
    eligible_niche
  ) |>
  pull(feature)

cat(
  "\nEligible niche states:",
  length(eligible),
  "\n"
)

# ============================================================
# ============================================================
# ============================================================
# 5. PREPARE PATIENT-LEVEL OUTCOMES
# ============================================================

# IMPORTANT:
# 20C contains only evaluable within-L1 observations.
# A missing patient-L1 combination must NOT automatically
# be interpreted as frequency=0.
#
# Therefore:
#   - existing patient-L1 rows are retained;
#   - absent L2 within an existing L1 would be zero;
#   - completely absent/non-evaluable L1 remains missing;
#   - models use available patient-level observations;
#   - minimum valid n is 12.

eligible_lookup <- support |>
  filter(
    eligible_niche
  ) |>
  dplyr::select(
    L1,
    L2,
    feature
  )

# Patient × L1 combinations that are actually evaluable
patient_L1 <- raw |>
  dplyr::select(
    patient_id,
    disease,
    L1,
    denominator_L1
  ) |>
  distinct() |>
  filter(
    is.finite(denominator_L1),
    denominator_L1 > 0
  )

# Complete L2 states ONLY inside evaluable patient-L1 pairs.
niche <- eligible_lookup |>
  inner_join(
    patient_L1,
    by="L1",
    relationship="many-to-many"
  ) |>
  left_join(
    raw |>
      dplyr::select(
        patient_id,
        disease,
        L1,
        L2,
        n,
        denominator_L1
      ) |>
      rename(
        denominator_observed=
          denominator_L1
      ),
    by=c(
      "patient_id",
      "disease",
      "L1",
      "L2"
    )
  ) |>
  mutate(
    n=replace_na(
      n,
      0
    ),

    # denominator comes from the evaluable patient-L1 pair
    frequency=
      n/denominator_L1,

    p_corrected=
      (n+0.5)/
      (denominator_L1+1),

    logit_frequency=
      qlogis(
        p_corrected
      )
  )

# ------------------------------------------------------------
# Model-level support after evaluability completion
# ------------------------------------------------------------

model_support <- niche |>
  group_by(
    feature,
    L1,
    L2
  ) |>
  summarise(
    n_valid=n_distinct(patient_id),

    n_HC=n_distinct(
      patient_id[disease=="HC"]
    ),

    n_FFA=n_distinct(
      patient_id[disease=="FFA"]
    ),

    n_AA=n_distinct(
      patient_id[disease=="AA"]
    ),

    .groups="drop"
  ) |>
  mutate(
    model_eligible=
      n_valid >= 12 &
      n_HC >= 3 &
      n_FFA >= 3 &
      n_AA >= 3
  )

write.csv(
  model_support,
  file.path(
    TAB,
    "22G3A_model_level_support.csv"
  ),
  row.names=FALSE
)

cat("\n===== MODEL-LEVEL SUPPORT =====\n")

print(
  model_support,
  n=Inf
)

eligible_model_features <- model_support |>
  filter(
    model_eligible
  ) |>
  pull(feature)

cat(
  "\nEligible states after model-level support:",
  length(eligible_model_features),
  "\n"
)

stopifnot(
  length(eligible_model_features)>0,
  all(
    model_support$n_valid[
      model_support$model_eligible
    ] >= 12
  )
)

# 6. PRIMARY TRAJECTORY MODELS
# ============================================================

cat("\n===== PRIMARY TRAJECTORY MODELS =====\n")

primary_list <- list()
model_store <- list()
k <- 1L

for(feature_i in eligible_model_features) {

  dd <- niche |>
    filter(
      feature==feature_i
    ) |>
    inner_join(
      bulge,
      by=c(
        "patient_id",
        "disease"
      )
    )

  stopifnot(
    nrow(dd) >= 12
  )

  # Global non-parametric association
  sp <- suppressWarnings(
    cor.test(
      dd$median_position,
      dd$frequency,
      method="spearman",
      exact=FALSE
    )
  )

  # M0 — trajectory only
  m0 <- lm(
    logit_frequency ~
      median_position,
    data=dd
  )

  # M1 — trajectory adjusted for disease
  m1 <- lm(
    logit_frequency ~
      median_position +
      disease,
    data=dd
  )

  # M2 — exploratory disease-specific coupling
  m2 <- lm(
    logit_frequency ~
      median_position *
      disease,
    data=dd
  )

  t0 <- broom::tidy(m0) |>
    filter(
      term=="median_position"
    )

  t1 <- broom::tidy(m1) |>
    filter(
      term=="median_position"
    )

  cmp <- anova(
    m1,
    m2
  )

  interaction_p <-
    cmp$`Pr(>F)`[2]

  primary_list[[k]] <- tibble(
    feature=feature_i,
    L1=dd$L1[1],
    L2=dd$L2[1],

    n=nrow(dd),

    spearman_rho=
      unname(
        sp$estimate
      ),

    spearman_p=
      sp$p.value,

    M0_beta=
      t0$estimate,

    M0_SE=
      t0$std.error,

    M0_p=
      t0$p.value,

    M1_adjusted_beta=
      t1$estimate,

    M1_adjusted_SE=
      t1$std.error,

    M1_adjusted_p=
      t1$p.value,

    interaction_p=
      interaction_p,

    R2_M0=
      summary(m0)$r.squared,

    R2_M1=
      summary(m1)$r.squared,

    delta_R2_disease=
      summary(m1)$r.squared -
      summary(m0)$r.squared
  )

  model_store[[feature_i]] <- list(
    data=dd,
    M0=m0,
    M1=m1,
    M2=m2
  )

  k <- k+1L
}

primary <- bind_rows(
  primary_list
) |>
  mutate(
    spearman_FDR=
      p.adjust(
        spearman_p,
        method="BH"
      ),

    M0_FDR=
      p.adjust(
        M0_p,
        method="BH"
      ),

    M1_adjusted_FDR=
      p.adjust(
        M1_adjusted_p,
        method="BH"
      ),

    interaction_FDR=
      p.adjust(
        interaction_p,
        method="BH"
      ),

    adjusted_direction=
      case_when(
        M1_adjusted_beta > 0 ~
          "later Bulge trajectory",

        M1_adjusted_beta < 0 ~
          "earlier Bulge trajectory",

        TRUE ~
          "neutral"
      )
  ) |>
  arrange(
    M1_adjusted_FDR,
    desc(abs(spearman_rho))
  )

write.csv(
  primary,
  file.path(
    TAB,
    "22G3A_primary_niche_trajectory_models.csv"
  ),
  row.names=FALSE
)

print(
  primary,
  n=Inf
)

# ============================================================
# 7. WITHIN-DISEASE DIRECTIONAL COUPLING
# ============================================================

cat("\n===== WITHIN-DISEASE DIRECTIONAL COUPLING =====\n")

within_list <- list()
k <- 1L

for(feature_i in eligible_model_features) {

  dd0 <- niche |>
    filter(
      feature==feature_i
    ) |>
    inner_join(
      bulge,
      by=c(
        "patient_id",
        "disease"
      )
    )

  for(dis_i in c(
    "HC",
    "FFA",
    "AA"
  )) {

    dd <- dd0 |>
      filter(
        disease==dis_i
      )

    rho <- NA_real_

    if(
      nrow(dd)>=4 &&
      sd(dd$median_position)>0 &&
      sd(dd$frequency)>0
    ) {
      rho <- suppressWarnings(
        cor(
          dd$median_position,
          dd$frequency,
          method="spearman"
        )
      )
    }

    within_list[[k]] <- tibble(
      feature=feature_i,
      L1=dd$L1[1],
      L2=dd$L2[1],
      disease=dis_i,
      n=nrow(dd),
      rho=rho
    )

    k <- k+1L
  }
}

within <- bind_rows(
  within_list
)

write.csv(
  within,
  file.path(
    TAB,
    "22G3A_within_disease_directional_coupling.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 8. BULGE FUNCTIONAL AXES
# ============================================================

cat("\n===== BULGE FUNCTIONAL AXES =====\n")

functional_axes <- c(
  "stem_niche_mean",
  "MHCII_mean",
  "cycling_mean",
  "CD200_z_mean"
)

functional_list <- list()
k <- 1L

for(feature_i in eligible_model_features) {

  dd <- niche |>
    filter(
      feature==feature_i
    ) |>
    inner_join(
      bulge,
      by=c(
        "patient_id",
        "disease"
      )
    )

  for(axis_i in functional_axes) {

    sp <- suppressWarnings(
      cor.test(
        dd[[axis_i]],
        dd$frequency,
        method="spearman",
        exact=FALSE
      )
    )

    f <- as.formula(
      paste(
        "logit_frequency ~",
        axis_i,
        "+ disease"
      )
    )

    fit <- lm(
      f,
      data=dd
    )

    tt <- broom::tidy(
      fit
    ) |>
      filter(
        term==axis_i
      )

    functional_list[[k]] <- tibble(
      feature=feature_i,
      L1=dd$L1[1],
      L2=dd$L2[1],
      bulge_axis=axis_i,

      rho=
        unname(
          sp$estimate
        ),

      spearman_p=
        sp$p.value,

      adjusted_beta=
        tt$estimate,

      adjusted_SE=
        tt$std.error,

      adjusted_p=
        tt$p.value
    )

    k <- k+1L
  }
}

functional <- bind_rows(
  functional_list
) |>
  group_by(
    bulge_axis
  ) |>
  mutate(
    spearman_FDR=
      p.adjust(
        spearman_p,
        method="BH"
      ),

    adjusted_FDR=
      p.adjust(
        adjusted_p,
        method="BH"
      )
  ) |>
  ungroup() |>
  arrange(
    adjusted_FDR,
    desc(abs(rho))
  )

write.csv(
  functional,
  file.path(
    TAB,
    "22G3A_functional_Bulge_niche_coupling.csv"
  ),
  row.names=FALSE
)

cat("\nTOP FUNCTIONAL COUPLINGS:\n")

print(
  functional |>
    slice_head(
      n=40
    ),
  n=Inf
)

# ============================================================
# 9. DISEASE EFFECTS
# ============================================================

cat("\n===== DISEASE EFFECTS =====\n")

disease_list <- list()
k <- 1L

for(feature_i in eligible_model_features) {

  dd <- niche |>
    filter(
      feature==feature_i
    )

  med <- dd |>
    group_by(
      disease
    ) |>
    summarise(
      median_frequency=
        median(
          frequency,
          na.rm=TRUE
        ),
      .groups="drop"
    )

  getmed <- function(x) {
    z <- med$median_frequency[
      med$disease==x
    ]
    if(length(z)) z else NA_real_
  }

  aa <- dd |>
    filter(
      disease=="AA"
    ) |>
    pull(frequency)

  ffa <- dd |>
    filter(
      disease=="FFA"
    ) |>
    pull(frequency)

  wt <- suppressWarnings(
    wilcox.test(
      aa,
      ffa,
      exact=FALSE
    )
  )

  kw <- kruskal.test(
    frequency ~ disease,
    data=dd
  )

  disease_list[[k]] <- tibble(
    feature=feature_i,
    L1=dd$L1[1],
    L2=dd$L2[1],

    median_HC=getmed("HC"),
    median_FFA=getmed("FFA"),
    median_AA=getmed("AA"),

    delta_AA_FFA=
      getmed("AA")-
      getmed("FFA"),

    KW_p=
      kw$p.value,

    AA_FFA_p=
      wt$p.value
  )

  k <- k+1L
}

disease_effects <- bind_rows(
  disease_list
) |>
  mutate(
    KW_FDR=
      p.adjust(
        KW_p,
        method="BH"
      ),

    AA_FFA_FDR=
      p.adjust(
        AA_FFA_p,
        method="BH"
      )
  )

write.csv(
  disease_effects,
  file.path(
    TAB,
    "22G3A_niche_disease_effects.csv"
  ),
  row.names=FALSE
)
# ============================================================
# 10. INTEGRATE WITHIN-DISEASE DIRECTIONALITY
# ============================================================

within_wide <- within |>
  dplyr::select(
    feature,
    disease,
    rho
  ) |>
  pivot_wider(
    names_from=disease,
    values_from=rho,
    names_prefix="rho_"
  )

# ============================================================
# 11. PRIOR BIOLOGICAL EVIDENCE
# ============================================================

# These flags DO NOT alter statistical significance.
# They are retained only for biological prioritization after
# the statistical analysis has been completed.

prior_evidence <- support |>
  filter(
    eligible_niche
  ) |>
  dplyr::select(
    feature,
    L1,
    L2
  ) |>
  mutate(
    prior_biological_evidence=case_when(

      L1=="Fibroblast" &
        L2=="PTGDS/IGF1 fibroblast"
        ~ "IGF trophic-support candidate",

      L1=="Smooth_muscle" &
        L2=="APM-bulge interface smooth muscle"
        ~ "Direct anatomical Bulge-interface candidate",

      L1=="Fibroblast" &
        L2=="COCH/AQP4 dermal papilla"
        ~ "Specialized follicular mesenchymal niche",

      TRUE
        ~ "None pre-specified"
    ),

    prior_flag=
      prior_biological_evidence !=
        "None pre-specified"
  )

# ============================================================
# 12. INTEGRATED NICHE CLASSIFICATION
# ============================================================

cat("\n===== INTEGRATED NICHE CLASSIFICATION =====\n")

integrated <- primary |>
  left_join(
    within_wide,
    by="feature"
  ) |>
  left_join(
    disease_effects,
    by=c(
      "feature",
      "L1",
      "L2"
    )
  ) |>
  left_join(
    prior_evidence |>
      dplyr::select(
        feature,
        prior_biological_evidence,
        prior_flag
      ),
    by="feature"
  ) |>
  mutate(
    global_strength=case_when(
      spearman_FDR < 0.05 &
        abs(spearman_rho) >= 0.50
        ~ "STRONG",

      spearman_FDR < 0.10 &
        abs(spearman_rho) >= 0.50
        ~ "SUPPORTED",

      abs(spearman_rho) >= 0.50
        ~ "EXPLORATORY",

      TRUE
        ~ "WEAK"
    ),

    adjusted_strength=case_when(
      M1_adjusted_FDR < 0.05
        ~ "STRONG",

      M1_adjusted_FDR < 0.10
        ~ "SUPPORTED",

      M1_adjusted_p < 0.05
        ~ "NOMINAL",

      TRUE
        ~ "NOT_SUPPORTED"
    ),

    interaction_strength=case_when(
      interaction_FDR < 0.10
        ~ "SUPPORTED",

      interaction_p < 0.05
        ~ "NOMINAL",

      TRUE
        ~ "NOT_SUPPORTED"
    ),

    direction_support_n=
      rowSums(
        cbind(
          !is.na(rho_HC) &
            sign(rho_HC)==
              sign(spearman_rho),

          !is.na(rho_FFA) &
            sign(rho_FFA)==
              sign(spearman_rho),

          !is.na(rho_AA) &
            sign(rho_AA)==
              sign(spearman_rho)
        ),
        na.rm=TRUE
      ),

    niche_class=case_when(

      M1_adjusted_FDR < 0.10
        ~ "TRAJECTORY_ASSOCIATED",

      M1_adjusted_p < 0.05
        ~ "TRAJECTORY_ASSOCIATED_NOMINAL",

      interaction_FDR < 0.10
        ~ "DISEASE_SPECIFIC_COUPLING",

      interaction_p < 0.05
        ~ "DISEASE_SPECIFIC_COUPLING_NOMINAL",

      spearman_FDR < 0.10 &
        abs(spearman_rho) >= 0.50
        ~ "DISEASE_STRUCTURED_CO_REMODELING",

      abs(spearman_rho) >= 0.50
        ~ "EXPLORATORY_CO_REMODELING",

      TRUE
        ~ "UNSUPPORTED"
    ),

    priority_score=
      abs(spearman_rho) +

      0.35*(
        M1_adjusted_p < 0.05
      ) +

      0.30*(
        M1_adjusted_FDR < 0.10
      ) +

      0.20*(
        interaction_p < 0.05
      ) +

      0.15*(
        direction_support_n >= 2
      ) +

      0.10*(
        AA_FFA_FDR < 0.10
      )
  ) |>
  arrange(
    factor(
      niche_class,
      levels=c(
        "TRAJECTORY_ASSOCIATED",
        "TRAJECTORY_ASSOCIATED_NOMINAL",
        "DISEASE_SPECIFIC_COUPLING",
        "DISEASE_SPECIFIC_COUPLING_NOMINAL",
        "DISEASE_STRUCTURED_CO_REMODELING",
        "EXPLORATORY_CO_REMODELING",
        "UNSUPPORTED"
      )
    ),
    desc(priority_score)
  )

write.csv(
  integrated,
  file.path(
    TAB,
    "22G3A_integrated_niche_classification.csv"
  ),
  row.names=FALSE
)

cat("\n===== PRIORITIZED NICHE RESULTS =====\n")

print(
  integrated |>
    dplyr::select(
      L1,
      L2,
      spearman_rho,
      spearman_FDR,
      M1_adjusted_beta,
      M1_adjusted_p,
      M1_adjusted_FDR,
      interaction_p,
      interaction_FDR,
      rho_HC,
      rho_FFA,
      rho_AA,
      delta_AA_FFA,
      AA_FFA_FDR,
      niche_class,
      prior_biological_evidence
    ),
  n=Inf
)

# ============================================================
# 13. TOP CANDIDATE LEDGER
# ============================================================

candidate_ledger <- integrated |>
  filter(
    niche_class !=
      "UNSUPPORTED" |
      prior_flag
  ) |>
  arrange(
    desc(priority_score)
  )

write.csv(
  candidate_ledger,
  file.path(
    TAB,
    "22G3A_candidate_niche_ledger.csv"
  ),
  row.names=FALSE
)

cat("\n===== CANDIDATE NICHE LEDGER =====\n")

print(
  candidate_ledger |>
    dplyr::select(
      L1,
      L2,
      spearman_rho,
      spearman_FDR,
      M1_adjusted_beta,
      M1_adjusted_p,
      interaction_p,
      rho_HC,
      rho_FFA,
      rho_AA,
      niche_class,
      prior_biological_evidence
    ),
  n=Inf
)

# ============================================================
# 14. PRE-SPECIFIED CANDIDATES
# ============================================================

prespecified <- integrated |>
  filter(
    prior_flag
  )

cat("\n===== PRE-SPECIFIED BIOLOGICAL CANDIDATES =====\n")

print(
  prespecified |>
    dplyr::select(
      L1,
      L2,
      spearman_rho,
      spearman_FDR,
      M1_adjusted_beta,
      M1_adjusted_p,
      M1_adjusted_FDR,
      interaction_p,
      rho_HC,
      rho_FFA,
      rho_AA,
      niche_class,
      prior_biological_evidence
    ),
  n=Inf
)

write.csv(
  prespecified,
  file.path(
    TAB,
    "22G3A_prespecified_candidate_results.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 15. FREEZE
# ============================================================

freeze <- list(
  module="22G3A",

  purpose=paste(
    "Patient-level analysis of mesenchymal, mural,",
    "vascular and lymphatic niche-state coupling to",
    "the frozen Bulge T3-to-T2-to-T1 trajectory."
  ),

  statistical_unit="patient",

  n_patients=
    n_distinct(
      bulge$patient_id
    ),

  disease_counts=
    bulge |>
      count(
        disease
      ),

  primary_predictor=
    "Bulge median_position",

  outcome_definition=
    paste(
      "continuity-corrected logit of frozen L2",
      "frequency within its corresponding L1 compartment"
    ),

  support_rule=paste(
    ">=40 total cells;",
    ">=12 positive patients;",
    ">=3 positive patients in each disease"
  ),

  niche_L1=NICHE_L1,

  support=support,

  eligible_features_pre_support=eligible,
  eligible_model_features=eligible_model_features,
  model_support=model_support,

  bulge_patient_data=bulge,

  primary=primary,

  within_disease=within,

  functional=functional,

  disease_effects=disease_effects,

  prior_evidence=prior_evidence,

  integrated=integrated,

  candidate_ledger=candidate_ledger,

  model_store=model_store
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G3A_Bulge_mesenchymal_vascular_niche_freeze.rds"
  )
)

# ============================================================
# 16. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

checks <- tibble(
  check=c(
    "16 patients",
    "Five niche L1 included",
    "Eligibility table generated",
    "At least one eligible niche state",
    "All model-eligible states have >=12 patient rows",
    "Primary model table complete",
    "Within-disease table complete",
    "Functional table complete",
    "Disease-effect table complete",
    "Integrated table complete",
    "PTGDS/IGF1 audited",
    "APM-bulge interface audited if eligible",
    "COCH/AQP4 audited if eligible",
    "Freeze exists"
  ),

  pass=c(
    n_distinct(
      bulge$patient_id
    )==16,

    all(
      NICHE_L1 %in%
        unique(raw$L1)
    ),

    nrow(support)>0,

    length(eligible)>0,

    all(
      model_support$n_valid[
        model_support$model_eligible
      ] >= 12
    ),

    nrow(primary)==
      length(eligible_model_features),

    nrow(within)==
      length(eligible_model_features)*3,

    nrow(functional)==
      length(eligible_model_features)*4,

    nrow(disease_effects)==
      length(eligible_model_features),

    nrow(integrated)==
      length(eligible_model_features),

    any(
      integrated$L2==
        "PTGDS/IGF1 fibroblast"
    ),

    if(
      "Smooth_muscle :: APM-bulge interface smooth muscle" %in% eligible
    ) {
      any(
        integrated$feature==
          "Smooth_muscle :: APM-bulge interface smooth muscle"
      )
    } else {
      TRUE
    },

    if(
      "Fibroblast :: COCH/AQP4 dermal papilla" %in% eligible
    ) {
      any(
        integrated$feature==
          "Fibroblast :: COCH/AQP4 dermal papilla"
      )
    } else {
      TRUE
    },

    file.exists(
      file.path(
        OBJ,
        "22G3A_Bulge_mesenchymal_vascular_niche_freeze.rds"
      )
    )
  )
)

print(
  checks,
  n=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G3A_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 17. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G3A — BULGE-CENTERED MESENCHYMAL/VASCULAR NICHE",
  "============================================================",
  "",
  paste0(
    "Patients: ",
    n_distinct(
      bulge$patient_id
    )
  ),

  paste0(
    "Raw niche states: ",
    nrow(support)
  ),

  paste0(
    "Eligible niche states for modelling: ",
    length(eligible_model_features)
  ),

  paste0(
    "Trajectory-associated: ",
    sum(
      integrated$niche_class %in% c(
          "TRAJECTORY_ASSOCIATED",
          "TRAJECTORY_ASSOCIATED_NOMINAL"
        )
    )
  ),

  paste0(
    "Disease-specific coupling: ",
    sum(
      integrated$niche_class %in% c(
          "DISEASE_SPECIFIC_COUPLING",
          "DISEASE_SPECIFIC_COUPLING_NOMINAL"
        )
    )
  ),

  paste0(
    "Disease-structured co-remodeling: ",
    sum(
      integrated$niche_class==
        "DISEASE_STRUCTURED_CO_REMODELING"
    )
  ),

  paste0(
    "Exploratory co-remodeling: ",
    sum(
      integrated$niche_class==
        "EXPLORATORY_CO_REMODELING"
    )
  ),

  "",
  "Interpretation constraints:",
  "- Patient is the statistical unit.",
  "- Frequencies are calculated within each frozen L1 compartment.",
  "- Disease-adjusted models are primary for trajectory-associated niche inference.",
  "- Interaction models are exploratory because n=16.",
  "- Prior biological evidence does not modify statistical significance.",
  "- Association does not establish causal direction.",
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
    "22G3A_audit.txt"
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
    "22G3A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G3A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

