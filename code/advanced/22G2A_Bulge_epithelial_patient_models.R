
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
})

cat("\n============================================================\n")
cat("22G2A — BULGE-TO-EPITHELIAL PROPAGATION\n")
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
  "22G2_EPITHELIAL_PROPAGATION",
  "22G2A_PATIENT_MODELLING"
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
# 2. LOAD FROZEN INPUTS
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

epi <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  filter(
    L1=="Epithelial",
    L2!="Bulge_Isthmus",
    eligible
  )

cat("\nPatients:",n_distinct(position$patient_id),"\n")
cat("External epithelial L2:",n_distinct(epi$L2),"\n")

# ============================================================
# 3. MASTER PATIENT DATA
# ============================================================

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
# 4. EPITHELIAL OUTCOMES
# ============================================================

# Continuity-corrected logit:
# (n + 0.5)/(denominator + 1)

epi <- epi |>
  mutate(
    p_corrected=
      (n+0.5)/
      (denominator_L1+1),

    logit_frequency=
      qlogis(
        p_corrected
      )
  )

support <- epi |>
  group_by(L2) |>
  summarise(
    n_patients=n_distinct(patient_id),
    total_cells=sum(n),
    zero_patients=sum(n==0),
    median_frequency=median(frequency),
    .groups="drop"
  )

write.csv(
  support,
  file.path(
    TAB,
    "22G2A_epithelial_L2_support.csv"
  ),
  row.names=FALSE
)

print(
  support,
  n=Inf
)

# ============================================================
# 5. PRIMARY MODELS
# ============================================================

results <- list()
k <- 1L

for(state in sort(unique(epi$L2))) {

  dd <- epi |>
    filter(L2==state) |>
    inner_join(
      bulge,
      by=c(
        "patient_id",
        "disease"
      )
    )

  # ----------------------------------------------------------
  # Spearman
  # ----------------------------------------------------------

  sp <- suppressWarnings(
    cor.test(
      dd$median_position,
      dd$frequency,
      method="spearman",
      exact=FALSE
    )
  )

  # ----------------------------------------------------------
  # M0: trajectory only
  # ----------------------------------------------------------

  m0 <- lm(
    logit_frequency ~
      median_position,
    data=dd
  )

  t0 <- broom::tidy(m0) |>
    filter(
      term=="median_position"
    )

  # ----------------------------------------------------------
  # M1: trajectory + disease
  # ----------------------------------------------------------

  m1 <- lm(
    logit_frequency ~
      median_position +
      disease,
    data=dd
  )

  t1 <- broom::tidy(m1) |>
    filter(
      term=="median_position"
    )

  # ----------------------------------------------------------
  # M2: trajectory × disease
  # Exploratory because n=16.
  # ----------------------------------------------------------

  m2 <- lm(
    logit_frequency ~
      median_position *
      disease,
    data=dd
  )

  a1 <- anova(
    m1,
    m2
  )

  interaction_p <-
    a1$`Pr(>F)`[2]

  results[[k]] <- tibble(
    L2=state,

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
      summary(m1)$r.squared
  )

  k <- k+1L
}

primary <- bind_rows(
  results
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
          "later trajectory",

        M1_adjusted_beta < 0 ~
          "earlier trajectory",

        TRUE ~
          "neutral"
      )
  ) |>
  arrange(
    M1_adjusted_FDR,
    desc(
      abs(
        M1_adjusted_beta
      )
    )
  )

write.csv(
  primary,
  file.path(
    TAB,
    "22G2A_primary_epithelial_propagation_models.csv"
  ),
  row.names=FALSE
)

cat("\n===== PRIMARY MODELS =====\n")

print(
  primary,
  n=Inf
)

# ============================================================
# 6. FUNCTIONAL BULGE AXES
# ============================================================

functional_axes <- c(
  "stem_niche_mean",
  "MHCII_mean",
  "cycling_mean",
  "CD200_z_mean"
)

functional <- list()
k <- 1L

for(state in sort(unique(epi$L2))) {

  dd <- epi |>
    filter(
      L2==state
    ) |>
    inner_join(
      bulge,
      by=c(
        "patient_id",
        "disease"
      )
    )

  for(axis in functional_axes) {

    x <- dd[[axis]]
    y <- dd$frequency

    sp <- suppressWarnings(
      cor.test(
        x,
        y,
        method="spearman",
        exact=FALSE
      )
    )

    # Adjusted model
    formula_i <- as.formula(
      paste(
        "logit_frequency ~",
        axis,
        "+ disease"
      )
    )

    fit <- lm(
      formula_i,
      data=dd
    )

    tt <- broom::tidy(
      fit
    ) |>
      filter(
        term==axis
      )

    functional[[k]] <- tibble(
      L2=state,
      bulge_axis=axis,
      rho=unname(
        sp$estimate
      ),
      spearman_p=sp$p.value,
      adjusted_beta=tt$estimate,
      adjusted_SE=tt$std.error,
      adjusted_p=tt$p.value
    )

    k <- k+1L
  }
}

functional <- bind_rows(
  functional
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
    "22G2A_functional_Bulge_epithelial_coupling.csv"
  ),
  row.names=FALSE
)

cat("\n===== FUNCTIONAL BULGE AXES =====\n")

print(
  functional,
  n=Inf
)

# ============================================================
# 7. DISEASE EFFECTS FOR ANATOMICAL INTERPRETATION
# ============================================================

disease_effects <- list()
k <- 1L

for(state in sort(unique(epi$L2))) {

  dd <- epi |>
    filter(
      L2==state
    )

  meds <- dd |>
    group_by(disease) |>
    summarise(
      median_frequency=
        median(frequency),
      .groups="drop"
    )

  getmed <- function(x) {
    z <- meds$median_frequency[
      meds$disease==x
    ]
    if(length(z)) z else NA_real_
  }

  aa <- dd |>
    filter(disease=="AA") |>
    pull(frequency)

  ffa <- dd |>
    filter(disease=="FFA") |>
    pull(frequency)

  wt <- suppressWarnings(
    wilcox.test(
      aa,
      ffa,
      exact=FALSE
    )
  )

  disease_effects[[k]] <- tibble(
    L2=state,
    median_HC=getmed("HC"),
    median_FFA=getmed("FFA"),
    median_AA=getmed("AA"),
    delta_AA_FFA=
      getmed("AA")-
      getmed("FFA"),
    AA_FFA_p=
      wt$p.value
  )

  k <- k+1L
}

disease_effects <- bind_rows(
  disease_effects
) |>
  mutate(
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
    "22G2A_epithelial_AA_vs_FFA_effects.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 8. INTEGRATED INTERPRETATION TABLE
# ============================================================

integrated <- primary |>
  left_join(
    disease_effects,
    by="L2"
  ) |>
  mutate(
    global_coupling=
      case_when(
        spearman_FDR < 0.05 ~
          "STRONG",
        spearman_FDR < 0.10 ~
          "SUPPORTED",
        abs(spearman_rho) >= 0.5 ~
          "EXPLORATORY",
        TRUE ~
          "WEAK"
      ),

    disease_adjusted_coupling=
      case_when(
        M1_adjusted_FDR < 0.05 ~
          "STRONG",
        M1_adjusted_FDR < 0.10 ~
          "SUPPORTED",
        M1_adjusted_p < 0.05 ~
          "NOMINAL",
        TRUE ~
          "NOT_SUPPORTED"
      ),

    interpretation=
      case_when(
        M1_adjusted_p < 0.05 ~
          "trajectory-associated beyond disease",

        spearman_p < 0.05 &
          M1_adjusted_p >= 0.05 ~
          "predominantly disease-structured co-remodeling",

        abs(spearman_rho) >= 0.5 ~
          "exploratory co-remodeling",

        TRUE ~
          "limited evidence"
      )
  )

write.csv(
  integrated,
  file.path(
    TAB,
    "22G2A_integrated_epithelial_propagation.csv"
  ),
  row.names=FALSE
)

cat("\n===== INTEGRATED EPITHELIAL PROPAGATION =====\n")

print(
  integrated |>
    dplyr::select(
      L2,
      spearman_rho,
      spearman_FDR,
      M1_adjusted_beta,
      M1_adjusted_p,
      M1_adjusted_FDR,
      interaction_p,
      delta_AA_FFA,
      AA_FFA_FDR,
      interpretation
    ),
  n=Inf
)

# ============================================================
# 9. FREEZE
# ============================================================

freeze <- list(
  module="22G2A",

  statistical_unit="patient",

  n_patients=16,

  primary_predictor=
    "Bulge median_position",

  outcome=
    "continuity-corrected logit within-epithelial L2 frequency",

  primary=primary,

  functional=functional,

  disease_effects=disease_effects,

  integrated=integrated,

  support=support,

  bulge_patient_data=bulge
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G2A_Bulge_epithelial_propagation_freeze.rds"
  )
)

# ============================================================
# 10. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "16 patients",
    "8 external epithelial L2",
    "Bulge excluded as outcome",
    "Primary table = 8 rows",
    "Functional table = 32 rows",
    "All primary models complete",
    "Freeze exists"
  ),

  pass=c(
    n_distinct(
      bulge$patient_id
    )==16,

    n_distinct(
      epi$L2
    )==8,

    !any(
      epi$L2==
        "Bulge_Isthmus"
    ),

    nrow(primary)==8,

    nrow(functional)==32,

    all(
      is.finite(
        primary$M1_adjusted_beta
      )
    ),

    file.exists(
      file.path(
        OBJ,
        "22G2A_Bulge_epithelial_propagation_freeze.rds"
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
    "22G2A_formal_checks.csv"
  ),
  row.names=FALSE
)

if(!all(checks$pass)) {
  stop(
    "22G2A formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G2A COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

