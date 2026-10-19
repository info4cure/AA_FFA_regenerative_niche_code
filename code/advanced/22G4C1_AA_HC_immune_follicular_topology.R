
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
})

cat("\n============================================================\n")
cat("22G4C1 — AA vs HC IMMUNE–FOLLICULAR TOPOLOGY\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJ,recursive=TRUE,showWarnings=FALSE)
dir.create(AUD,recursive=TRUE,showWarnings=FALSE)

L2_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

IMMUNE_ELIGIBLE <- file.path(
  TAB,
  "22G4C_AA_HC_eligible_immune_states.csv"
)

TARGET_PANEL <- file.path(
  TAB,
  "22G4C_lower_follicle_target_panel.csv"
)

BULGE_FILE <- file.path(
  TAB,
  "22G4C_AA_HC_Bulge_comparator.csv"
)

stopifnot(
  file.exists(L2_FILE),
  file.exists(IMMUNE_ELIGIBLE),
  file.exists(TARGET_PANEL),
  file.exists(BULGE_FILE)
)

# ============================================================
# 1. LOAD AA + HC
# ============================================================

l2 <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  filter(
    disease %in% c("HC","AA")
  ) |>
  mutate(
    disease=factor(
      disease,
      levels=c("HC","AA")
    )
  )

immune_def <- read.csv(
  IMMUNE_ELIGIBLE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

targets <- read.csv(
  TARGET_PANEL,
  stringsAsFactors=FALSE
) |>
  as_tibble()

bulge <- read.csv(
  BULGE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
  mutate(
    disease=factor(
      disease,
      levels=c("HC","AA")
    )
  )

cat("\nPatients:",n_distinct(l2$patient_id),"\n")
cat("Eligible immune states:",nrow(immune_def),"\n")
cat("Structural targets:",nrow(targets),"\n")

# ============================================================
# 2. IMMUNE MATRIX
# ============================================================

immune <- l2 |>
  inner_join(
    immune_def |>
      dplyr::select(
        L1,
        L2,
        feature
      ),
    by=c("L1","L2"),
    relationship="many-to-many"
  ) |>
  mutate(
    immune_p=
      (n+0.5)/
      (denominator_L1+1),

    immune_logit=
      qlogis(immune_p)
  ) |>
  dplyr::select(
    patient_id,
    disease,
    immune_L1=L1,
    immune_L2=L2,
    immune_feature=feature,
    immune_n=n,
    immune_frequency=frequency,
    immune_logit
  )

# ============================================================
# 3. STRUCTURAL TARGET MATRIX
# ============================================================

structure_long <- targets |>
  inner_join(
    l2,
    by=c("L1","L2"),
    relationship="many-to-many"
  ) |>
  mutate(
    target_p=
      (n+0.5)/
      (denominator_L1+1),

    target_logit=
      qlogis(target_p)
  ) |>
  dplyr::select(
    patient_id,
    disease,
    anatomical_axis,
    role,
    target_L1=L1,
    target_L2=L2,
    target_variable=variable,
    target_n=n,
    target_frequency=frequency,
    target_logit
  )

stopifnot(
  n_distinct(
    structure_long$target_variable
  )==9
)

# ============================================================
# 4. ALL 81 IMMUNE × STRUCTURAL PAIRS
# ============================================================

pair_def <- tidyr::crossing(
  immune_feature=
    unique(immune$immune_feature),

  target_variable=
    unique(structure_long$target_variable)
)

stopifnot(
  nrow(pair_def)==81
)

cat("\nStructural immune-target pairs:",nrow(pair_def),"\n")

# ============================================================
# 5. MODEL FUNCTION
# ============================================================

fit_pair <- function(
  immune_name,
  target_name
) {

  ii <- immune |>
    filter(
      immune_feature==
        immune_name
    )

  tt <- structure_long |>
    filter(
      target_variable==
        target_name
    )

  dd <- ii |>
    inner_join(
      tt,
      by=c(
        "patient_id",
        "disease"
      )
    ) |>
    filter(
      is.finite(immune_logit),
      is.finite(target_logit)
    )

  stopifnot(
    nrow(dd)==11
  )

  # ----------------------------------------------------------
  # Disease-specific Spearman correlations
  # ----------------------------------------------------------

  rho_HC <- suppressWarnings(
    cor(
      dd$immune_frequency[
        dd$disease=="HC"
      ],
      dd$target_frequency[
        dd$disease=="HC"
      ],
      method="spearman"
    )
  )

  rho_AA <- suppressWarnings(
    cor(
      dd$immune_frequency[
        dd$disease=="AA"
      ],
      dd$target_frequency[
        dd$disease=="AA"
      ],
      method="spearman"
    )
  )

  # ----------------------------------------------------------
  # Pooled descriptive correlation
  # ----------------------------------------------------------

  sp <- suppressWarnings(
    cor.test(
      dd$immune_frequency,
      dd$target_frequency,
      method="spearman",
      exact=FALSE
    )
  )

  # ----------------------------------------------------------
  # Shared coupling model
  # target ~ immune + disease
  # ----------------------------------------------------------

  m_shared <- lm(
    target_logit ~
      immune_logit +
      disease,
    data=dd
  )

  shared_term <- broom::tidy(
    m_shared
  ) |>
    filter(
      term=="immune_logit"
    )

  # ----------------------------------------------------------
  # AA-specific interaction model
  #
  # HC is reference:
  # immune_logit coefficient = HC slope
  # immune_logit:diseaseAA = change in slope in AA
  # ----------------------------------------------------------

  m_interaction <- lm(
    target_logit ~
      immune_logit *
      disease,
    data=dd
  )

  tt_int <- broom::tidy(
    m_interaction
  )

  HC_term <- tt_int |>
    filter(
      term=="immune_logit"
    )

  interaction_term <- tt_int |>
    filter(
      term=="immune_logit:diseaseAA"
    )

  beta_HC <-
    HC_term$estimate

  beta_interaction <-
    interaction_term$estimate

  beta_AA <-
    beta_HC +
    beta_interaction

  # AA slope significance directly by releveling disease.
  dd_AAref <- dd |>
    mutate(
      disease_AAref=
        relevel(
          disease,
          ref="AA"
        )
    )

  m_AAref <- lm(
    target_logit ~
      immune_logit *
      disease_AAref,
    data=dd_AAref
  )

  AA_term <- broom::tidy(
    m_AAref
  ) |>
    filter(
      term=="immune_logit"
    )

  meta <- tt |>
    distinct(
      anatomical_axis,
      role,
      target_L1,
      target_L2,
      target_variable
    )

  imeta <- ii |>
    distinct(
      immune_L1,
      immune_L2,
      immune_feature
    )

  tibble(
    immune_L1=
      imeta$immune_L1,

    immune_L2=
      imeta$immune_L2,

    immune_feature=
      immune_name,

    anatomical_axis=
      meta$anatomical_axis,

    role=
      meta$role,

    target_L1=
      meta$target_L1,

    target_L2=
      meta$target_L2,

    target_variable=
      target_name,

    n=11,

    rho_pooled=
      unname(sp$estimate),

    pooled_p=
      sp$p.value,

    rho_HC=
      rho_HC,

    rho_AA=
      rho_AA,

    delta_rho_AA_HC=
      rho_AA-rho_HC,

    shared_beta=
      shared_term$estimate,

    shared_SE=
      shared_term$std.error,

    shared_p=
      shared_term$p.value,

    beta_HC=
      beta_HC,

    beta_AA=
      beta_AA,

    AA_slope_SE=
      AA_term$std.error,

    AA_slope_p=
      AA_term$p.value,

    interaction_beta=
      beta_interaction,

    interaction_SE=
      interaction_term$std.error,

    interaction_p=
      interaction_term$p.value,

    R2_shared=
      summary(m_shared)$r.squared,

    R2_interaction=
      summary(m_interaction)$r.squared
  )
}

# ============================================================
# 6. FIT STRUCTURAL NETWORK
# ============================================================

cat("\n===== FIT 81 STRUCTURAL CONNECTIONS =====\n")

res <- vector(
  "list",
  nrow(pair_def)
)

for(i in seq_len(nrow(pair_def))) {

  res[[i]] <- fit_pair(
    pair_def$immune_feature[i],
    pair_def$target_variable[i]
  )
}

structural <- bind_rows(res)

stopifnot(
  nrow(structural)==81
)

# ============================================================
# 7. FDR BY ANATOMICAL HYPOTHESIS FAMILY
# ============================================================

structural <- structural |>
  group_by(
    anatomical_axis
  ) |>
  mutate(
    pooled_FDR=
      p.adjust(
        pooled_p,
        method="BH"
      ),

    shared_FDR=
      p.adjust(
        shared_p,
        method="BH"
      ),

    AA_slope_FDR=
      p.adjust(
        AA_slope_p,
        method="BH"
      ),

    interaction_FDR=
      p.adjust(
        interaction_p,
        method="BH"
      )
  ) |>
  ungroup()

# ============================================================
# 8. CLASSIFY TOPOLOGY
# ============================================================

structural <- structural |>
  mutate(
    topology_class=case_when(

      interaction_FDR < 0.10 &
        AA_slope_p < 0.05
        ~ "AA_SPECIFIC_SUPPORTED",

      interaction_p < 0.05 &
        AA_slope_p < 0.05
        ~ "AA_SPECIFIC_NOMINAL",

      AA_slope_FDR < 0.10 &
        abs(rho_AA) >= 0.50
        ~ "AA_COUPLING_SUPPORTED",

      AA_slope_p < 0.05 &
        abs(rho_AA) >= 0.50
        ~ "AA_COUPLING_NOMINAL",

      shared_FDR < 0.10
        ~ "SHARED_SUPPORTED",

      shared_p < 0.05
        ~ "SHARED_NOMINAL",

      abs(rho_AA) >= 0.70
        ~ "AA_EXPLORATORY_STRONG",

      TRUE
        ~ "UNSUPPORTED"
    ),

    AA_direction=case_when(
      beta_AA > 0
        ~ "POSITIVE",
      beta_AA < 0
        ~ "NEGATIVE",
      TRUE
        ~ "NEUTRAL"
    )
  )

# ============================================================
# 9. BULGE COMPARATOR — 45 CONNECTIONS
# ============================================================

bulge_features <- c(
  "median_position",
  "stem_niche_mean",
  "MHCII_mean",
  "cycling_mean",
  "CD200_z_mean"
)

bulge_results <- list()

k <- 1L

for(im in unique(immune$immune_feature)) {

  ii <- immune |>
    filter(
      immune_feature==im
    )

  for(v in bulge_features) {

    dd <- ii |>
      inner_join(
        bulge |>
          dplyr::select(
            patient_id,
            disease,
            all_of(v)
          ),
        by=c(
          "patient_id",
          "disease"
        )
      ) |>
      rename(
        outcome=
          all_of(v)
      )

    rho_HC <- suppressWarnings(
      cor(
        dd$immune_frequency[
          dd$disease=="HC"
        ],
        dd$outcome[
          dd$disease=="HC"
        ],
        method="spearman"
      )
    )

    rho_AA <- suppressWarnings(
      cor(
        dd$immune_frequency[
          dd$disease=="AA"
        ],
        dd$outcome[
          dd$disease=="AA"
        ],
        method="spearman"
      )
    )

    m1 <- lm(
      outcome ~
        immune_logit +
        disease,
      data=dd
    )

    m2 <- lm(
      outcome ~
        immune_logit *
        disease,
      data=dd
    )

    t1 <- broom::tidy(m1) |>
      filter(
        term=="immune_logit"
      )

    t2 <- broom::tidy(m2)

    hc <- t2 |>
      filter(
        term=="immune_logit"
      )

    it <- t2 |>
      filter(
        term=="immune_logit:diseaseAA"
      )

    dd2 <- dd |>
      mutate(
        disease_AAref=
          relevel(
            disease,
            ref="AA"
          )
      )

    m3 <- lm(
      outcome ~
        immune_logit *
        disease_AAref,
      data=dd2
    )

    aa <- broom::tidy(m3) |>
      filter(
        term=="immune_logit"
      )

    bulge_results[[k]] <- tibble(
      immune_feature=im,
      bulge_feature=v,
      rho_HC=rho_HC,
      rho_AA=rho_AA,
      shared_beta=t1$estimate,
      shared_p=t1$p.value,
      beta_HC=hc$estimate,
      beta_AA=aa$estimate,
      AA_slope_p=aa$p.value,
      interaction_beta=it$estimate,
      interaction_p=it$p.value
    )

    k <- k+1L
  }
}

bulge_coupling <- bind_rows(
  bulge_results
) |>
  group_by(
    bulge_feature
  ) |>
  mutate(
    shared_FDR=
      p.adjust(
        shared_p,
        method="BH"
      ),

    AA_slope_FDR=
      p.adjust(
        AA_slope_p,
        method="BH"
      ),

    interaction_FDR=
      p.adjust(
        interaction_p,
        method="BH"
      )
  ) |>
  ungroup()

stopifnot(
  nrow(bulge_coupling)==45
)

# ============================================================
# 10. SAVE
# ============================================================

write.csv(
  structural,
  file.path(
    TAB,
    "22G4C1_ALL_AA_HC_immune_follicular_connections.csv"
  ),
  row.names=FALSE
)

prioritized <- structural |>
  filter(
    topology_class !=
      "UNSUPPORTED"
  ) |>
  arrange(
    factor(
      topology_class,
      levels=c(
        "AA_SPECIFIC_SUPPORTED",
        "AA_COUPLING_SUPPORTED",
        "AA_SPECIFIC_NOMINAL",
        "AA_COUPLING_NOMINAL",
        "SHARED_SUPPORTED",
        "SHARED_NOMINAL",
        "AA_EXPLORATORY_STRONG"
      )
    ),
    interaction_FDR,
    AA_slope_FDR
  )

write.csv(
  prioritized,
  file.path(
    TAB,
    "22G4C1_prioritized_AA_HC_connections.csv"
  ),
  row.names=FALSE
)

write.csv(
  bulge_coupling,
  file.path(
    TAB,
    "22G4C1_AA_HC_immune_Bulge_comparator.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 11. REPORT
# ============================================================

cat("\n===== TOPOLOGY CLASS COUNTS =====\n")

print(
  structural |>
    count(
      anatomical_axis,
      topology_class
    ) |>
    arrange(
      anatomical_axis,
      topology_class
    ),
  n=Inf
)

cat("\n===== PRIORITIZED AA/HC CONNECTIONS =====\n")

print(
  prioritized |>
    dplyr::select(
      immune_L1,
      immune_L2,
      anatomical_axis,
      target_L2,
      rho_HC,
      rho_AA,
      delta_rho_AA_HC,
      shared_beta,
      shared_p,
      shared_FDR,
      beta_HC,
      beta_AA,
      AA_slope_p,
      AA_slope_FDR,
      interaction_beta,
      interaction_p,
      interaction_FDR,
      topology_class,
      AA_direction
    ),
  n=Inf
)

cat("\n===== AA-SPECIFIC / AA-COUPLED LOWER-FOLLICLE EDGES =====\n")

lower_AA <- prioritized |>
  filter(
    anatomical_axis %in%
      c(
        "LOWER_FOLLICLE",
        "LOWER_MESENCHYME"
      ),
    topology_class %in%
      c(
        "AA_SPECIFIC_SUPPORTED",
        "AA_SPECIFIC_NOMINAL",
        "AA_COUPLING_SUPPORTED",
        "AA_COUPLING_NOMINAL",
        "AA_EXPLORATORY_STRONG"
      )
  ) |>
  arrange(
    factor(
      topology_class,
      levels=c(
        "AA_SPECIFIC_SUPPORTED",
        "AA_COUPLING_SUPPORTED",
        "AA_SPECIFIC_NOMINAL",
        "AA_COUPLING_NOMINAL",
        "AA_EXPLORATORY_STRONG"
      )
    ),
    interaction_FDR,
    AA_slope_FDR,
    desc(
      abs(rho_AA)
    )
  )

print(
  lower_AA |>
    dplyr::select(
      immune_L1,
      immune_L2,
      anatomical_axis,
      target_L2,
      rho_HC,
      rho_AA,
      delta_rho_AA_HC,
      beta_HC,
      beta_AA,
      AA_slope_p,
      AA_slope_FDR,
      interaction_p,
      interaction_FDR,
      topology_class
    ),
  n=Inf
)

write.csv(
  lower_AA,
  file.path(
    TAB,
    "22G4C1_AA_lower_follicle_candidate_edges.csv"
  ),
  row.names=FALSE
)

cat("\n===== UPPER-COMPARATOR EDGES =====\n")

upper_AA <- prioritized |>
  filter(
    anatomical_axis==
      "UPPER_COMPARATOR"
  )

print(
  upper_AA |>
    dplyr::select(
      immune_L1,
      immune_L2,
      target_L2,
      rho_HC,
      rho_AA,
      delta_rho_AA_HC,
      beta_HC,
      beta_AA,
      AA_slope_p,
      AA_slope_FDR,
      interaction_p,
      interaction_FDR,
      topology_class
    ),
  n=Inf
)

write.csv(
  upper_AA,
  file.path(
    TAB,
    "22G4C1_AA_upper_comparator_edges.csv"
  ),
  row.names=FALSE
)

cat("\n===== BULGE COMPARATOR — TOP RESULTS =====\n")

bulge_prioritized <- bulge_coupling |>
  mutate(
    bulge_class=case_when(

      interaction_FDR < 0.10 &
        AA_slope_p < 0.05
        ~ "AA_SPECIFIC_SUPPORTED",

      interaction_p < 0.05 &
        AA_slope_p < 0.05
        ~ "AA_SPECIFIC_NOMINAL",

      AA_slope_FDR < 0.10 &
        abs(rho_AA) >= 0.50
        ~ "AA_COUPLING_SUPPORTED",

      AA_slope_p < 0.05 &
        abs(rho_AA) >= 0.50
        ~ "AA_COUPLING_NOMINAL",

      shared_FDR < 0.10
        ~ "SHARED_SUPPORTED",

      shared_p < 0.05
        ~ "SHARED_NOMINAL",

      abs(rho_AA) >= 0.70
        ~ "AA_EXPLORATORY_STRONG",

      TRUE
        ~ "UNSUPPORTED"
    )
  ) |>
  arrange(
    factor(
      bulge_class,
      levels=c(
        "AA_SPECIFIC_SUPPORTED",
        "AA_COUPLING_SUPPORTED",
        "AA_SPECIFIC_NOMINAL",
        "AA_COUPLING_NOMINAL",
        "SHARED_SUPPORTED",
        "SHARED_NOMINAL",
        "AA_EXPLORATORY_STRONG",
        "UNSUPPORTED"
      )
    ),
    interaction_FDR,
    AA_slope_FDR
  )

print(
  bulge_prioritized |>
    filter(
      bulge_class !=
        "UNSUPPORTED"
    ),
  n=Inf
)

write.csv(
  bulge_prioritized,
  file.path(
    TAB,
    "22G4C1_AA_HC_immune_Bulge_comparator_classified.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 12. DISEASE EFFECTS ON TARGETS
# ============================================================

cat("\n===== AA vs HC TARGET DISEASE EFFECTS =====\n")

target_disease_effects <- structure_long |>
  group_by(
    anatomical_axis,
    role,
    target_L1,
    target_L2,
    target_variable
  ) |>
  group_modify(
    ~{
      dd <- .x

      wt <- suppressWarnings(
        wilcox.test(
          target_frequency ~
            disease,
          data=dd,
          exact=FALSE
        )
      )

      tibble(
        median_HC=
          median(
            dd$target_frequency[
              dd$disease=="HC"
            ],
            na.rm=TRUE
          ),

        median_AA=
          median(
            dd$target_frequency[
              dd$disease=="AA"
            ],
            na.rm=TRUE
          ),

        delta_AA_HC=
          median_AA-
          median_HC,

        AA_HC_p=
          wt$p.value
      )
    }
  ) |>
  ungroup() |>
  group_by(
    anatomical_axis
  ) |>
  mutate(
    AA_HC_FDR=
      p.adjust(
        AA_HC_p,
        method="BH"
      )
  ) |>
  ungroup() |>
  arrange(
    anatomical_axis,
    AA_HC_FDR
  )

print(
  target_disease_effects,
  n=Inf
)

write.csv(
  target_disease_effects,
  file.path(
    TAB,
    "22G4C1_AA_HC_target_disease_effects.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 13. IMMUNE DISEASE EFFECTS
# ============================================================

cat("\n===== AA vs HC IMMUNE DISEASE EFFECTS =====\n")

immune_disease_effects <- immune |>
  group_by(
    immune_L1,
    immune_L2,
    immune_feature
  ) |>
  group_modify(
    ~{
      dd <- .x

      wt <- suppressWarnings(
        wilcox.test(
          immune_frequency ~
            disease,
          data=dd,
          exact=FALSE
        )
      )

      tibble(
        median_HC=
          median(
            dd$immune_frequency[
              dd$disease=="HC"
            ],
            na.rm=TRUE
          ),

        median_AA=
          median(
            dd$immune_frequency[
              dd$disease=="AA"
            ],
            na.rm=TRUE
          ),

        delta_AA_HC=
          median_AA-
          median_HC,

        AA_HC_p=
          wt$p.value
      )
    }
  ) |>
  ungroup() |>
  mutate(
    AA_HC_FDR=
      p.adjust(
        AA_HC_p,
        method="BH"
      )
  ) |>
  arrange(
    AA_HC_FDR
  )

print(
  immune_disease_effects,
  n=Inf
)

write.csv(
  immune_disease_effects,
  file.path(
    TAB,
    "22G4C1_AA_HC_immune_disease_effects.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 14. FREEZE
# ============================================================

freeze <- list(
  module="22G4C1",

  purpose=paste(
    "AA-versus-HC patient-level immune-follicular topology",
    "with explicit lower-follicle, lower-mesenchymal,",
    "upper-follicular comparator and Bulge analyses."
  ),

  statistical_unit="patient",

  disease_universe=c(
    "HC",
    "AA"
  ),

  n_HC=5,
  n_AA=6,

  n_immune_states=
    n_distinct(
      immune$immune_feature
    ),

  n_structural_targets=
    n_distinct(
      structure_long$target_variable
    ),

  n_structural_connections=
    nrow(
      structural
    ),

  n_Bulge_connections=
    nrow(
      bulge_coupling
    ),

  structural_connections=
    structural,

  prioritized_connections=
    prioritized,

  lower_follicle_candidates=
    lower_AA,

  upper_comparator_candidates=
    upper_AA,

  Bulge_comparator=
    bulge_prioritized,

  target_disease_effects=
    target_disease_effects,

  immune_disease_effects=
    immune_disease_effects,

  interpretation_constraint=paste(
    "Patient-level association does not establish anatomical",
    "proximity, direct cell-cell interaction, mediation or",
    "causal direction."
  )
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4C1_AA_HC_immune_follicular_topology_freeze.rds"
  )
)

# ============================================================
# 15. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "11 AA+HC patients",
    "9 eligible immune states",
    "9 structural targets",
    "81 structural connections",
    "45 Bulge comparator connections",
    "Lower follicle represented",
    "Lower mesenchyme represented",
    "Upper comparator represented",
    "Target disease effects generated",
    "Immune disease effects generated",
    "Full structural table exists",
    "Prioritized table exists",
    "Freeze exists"
  ),

  pass=c(
    n_distinct(
      l2$patient_id
    )==11,

    n_distinct(
      immune$immune_feature
    )==9,

    n_distinct(
      structure_long$target_variable
    )==9,

    nrow(
      structural
    )==81,

    nrow(
      bulge_coupling
    )==45,

    "LOWER_FOLLICLE" %in%
      structural$anatomical_axis,

    "LOWER_MESENCHYME" %in%
      structural$anatomical_axis,

    "UPPER_COMPARATOR" %in%
      structural$anatomical_axis,

    nrow(
      target_disease_effects
    )==9,

    nrow(
      immune_disease_effects
    )==9,

    file.exists(
      file.path(
        TAB,
        "22G4C1_ALL_AA_HC_immune_follicular_connections.csv"
      )
    ),

    file.exists(
      file.path(
        TAB,
        "22G4C1_prioritized_AA_HC_connections.csv"
      )
    ),

    file.exists(
      file.path(
        OBJ,
        "22G4C1_AA_HC_immune_follicular_topology_freeze.rds"
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
    "22G4C1_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 16. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C1 — AA vs HC IMMUNE–FOLLICULAR TOPOLOGY",
  "============================================================",
  "",
  "FFA excluded from all analyses.",
  "HC patients: 5",
  "AA patients: 6",
  "",
  paste0(
    "Immune states: ",
    n_distinct(
      immune$immune_feature
    )
  ),
  paste0(
    "Structural targets: ",
    n_distinct(
      structure_long$target_variable
    )
  ),
  paste0(
    "Structural connections: ",
    nrow(
      structural
    )
  ),
  paste0(
    "Bulge comparator connections: ",
    nrow(
      bulge_coupling
    )
  ),
  "",
  "Primary anatomical hypotheses:",
  "- LOWER_FOLLICLE: ORS, ORS_basal, IRS_companion, Hair_shaft",
  "- LOWER_MESENCHYME: COCH/AQP4 dermal papilla, COL11A1/DPEP1 dermal sheath",
  "- UPPER_COMPARATOR: three IFE/Infundibulum states",
  "- BULGE: comparator rather than assumed AA inflammatory epicentre",
  "",
  "Primary topology test:",
  "immune x disease interaction with HC as reference.",
  "",
  "Interpretation:",
  "AA-specific coupling requires evidence that the immune-target",
  "relationship differs between AA and HC; strong AA-only",
  "correlations without interaction support remain exploratory.",
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
    "22G4C1_audit.txt"
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
# 17. FINAL STATUS
# ============================================================

if(!all(checks$pass)) {
  stop(
    "22G4C1 formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C1 COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

