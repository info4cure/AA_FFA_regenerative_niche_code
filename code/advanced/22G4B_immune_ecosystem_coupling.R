
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(broom)
})

cat("\n============================================================\n")
cat("22G4B — IMMUNE × ECOSYSTEM COUPLING ATLAS\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4B_ECOSYSTEM_COUPLING"
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

MASTER_FILE <- file.path(
  TAB,
  "22G4B_ecosystem_patient_master.csv"
)

IMMUNE_FILE <- file.path(
  TAB,
  "22G4B_immune_patient_source.csv"
)

MANIFEST_FILE <- file.path(
  TAB,
  "22G4B_ecosystem_feature_manifest.csv"
)

stopifnot(
  file.exists(MASTER_FILE),
  file.exists(IMMUNE_FILE),
  file.exists(MANIFEST_FILE)
)

# ============================================================
# 2. LOAD FROZEN PREFLIGHT INPUTS
# ============================================================

master <- read.csv(
  MASTER_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
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

immune <- read.csv(
  IMMUNE_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble() |>
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

manifest <- read.csv(
  MANIFEST_FILE,
  stringsAsFactors=FALSE
) |>
  as_tibble()

cat("\nPatients:",n_distinct(master$patient_id),"\n")
cat("Ecosystem features:",nrow(manifest),"\n")
cat("Immune states:",n_distinct(immune$feature),"\n")

# ============================================================
# 3. DEFINE OUTCOME TRANSFORMATION
# ============================================================

# Bulge features are continuous frozen scores and remain
# on their original scales.
#
# Epithelial / niche features are within-L1 proportions.
# They are continuity-corrected before logit transformation.

proportion_features <- manifest |>
  filter(
    layer != "BULGE"
  ) |>
  pull(variable)

analysis_master <- master

for(v in proportion_features) {

  # Frequencies are already proportions in [0,1].
  # Use a small sample-size-independent correction only
  # for exact boundary values.

  x <- analysis_master[[v]]

  eps <- 0.5 / 10000

  x2 <- pmin(
    pmax(
      x,
      eps
    ),
    1-eps
  )

  analysis_master[[paste0(v,"__analysis")]] <-
    qlogis(x2)
}

for(v in manifest$variable[
  manifest$layer=="BULGE"
]) {

  analysis_master[[paste0(v,"__analysis")]] <-
    analysis_master[[v]
    ]
}

# ============================================================
# 4. IMMUNE PREDICTOR TRANSFORMATION
# ============================================================

# Use continuity-corrected logit of within-L1 immune frequency.

immune <- immune |>
  mutate(
    immune_p=
      (n+0.5)/
      (denominator_L1+1),

    immune_logit=
      qlogis(
        immune_p
      )
  )

# ============================================================
# 5. BUILD ALL 135 PAIRS
# ============================================================

immune_states <- immune |>
  distinct(
    L1,
    L2,
    feature
  ) |>
  arrange(
    L1,
    L2
  )

pair_manifest <- tidyr::crossing(
  immune_feature=
    immune_states$feature,

  ecosystem_variable=
    manifest$variable
) |>
  left_join(
    immune_states |>
      rename(
        immune_feature=feature,
        immune_L1=L1,
        immune_L2=L2
      ),
    by="immune_feature"
  ) |>
  left_join(
    manifest |>
      rename(
        ecosystem_variable=variable
      ),
    by="ecosystem_variable"
  )

stopifnot(
  nrow(pair_manifest)==135
)

cat("\nTotal immune × ecosystem pairs:",nrow(pair_manifest),"\n")

# ============================================================
# 6. GLOBAL + DISEASE-ADJUSTED MODELS
# ============================================================

cat("\n===== FIT 135 COUPLING MODELS =====\n")

results <- vector(
  "list",
  nrow(pair_manifest)
)

for(i in seq_len(nrow(pair_manifest))) {

  pm <- pair_manifest[i,]

  immune_i <- immune |>
    filter(
      feature==
        pm$immune_feature
    ) |>
    dplyr::select(
      patient_id,
      disease,
      immune_frequency=frequency,
      immune_logit
    )

  outcome_col <-
    paste0(
      pm$ecosystem_variable,
      "__analysis"
    )

  dd <- analysis_master |>
    dplyr::select(
      patient_id,
      disease,
      all_of(
        pm$ecosystem_variable
      ),
      all_of(
        outcome_col
      )
    ) |>
    rename(
      ecosystem_raw=
        all_of(
          pm$ecosystem_variable
        ),

      ecosystem_analysis=
        all_of(
          outcome_col
        )
    ) |>
    inner_join(
      immune_i,
      by=c(
        "patient_id",
        "disease"
      )
    ) |>
    filter(
      is.finite(ecosystem_analysis),
      is.finite(immune_logit)
    )

  # Spearman uses the untransformed biological quantities.
  sp <- suppressWarnings(
    cor.test(
      dd$immune_frequency,
      dd$ecosystem_raw,
      method="spearman",
      exact=FALSE
    )
  )

  # Unadjusted
  m0 <- lm(
    ecosystem_analysis ~
      immune_logit,
    data=dd
  )

  # Disease-adjusted primary model
  m1 <- lm(
    ecosystem_analysis ~
      immune_logit +
      disease,
    data=dd
  )

  # Exploratory disease interaction
  m2 <- lm(
    ecosystem_analysis ~
      immune_logit *
      disease,
    data=dd
  )

  t0 <- broom::tidy(m0) |>
    filter(
      term=="immune_logit"
    )

  t1 <- broom::tidy(m1) |>
    filter(
      term=="immune_logit"
    )

  cmp <- anova(
    m1,
    m2
  )

  interaction_p <-
    cmp$`Pr(>F)`[2]

  # Within-disease directional robustness
  rho_by_disease <- c(
    HC=NA_real_,
    FFA=NA_real_,
    AA=NA_real_
  )

  for(dis_i in names(rho_by_disease)) {

    ddis <- dd |>
      filter(
        disease==dis_i
      )

    if(
      nrow(ddis)>=4 &&
      sd(ddis$immune_frequency)>0 &&
      sd(ddis$ecosystem_raw)>0
    ) {

      rho_by_disease[dis_i] <-
        suppressWarnings(
          cor(
            ddis$immune_frequency,
            ddis$ecosystem_raw,
            method="spearman"
          )
        )
    }
  }

  results[[i]] <- tibble(
    immune_L1=
      pm$immune_L1,

    immune_L2=
      pm$immune_L2,

    immune_feature=
      pm$immune_feature,

    layer=
      pm$layer,

    ecosystem_variable=
      pm$ecosystem_variable,

    n=
      nrow(dd),

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

    adjusted_beta=
      t1$estimate,

    adjusted_SE=
      t1$std.error,

    adjusted_p=
      t1$p.value,

    interaction_p=
      interaction_p,

    rho_HC=
      rho_by_disease["HC"],

    rho_FFA=
      rho_by_disease["FFA"],

    rho_AA=
      rho_by_disease["AA"],

    R2_unadjusted=
      summary(m0)$r.squared,

    R2_adjusted=
      summary(m1)$r.squared
  )
}

coupling <- bind_rows(
  results
)

cat("\nRaw coupling rows:",nrow(coupling),"\n")



# ============================================================
# 7. MULTIPLE-TESTING CORRECTION BY BIOLOGICAL LAYER
# ============================================================

cat("\n===== MULTIPLE-TESTING CORRECTION BY LAYER =====\n")

coupling <- coupling |>
  group_by(
    layer
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

    adjusted_FDR=
      p.adjust(
        adjusted_p,
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
# 8. DIRECTIONAL ROBUSTNESS
# ============================================================

cat("\n===== DIRECTIONAL ROBUSTNESS =====\n")

coupling <- coupling |>
  rowwise() |>
  mutate(
    directional_support_n=
      sum(
        c(
          !is.na(rho_HC) &&
            sign(rho_HC)==sign(adjusted_beta),

          !is.na(rho_FFA) &&
            sign(rho_FFA)==sign(adjusted_beta),

          !is.na(rho_AA) &&
            sign(rho_AA)==sign(adjusted_beta)
        )
      )
  ) |>
  ungroup()

# ============================================================
# 9. CONNECTION CLASSIFICATION
# ============================================================

cat("\n===== CONNECTION CLASSIFICATION =====\n")

coupling <- coupling |>
  mutate(
    connection_class=case_when(

      adjusted_FDR < 0.05
        ~ "CORE",

      adjusted_FDR < 0.10
        ~ "SUPPORTED",

      adjusted_p < 0.05
        ~ "NOMINAL",

      spearman_FDR < 0.10 &
        abs(spearman_rho) >= 0.50
        ~ "GLOBAL_SUPPORTED",

      abs(spearman_rho) >= 0.50
        ~ "EXPLORATORY",

      TRUE
        ~ "UNSUPPORTED"
    ),

    direction=case_when(
      adjusted_beta > 0
        ~ "POSITIVE",

      adjusted_beta < 0
        ~ "NEGATIVE",

      TRUE
        ~ "NEUTRAL"
    ),

    interaction_class=case_when(
      interaction_FDR < 0.10
        ~ "SUPPORTED",

      interaction_p < 0.05
        ~ "NOMINAL",

      TRUE
        ~ "NONE"
    ),

    evidence_score=
      abs(spearman_rho) +

      0.50*(
        adjusted_FDR < 0.05
      ) +

      0.35*(
        adjusted_FDR < 0.10
      ) +

      0.20*(
        adjusted_p < 0.05
      ) +

      0.15*(
        directional_support_n >= 2
      )
  ) |>
  arrange(
    factor(
      connection_class,
      levels=c(
        "CORE",
        "SUPPORTED",
        "NOMINAL",
        "GLOBAL_SUPPORTED",
        "EXPLORATORY",
        "UNSUPPORTED"
      )
    ),
    desc(
      evidence_score
    )
  )

# ============================================================
# 10. SAVE FULL ATLAS
# ============================================================

cat("\n===== SAVE FULL COUPLING ATLAS =====\n")

write.csv(
  coupling,
  file.path(
    TAB,
    "22G4B_ALL_immune_ecosystem_connections.csv"
  ),
  row.names=FALSE
)

supported <- coupling |>
  filter(
    connection_class !=
      "UNSUPPORTED"
  )

write.csv(
  supported,
  file.path(
    TAB,
    "22G4B_supported_immune_ecosystem_connections.csv"
  ),
  row.names=FALSE
)

cat("\n===== CONNECTION CLASS COUNTS =====\n")

print(
  coupling |>
    count(
      layer,
      connection_class
    ) |>
    arrange(
      layer,
      factor(
        connection_class,
        levels=c(
          "CORE",
          "SUPPORTED",
          "NOMINAL",
          "GLOBAL_SUPPORTED",
          "EXPLORATORY",
          "UNSUPPORTED"
        )
      )
    ),
  n=Inf
)

cat("\n===== SUPPORTED / PRIORITIZED CONNECTIONS =====\n")

print(
  supported |>
    dplyr::select(
      immune_L1,
      immune_L2,
      layer,
      ecosystem_variable,
      spearman_rho,
      spearman_FDR,
      adjusted_beta,
      adjusted_p,
      adjusted_FDR,
      interaction_p,
      interaction_FDR,
      rho_HC,
      rho_FFA,
      rho_AA,
      directional_support_n,
      connection_class
    ),
  n=Inf
)


# ============================================================
# 11. LAYER-SPECIFIC LEDGERS
# ============================================================

for(layer_i in unique(coupling$layer)) {

  z <- coupling |>
    filter(
      layer==layer_i
    ) |>
    arrange(
      adjusted_FDR,
      desc(abs(spearman_rho))
    )

  stem <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    layer_i
  )

  write.csv(
    z,
    file.path(
      TAB,
      paste0(
        "22G4B_",
        stem,
        "_connections.csv"
      )
    ),
    row.names=FALSE
  )
}

# ============================================================
# 12. CONNECTIONS TO ALREADY-FROZEN OUTWARD MODEL
# ============================================================

# These are the biologically central nodes established in
# previous frozen modules. No statistical preference is added;
# this is simply a convenient paper-oriented extraction.

central_nodes <- c(
  "median_position",
  "MHCII_mean",
  "stem_niche_mean",
  "cycling_mean",
  "CD200_z_mean",
  "ORS",
  "IRS_companion",
  "Hair_shaft",
  "PTGDS_IGF1",
  "COL11A1_DPEP1",
  "APOE_CCL19",
  "Vascular_epithelial_like"
)

paper_network <- coupling |>
  filter(
    ecosystem_variable %in%
      central_nodes,
    connection_class !=
      "UNSUPPORTED"
  ) |>
  arrange(
    factor(
      connection_class,
      levels=c(
        "CORE",
        "SUPPORTED",
        "NOMINAL",
        "GLOBAL_SUPPORTED",
        "EXPLORATORY"
      )
    ),
    desc(evidence_score)
  )

write.csv(
  paper_network,
  file.path(
    TAB,
    "22G4B_paper_network_candidate_edges.csv"
  ),
  row.names=FALSE
)

cat("\n===== PAPER NETWORK CANDIDATE EDGES =====\n")

print(
  paper_network |>
    dplyr::select(
      immune_L2,
      layer,
      ecosystem_variable,
      adjusted_beta,
      adjusted_FDR,
      rho_HC,
      rho_FFA,
      rho_AA,
      connection_class
    ),
  n=Inf
)

# ============================================================
# 13. FREEZE
# ============================================================

freeze <- list(
  module="22G4B",

  purpose=paste(
    "Patient-level immune-centered coupling atlas spanning",
    "Bulge state, follicular epithelium, superficial epithelium",
    "and the mesenchymal/vascular niche."
  ),

  statistical_unit="patient",

  n_patients=
    n_distinct(
      master$patient_id
    ),

  n_immune_states=
    n_distinct(
      immune$feature
    ),

  n_ecosystem_features=
    nrow(
      manifest
    ),

  n_connections=
    nrow(
      coupling
    ),

  FDR_strategy=
    "Benjamini-Hochberg correction within biological layer",

  interpretation_constraint=paste(
    "Associations represent patient-level co-variation.",
    "Directness, mediation and causal direction are not inferred."
  ),

  pair_manifest=
    pair_manifest,

  coupling=
    coupling,

  supported=
    supported,

  paper_network=
    paper_network
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4B_immune_ecosystem_coupling_freeze.rds"
  )
)

# ============================================================
# 14. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "16 patients",
    "9 immune states",
    "15 ecosystem features",
    "135 immune-ecosystem pairs",
    "135 fitted models",
    "Four biological layers",
    "FDR available for all models",
    "Freeze exists",
    "Full connection table exists",
    "Supported connection table exists",
    "Paper-network edge table exists"
  ),

  pass=c(
    n_distinct(
      master$patient_id
    ) == 16,

    n_distinct(
      immune$feature
    ) == 9,

    nrow(
      manifest
    ) == 15,

    nrow(
      pair_manifest
    ) == 135,

    nrow(
      coupling
    ) == 135,

    setequal(
      unique(
        coupling$layer
      ),
      c(
        "BULGE",
        "FOLLICULAR_EPITHELIUM",
        "SUPERFICIAL_EPITHELIUM",
        "MESENCHYMAL_NICHE"
      )
    ),

    all(
      is.finite(
        coupling$spearman_FDR
      )
    ) &&
      all(
        is.finite(
          coupling$adjusted_FDR
        )
      ) &&
      all(
        is.finite(
          coupling$interaction_FDR
        )
      ),

    file.exists(
      file.path(
        OBJ,
        "22G4B_immune_ecosystem_coupling_freeze.rds"
      )
    ),

    file.exists(
      file.path(
        TAB,
        "22G4B_ALL_immune_ecosystem_connections.csv"
      )
    ),

    file.exists(
      file.path(
        TAB,
        "22G4B_supported_immune_ecosystem_connections.csv"
      )
    ),

    file.exists(
      file.path(
        TAB,
        "22G4B_paper_network_candidate_edges.csv"
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
    "22G4B_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 15. SUMMARY LEDGER
# ============================================================

class_summary <- coupling |>
  count(
    layer,
    connection_class,
    name="n_connections"
  ) |>
  arrange(
    layer,
    factor(
      connection_class,
      levels=c(
        "CORE",
        "SUPPORTED",
        "NOMINAL",
        "GLOBAL_SUPPORTED",
        "EXPLORATORY",
        "UNSUPPORTED"
      )
    )
  )

write.csv(
  class_summary,
  file.path(
    TAB,
    "22G4B_connection_class_summary.csv"
  ),
  row.names=FALSE
)

cat("\n===== CONNECTION SUMMARY =====\n")

print(
  class_summary,
  n=Inf
)

# ============================================================
# 16. IMMUNE-STATE SUMMARY
# ============================================================

immune_summary <- coupling |>
  group_by(
    immune_L1,
    immune_L2,
    immune_feature
  ) |>
  summarise(
    n_connections=n(),

    n_CORE=
      sum(
        connection_class=="CORE"
      ),

    n_SUPPORTED=
      sum(
        connection_class=="SUPPORTED"
      ),

    n_NOMINAL=
      sum(
        connection_class=="NOMINAL"
      ),

    n_GLOBAL_SUPPORTED=
      sum(
        connection_class=="GLOBAL_SUPPORTED"
      ),

    n_EXPLORATORY=
      sum(
        connection_class=="EXPLORATORY"
      ),

    best_adjusted_FDR=
      min(
        adjusted_FDR,
        na.rm=TRUE
      ),

    max_abs_adjusted_beta=
      max(
        abs(
          adjusted_beta
        ),
        na.rm=TRUE
      ),

    .groups="drop"
  ) |>
  mutate(
    n_supported_total=
      n_CORE +
      n_SUPPORTED +
      n_NOMINAL +
      n_GLOBAL_SUPPORTED +
      n_EXPLORATORY
  ) |>
  arrange(
    desc(n_CORE),
    desc(n_SUPPORTED),
    desc(n_NOMINAL),
    best_adjusted_FDR
  )

write.csv(
  immune_summary,
  file.path(
    TAB,
    "22G4B_immune_state_connection_summary.csv"
  ),
  row.names=FALSE
)

cat("\n===== IMMUNE STATE SUMMARY =====\n")

print(
  immune_summary,
  n=Inf
)

# ============================================================
# 17. LAYER × IMMUNE SUMMARY
# ============================================================

layer_immune_summary <- coupling |>
  group_by(
    immune_L1,
    immune_L2,
    layer
  ) |>
  summarise(
    n_tested=n(),

    n_CORE=
      sum(
        connection_class=="CORE"
      ),

    n_SUPPORTED=
      sum(
        connection_class=="SUPPORTED"
      ),

    n_NOMINAL=
      sum(
        connection_class=="NOMINAL"
      ),

    best_FDR=
      min(
        adjusted_FDR,
        na.rm=TRUE
      ),

    strongest_beta=
      adjusted_beta[
        which.max(
          abs(
            adjusted_beta
          )
        )
      ],

    strongest_target=
      ecosystem_variable[
        which.max(
          abs(
            adjusted_beta
          )
        )
      ],

    .groups="drop"
  ) |>
  arrange(
    layer,
    best_FDR
  )

write.csv(
  layer_immune_summary,
  file.path(
    TAB,
    "22G4B_layer_by_immune_summary.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 18. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4B — IMMUNE × ECOSYSTEM COUPLING ATLAS",
  "============================================================",
  "",

  paste0(
    "Patients: ",
    n_distinct(
      master$patient_id
    )
  ),

  paste0(
    "Eligible immune states: ",
    n_distinct(
      immune$feature
    )
  ),

  paste0(
    "Ecosystem features: ",
    nrow(
      manifest
    )
  ),

  paste0(
    "Connections tested: ",
    nrow(
      coupling
    )
  ),

  "",

  "Biological layers:",
  "- BULGE",
  "- FOLLICULAR_EPITHELIUM",
  "- SUPERFICIAL_EPITHELIUM",
  "- MESENCHYMAL_NICHE",

  "",

  paste0(
    "CORE connections: ",
    sum(
      coupling$connection_class=="CORE"
    )
  ),

  paste0(
    "SUPPORTED connections: ",
    sum(
      coupling$connection_class=="SUPPORTED"
    )
  ),

  paste0(
    "NOMINAL connections: ",
    sum(
      coupling$connection_class=="NOMINAL"
    )
  ),

  paste0(
    "GLOBAL_SUPPORTED connections: ",
    sum(
      coupling$connection_class=="GLOBAL_SUPPORTED"
    )
  ),

  paste0(
    "EXPLORATORY connections: ",
    sum(
      coupling$connection_class=="EXPLORATORY"
    )
  ),

  "",

  "Statistical framework:",
  "- Patient is the statistical unit.",
  "- Spearman correlation describes global monotonic coupling.",
  "- Disease-adjusted linear models are the primary inferential models.",
  "- Disease interactions are exploratory.",
  "- BH correction is performed independently within each biological layer.",

  "",

  "Interpretation constraints:",
  "- Coupling represents patient-level co-variation.",
  "- A connection does not establish direct cell-cell communication.",
  "- A connection does not establish causal direction.",
  "- Indirect paths through fibroblast, vascular or epithelial states remain possible.",
  "- Direct versus indirect topology will be assessed only after conditional modelling and LR integration.",

  "",

  paste0(
    "FORMAL CHECKS: ",
    sum(
      checks$pass
    ),
    "/",
    nrow(
      checks
    ),
    " PASS"
  ),

  "",

  if(
    all(
      checks$pass
    )
  ) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUD,
    "22G4B_audit.txt"
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
# 19. FINAL STATUS
# ============================================================

if(
  !all(
    checks$pass
  )
) {
  stop(
    "22G4B formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4B COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

