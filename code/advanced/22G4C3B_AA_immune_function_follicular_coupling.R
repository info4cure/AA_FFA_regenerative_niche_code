
# ============================================================
# 22G4C3B — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY COUPLING
#
# AA ONLY
# Patient = statistical unit
#
# Frozen universe:
# 7 pre-specified immune functional axes
# 9 pre-specified structural targets
# 63 pre-specified connections
#
# Analyses:
# 1. Spearman coupling
# 2. >=5-cell sensitivity
# 3. leave-one-patient-out robustness
# 4. lower-versus-upper anatomical preference
#
# No feature selection from observed results.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

cat("\n============================================================\n")
cat("22G4C3B — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY COUPLING\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

BASE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE"
)

PREF <- file.path(
  BASE,
  "22G4C3_IMMUNE_FUNCTION_FOLLICULAR_TOPOLOGY",
  "22G4C3A_PREFLIGHT"
)

OUT <- file.path(
  BASE,
  "22G4C3_IMMUNE_FUNCTION_FOLLICULAR_TOPOLOGY",
  "22G4C3B_AA_COUPLING"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJ,recursive=TRUE,showWarnings=FALSE)
dir.create(AUD,recursive=TRUE,showWarnings=FALSE)

PAIR_FILE <- file.path(
  PREF,
  "tables",
  "22G4C3A_AA_patient_pair_master.csv"
)

SUPPORT_FILE <- file.path(
  PREF,
  "tables",
  "22G4C3A_connection_support.csv"
)

stopifnot(
  file.exists(PAIR_FILE),
  file.exists(SUPPORT_FILE)
)

pair_master <- read.csv(
  PAIR_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  as_tibble()

connection_support <- read.csv(
  SUPPORT_FILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
) |>
  as_tibble()

cat("\n===== FROZEN INPUT =====\n")
cat("Patient-level pair rows:",nrow(pair_master),"\n")
cat("Pre-specified connections:",nrow(connection_support),"\n")
cat("AA patients:",n_distinct(pair_master$patient_id),"\n")

stopifnot(
  nrow(connection_support)==63,
  all(pair_master$disease=="AA")
)

# ============================================================
# 1. CONNECTION IDENTIFIER
# ============================================================

pair_master <- pair_master |>
  mutate(
    connection_id=paste(
      immune_axis,
      target_variable,
      sep="__"
    )
  )

connection_support <- connection_support |>
  mutate(
    connection_id=paste(
      immune_axis,
      target_variable,
      sep="__"
    )
  )

stopifnot(
  n_distinct(connection_support$connection_id)==63
)

# ============================================================
# 2. SAFE SPEARMAN FUNCTION
# ============================================================

safe_spearman <- function(d) {

  d <- d |>
    filter(
      is.finite(immune_score),
      is.finite(structural_frequency)
    )

  n <- nrow(d)

  if(
    n < 3 ||
    n_distinct(d$immune_score) < 3 ||
    n_distinct(d$structural_frequency) < 3
  ) {
    return(
      tibble(
        n=n,
        rho=NA_real_,
        p=NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      d$immune_score,
      d$structural_frequency,
      method="spearman",
      exact=FALSE
    )
  )

  tibble(
    n=n,
    rho=unname(ct$estimate),
    p=ct$p.value
  )
}

# ============================================================
# 3. PRIMARY ALL-PATIENT COUPLING
# ============================================================

cat("\n===== PRIMARY ALL-PATIENT SPEARMAN COUPLING =====\n")

primary <- pair_master |>
  group_by(
    connection_id,
    immune_L2,
    immune_axis,
    priority,
    anatomical_axis,
    target_L1,
    target_L2,
    target_variable,
    role
  ) |>
  group_modify(
    ~safe_spearman(.x)
  ) |>
  ungroup()

stopifnot(
  nrow(primary)==63
)

# ============================================================
# 4. >=5-CELL IMMUNE SUPPORT SENSITIVITY
# ============================================================

cat("\n===== >=5-CELL SENSITIVITY =====\n")

ge5 <- pair_master |>
  filter(
    immune_n_cells >= 5
  ) |>
  group_by(
    connection_id,
    immune_L2,
    immune_axis,
    priority,
    anatomical_axis,
    target_L1,
    target_L2,
    target_variable,
    role
  ) |>
  group_modify(
    ~safe_spearman(.x)
  ) |>
  ungroup() |>
  rename(
    n_GE5=n,
    rho_GE5=rho,
    p_GE5=p
  )

# ============================================================
# 5. LEAVE-ONE-PATIENT-OUT ROBUSTNESS
# ============================================================

cat("\n===== LEAVE-ONE-PATIENT-OUT ROBUSTNESS =====\n")

connection_ids <- unique(
  pair_master$connection_id
)

lopo_list <- vector(
  "list",
  length(connection_ids)
)

for(i in seq_along(connection_ids)) {

  cid <- connection_ids[i]

  d <- pair_master |>
    filter(
      connection_id==cid
    )

  patients <- sort(
    unique(
      d$patient_id
    )
  )

  one_connection <- list()

  for(j in seq_along(patients)) {

    omitted <- patients[j]

    dd <- d |>
      filter(
        patient_id != omitted
      )

    s <- safe_spearman(dd)

    one_connection[[j]] <- tibble(
      connection_id=cid,
      omitted_patient=omitted,
      n=s$n,
      rho=s$rho,
      p=s$p
    )
  }

  lopo_list[[i]] <- bind_rows(
    one_connection
  )
}

lopo <- bind_rows(
  lopo_list
)

# ============================================================
# 6. LOPO SUMMARY
# ============================================================

primary_sign <- primary |>
  dplyr::select(
    connection_id,
    primary_rho=rho
  )

lopo_summary <- lopo |>
  left_join(
    primary_sign,
    by="connection_id"
  ) |>
  group_by(
    connection_id
  ) |>
  summarise(
    n_LOPO=n(),

    n_LOPO_finite=
      sum(
        is.finite(rho)
      ),

    median_LOPO_rho=
      ifelse(
        any(
          is.finite(rho)
        ),
        median(
          rho,
          na.rm=TRUE
        ),
        NA_real_
      ),

    min_LOPO_rho=
      ifelse(
        any(
          is.finite(rho)
        ),
        min(
          rho,
          na.rm=TRUE
        ),
        NA_real_
      ),

    max_LOPO_rho=
      ifelse(
        any(
          is.finite(rho)
        ),
        max(
          rho,
          na.rm=TRUE
        ),
        NA_real_
      ),

    same_sign_n=
      sum(
        is.finite(rho) &
        is.finite(primary_rho) &
        sign(rho)==sign(primary_rho)
      ),

    same_sign_fraction=
      ifelse(
        sum(
          is.finite(rho)
        ) > 0,
        sum(
          is.finite(rho) &
          is.finite(primary_rho) &
          sign(rho)==sign(primary_rho)
        ) /
          sum(
            is.finite(rho)
          ),
        NA_real_
      ),

    .groups="drop"
  )

# ============================================================
# 7. INTEGRATE PRIMARY + GE5 + LOPO
# ============================================================

integrated <- primary |>
  left_join(
    ge5 |>
      dplyr::select(
        connection_id,
        n_GE5,
        rho_GE5,
        p_GE5
      ),
    by="connection_id"
  ) |>
  left_join(
    lopo_summary,
    by="connection_id"
  ) |>
  left_join(
    connection_support |>
      dplyr::select(
        connection_id,
        immune_total_cells,
        immune_min_cells,
        immune_median_cells,
        immune_n_ge5,
        immune_n_ge10,
        structural_total_cells,
        structural_zero_patients,
        structural_unique_values,
        analysis_tier
      ),
    by="connection_id"
  ) |>
  mutate(
    abs_rho=abs(rho),

    GE5_direction_concordant=
      case_when(
        !is.finite(rho_GE5)
          ~ NA,
        sign(rho_GE5)==sign(rho)
          ~ TRUE,
        TRUE
          ~ FALSE
      ),

    sparse_structural=
      structural_zero_patients > 0 |
      structural_total_cells < 50,

    lower_domain=
      anatomical_axis %in%c(
        "LOWER_FOLLICLE",
        "LOWER_MESENCHYME"
      )
  )

# ============================================================
# 8. MULTIPLE TESTING — PRE-SPECIFIED FAMILIES
# ============================================================

cat("\n===== MULTIPLE-TESTING FAMILIES =====\n")

integrated <- integrated |>
  mutate(
    testing_family=case_when(

      analysis_tier=="PRIMARY_CD8" &
        lower_domain
        ~ "PRIMARY_CD8_LOWER",

      analysis_tier=="PRIMARY_CD8" &
        anatomical_axis=="UPPER_COMPARATOR"
        ~ "PRIMARY_CD8_UPPER",

      analysis_tier=="PRIMARY_NK_LOW_SUPPORT" &
        lower_domain
        ~ "NK_LOWER_EXPLORATORY",

      analysis_tier=="PRIMARY_NK_LOW_SUPPORT" &
        anatomical_axis=="UPPER_COMPARATOR"
        ~ "NK_UPPER_EXPLORATORY",

      analysis_tier=="SECONDARY" &
        lower_domain
        ~ "MAC_LOWER_SECONDARY",

      analysis_tier=="SECONDARY" &
        anatomical_axis=="UPPER_COMPARATOR"
        ~ "MAC_UPPER_SECONDARY",

      TRUE
        ~ "OTHER"
    )
  ) |>
  group_by(
    testing_family
  ) |>
  mutate(
    FDR_family=
      p.adjust(
        p,
        method="BH"
      )
  ) |>
  ungroup()

family_counts <- integrated |>
  count(
    testing_family,
    name="n_tests"
  ) |>
  arrange(
    testing_family
  )

print(
  family_counts,
  n=Inf,
  width=Inf
)

# ============================================================
# 9. PRE-SPECIFIED CLASSIFICATION
# ============================================================

integrated <- integrated |>
  mutate(
    coupling_class=case_when(

      # ------------------------------------------------------
      # CD8 — confirmatory / primary
      # ------------------------------------------------------

      analysis_tier=="PRIMARY_CD8" &
        is.finite(rho) &
        abs_rho >= 0.70 &
        p < 0.05 &
        same_sign_fraction >= 0.80 &
        GE5_direction_concordant %in% TRUE
        ~ "PRIMARY_ROBUST",

      analysis_tier=="PRIMARY_CD8" &
        is.finite(rho) &
        abs_rho >= 0.60 &
        same_sign_fraction >= 0.80 &
        GE5_direction_concordant %in% TRUE
        ~ "PRIMARY_SUGGESTIVE",

      # ------------------------------------------------------
      # NK — exploratory because n=4
      # ------------------------------------------------------

      analysis_tier=="PRIMARY_NK_LOW_SUPPORT" &
        is.finite(rho) &
        abs_rho >= 0.80 &
        same_sign_fraction >= 0.75
        ~ "NK_EXPLORATORY_STRONG",

      analysis_tier=="PRIMARY_NK_LOW_SUPPORT" &
        is.finite(rho) &
        abs_rho >= 0.60 &
        same_sign_fraction >= 0.75
        ~ "NK_EXPLORATORY",

      # ------------------------------------------------------
      # Macrophage — secondary
      # ------------------------------------------------------

      analysis_tier=="SECONDARY" &
        is.finite(rho) &
        abs_rho >= 0.70 &
        p < 0.05 &
        same_sign_fraction >= 0.80 &
        GE5_direction_concordant %in% TRUE
        ~ "SECONDARY_ROBUST",

      analysis_tier=="SECONDARY" &
        is.finite(rho) &
        abs_rho >= 0.60 &
        same_sign_fraction >= 0.80 &
        GE5_direction_concordant %in% TRUE
        ~ "SECONDARY_SUGGESTIVE",

      TRUE
        ~ "UNSUPPORTED"
    ),

    FDR_supported=
      is.finite(FDR_family) &
      FDR_family < 0.05
  )

# ============================================================
# 10. NK PROGRAMME COHERENCE
# ============================================================

cat("\n===== NK PROGRAMME COHERENCE =====\n")

nk_coherence <- integrated |>
  filter(
    immune_axis %in%c(
      "NK_CYTOTOXICITY",
      "NK_IFNG"
    )
  ) |>
  dplyr::select(
    anatomical_axis,
    target_variable,
    immune_axis,
    rho,
    abs_rho,
    coupling_class
  ) |>
  pivot_wider(
    names_from=immune_axis,
    values_from=c(
      rho,
      abs_rho,
      coupling_class
    )
  ) |>
  mutate(
    NK_same_direction=
      is.finite(
        rho_NK_CYTOTOXICITY
      ) &
      is.finite(
        rho_NK_IFNG
      ) &
      sign(
        rho_NK_CYTOTOXICITY
      )==
      sign(
        rho_NK_IFNG
      ),

    NK_both_abs_ge_06=
      abs_rho_NK_CYTOTOXICITY >= 0.60 &
      abs_rho_NK_IFNG >= 0.60,

    NK_both_abs_ge_08=
      abs_rho_NK_CYTOTOXICITY >= 0.80 &
      abs_rho_NK_IFNG >= 0.80,

    NK_coherence_class=case_when(
      NK_same_direction &
        NK_both_abs_ge_08
        ~ "NK_COHERENT_STRONG",

      NK_same_direction &
        NK_both_abs_ge_06
        ~ "NK_COHERENT",

      NK_same_direction
        ~ "NK_DIRECTION_CONCORDANT",

      TRUE
        ~ "NK_DISCORDANT"
    )
  )

print(
  nk_coherence,
  n=Inf,
  width=Inf
)

# ============================================================
# 11. LOWER-versus-UPPER TOPOLOGICAL PREFERENCE
# ============================================================

cat("\n===== LOWER-vs-UPPER TOPOLOGICAL PREFERENCE =====\n")

# Descriptive topology metric.
#
# We compare absolute coupling strength rather than signed rho,
# because different structural outcomes can legitimately move
# in opposite directions during follicular injury.
#
# lower_primary = LOWER_FOLLICLE + LOWER_MESENCHYME
# upper comparator = IFE/Infundibulum states.

topology_preference <- integrated |>
  group_by(
    immune_L2,
    immune_axis,
    priority,
    analysis_tier
  ) |>
  summarise(
    n_lower=
      sum(
        lower_domain &
        is.finite(rho)
      ),

    n_upper=
      sum(
        anatomical_axis=="UPPER_COMPARATOR" &
        is.finite(rho)
      ),

    median_abs_rho_lower=
      ifelse(
        any(
          lower_domain &
          is.finite(rho)
        ),
        median(
          abs_rho[
            lower_domain &
            is.finite(rho)
          ],
          na.rm=TRUE
        ),
        NA_real_
      ),

    median_abs_rho_upper=
      ifelse(
        any(
          anatomical_axis=="UPPER_COMPARATOR" &
          is.finite(rho)
        ),
        median(
          abs_rho[
            anatomical_axis=="UPPER_COMPARATOR" &
            is.finite(rho)
          ],
          na.rm=TRUE
        ),
        NA_real_
      ),

    max_abs_rho_lower=
      ifelse(
        any(
          lower_domain &
          is.finite(rho)
        ),
        max(
          abs_rho[
            lower_domain &
            is.finite(rho)
          ],
          na.rm=TRUE
        ),
        NA_real_
      ),

    max_abs_rho_upper=
      ifelse(
        any(
          anatomical_axis=="UPPER_COMPARATOR" &
          is.finite(rho)
        ),
        max(
          abs_rho[
            anatomical_axis=="UPPER_COMPARATOR" &
            is.finite(rho)
          ],
          na.rm=TRUE
        ),
        NA_real_
      ),

    n_lower_abs_ge_06=
      sum(
        lower_domain &
        is.finite(rho) &
        abs_rho >= 0.60
      ),

    n_upper_abs_ge_06=
      sum(
        anatomical_axis=="UPPER_COMPARATOR" &
        is.finite(rho) &
        abs_rho >= 0.60
      ),

    n_lower_abs_ge_07=
      sum(
        lower_domain &
        is.finite(rho) &
        abs_rho >= 0.70
      ),

    n_upper_abs_ge_07=
      sum(
        anatomical_axis=="UPPER_COMPARATOR" &
        is.finite(rho) &
        abs_rho >= 0.70
      ),

    .groups="drop"
  ) |>
  mutate(
    delta_median_abs_rho_lower_upper=
      median_abs_rho_lower -
      median_abs_rho_upper,

    topology_preference=case_when(
      !is.finite(
        delta_median_abs_rho_lower_upper
      )
        ~ "NOT_EVALUABLE",

      delta_median_abs_rho_lower_upper >= 0.20
        ~ "LOWER_PREFERRED",

      delta_median_abs_rho_lower_upper <= -0.20
        ~ "UPPER_PREFERRED",

      TRUE
        ~ "NO_CLEAR_PREFERENCE"
    )
  )

print(
  topology_preference,
  n=Inf,
  width=Inf
)

# ============================================================
# 12. PRIORITIZED CONNECTION LEDGER
# ============================================================

cat("\n===== PRIORITIZED CONNECTION LEDGER =====\n")

class_order <- c(
  "PRIMARY_ROBUST",
  "PRIMARY_SUGGESTIVE",
  "NK_EXPLORATORY_STRONG",
  "NK_EXPLORATORY",
  "SECONDARY_ROBUST",
  "SECONDARY_SUGGESTIVE",
  "UNSUPPORTED"
)

prioritized <- integrated |>
  filter(
    coupling_class != "UNSUPPORTED"
  ) |>
  mutate(
    class_rank=
      match(
        coupling_class,
        class_order
      )
  ) |>
  arrange(
    class_rank,
    FDR_family,
    p,
    desc(abs_rho)
  ) |>
  dplyr::select(
    -class_rank
  )

print(
  prioritized |>
    dplyr::select(
      immune_L2,
      immune_axis,
      anatomical_axis,
      target_variable,
      n,
      rho,
      p,
      FDR_family,
      n_GE5,
      rho_GE5,
      p_GE5,
      same_sign_fraction,
      median_LOPO_rho,
      GE5_direction_concordant,
      sparse_structural,
      coupling_class,
      FDR_supported
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 13. LOWER-FOLLICLE / MESENCHYME PRIORITY LEDGER
# ============================================================

cat("\n===== LOWER-UNIT PRIORITY LEDGER =====\n")

lower_priority <- integrated |>
  filter(
    lower_domain
  ) |>
  mutate(
    class_rank=
      match(
        coupling_class,
        class_order
      )
  ) |>
  arrange(
    class_rank,
    FDR_family,
    p,
    desc(abs_rho)
  ) |>
  dplyr::select(
    -class_rank
  )

print(
  lower_priority |>
    dplyr::select(
      immune_axis,
      anatomical_axis,
      target_variable,
      n,
      rho,
      p,
      FDR_family,
      rho_GE5,
      same_sign_fraction,
      sparse_structural,
      coupling_class
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 14. UPPER COMPARATOR LEDGER
# ============================================================

cat("\n===== UPPER COMPARATOR LEDGER =====\n")

upper_comparator <- integrated |>
  filter(
    anatomical_axis=="UPPER_COMPARATOR"
  ) |>
  arrange(
    immune_axis,
    desc(abs_rho)
  )

print(
  upper_comparator |>
    dplyr::select(
      immune_axis,
      target_variable,
      n,
      rho,
      p,
      FDR_family,
      rho_GE5,
      same_sign_fraction,
      coupling_class
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 15. CONNECTION CLASS COUNTS
# ============================================================

cat("\n===== CONNECTION CLASS COUNTS =====\n")

class_counts <- integrated |>
  count(
    analysis_tier,
    coupling_class,
    name="n_connections"
  ) |>
  arrange(
    analysis_tier,
    coupling_class
  )

print(
  class_counts,
  n=Inf,
  width=Inf
)

# ============================================================
# 16. SAVE TABLES
# ============================================================

write.csv(
  primary,
  file.path(
    TAB,
    "22G4C3B_primary_AA_Spearman_coupling.csv"
  ),
  row.names=FALSE
)

write.csv(
  ge5,
  file.path(
    TAB,
    "22G4C3B_GE5_sensitivity.csv"
  ),
  row.names=FALSE
)

write.csv(
  lopo,
  file.path(
    TAB,
    "22G4C3B_leave_one_patient_out.csv"
  ),
  row.names=FALSE
)

write.csv(
  lopo_summary,
  file.path(
    TAB,
    "22G4C3B_LOPO_summary.csv"
  ),
  row.names=FALSE
)

write.csv(
  integrated,
  file.path(
    TAB,
    "22G4C3B_integrated_AA_coupling.csv"
  ),
  row.names=FALSE
)

write.csv(
  prioritized,
  file.path(
    TAB,
    "22G4C3B_prioritized_connections.csv"
  ),
  row.names=FALSE
)

write.csv(
  lower_priority,
  file.path(
    TAB,
    "22G4C3B_lower_unit_ledger.csv"
  ),
  row.names=FALSE
)

write.csv(
  upper_comparator,
  file.path(
    TAB,
    "22G4C3B_upper_comparator_ledger.csv"
  ),
  row.names=FALSE
)

write.csv(
  nk_coherence,
  file.path(
    TAB,
    "22G4C3B_NK_programme_coherence.csv"
  ),
  row.names=FALSE
)

write.csv(
  topology_preference,
  file.path(
    TAB,
    "22G4C3B_lower_upper_topology_preference.csv"
  ),
  row.names=FALSE
)

write.csv(
  family_counts,
  file.path(
    TAB,
    "22G4C3B_testing_family_counts.csv"
  ),
  row.names=FALSE
)

write.csv(
  class_counts,
  file.path(
    TAB,
    "22G4C3B_connection_class_counts.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 17. FREEZE
# ============================================================

freeze <- list(
  module="22G4C3B",

  purpose=paste(
    "AA-only patient-level coupling between pre-specified",
    "immune functional programmes and frozen follicular",
    "structural states."
  ),

  disease_universe="AA only",

  statistical_unit="patient",

  n_prespecified_connections=
    nrow(
      integrated
    ),

  primary_coupling=
    primary,

  GE5_sensitivity=
    ge5,

  LOPO=
    lopo,

  LOPO_summary=
    lopo_summary,

  integrated=
    integrated,

  prioritized=
    prioritized,

  NK_coherence=
    nk_coherence,

  topology_preference=
    topology_preference,

  classification_rules=list(
    PRIMARY_ROBUST=paste(
      "|rho| >=0.70; nominal p<0.05;",
      "LOPO same-sign fraction >=0.80;",
      "GE5 same direction."
    ),

    PRIMARY_SUGGESTIVE=paste(
      "|rho| >=0.60;",
      "LOPO same-sign fraction >=0.80;",
      "GE5 same direction."
    ),

    NK_EXPLORATORY_STRONG=paste(
      "|rho| >=0.80;",
      "LOPO same-sign fraction >=0.75."
    ),

    NK_EXPLORATORY=paste(
      "|rho| >=0.60;",
      "LOPO same-sign fraction >=0.75."
    )
  ),

  interpretation_constraints=c(
    "AA-only analysis.",
    "Patient is the statistical unit.",
    "Spearman coupling is associative, not causal.",
    "NK n=4 is exploratory regardless of nominal P value.",
    "Hair shaft is flagged as sparse because two AA patients have zero cells.",
    "FDR is calculated within pre-specified biological testing families.",
    "Lower-versus-upper preference is descriptive and based on absolute rho.",
    "No post hoc feature selection was used to define the tested universe."
  )
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4C3B_AA_immune_function_follicular_coupling_freeze.rds"
  )
)

# ============================================================
# 18. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

checks <- tibble(
  check=c(
    "AA-only analysis",
    "63 pre-specified connections",
    "27 primary CD8 connections",
    "18 NK exploratory connections",
    "18 macrophage secondary connections",
    "Primary Spearman table complete",
    "GE5 sensitivity generated",
    "LOPO generated",
    "LOPO summary complete",
    "Six testing families represented",
    "Family FDR available",
    "NK coherence table generated",
    "Topology preference table has seven immune axes",
    "Lower-unit ledger generated",
    "Upper-comparator ledger generated",
    "Freeze exists"
  ),

  pass=c(
    all(
      pair_master$disease=="AA"
    ),

    nrow(
      integrated
    )==63,

    sum(
      integrated$analysis_tier==
        "PRIMARY_CD8"
    )==27,

    sum(
      integrated$analysis_tier==
        "PRIMARY_NK_LOW_SUPPORT"
    )==18,

    sum(
      integrated$analysis_tier==
        "SECONDARY"
    )==18,

    nrow(
      primary
    )==63,

    nrow(
      ge5
    )>0,

    nrow(
      lopo
    )>0,

    nrow(
      lopo_summary
    )==63,

    n_distinct(
      integrated$testing_family
    )==6,

    all(
      is.finite(
        integrated$FDR_family
      ) |
      !is.finite(
        integrated$p
      )
    ),

    nrow(
      nk_coherence
    )==9,

    nrow(
      topology_preference
    )==7,

    nrow(
      lower_priority
    )==42,

    nrow(
      upper_comparator
    )==21,

    file.exists(
      file.path(
        OBJ,
        "22G4C3B_AA_immune_function_follicular_coupling_freeze.rds"
      )
    )
  )
)

print(
  checks,
  n=Inf,
  width=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G4C3B_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 19. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C3B — AA IMMUNE FUNCTION × FOLLICULAR TOPOLOGY COUPLING",
  "============================================================",
  "",
  "Disease universe: AA only",
  "Statistical unit: patient",
  "",
  paste0(
    "Pre-specified connections tested: ",
    nrow(integrated)
  ),
  "",
  "Primary CD8 programmes:",
  "- cytotoxicity",
  "- IFN-gamma response",
  "- IFN-alpha response",
  "",
  "Exploratory NK programmes:",
  "- cytotoxicity",
  "- IFN-gamma response",
  "",
  "Secondary macrophage programmes:",
  "- IFN-alpha response",
  "- TNFA/NFKB signaling",
  "",
  "Structural domains:",
  "- LOWER_FOLLICLE: ORS, ORS basal, IRS/companion, Hair shaft",
  "- LOWER_MESENCHYME: dermal papilla and dermal sheath",
  "- UPPER_COMPARATOR: three IFE/Infundibulum states",
  "",
  "Robustness framework:",
  "- primary Spearman coupling",
  "- >=5 immune-cell sensitivity",
  "- leave-one-patient-out directional robustness",
  "- biological-family BH correction",
  "- descriptive lower-versus-upper topological preference",
  "",
  paste0(
    "PRIMARY_ROBUST: ",
    sum(
      integrated$coupling_class==
        "PRIMARY_ROBUST",
      na.rm=TRUE
    )
  ),
  paste0(
    "PRIMARY_SUGGESTIVE: ",
    sum(
      integrated$coupling_class==
        "PRIMARY_SUGGESTIVE",
      na.rm=TRUE
    )
  ),
  paste0(
    "NK_EXPLORATORY_STRONG: ",
    sum(
      integrated$coupling_class==
        "NK_EXPLORATORY_STRONG",
      na.rm=TRUE
    )
  ),
  paste0(
    "NK_EXPLORATORY: ",
    sum(
      integrated$coupling_class==
        "NK_EXPLORATORY",
      na.rm=TRUE
    )
  ),
  paste0(
    "SECONDARY_ROBUST: ",
    sum(
      integrated$coupling_class==
        "SECONDARY_ROBUST",
      na.rm=TRUE
    )
  ),
  paste0(
    "SECONDARY_SUGGESTIVE: ",
    sum(
      integrated$coupling_class==
        "SECONDARY_SUGGESTIVE",
      na.rm=TRUE
    )
  ),
  paste0(
    "Family-FDR <0.05 connections: ",
    sum(
      integrated$FDR_supported,
      na.rm=TRUE
    )
  ),
  "",
  paste(
    "NK results remain exploratory because only four AA",
    "patients have NK pseudobulks."
  ),
  paste(
    "Hair-shaft associations are flagged as structurally sparse",
    "because two AA patients have zero Hair-shaft cells."
  ),
  paste(
    "Lower-versus-upper preference is descriptive and must not",
    "be interpreted as formal evidence of anatomical causality."
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
    "22G4C3B_audit.txt"
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
# 20. FINAL STATUS
# ============================================================

if(!all(checks$pass)) {
  stop(
    "22G4C3B formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C3B COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

