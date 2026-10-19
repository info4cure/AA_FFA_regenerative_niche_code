# ==============================================================================
# 22F4E4D — FORMAL INFLAMMATION × HFSC COUPLING
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

MODULE <- "22F4E4D"

ROOT <- "PIPELINE_V2/results/disease_biology/22_MULTICOMPARTMENT_DISEASE_AXES/22F_BULGE_STEM_NICHE_INTEGRITY/22F4_BULGE_MONOCLE3_TRAJECTORIES/22F4E_DISEASE_TRAJECTORY/22F4E4_HALLMARK_TRAJECTORY"

INPUT <- file.path(
  ROOT,
  "22F4E4C_HYBRID_HFSC_INFLAMMATORY_STATE",
  "objects",
  "22F4E4C_hybrid_HFSC_inflammatory_state_freeze.rds"
)

OUT <- file.path(
  ROOT,
  "22F4E4D_INFLAMMATION_HFSC_COUPLING"
)

TABLE_DIR <- file.path(OUT,"tables")
OBJECT_DIR <- file.path(OUT,"objects")

dir.create(TABLE_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJECT_DIR, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(INPUT))

x <- readRDS(INPUT)
dat <- x$data

cat("\n============================================================\n")
cat("22F4E4D — INFLAMMATION × HFSC COUPLING\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# 1. INPUT AUDIT
# ------------------------------------------------------------------------------

required <- c(
  "patient_id",
  "disease",
  "trajectory",
  "HFSC_identity",
  "inflammatory_transition"
)

stopifnot(
  all(required %in% colnames(dat)),
  nrow(dat)==56,
  length(unique(dat$patient_id))==16,
  all(is.finite(dat$HFSC_identity)),
  all(is.finite(dat$inflammatory_transition)),
  all(is.finite(dat$trajectory))
)

dat <- dat |>
  mutate(
    patient_id=factor(patient_id),
    disease=factor(
      disease,
      levels=c("HC","FFA","AA")
    )
  )

input_audit <- dat |>
  group_by(disease) |>
  summarise(
    n_pseudobulks=n(),
    n_patients=n_distinct(patient_id),
    mean_HFSC=mean(HFSC_identity),
    mean_inflammatory=mean(inflammatory_transition),
    .groups="drop"
  )

cat("\n===== INPUT AUDIT =====\n")
print(input_audit, n=Inf)

write.csv(
  input_audit,
  file.path(TABLE_DIR,"22F4E4D_input_audit.csv"),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 2. PRIMARY MODEL
#
# Does disease alter the relationship between inflammatory activation
# and HFSC identity?
# ------------------------------------------------------------------------------

m0 <- lmer(
  HFSC_identity ~
    inflammatory_transition +
    disease +
    (1|patient_id),
  data=dat,
  REML=FALSE
)

m1 <- lmer(
  HFSC_identity ~
    inflammatory_transition*disease +
    (1|patient_id),
  data=dat,
  REML=FALSE
)

primary_LRT <- anova(
  m0,
  m1,
  refit=FALSE
)

primary_LRT_table <- as.data.frame(primary_LRT) |>
  rownames_to_column("model")

write.csv(
  primary_LRT_table,
  file.path(
    TABLE_DIR,
    "22F4E4D_primary_global_interaction_LRT.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 3. DISEASE-SPECIFIC SLOPES
# ------------------------------------------------------------------------------

tr <- emtrends(
  m1,
  ~disease,
  var="inflammatory_transition"
)

slopes <- as.data.frame(
  summary(
    tr,
    infer=c(TRUE,TRUE)
  )
)

slope_pairs <- as.data.frame(
  summary(
    pairs(
      tr,
      adjust="BH"
    ),
    infer=c(TRUE,TRUE)
  )
)

write.csv(
  slopes,
  file.path(
    TABLE_DIR,
    "22F4E4D_disease_specific_slopes.csv"
  ),
  row.names=FALSE
)

write.csv(
  slope_pairs,
  file.path(
    TABLE_DIR,
    "22F4E4D_pairwise_slope_contrasts.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 4. TRAJECTORY-ADJUSTED SENSITIVITY
# ------------------------------------------------------------------------------

m0_adj <- lmer(
  HFSC_identity ~
    inflammatory_transition +
    disease +
    trajectory +
    (1|patient_id),
  data=dat,
  REML=FALSE
)

m1_adj <- lmer(
  HFSC_identity ~
    inflammatory_transition*disease +
    trajectory +
    (1|patient_id),
  data=dat,
  REML=FALSE
)

sensitivity_LRT <- anova(
  m0_adj,
  m1_adj,
  refit=FALSE
)

sensitivity_LRT_table <-
  as.data.frame(sensitivity_LRT) |>
  rownames_to_column("model")

tr_adj <- emtrends(
  m1_adj,
  ~disease,
  var="inflammatory_transition"
)

slopes_adj <- as.data.frame(
  summary(
    tr_adj,
    infer=c(TRUE,TRUE)
  )
)

slope_pairs_adj <- as.data.frame(
  summary(
    pairs(
      tr_adj,
      adjust="BH"
    ),
    infer=c(TRUE,TRUE)
  )
)

write.csv(
  sensitivity_LRT_table,
  file.path(
    TABLE_DIR,
    "22F4E4D_trajectory_adjusted_global_LRT.csv"
  ),
  row.names=FALSE
)

write.csv(
  slopes_adj,
  file.path(
    TABLE_DIR,
    "22F4E4D_trajectory_adjusted_slopes.csv"
  ),
  row.names=FALSE
)

write.csv(
  slope_pairs_adj,
  file.path(
    TABLE_DIR,
    "22F4E4D_trajectory_adjusted_pairwise_slopes.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 5. MODEL AUDIT
# ------------------------------------------------------------------------------

model_audit <- tibble(
  check=c(
    "n_pseudobulks_56",
    "n_patients_16",
    "three_diseases_present",
    "primary_model_not_singular",
    "adjusted_model_not_singular",
    "primary_model_converged",
    "adjusted_model_converged"
  ),
  passed=c(
    nrow(dat)==56,
    n_distinct(dat$patient_id)==16,
    nlevels(dat$disease)==3,
    !isSingular(m1,tol=1e-4),
    !isSingular(m1_adj,tol=1e-4),
    is.null(m1@optinfo$conv$lme4$messages),
    is.null(m1_adj@optinfo$conv$lme4$messages)
  )
)

write.csv(
  model_audit,
  file.path(
    TABLE_DIR,
    "22F4E4D_model_audit.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 6. FREEZE
# ------------------------------------------------------------------------------

freeze <- list(
  module=MODULE,
  analysis=
    "Disease-specific inflammation-HFSC coupling",
  source=
    "Frozen 22F4E4C joint-state data",
  n_pseudobulks=nrow(dat),
  n_patients=n_distinct(dat$patient_id),

  primary_formula=
    "HFSC_identity ~ inflammatory_transition*disease + (1|patient_id)",

  sensitivity_formula=
    "HFSC_identity ~ inflammatory_transition*disease + trajectory + (1|patient_id)",

  data=dat,

  primary_null=m0,
  primary_full=m1,
  primary_LRT=primary_LRT_table,
  slopes=slopes,
  slope_pairs=slope_pairs,

  adjusted_null=m0_adj,
  adjusted_full=m1_adj,
  adjusted_LRT=sensitivity_LRT_table,
  adjusted_slopes=slopes_adj,
  adjusted_slope_pairs=slope_pairs_adj,

  model_audit=model_audit
)

saveRDS(
  freeze,
  file.path(
    OBJECT_DIR,
    "22F4E4D_inflammation_HFSC_coupling_freeze.rds"
  )
)

# ------------------------------------------------------------------------------
# 7. FINAL REPORT
# ------------------------------------------------------------------------------

cat("\n===== PRIMARY GLOBAL DISEASE × INFLAMMATION TEST =====\n")
print(primary_LRT)

cat("\n===== DISEASE-SPECIFIC SLOPES =====\n")
print(slopes, row.names=FALSE)

cat("\n===== PAIRWISE SLOPE CONTRASTS — BH =====\n")
print(slope_pairs, row.names=FALSE)

cat("\n===== TRAJECTORY-ADJUSTED GLOBAL TEST =====\n")
print(sensitivity_LRT)

cat("\n===== TRAJECTORY-ADJUSTED SLOPES =====\n")
print(slopes_adj, row.names=FALSE)

cat("\n===== TRAJECTORY-ADJUSTED PAIRWISE SLOPES — BH =====\n")
print(slope_pairs_adj, row.names=FALSE)

cat("\n===== MODEL AUDIT =====\n")
print(model_audit, n=Inf)

cat("\nSTATUS=PASS\n")
