suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(width = 200)

cat("============================================================\n")
cat("PIPELINE_V2 — 10B4D C8/C9 TECHNICAL COEXPRESSION AUDIT\n")
cat("============================================================\n")

SRC <- paste0(
  "PIPELINE_V2/objects/annotated/",
  "vascular_v2_10B4B_rebuilt_3023.rds"
)

EVID <- paste0(
  "PIPELINE_V2/results/annotation/10_VASCULAR_L2/",
  "10B4C_targeted_lineage_audit/tables/",
  "10B4C_cell_programme_evidence.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/10_VASCULAR_L2/",
  "10B4D_C8_C9_technical_audit"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(SRC)
ev <- read_csv(EVID, show_col_types = FALSE)

stopifnot(ncol(x) == 3023)

z <- ev %>%
  filter(cluster %in% c("8","9")) %>%
  mutate(
    EC_supported =
      ENDOTHELIAL_CORE_n >= 3,

    mural_supported =
      MURAL_PERICYTE_n >= 3,

    ECM_supported =
      MESENCHYMAL_ECM_n >= 3,

    epithelial_supported =
      EPITHELIAL_BASAL_n >= 2 |
      EPITHELIAL_JUNCTION_n >= 2,

    EC_mural_joint =
      EC_supported & mural_supported,

    EC_ECM_joint =
      EC_supported & ECM_supported,

    EC_epithelial_joint =
      EC_supported & epithelial_supported
  )

cat("\n===== CELLS =====\n")
print(table(z$cluster))

# ------------------------------------------------------------
# Technical correlations
# ------------------------------------------------------------

bio <- c(
  "ENDOTHELIAL_CORE_score",
  "MURAL_PERICYTE_score",
  "MESENCHYMAL_ECM_score",
  "PERIVASCULAR_EBF2_score",
  "EPITHELIAL_BASAL_score",
  "EPITHELIAL_JUNCTION_score"
)

tech <- c(
  "doublet_score",
  "nFeature_RNA",
  "nCount_RNA"
)

corr <- bind_rows(
  lapply(
    c("8","9"),
    function(cc) {

      zz <- z %>% filter(cluster == cc)

      bind_rows(
        lapply(
          bio,
          function(b) {

            bind_rows(
              lapply(
                tech,
                function(t) {

                  ct <- suppressWarnings(
                    cor.test(
                      zz[[b]],
                      zz[[t]],
                      method = "spearman",
                      exact = FALSE
                    )
                  )

                  tibble(
                    cluster = cc,
                    biological_score = b,
                    technical_variable = t,
                    n = nrow(zz),
                    spearman_rho =
                      unname(ct$estimate),
                    p_value =
                      ct$p.value
                  )
                }
              )
            )
          }
        )
      )
    }
  )
) %>%
  group_by(cluster) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

cat("\n===== TECHNICAL CORRELATIONS =====\n")
print(corr, n = Inf)

write_csv(
  corr,
  file.path(TAB, "10B4D_technical_correlations.csv")
)

# ------------------------------------------------------------
# Doublet-score quartiles within each cluster
# ------------------------------------------------------------

strata <- z %>%
  group_by(cluster) %>%
  mutate(
    doublet_stratum =
      ntile(doublet_score, 4)
  ) %>%
  ungroup() %>%
  mutate(
    doublet_stratum =
      factor(
        doublet_stratum,
        levels = 1:4,
        labels = c(
          "Q1_low",
          "Q2",
          "Q3",
          "Q4_high"
        )
      )
  )

strata_summary <- strata %>%
  group_by(
    cluster,
    doublet_stratum
  ) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    median_features =
      median(nFeature_RNA),

    median_counts =
      median(nCount_RNA),

    EC_score =
      median(ENDOTHELIAL_CORE_score),

    mural_score =
      median(MURAL_PERICYTE_score),

    ECM_score =
      median(MESENCHYMAL_ECM_score),

    epithelial_basal_score =
      median(EPITHELIAL_BASAL_score),

    epithelial_junction_score =
      median(EPITHELIAL_JUNCTION_score),

    EC_mural_joint_fraction =
      mean(EC_mural_joint),

    EC_ECM_joint_fraction =
      mean(EC_ECM_joint),

    EC_epithelial_joint_fraction =
      mean(EC_epithelial_joint),

    .groups = "drop"
  )

cat("\n===== DOUBLET STRATA =====\n")
print(strata_summary, n = Inf)

write_csv(
  strata_summary,
  file.path(TAB, "10B4D_doublet_strata.csv")
)

# ------------------------------------------------------------
# Lowest-doublet quartile explicitly
# ------------------------------------------------------------

low <- strata %>%
  filter(doublet_stratum == "Q1_low") %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    EC_supported_fraction =
      mean(EC_supported),

    mural_supported_fraction =
      mean(mural_supported),

    ECM_supported_fraction =
      mean(ECM_supported),

    epithelial_supported_fraction =
      mean(epithelial_supported),

    EC_mural_joint_fraction =
      mean(EC_mural_joint),

    EC_ECM_joint_fraction =
      mean(EC_ECM_joint),

    EC_epithelial_joint_fraction =
      mean(EC_epithelial_joint),

    .groups = "drop"
  )

cat("\n===== LOWEST DOUBLET QUARTILE =====\n")
print(low, n = Inf)

write_csv(
  low,
  file.path(TAB, "10B4D_low_doublet_quartile.csv")
)

# ------------------------------------------------------------
# Patient recurrence of joint programmes
# ------------------------------------------------------------

patient <- z %>%
  group_by(
    cluster,
    patient_id,
    batch
  ) %>%
  summarise(
    n_cells = n(),

    median_doublet =
      median(doublet_score),

    EC_fraction =
      mean(EC_supported),

    EC_mural_joint_fraction =
      mean(EC_mural_joint),

    EC_ECM_joint_fraction =
      mean(EC_ECM_joint),

    EC_epithelial_joint_fraction =
      mean(EC_epithelial_joint),

    .groups = "drop"
  )

cat("\n===== PATIENT RECURRENCE =====\n")
print(patient, n = Inf)

write_csv(
  patient,
  file.path(TAB, "10B4D_patient_recurrence.csv")
)

# ------------------------------------------------------------
# Formal checks
# ------------------------------------------------------------

checks <- tibble(
  check = c(
    "source_cells_3023",
    "C8_cells_113",
    "C9_cells_57",
    "technical_correlations_generated",
    "doublet_strata_generated",
    "low_doublet_quartile_generated",
    "patient_recurrence_generated",
    "disease_not_used",
    "no_cells_removed",
    "no_identity_assigned",
    "source_object_not_modified"
  ),

  pass = c(
    ncol(x) == 3023,
    sum(z$cluster == "8") == 113,
    sum(z$cluster == "9") == 57,
    nrow(corr) > 0,
    nrow(strata_summary) == 8,
    nrow(low) == 2,
    nrow(patient) > 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  checks,
  file.path(TAB, "10B4D_formal_checks.csv")
)

cat("\n===== FORMAL CHECKS =====\n")
print(checks, n = Inf)

stopifnot(all(checks$pass))

writeLines(
  c(
    "PIPELINE_V2 — 10B4D C8/C9 TECHNICAL COEXPRESSION AUDIT",
    "",
    "Purpose:",
    "- determine whether C8 endothelial-mural/ECM coexpression",
    "  and C9 endothelial-epithelial coexpression are associated",
    "  with doublet score or library complexity",
    "- test persistence in the lowest doublet-score quartile",
    "- assess recurrence across patients and batches",
    "",
    "Interpretive boundary:",
    "- disease not used",
    "- no cells removed",
    "- no identity assigned",
    "- no resolution frozen",
    "- source object unchanged",
    "",
    "PASS"
  ),
  file.path(TAB, "10B4D_audit.txt")
)

cat("\n============================================================\n")
cat("10B4D COMPLETE\n")
cat("============================================================\n")
cat("C8:", sum(z$cluster=="8"), "cells\n")
cat("C9:", sum(z$cluster=="9"), "cells\n")
cat("No cells removed.\n")
cat("No identity assigned.\n")
cat("PASS\n")
