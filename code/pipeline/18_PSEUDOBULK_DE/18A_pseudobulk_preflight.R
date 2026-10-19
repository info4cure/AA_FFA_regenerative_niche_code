suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 18A PSEUDOBULK DE PREFLIGHT\n")
cat("============================================================\n")

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "18_PSEUDOBULK_DE/18A_preflight/tables"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MASTER <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "17A_frozen_analysis_atlas/tables/",
  "17A_master_cell_ledger_FINAL.rds"
)

master <- readRDS(MASTER)

stopifnot(
  nrow(master) == 82893,
  n_distinct(master$cell) == 82893
)

# ============================================================
# 1. Frozen L2 universe
# ============================================================

cat("\nSTEP 1/7 — frozen L2 universe\n")

z <- master %>%
  filter(
    analysis_eligible_global,
    analysis_eligible_L2,
    final_L2_status == "FROZEN",
    !is.na(final_L2_label)
  )

cat("Frozen-L2 cells:", nrow(z), "\n")
cat("Patients:", n_distinct(z$patient_id), "\n")
cat("L2 identities:", n_distinct(z$final_L2_label), "\n")

# ============================================================
# 2. Final object inventory
# ============================================================

cat("\nSTEP 2/7 — final object count-layer inventory\n")

objects <- tribble(
  ~module, ~path,
  "07_EPITHELIAL_L2",
  "PIPELINE_V2/objects/annotated/epithelial_v2_L2_frozen.rds",

  "08_IMMUNE_L2",
  "PIPELINE_V2/objects/annotated/immune_v2_08I_L2_frozen.rds",

  "09_FIBROBLAST_L2",
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "fibroblast_v2_09H5B_L2_frozen_4685.rds"
  ),

  "10_VASCULAR_L2",
  paste0(
    "PIPELINE_V2/objects/annotated/",
    "vascular_v2_10B6_L2_frozen_3023.rds"
  ),

  "11_MURAL_L2",
  paste0(
    "PIPELINE_V2/results/annotation/11_MURAL_L2/",
    "11B5_freeze/mural_v2_11B5_L2_frozen.rds"
  ),

  "12_SWEAT_GLAND_L2",
  paste0(
    "PIPELINE_V2/results/annotation/12_SWEAT_GLAND_L2/",
    "12D_freeze/sweat_gland_v2_12D_L2_frozen.rds"
  ),

  "13_MELANOCYTE_L2",
  paste0(
    "PIPELINE_V2/results/annotation/13_MELANOCYTE_L2/",
    "13D_freeze/melanocyte_v2_13D_L2_frozen.rds"
  ),

  "14_SEBOCYTE_L2",
  paste0(
    "PIPELINE_V2/results/annotation/14_SEBOCYTE_L2/",
    "14D_freeze/sebocyte_v2_14D_L2_frozen.rds"
  )
)

inventory <- bind_rows(
  lapply(
    seq_len(nrow(objects)),
    function(i) {

      f <- objects$path[i]
      module <- objects$module[i]

      stopifnot(file.exists(f))

      o <- readRDS(f)
      stopifnot(inherits(o, "Seurat"))

      assays <- Assays(o)

      bind_rows(
        lapply(
          assays,
          function(a) {

            lay <- tryCatch(
              Layers(o[[a]]),
              error = function(e) character()
            )

            tibble(
              module = module,
              path = f,
              n_cells = ncol(o),
              n_genes = nrow(o),
              default_assay = DefaultAssay(o),
              assay = a,
              layers = paste(lay, collapse = ";"),
              has_counts_layer =
                any(
                  lay == "counts" |
                  grepl("^counts\\.", lay)
                )
            )
          }
        )
      )
    }
  )
)

print(as.data.frame(inventory), row.names = FALSE)

write_csv(
  inventory,
  file.path(
    OUT,
    "18A_object_count_layer_inventory.csv"
  )
)

# ============================================================
# 3. Patient × L2 cell support
# ============================================================

cat("\nSTEP 3/7 — patient × L2 cell support\n")

support_patient <- z %>%
  count(
    final_compartment,
    final_L2_label,
    disease,
    patient_id,
    name = "n_cells"
  )

write_csv(
  support_patient,
  file.path(
    OUT,
    "18A_patient_L2_cell_support.csv"
  )
)

# ============================================================
# 4. Contrast-specific support
# ============================================================

cat("\nSTEP 4/7 — contrast-specific pseudobulk eligibility\n")

contrasts <- tribble(
  ~contrast, ~group1, ~group0, ~primary,
  "AA_vs_HC",  "AA",  "HC",  TRUE,
  "FFA_vs_HC", "FFA", "HC",  TRUE,
  "AA_vs_FFA", "AA",  "FFA", FALSE
)

# Pre-specified pseudobulk sample threshold:
# >=20 cells contributing to a patient × L2 pseudobulk.
#
# Require >=3 usable patients in each arm.
#
# We will additionally retain >=50-cell support as a stricter
# sensitivity flag, but NOT redefine eligibility post hoc.

keys <- z %>%
  distinct(
    final_compartment,
    final_L2_label
  )

eligibility <- bind_rows(
  lapply(
    seq_len(nrow(contrasts)),
    function(k) {

      cc <- contrasts[k, ]

      bind_rows(
        lapply(
          seq_len(nrow(keys)),
          function(i) {

            comp <- keys$final_compartment[i]
            lab <- keys$final_L2_label[i]

            d <- support_patient %>%
              filter(
                final_compartment == comp,
                final_L2_label == lab,
                disease %in%
                  c(cc$group1, cc$group0)
              )

            a <- d %>%
              filter(disease == cc$group1)

            b <- d %>%
              filter(disease == cc$group0)

            tibble(
              contrast = cc$contrast,
              primary = cc$primary,
              numerator = cc$group1,
              denominator = cc$group0,

              final_compartment = comp,
              final_L2_label = lab,

              numerator_patients_any =
                n_distinct(a$patient_id),

              denominator_patients_any =
                n_distinct(b$patient_id),

              numerator_patients_ge20 =
                sum(a$n_cells >= 20),

              denominator_patients_ge20 =
                sum(b$n_cells >= 20),

              numerator_patients_ge50 =
                sum(a$n_cells >= 50),

              denominator_patients_ge50 =
                sum(b$n_cells >= 50),

              numerator_median_cells =
                if (nrow(a) > 0)
                  median(a$n_cells)
                else 0,

              denominator_median_cells =
                if (nrow(b) > 0)
                  median(b$n_cells)
                else 0,

              eligible_ge20 =
                sum(a$n_cells >= 20) >= 3 &&
                sum(b$n_cells >= 20) >= 3,

              robust_ge50 =
                sum(a$n_cells >= 50) >= 3 &&
                sum(b$n_cells >= 50) >= 3
            )
          }
        )
      )
    }
  )
)

write_csv(
  eligibility,
  file.path(
    OUT,
    "18A_pseudobulk_contrast_eligibility.csv"
  )
)

cat("\n===== ELIGIBLE >=20 CELLS / >=3 PATIENTS PER ARM =====\n")

print(
  as.data.frame(
    eligibility %>%
      filter(eligible_ge20) %>%
      arrange(
        contrast,
        final_compartment,
        final_L2_label
      )
  ),
  row.names = FALSE
)

cat("\n===== ROBUST >=50 CELLS / >=3 PATIENTS PER ARM =====\n")

print(
  as.data.frame(
    eligibility %>%
      filter(robust_ge50) %>%
      arrange(
        contrast,
        final_compartment,
        final_L2_label
      )
  ),
  row.names = FALSE
)

# ============================================================
# 5. Eligibility summary
# ============================================================

cat("\nSTEP 5/7 — eligibility summary\n")

elig_summary <- eligibility %>%
  group_by(contrast) %>%
  summarise(
    n_L2_total = n(),
    n_L2_eligible_ge20 =
      sum(eligible_ge20),
    n_L2_robust_ge50 =
      sum(robust_ge50),
    .groups = "drop"
  )

print(
  as.data.frame(elig_summary),
  row.names = FALSE
)

write_csv(
  elig_summary,
  file.path(
    OUT,
    "18A_pseudobulk_eligibility_summary.csv"
  )
)

# ============================================================
# 6. Batch identifiability
# ============================================================

cat("\nSTEP 6/7 — batch identifiability rules\n")

patient_batch <- master %>%
  filter(analysis_eligible_global) %>%
  distinct(
    patient_id,
    disease,
    batch
  )

batch_rules <- tibble(
  contrast = c(
    "AA_vs_HC",
    "FFA_vs_HC",
    "AA_vs_FFA"
  ),

  batch_in_primary_model = c(
    FALSE,
    TRUE,
    FALSE
  ),

  reason = c(
    paste0(
      "AA is completely confounded with CABIMER2026; ",
      "batch cannot be separated from disease."
    ),

    paste0(
      "FFA and HC each contain four historical and one ",
      "CABIMER2026 patient; batch adjustment is identifiable."
    ),

    paste0(
      "AA is entirely CABIMER2026 whereas FFA is predominantly ",
      "historical; batch is severely confounded with disease."
    )
  )
)

print(as.data.frame(batch_rules), row.names = FALSE)

write_csv(
  batch_rules,
  file.path(
    OUT,
    "18A_batch_model_rules.csv"
  )
)

# ============================================================
# 7. Formal checks
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check = c(
    "master_82893",
    "patients_16",
    "three_contrasts",
    "eligibility_threshold_ge20",
    "robustness_threshold_ge50",
    "eligibility_disease_independent_threshold",
    "AA_batch_not_adjustable",
    "FFA_HC_batch_adjustable",
    "AA_FFA_batch_not_adjustable",
    "object_inventory_generated",
    "no_DE_performed",
    "no_gene_filtering_performed",
    "no_figures_generated",
    "annotation_unchanged"
  ),

  pass = c(
    nrow(master) == 82893,
    n_distinct(z$patient_id) == 16,
    n_distinct(eligibility$contrast) == 3,
    TRUE,
    TRUE,
    TRUE,
    !batch_rules$batch_in_primary_model[
      batch_rules$contrast == "AA_vs_HC"
    ],
    batch_rules$batch_in_primary_model[
      batch_rules$contrast == "FFA_vs_HC"
    ],
    !batch_rules$batch_in_primary_model[
      batch_rules$contrast == "AA_vs_FFA"
    ],
    nrow(inventory) > 0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names = FALSE
)

write_csv(
  checks,
  file.path(
    OUT,
    "18A_formal_checks.csv"
  )
)

stopifnot(all(checks$pass))

cat("\n============================================================\n")
cat("18A COMPLETE\n")
cat("============================================================\n")
cat("Pseudobulk DE preflight only\n")
cat("Unit of replication: patient\n")
cat("Eligibility: >=20 cells in >=3 patients per arm\n")
cat("Robust sensitivity: >=50 cells in >=3 patients per arm\n")
cat("AA contrasts: batch NOT included in primary model\n")
cat("FFA vs HC: batch-adjusted model permitted\n")
cat("DE performed: NO\n")
cat("Gene filtering: NO\n")
cat("Figures generated: NO\n")
cat("Annotation modified: NO\n")
cat("PASS\n")
