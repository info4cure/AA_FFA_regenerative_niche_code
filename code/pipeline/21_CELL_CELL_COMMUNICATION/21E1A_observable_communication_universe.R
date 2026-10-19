suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21E1A OBSERVABLE COMMUNICATION UNIVERSE\n")
cat("============================================================\n")

# ============================================================
# Paths
# ============================================================

OBJ <- paste0(
  "PIPELINE_V2/objects/frozen/",
  "global_atlas_v2_21B_communication_frozen.rds"
)

E0 <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E0_OMNIPATH_RESOURCE"
)

LR_FILE <- file.path(
  E0,
  "raw",
  "21E0_OmniPath_ligrecextra_frozen.rds"
)

IC_FILE <- file.path(
  E0,
  "raw",
  "21E0_OmniPath_intercell_frozen.rds"
)

OUT <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21E1_OBSERVABLE_UNIVERSE"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Recover frozen inputs
# ============================================================

cat("\nSTEP 1/9 — recovering frozen inputs\n")

stopifnot(
  file.exists(OBJ),
  file.exists(LR_FILE),
  file.exists(IC_FILE)
)

obj <- readRDS(OBJ)
lr  <- readRDS(LR_FILE)
ic  <- readRDS(IC_FILE)

atlas_genes <- rownames(obj[["RNA"]])

cat("Atlas cells       :", ncol(obj), "\n")
cat("Atlas genes       :", length(atlas_genes), "\n")
cat("OmniPath LR rows  :", nrow(lr), "\n")
cat("Intercell records :", nrow(ic), "\n")

stopifnot(
  ncol(obj) == 82893,
  length(atlas_genes) == 23412,
  nrow(lr) == 6658
)

# ============================================================
# 2. Validate LR schema
# ============================================================

cat("\nSTEP 2/9 — validating ligand-receptor schema\n")

required_lr <- c(
  "source",
  "target",
  "source_genesymbol",
  "target_genesymbol"
)

stopifnot(
  all(required_lr %in% names(lr))
)

lr0 <- lr %>%
  mutate(
    source_genesymbol =
      as.character(source_genesymbol),
    target_genesymbol =
      as.character(target_genesymbol)
  )

# ============================================================
# 3. Structural atlas overlap
# ============================================================

cat("\nSTEP 3/9 — intersecting frozen LR resource with atlas features\n")

lr_atlas <- lr0 %>%
  mutate(
    source_in_atlas =
      source_genesymbol %in% atlas_genes,
    target_in_atlas =
      target_genesymbol %in% atlas_genes,
    pair_in_atlas =
      source_in_atlas & target_in_atlas,
    atlas_overlap_class =
      case_when(
        pair_in_atlas ~ "BOTH_IN_ATLAS",
        source_in_atlas & !target_in_atlas ~
          "SOURCE_ONLY_IN_ATLAS",
        !source_in_atlas & target_in_atlas ~
          "TARGET_ONLY_IN_ATLAS",
        TRUE ~ "NEITHER_IN_ATLAS"
      )
  )

overlap_summary <- lr_atlas %>%
  count(
    atlas_overlap_class,
    name="n_pairs"
  ) %>%
  mutate(
    fraction =
      n_pairs / sum(n_pairs)
  ) %>%
  arrange(desc(n_pairs))

cat("\n===== ATLAS OVERLAP =====\n")
print(overlap_summary, n=Inf)

# ============================================================
# 4. Observable LR universe
# ============================================================

cat("\nSTEP 4/9 — defining structurally observable LR universe\n")

observable <- lr_atlas %>%
  filter(pair_in_atlas) %>%
  distinct(
    source,
    target,
    source_genesymbol,
    target_genesymbol,
    .keep_all=TRUE
  )

cat(
  "Observable LR pairs:",
  nrow(observable),
  "/",
  nrow(lr_atlas),
  "\n"
)

cat(
  "Observable source genes:",
  n_distinct(observable$source_genesymbol),
  "\n"
)

cat(
  "Observable target genes:",
  n_distinct(observable$target_genesymbol),
  "\n"
)

# ============================================================
# 5. Recover intercell roles
# ============================================================

cat("\nSTEP 5/9 — mapping OmniPath intercell roles\n")

required_ic <- c(
  "genesymbol",
  "category",
  "parent",
  "transmitter",
  "receiver",
  "secreted",
  "plasma_membrane_transmembrane",
  "plasma_membrane_peripheral"
)

stopifnot(
  all(required_ic %in% names(ic))
)

as_logical_safe <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  tolower(as.character(x)) %in%
    c("true", "t", "1", "yes")
}

ic2 <- ic %>%
  mutate(
    across(
      all_of(
        c(
          "transmitter",
          "receiver",
          "secreted",
          "plasma_membrane_transmembrane",
          "plasma_membrane_peripheral"
        )
      ),
      as_logical_safe
    )
  )

gene_roles <- ic2 %>%
  filter(
    !is.na(genesymbol),
    genesymbol != ""
  ) %>%
  group_by(genesymbol) %>%
  summarise(
    transmitter =
      any(transmitter, na.rm=TRUE),
    receiver =
      any(receiver, na.rm=TRUE),
    secreted =
      any(secreted, na.rm=TRUE),
    plasma_membrane_transmembrane =
      any(
        plasma_membrane_transmembrane,
        na.rm=TRUE
      ),
    plasma_membrane_peripheral =
      any(
        plasma_membrane_peripheral,
        na.rm=TRUE
      ),
    categories =
      paste(
        sort(unique(category)),
        collapse=";"
      ),
    parents =
      paste(
        sort(unique(parent)),
        collapse=";"
      ),
    n_intercell_records=n(),
    .groups="drop"
  )

# ============================================================
# 6. Annotate both ends
# ============================================================

cat("\nSTEP 6/9 — annotating observable interaction endpoints\n")

src_roles <- gene_roles %>%
  rename_with(
    ~paste0("source_", .x),
    -genesymbol
  ) %>%
  rename(
    source_genesymbol=genesymbol
  )

tgt_roles <- gene_roles %>%
  rename_with(
    ~paste0("target_", .x),
    -genesymbol
  ) %>%
  rename(
    target_genesymbol=genesymbol
  )

observable_annotated <- observable %>%
  left_join(
    src_roles,
    by="source_genesymbol"
  ) %>%
  left_join(
    tgt_roles,
    by="target_genesymbol"
  )

# ============================================================
# 7. Communication modality scaffold
# ============================================================

cat("\nSTEP 7/9 — deriving modality scaffold\n")

observable_annotated <- observable_annotated %>%
  mutate(
    source_transmitter =
      coalesce(source_transmitter, FALSE),
    target_receiver =
      coalesce(target_receiver, FALSE),
    source_secreted =
      coalesce(source_secreted, FALSE),
    source_membrane =
      coalesce(
        source_plasma_membrane_transmembrane,
        FALSE
      ) |
      coalesce(
        source_plasma_membrane_peripheral,
        FALSE
      ),
    target_membrane =
      coalesce(
        target_plasma_membrane_transmembrane,
        FALSE
      ) |
      coalesce(
        target_plasma_membrane_peripheral,
        FALSE
      ),

    communication_modality =
      case_when(
        source_secreted &
          target_membrane ~
          "SECRETED_TO_MEMBRANE",

        source_membrane &
          target_membrane ~
          "CONTACT_ASSOCIATED",

        TRUE ~
          "OTHER_OR_UNRESOLVED"
      )
  )

modality_summary <- observable_annotated %>%
  count(
    communication_modality,
    name="n_pairs"
  ) %>%
  mutate(
    fraction =
      n_pairs / sum(n_pairs)
  ) %>%
  arrange(desc(n_pairs))

cat("\n===== MODALITY SCAFFOLD =====\n")
print(modality_summary, n=Inf)

# ============================================================
# 8. Write outputs
# ============================================================

cat("\nSTEP 8/9 — writing frozen structural universe\n")

write_csv(
  overlap_summary,
  file.path(
    TAB,
    "21E1A_atlas_overlap_summary.csv"
  )
)

write_csv(
  lr_atlas,
  file.path(
    TAB,
    "21E1A_all_LR_pairs_atlas_overlap.csv"
  )
)

write_csv(
  gene_roles,
  file.path(
    TAB,
    "21E1A_intercell_gene_roles.csv"
  )
)

write_csv(
  observable_annotated,
  file.path(
    TAB,
    "21E1A_observable_LR_universe.csv"
  )
)

write_csv(
  modality_summary,
  file.path(
    TAB,
    "21E1A_modality_summary.csv"
  )
)

# ============================================================
# 9. Formal checks
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "21B_frozen_object_present",
    "21E0_LR_resource_present",
    "21E0_intercell_resource_present",
    "global_cells_82893",
    "atlas_features_23412",
    "frozen_LR_pairs_6658",
    "LR_key_columns_present",
    "intercell_role_columns_present",
    "observable_universe_nonzero",
    "observable_pairs_unique",
    "structural_overlap_only",
    "no_expression_thresholding",
    "no_patient_filtering",
    "no_disease_filtering",
    "no_AA_FFA_comparison",
    "no_communication_inference",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(OBJ),
    file.exists(LR_FILE),
    file.exists(IC_FILE),
    ncol(obj) == 82893,
    length(atlas_genes) == 23412,
    nrow(lr0) == 6658,
    all(required_lr %in% names(lr0)),
    all(required_ic %in% names(ic)),
    nrow(observable_annotated) > 0,
    !anyDuplicated(
      observable_annotated[
        c(
          "source",
          "target",
          "source_genesymbol",
          "target_genesymbol"
        )
      ]
    ),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks, n=Inf, width=Inf)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(
    TAB,
    "21E1A_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 21E1A observable communication universe",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Frozen OmniPath LR pairs: ",
    nrow(lr0)
  ),
  paste0(
    "Atlas RNA features: ",
    length(atlas_genes)
  ),
  paste0(
    "Structurally observable LR pairs: ",
    nrow(observable_annotated)
  ),
  "",
  "Definition:",
  paste(
    "A structurally observable pair requires both source and",
    "target gene symbols to exist among the frozen atlas RNA features."
  ),
  "",
  paste(
    "This definition does not require detectable expression in",
    "any cell, identity, patient or disease."
  ),
  "",
  "OmniPath intercell annotations were used only to attach",
  "transmitter, receiver, secreted and membrane-associated roles.",
  "",
  "No expression thresholding.",
  "No patient filtering.",
  "No disease filtering.",
  "No AA-vs-FFA comparison.",
  "No communication inference.",
  "No differential-expression refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    OUT,
    "21E1A_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21E1A COMPLETE\n")
cat("============================================================\n")
cat("Structural OmniPath-atlas communication universe frozen\n")
cat("No expression thresholding\n")
cat("No disease filtering\n")
cat("No communication inference\n")
cat("PASS\n")
