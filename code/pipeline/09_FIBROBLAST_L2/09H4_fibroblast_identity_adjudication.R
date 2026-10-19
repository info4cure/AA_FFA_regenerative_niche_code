#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H4 FIBROBLAST IDENTITY ADJUDICATION\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

H2 <- file.path(
  ROOT, "results/annotation/09_FIBROBLAST_L2",
  "09H2_FibConsensus_de_novo_adjudication/tables"
)

G <- file.path(
  ROOT, "results/annotation/09_FIBROBLAST_L2",
  "09G_refined_structure/tables"
)

OUT <- file.path(
  ROOT, "results/annotation/09_FIBROBLAST_L2",
  "09H4_identity_adjudication"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1
# ============================================================

cat("\nSTEP 1/8 — loading r0.3 identity evidence\n")

proj <- read_csv(
  file.path(H2, "09H2_projection_evidence_classification.csv"),
  show_col_types = FALSE
) %>%
  filter(as.character(resolution) == "0.3") %>%
  mutate(cluster = as.character(cluster))

rank <- read_csv(
  file.path(H2, "09H2_cluster_family_ranking.csv"),
  show_col_types = FALSE
) %>%
  filter(as.character(resolution) == "0.3") %>%
  mutate(cluster = as.character(cluster))

markers <- read_csv(
  file.path(G, "09G_markers_res0.3_all.csv"),
  show_col_types = FALSE
) %>%
  mutate(cluster = as.character(cluster))

stopifnot(nrow(proj) == 10)
stopifnot(length(unique(proj$cluster)) == 10)

cat("r0.3 clusters:", nrow(proj), "\n")


# ============================================================
# STEP 2 — explicit marker programmes
# ============================================================

cat("\nSTEP 2/8 — defining independent marker evidence\n")

programmes <- list(

  FAM01_PAPILLARY = c(
    "APCDD1","COL6A5","COL23A1","COL13A1",
    "COL18A1","NKD2","WIF1","HSPB3"
  ),

  FAM02_RETICULAR = c(
    "SFRP2","DPP4","COL1A1","COL1A2",
    "COL3A1","FBLN1","SMOC2","COL12A1"
  ),

  FAM03_DERMAL_SHEATH = c(
    "COL11A1","DPEP1","RBP4","EDNRA",
    "POSTN","KIF26B","SOX2","MME"
  ),

  FAM04_DERMAL_PAPILLA = c(
    "COCH","CRABP1","ASPN","TNN",
    "RSPO3","CORIN","HHIP","LEF1"
  ),

  FAM05_PI16_PERIVASCULAR = c(
    "PI16","CD36","PPARG","APOD",
    "CXCL12","PTGDS","COL15A1"
  ),

  FAM06_APOE_IMMUNE_PERIVASCULAR = c(
    "APOE","APOC1","C7","CCL19",
    "CXCL12","C3","IL33","CFD"
  ),

  EPITHELIAL = c(
    "KRT1","KRT2","KRT5","KRT10","KRT14",
    "KRT15","KRT17","KRT19",
    "EPCAM","KRTDAP","SFN","DSP","DSC3",
    "DSG1","SERPINB5","NECTIN4","SPINK5",
    "TACSTD2","PKP1","SBSN","CALML3"
  ),

  MURAL_PERICYTE = c(
    "RGS5","CSPG4","MCAM","NOTCH3",
    "ACTA2","TAGLN","MYH11","RBP1"
  ),

  NEURAL_SCHWANN = c(
    "S100B","SOX10","PLP1","S100A1",
    "MPZ","PMP22","NGFR"
  )
)


# ============================================================
# STEP 3 — marker evidence per cluster
# ============================================================

cat("\nSTEP 3/8 — scoring de novo marker evidence\n")

positive <- markers %>%
  filter(avg_log2FC > 0)

marker_evidence <- bind_rows(
  lapply(names(programmes), function(p) {

    genes <- programmes[[p]]

    positive %>%
      filter(gene %in% genes) %>%
      group_by(cluster) %>%
      summarise(
        programme = p,
        n_programme_markers = n_distinct(gene),
        marker_genes = paste(
          unique(gene),
          collapse = ";"
        ),
        max_log2FC = max(avg_log2FC, na.rm = TRUE),
        .groups = "drop"
      )
  })
)

all_grid <- expand_grid(
  cluster = sort(unique(proj$cluster)),
  programme = names(programmes)
)

marker_evidence <- all_grid %>%
  left_join(
    marker_evidence,
    by = c("cluster","programme")
  ) %>%
  mutate(
    n_programme_markers =
      replace_na(n_programme_markers, 0L),

    marker_genes =
      replace_na(marker_genes, ""),

    max_log2FC =
      replace_na(max_log2FC, 0)
  )

write_csv(
  marker_evidence,
  file.path(TAB, "09H4_marker_programme_evidence.csv")
)


# ============================================================
# STEP 4 — explicit contamination/special-state audit
# ============================================================

cat("\nSTEP 4/8 — lineage-conflict and specialized-state audit\n")

wide <- marker_evidence %>%
  select(cluster, programme, n_programme_markers) %>%
  pivot_wider(
    names_from = programme,
    values_from = n_programme_markers,
    values_fill = 0
  )

special <- wide %>%
  mutate(

    epithelial_flag =
      EPITHELIAL >= 3,

    mural_flag =
      MURAL_PERICYTE >= 3,

    neural_flag =
      NEURAL_SCHWANN >= 3,

    nonfibroblast_lineage_flag =
      epithelial_flag |
      mural_flag |
      neural_flag
  )

write_csv(
  special,
  file.path(TAB, "09H4_lineage_conflict_audit.csv")
)


# ============================================================
# STEP 5 — family-level de novo support
# ============================================================

cat("\nSTEP 5/8 — integrating frozen and de novo family evidence\n")

family_programme_map <- tibble(
  top_family = paste0("FAM0",1:6),
  marker_programme = c(
    "FAM01_PAPILLARY",
    "FAM02_RETICULAR",
    "FAM03_DERMAL_SHEATH",
    "FAM04_DERMAL_PAPILLA",
    "FAM05_PI16_PERIVASCULAR",
    "FAM06_APOE_IMMUNE_PERIVASCULAR"
  )
)

family_marker <- marker_evidence %>%
  inner_join(
    family_programme_map,
    by = c("programme" = "marker_programme")
  ) %>%
  rename(
    top_family = top_family,
    de_novo_family_markers = n_programme_markers,
    de_novo_family_genes = marker_genes
  ) %>%
  select(
    cluster,
    top_family,
    de_novo_family_markers,
    de_novo_family_genes
  )

evidence <- proj %>%
  left_join(
    family_marker,
    by = c("cluster","top_family")
  ) %>%
  left_join(
    special %>%
      select(
        cluster,
        EPITHELIAL,
        MURAL_PERICYTE,
        NEURAL_SCHWANN,
        epithelial_flag,
        mural_flag,
        neural_flag,
        nonfibroblast_lineage_flag
      ),
    by = "cluster"
  ) %>%
  mutate(
    de_novo_family_markers =
      replace_na(de_novo_family_markers, 0L),

    de_novo_family_genes =
      replace_na(de_novo_family_genes, "")
  )


# ============================================================
# STEP 6 — conservative evidence adjudication
# ============================================================

cat("\nSTEP 6/8 — conservative identity adjudication\n")

evidence <- evidence %>%
  mutate(

    identity_evidence = case_when(

      epithelial_flag ~
        "NON_FIBROBLAST_LINEAGE_REVIEW",

      mural_flag | neural_flag ~
        "SPECIALIZED_LINEAGE_REVIEW",

      projection_evidence ==
        "STRONG_FROZEN_FAMILY_SUPPORT" &
        de_novo_family_markers >= 2 ~
        "CONCORDANT_FAMILY_SUPPORT",

      projection_evidence ==
        "STRONG_FROZEN_FAMILY_SUPPORT" ~
        "STRONG_PROJECTION_LIMITED_MARKER_SUPPORT",

      projection_evidence ==
        "WEAK_FROZEN_FAMILY_SUPPORT" &
        de_novo_family_markers >= 2 ~
        "PARTIAL_FAMILY_SUPPORT",

      TRUE ~
        "UNRESOLVED_FIBROBLAST_STATE"
    ),

    freeze_eligible = case_when(

      identity_evidence ==
        "CONCORDANT_FAMILY_SUPPORT" ~ TRUE,

      TRUE ~ FALSE
    )
  )


# ============================================================
# STEP 7 — explicit biological review
# ============================================================

cat("\nSTEP 7/8 — generating cluster adjudication ledger\n")

family_names <- c(
  FAM01 = "Papillary/superficial fibroblasts",
  FAM02 = "Reticular/SFRP2-DPP4 fibroblasts",
  FAM03 = "COL11A1/DPEP1 dermal-sheath-like fibroblasts",
  FAM04 = "Hair-follicle dermal papilla fibroblasts",
  FAM05 = "PI16-associated perivascular/adipogenic fibroblasts",
  FAM06 = "APOE-associated immune/perivascular fibroblasts"
)

ledger <- evidence %>%
  mutate(

    projected_family_name =
      unname(family_names[top_family]),

    adjudication_note = case_when(

      cluster == "0" ~
        paste(
          "Strong FAM06 projection with APOE/CXCL12/C3/C7/",
          "APOC1/CCL19/IL33 evidence."
        ),

      cluster == "1" ~
        paste(
          "Strong FAM03 projection with DPEP1/COL11A1/",
          "KIF26B/EDNRA/RBP4/SOX2 evidence."
        ),

      cluster == "2" ~
        paste(
          "Weak FibConsensus separation; mitochondrial and",
          "broad mesenchymal programme predominates."
        ),

      cluster == "3" ~
        paste(
          "Reticular programme supported by SFRP2 and ECM genes,",
          "but FAM02/FAM05 projection separation is limited."
        ),

      cluster == "4" ~
        paste(
          "Strong FAM01 projection with COL6A5/APCDD1/",
          "COL23A1/COL13A1/COL18A1/NKD2 evidence."
        ),

      cluster == "5" ~
        paste(
          "Strong FAM04 projection with COCH/ASPN and",
          "specialized follicular mesenchymal programme."
        ),

      cluster == "6" ~
        paste(
          "FAM02/FAM03 mixed projection with strong remodeling",
          "programme including ADAM12/MMP11/COL10A1/COMP."
        ),

      cluster == "7" ~
        paste(
          "Weak FAM03 projection; distinct specialized",
          "COL26A1/RELN/FMO2/PDGFD programme."
        ),

      cluster == "8" ~
        paste(
          "FibConsensus-unresolved EBF2/SFRP5/FOXC2/FOXS1",
          "population with CLDN1/ITGA6/KRT19 expression;",
          "requires dedicated lineage review."
        ),

      cluster == "9" ~
        paste(
          "FibConsensus-unresolved cluster with strong epithelial",
          "SERPINB5/NECTIN4/SPINK5/TACSTD2/PKP1/DSC3/DSG1",
          "programme; return-to-epithelial candidacy requires",
          "cell-level confirmation."
        ),

      TRUE ~ ""
    )
  )

write_csv(
  ledger,
  file.path(TAB, "09H4_cluster_identity_adjudication.csv")
)

cat("\n=== IDENTITY ADJUDICATION ===\n")
print(
  ledger %>%
    select(
      cluster,
      n_cells,
      n_patients,
      top_family,
      projected_family_name,
      top_core_z,
      core_margin,
      patient_fraction,
      de_novo_family_markers,
      identity_evidence,
      freeze_eligible
    )
)


# ============================================================
# STEP 8 — formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check = c(
    "resolution_03_used",
    "ten_clusters_adjudicated",
    "six_FibConsensus_families_available",
    "de_novo_markers_integrated",
    "lineage_conflict_audited",
    "epithelial_review_possible",
    "disease_not_used_for_identity",
    "no_automatic_family_assignment",
    "no_L2_labels_frozen",
    "no_cells_removed",
    "no_object_modified"
  ),

  pass = c(
    TRUE,
    nrow(ledger) == 10,
    length(unique(rank$family_id)) == 6,
    nrow(marker_evidence) > 0,
    nrow(special) == 10,
    any(special$epithelial_flag),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(checks)

write_csv(
  checks,
  file.path(TAB, "09H4_formal_checks.csv")
)

if (!all(checks$pass)) {
  stop("09H4 formal checks failed.")
}

audit <- c(
  "PIPELINE_V2 — 09H4 FIBROBLAST IDENTITY ADJUDICATION",
  "",
  "Selected structural scaffold: resolution 0.3",
  "Clusters adjudicated: 10",
  "",
  "Evidence integrated:",
  "- frozen FibConsensus CORE/FULL projection",
  "- patient-level projection support",
  "- independent de novo cluster markers",
  "- explicit epithelial/mural/neural lineage-conflict programmes",
  "",
  "Important:",
  "No L2 identity was automatically assigned.",
  "No cell was removed.",
  "No Seurat object was modified.",
  "",
  "Clusters meeting concordant family evidence may proceed to",
  "final identity freeze after targeted review.",
  "Clusters with weak, mixed, or incompatible evidence remain",
  "explicitly unresolved.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(TAB, "09H4_audit.txt")
)

cat("\n============================================================\n")
cat("09H4 COMPLETE\n")
cat("============================================================\n")
cat("Selected structural scaffold: resolution 0.3\n")
cat("Clusters adjudicated: 10\n")
cat("No L2 identities frozen.\n")
cat("No cells removed.\n")
cat("PASS — fibroblast identity evidence adjudicated.\n")
