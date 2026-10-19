
# ============================================================
# 22G4C2B — TARGETED IMMUNE FUNCTIONAL SCORING
#
# Frozen source:
# 18B patient x frozen-L2 pseudobulk counts
#
# Populations:
# CD4 T
# CD8 T
# NK
# Macrophage/monocyte
#
# Disease universe:
# AA + HC
#
# Functional panel:
# 7 pre-specified Hallmarks
# + targeted cytotoxicity score for CD8/NK
#
# IMPORTANT:
# - No L2 reannotation
# - No pseudobulk reconstruction
# - No disease testing in this module
# - n_cells retained for downstream sensitivity analyses
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

cat("\n============================================================\n")
cat("22G4C2B — TARGETED IMMUNE FUNCTIONAL SCORING\n")
cat("============================================================\n")

ROOT <- "PIPELINE_V2"

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G4_IMMUNE_ECOSYSTEM",
  "22G4C_AA_HC_LOWER_FOLLICLE",
  "22G4C2_IMMUNE_FUNCTION"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJ,recursive=TRUE,showWarnings=FALSE)
dir.create(AUD,recursive=TRUE,showWarnings=FALSE)

PB <- file.path(
  ROOT,
  "results/disease_biology",
  "18_PSEUDOBULK_DE",
  "18B_pseudobulk"
)

COUNT_FILE <- file.path(
  PB,
  "objects",
  "18B_08_IMMUNE_L2_pseudobulk_counts.rds"
)

META_FILE <- file.path(
  PB,
  "objects",
  "18B_pseudobulk_metadata_ALL.rds"
)

GSFILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_PATHWAY_ANALYSIS",
  "20A_HALLMARK_GSEA",
  "gene_sets",
  "20A_Hallmark_gene_sets_frozen.rds"
)

stopifnot(
  file.exists(COUNT_FILE),
  file.exists(META_FILE),
  file.exists(GSFILE)
)

# ============================================================
# 1. LOAD FROZEN SOURCES
# ============================================================

counts_all <- readRDS(
  COUNT_FILE
)

md_all <- readRDS(
  META_FILE
) |>
  as_tibble()

hallmarks_all <- readRDS(
  GSFILE
)

cat("\n===== FROZEN INPUT =====\n")
cat("Genes:",nrow(counts_all),"\n")
cat("Immune pseudobulks:",ncol(counts_all),"\n")
cat("Hallmarks:",length(hallmarks_all),"\n")

# ============================================================
# 2. TARGET DEFINITIONS
# ============================================================

TARGET_L2 <- c(
  "CD4 T",
  "CD8 T",
  "NK",
  "Macrophage/monocyte"
)

HALLMARK_PANEL <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_APOPTOSIS"
)

CYTOTOXICITY_GENES <- c(
  "NKG7",
  "GNLY",
  "PRF1",
  "GZMB",
  "GZMH",
  "GZMK",
  "CTSW",
  "CCL5"
)

stopifnot(
  all(
    HALLMARK_PANEL %in%
      names(hallmarks_all)
  )
)

hallmarks <- hallmarks_all[
  HALLMARK_PANEL
]

# ============================================================
# 3. TARGET METADATA — AA + HC ONLY
# ============================================================

md <- md_all |>
  filter(
    module=="08_IMMUNE_L2",
    disease %in%c("HC","AA"),
    final_L2_label %in%TARGET_L2
  ) |>
  arrange(
    final_L2_label,
    disease,
    patient_id
  )

stopifnot(
  all(
    md$pseudobulk_id %in%
      colnames(counts_all)
  )
)

cat("\n===== TARGET PSEUDOBULKS =====\n")

print(
  md |>
    count(
      final_L2_label,
      disease,
      name="n_pseudobulks"
    ),
  n=Inf,
  width=Inf
)

# ============================================================
# 4. HALLMARK GSVA — SEPARATELY WITHIN EACH L2
# ============================================================

cat("\n===== HALLMARK GSVA =====\n")

hallmark_list <- list()
logcpm_list <- list()

counter <- 0L

for(l2 in TARGET_L2) {

  cat("\n---",l2,"---\n")

  submd <- md |>
    filter(
      final_L2_label==l2
    )

  ids <- submd$pseudobulk_id

  counts <- counts_all[
    ,
    ids,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(counts),
      ids
    )
  )

  # Same principle as frozen 20C.
  keep <- Matrix::rowSums(counts) > 0

  counts <- counts[
    keep,
    ,
    drop=FALSE
  ]

  y <- edgeR::DGEList(
    counts=counts
  )

  y <- edgeR::calcNormFactors(
    y,
    method="TMM"
  )

  logcpm <- edgeR::cpm(
    y,
    log=TRUE,
    prior.count=2
  )

  stopifnot(
    !anyDuplicated(
      rownames(logcpm)
    )
  )

  param <- GSVA::gsvaParam(
    exprData=logcpm,
    geneSets=hallmarks,
    minSize=10,
    maxSize=500,
    kcdf="Gaussian"
  )

  scores <- GSVA::gsva(
    param,
    verbose=FALSE
  )

  stopifnot(
    ncol(scores)==nrow(submd),
    identical(
      colnames(scores),
      submd$pseudobulk_id
    )
  )

  counter <- counter+1L

  hallmark_list[[counter]] <-
    as.data.frame(scores) |>
    rownames_to_column(
      "programme"
    ) |>
    pivot_longer(
      cols=-programme,
      names_to="pseudobulk_id",
      values_to="score"
    ) |>
    left_join(
      submd |>
        dplyr::select(
          pseudobulk_id,
          patient_id,
          disease,
          batch,
          final_L2_label,
          n_cells
        ),
      by="pseudobulk_id"
    ) |>
    mutate(
      score_type="HALLMARK",
      .before=1
    )

  logcpm_list[[l2]] <- logcpm

  cat(
    "pseudobulks:",
    ncol(logcpm),
    "| genes:",
    nrow(logcpm),
    "\n"
  )
}

hallmark_scores <- bind_rows(
  hallmark_list
)

# ============================================================
# 5. CYTOTOXICITY SCORE — CD8 + NK
# ============================================================

cat("\n===== CYTOTOXICITY SCORE =====\n")

cyto_list <- list()

for(l2 in c("CD8 T","NK")) {

  logcpm <- logcpm_list[[l2]]

  genes_present <- intersect(
    CYTOTOXICITY_GENES,
    rownames(logcpm)
  )

  cat(
    l2,
    ":",
    length(genes_present),
    "/",
    length(CYTOTOXICITY_GENES),
    "genes present\n"
  )

  print(
    genes_present
  )

  stopifnot(
    length(genes_present)>=5
  )

  expr <- logcpm[
    genes_present,
    ,
    drop=FALSE
  ]

  # Gene-wise standardization across patient pseudobulks
  # within the same frozen immune L2.
  z <- t(
    scale(
      t(expr)
    )
  )

  # Protect against zero-variance genes.
  z[
    !is.finite(z)
  ] <- 0

  cyto_score <- colMeans(
    z,
    na.rm=TRUE
  )

  submd <- md |>
    filter(
      final_L2_label==l2
    )

  cyto_list[[l2]] <- tibble(
    score_type="CUSTOM",
    programme="CYTOTOXICITY",
    pseudobulk_id=
      names(cyto_score),
    score=
      as.numeric(
        cyto_score
      )
  ) |>
    left_join(
      submd |>
        dplyr::select(
          pseudobulk_id,
          patient_id,
          disease,
          batch,
          final_L2_label,
          n_cells
        ),
      by="pseudobulk_id"
    )
}

cytotoxicity_scores <- bind_rows(
  cyto_list
)

# ============================================================
# 6. INTEGRATED FUNCTIONAL SCORE TABLE
# ============================================================

functional_scores <- bind_rows(
  hallmark_scores,
  cytotoxicity_scores
) |>
  mutate(
    support_tier=case_when(
      n_cells >= 20
        ~ "HIGH",
      n_cells >= 10
        ~ "MODERATE",
      n_cells >= 5
        ~ "LOW",
      TRUE
        ~ "VERY_LOW"
    ),

    analysis_population=case_when(
      final_L2_label %in%
        c(
          "CD4 T",
          "Macrophage/monocyte"
        )
        ~ "TIER1",

      final_L2_label %in%
        c(
          "CD8 T",
          "NK"
        )
        ~ "TIER2",

      TRUE
        ~ "OTHER"
    )
  ) |>
  arrange(
    final_L2_label,
    programme,
    disease,
    patient_id
  )

# ============================================================
# 7. SUPPORT LEDGER
# ============================================================

support_ledger <- functional_scores |>
  distinct(
    patient_id,
    disease,
    batch,
    final_L2_label,
    n_cells,
    support_tier,
    analysis_population
  ) |>
  group_by(
    final_L2_label,
    disease
  ) |>
  summarise(
    n_patients=n(),
    total_cells=sum(n_cells),
    min_cells=min(n_cells),
    median_cells=median(n_cells),
    max_cells=max(n_cells),

    n_HIGH=
      sum(
        support_tier=="HIGH"
      ),

    n_MODERATE=
      sum(
        support_tier=="MODERATE"
      ),

    n_LOW=
      sum(
        support_tier=="LOW"
      ),

    n_VERY_LOW=
      sum(
        support_tier=="VERY_LOW"
      ),

    .groups="drop"
  )

cat("\n===== SUPPORT LEDGER =====\n")

print(
  support_ledger,
  n=Inf,
  width=Inf
)

# ============================================================
# 8. SCORE COMPLETENESS
# ============================================================

score_audit <- functional_scores |>
  group_by(
    final_L2_label,
    programme,
    score_type
  ) |>
  summarise(
    n_scores=n(),
    n_finite=
      sum(
        is.finite(score)
      ),

    n_HC=
      sum(
        disease=="HC"
      ),

    n_AA=
      sum(
        disease=="AA"
      ),

    .groups="drop"
  )

cat("\n===== SCORE COMPLETENESS =====\n")

print(
  score_audit,
  n=Inf,
  width=Inf
)

# ============================================================
# 9. SAVE TABLES
# ============================================================

write.csv(
  functional_scores,
  file.path(
    TAB,
    "22G4C2B_targeted_immune_functional_scores.csv"
  ),
  row.names=FALSE
)

write.csv(
  support_ledger,
  file.path(
    TAB,
    "22G4C2B_support_ledger.csv"
  ),
  row.names=FALSE
)

write.csv(
  score_audit,
  file.path(
    TAB,
    "22G4C2B_score_completeness.csv"
  ),
  row.names=FALSE
)

write.csv(
  tibble(
    programme="CYTOTOXICITY",
    gene=CYTOTOXICITY_GENES
  ),
  file.path(
    TAB,
    "22G4C2B_cytotoxicity_gene_set.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 10. FREEZE
# ============================================================

freeze <- list(
  module="22G4C2B",

  purpose=paste(
    "Targeted functional scoring of frozen AA/HC immune",
    "patient-level pseudobulks including NK."
  ),

  source_counts=
    COUNT_FILE,

  source_metadata=
    META_FILE,

  source_hallmarks=
    GSFILE,

  populations=
    TARGET_L2,

  Hallmark_panel=
    HALLMARK_PANEL,

  cytotoxicity_genes=
    CYTOTOXICITY_GENES,

  normalization=
    "TMM + logCPM(prior.count=2), separately within each L2",

  Hallmark_scoring=
    "GSVA gsvaParam; Gaussian; minSize=10; maxSize=500",

  cytotoxicity_scoring=paste(
    "Mean gene-wise z-scored logCPM within frozen L2;",
    "CD8 T and NK only."
  ),

  disease_testing=
    "NONE",

  functional_scores=
    functional_scores,

  support_ledger=
    support_ledger
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G4C2B_targeted_immune_functional_scoring_freeze.rds"
  )
)

# ============================================================
# 11. FORMAL CHECKS
# ============================================================

checks <- tibble(
  check=c(
    "Frozen 18B counts used",
    "Four target immune populations",
    "CD4 represented",
    "CD8 represented",
    "NK represented",
    "Macrophage represented",
    "Seven Hallmarks scored",
    "CD8 cytotoxicity scored",
    "NK cytotoxicity scored",
    "All scores finite",
    "n_cells retained",
    "No disease testing performed",
    "Freeze exists"
  ),

  pass=c(
    file.exists(
      COUNT_FILE
    ),

    n_distinct(
      functional_scores$final_L2_label
    )==4,

    "CD4 T" %in%
      functional_scores$final_L2_label,

    "CD8 T" %in%
      functional_scores$final_L2_label,

    "NK" %in%
      functional_scores$final_L2_label,

    "Macrophage/monocyte" %in%
      functional_scores$final_L2_label,

    n_distinct(
      hallmark_scores$programme
    )==7,

    any(
      cytotoxicity_scores$final_L2_label==
        "CD8 T"
    ),

    any(
      cytotoxicity_scores$final_L2_label==
        "NK"
    ),

    all(
      is.finite(
        functional_scores$score
      )
    ),

    "n_cells" %in%
      names(
        functional_scores
      ),

    TRUE,

    file.exists(
      file.path(
        OBJ,
        "22G4C2B_targeted_immune_functional_scoring_freeze.rds"
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
    "22G4C2B_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 12. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G4C2B — TARGETED IMMUNE FUNCTIONAL SCORING",
  "============================================================",
  "",
  "Disease universe: AA + HC",
  "",
  "Frozen immune populations:",
  "- CD4 T",
  "- CD8 T",
  "- NK",
  "- Macrophage/monocyte",
  "",
  "Seven pre-specified Hallmarks scored.",
  "Custom cytotoxicity scored for CD8 T and NK.",
  "",
  paste0(
    "Functional score rows: ",
    nrow(
      functional_scores
    )
  ),
  "",
  "Support tiers:",
  "HIGH: >=20 cells
",
  "MODERATE: 10-19 cells",
  "LOW: 5-9 cells",
  "VERY_LOW: <5 cells",
  "",
  "Primary statistical tiers:",
  "TIER1: CD4 T and Macrophage/monocyte",
  "TIER2: CD8 T and NK",
  "",
  paste(
    "NK is retained as a biologically prioritized population",
    "despite low patient-level cell counts."
  ),
  paste(
    "All NK analyses must preserve n_cells and include",
    "cell-support sensitivity analyses downstream."
  ),
  "",
  "Normalization:",
  "TMM + logCPM with prior.count=2, independently within each frozen L2.",
  "",
  "Hallmark scoring:",
  "GSVA gsvaParam; Gaussian kernel; minSize=10; maxSize=500.",
  "",
  "Cytotoxicity scoring:",
  paste(
    "Mean gene-wise standardized logCPM for",
    paste(
      CYTOTOXICITY_GENES,
      collapse=", "
    )
  ),
  "",
  paste(
    "No immune identities, pseudobulk counts or frozen",
    "patient x L2 assignments were re-estimated."
  ),
  "",
  "No AA-vs-HC disease testing was performed in this module.",
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
    "22G4C2B_audit.txt"
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
# 13. FINAL STATUS
# ============================================================

if(!all(checks$pass)) {
  stop(
    "22G4C2B formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G4C2B COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

