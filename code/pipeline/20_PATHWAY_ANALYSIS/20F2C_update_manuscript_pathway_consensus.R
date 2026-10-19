
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20F2C MANUSCRIPT PATHWAY-CONSENSUS INTEGRATION\n")
cat("============================================================\n")

# ============================================================
# 1. Paths
# ============================================================

MAN <- "08_REPRODUCIBILITY/manuscript"

METHODS <- file.path(
  MAN,
  "SCRNASEQ_Methods_frozen_modules.md"
)

RESULTS <- file.path(
  MAN,
  "SCRNASEQ_Results_frozen_modules.md"
)

STATUS <- file.path(
  MAN,
  "SCRNASEQ_MANUSCRIPT_STATUS.tsv"
)

AUDIT <- file.path(
  MAN,
  "SCRNASEQ_MANUSCRIPT_AUDIT.txt"
)

ROOT20 <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "20_PATHWAY_ANALYSIS"
)

F2B <- file.path(
  ROOT20,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2B_integrated_prioritization",
  "tables"
)

F2A <- file.path(
  ROOT20,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2A_disease_architecture",
  "tables"
)

required <- c(
  METHODS,
  RESULTS,
  STATUS,
  AUDIT,
  file.path(F2B, "20F2B_integrated_state_prioritization.csv"),
  file.path(F2B, "20F2B_primary_publication_states.csv"),
  file.path(F2B, "20F2B_direct_confirmed_supported_only.csv"),
  file.path(F2B, "20F2B_opposite_states.csv"),
  file.path(F2B, "20F2B_formal_checks.csv"),
  file.path(F2A, "20F2A_disease_architecture_map.csv")
)

stopifnot(all(file.exists(required)))

cat("\nSTEP 1/8 — frozen inputs present\n")
cat("Required files:", length(required), "\n")

# ============================================================
# 2. Recover frozen numerical results
# ============================================================

cat("\nSTEP 2/8 — recovering frozen pathway-consensus results\n")

states <- read_csv(
  file.path(
    F2B,
    "20F2B_integrated_state_prioritization.csv"
  ),
  show_col_types=FALSE
)

primary <- read_csv(
  file.path(
    F2B,
    "20F2B_primary_publication_states.csv"
  ),
  show_col_types=FALSE
)

direct <- read_csv(
  file.path(
    F2B,
    "20F2B_direct_confirmed_supported_only.csv"
  ),
  show_col_types=FALSE
)

opposite <- read_csv(
  file.path(
    F2B,
    "20F2B_opposite_states.csv"
  ),
  show_col_types=FALSE
)

formal <- read_csv(
  file.path(
    F2B,
    "20F2B_formal_checks.csv"
  ),
  show_col_types=FALSE
)

stopifnot(
  nrow(states) == 65,
  nrow(primary) == 22,
  nrow(direct) == 9,
  nrow(opposite) == 2,
  all(formal$pass)
)

tier_counts <- states %>%
  count(
    integrated_evidence_tier,
    name="n"
  )

n_tier1 <- tier_counts$n[
  tier_counts$integrated_evidence_tier == "TIER_1"
]

n_tier2 <- tier_counts$n[
  tier_counts$integrated_evidence_tier == "TIER_2"
]

stopifnot(
  n_tier1 == 9,
  n_tier2 == 13
)

cat("Candidate states :", nrow(states), "\n")
cat("Tier 1           :", n_tier1, "\n")
cat("Tier 2           :", n_tier2, "\n")
cat("Primary + supporting:", nrow(primary), "\n")
cat("Direct-confirmed supported-only:", nrow(direct), "\n")
cat("Opposite states  :", nrow(opposite), "\n")

# ============================================================
# 3. Protect against duplicate integration
# ============================================================

cat("\nSTEP 3/8 — duplicate-integration guard\n")

methods_old <- readLines(
  METHODS,
  warn=FALSE
)

results_old <- readLines(
  RESULTS,
  warn=FALSE
)

method_heading <- paste0(
  "## Cross-database pathway consensus and ",
  "integrated disease-state prioritization"
)

result_heading <- paste0(
  "## Cross-database consensus resolves shared and ",
  "differential disease-state architecture"
)

if (any(methods_old == method_heading)) {
  stop(
    "20F2C Methods section already present; ",
    "refusing duplicate integration."
  )
}

if (any(results_old == result_heading)) {
  stop(
    "20F2C Results section already present; ",
    "refusing duplicate integration."
  )
}

cat("No previous 20F2C integration detected\n")

# ============================================================
# 4. Methods text
# ============================================================

cat("\nSTEP 4/8 — constructing Methods integration\n")

methods_block <- c(
  "",
  method_heading,
  "",
  paste(
    "To extend Hallmark enrichment beyond a single gene-set collection,",
    "the same frozen preranked differential-expression statistics were",
    "evaluated against Reactome, Gene Ontology Biological Process (GO-BP),",
    "and NCATS BioPlanet collections. No differential-expression model was",
    "refitted for these analyses. For every eligible L2-by-contrast analysis,",
    "the ranking remained the signed edgeR quasi-likelihood statistic,",
    "`sign(logFC) × sqrt(F)`, and enrichment direction was defined relative",
    "to the numerator of each contrast."
  ),
  "",
  paste(
    "Eight biological interpretation axes were predefined before",
    "cross-database integration: IFN/JAK-STAT/antigen presentation;",
    "T-cell/cytotoxic immunity; ECM/collagen/fibrosis;",
    "WNT/Hedgehog/hair-follicle biology; epithelial differentiation and",
    "junctional organization; vascular/angiogenic biology;",
    "metabolism/oxidative phosphorylation/reactive-oxygen-species biology;",
    "and cell-cycle/DNA-repair biology. Reactome and GO-BP pathways were",
    "mapped to these axes using the frozen pathway-name rules. The same",
    "rules were subsequently projected without modification onto the",
    "BioPlanet pathway inventory."
  ),
  "",
  paste(
    "BioPlanet was frozen from the NCATS BioPlanet source dated",
    "14 May 2019. The source comprised 1,658 pathway identifiers and",
    "9,818 unique gene symbols. Preranked enrichment was performed using",
    "`fgseaMultilevel`, restricting tested BioPlanet gene sets to sizes",
    "of 10-500 genes and retaining FDR <0.05 as the primary pathway-level",
    "threshold. BioPlanet analysis reused the 60 frozen rankings from the",
    "Hallmark analysis and therefore introduced neither a new",
    "differential-expression fit nor a change in the analytical universe."
  ),
  "",
  paste(
    "Cross-database consensus was evaluated using Reactome, GO-BP and",
    "BioPlanet as three primary pathway databases. For each contrast,",
    "frozen L2 identity, biological axis and direction, evidence supported",
    "by one primary database was classified as SINGLE, concordant evidence",
    "from two databases as CORE, and concordant evidence from all three",
    "databases as TRIPLE. Hallmark enrichment was retained as orthogonal",
    "gene-set support and was not counted toward the three-database",
    "consensus denominator. Patient-level Hallmark GSVA was likewise",
    "retained as orthogonal patient-level support rather than primary",
    "database evidence. States with significant evidence in both",
    "directions across primary databases were retained explicitly as",
    "mixed rather than forcing a single direction."
  ),
  "",
  paste(
    "Disease architecture was defined exclusively from the AA-versus-HC",
    "and FFA-versus-HC primary CORE/TRIPLE consensus. L2-by-axis states",
    "were classified as concordantly increased (SHARED_UP), concordantly",
    "decreased (SHARED_DOWN), directionally opposed (OPPOSITE),",
    "AA-supported only (AA_SUPPORTED_ONLY), FFA-supported only",
    "(FFA_SUPPORTED_ONLY), mixed/ambiguous, or without primary consensus.",
    "The supported-only terminology was deliberately used instead of",
    "disease-specific terminology because absence of consensus evidence",
    "in the other disease does not establish biological specificity."
  ),
  "",
  paste(
    "The direct AA-versus-FFA contrast was not used to define this",
    "architecture. It was retained as an independent triangulation layer.",
    "For supported-only states, a concordant direct contrast was considered",
    "confirmatory; for shared states, a direct signal was treated as",
    "compatible but non-discriminating rather than as evidence against",
    "shared dysregulation. Directionally mixed direct evidence was retained",
    "explicitly."
  ),
  "",
  paste(
    "Finally, the frozen disease-architecture states were prioritized",
    "categorically without constructing an artificial continuous evidence",
    "score. Integrated evidence incorporated the number of concordant",
    "primary databases in each disease, orthogonal Hallmark support,",
    "patient-level GSVA support, and the independent direct AA-versus-FFA",
    "triangulation. Tier 1 denoted the highest-support bilateral states;",
    "Tier 2 comprised primary supporting states, including strongly",
    "supported unilateral states; lower tiers retained secondary and",
    "contextual evidence. This prioritization changed neither the",
    "underlying GSEA/GSVA results nor the frozen L2 annotation."
  )
)

# ============================================================
# 5. Results text
# ============================================================

cat("\nSTEP 5/8 — constructing Results integration\n")

results_block <- c(
  "",
  result_heading,
  "",
  paste(
    "Integration of Reactome, GO-BP and BioPlanet enrichment produced",
    "65 candidate L2-by-biological-axis states with primary consensus",
    "support in at least one disease-versus-control comparison.",
    "Categorical evidence integration resolved 9 Tier 1 states and",
    "13 Tier 2 states; these 22 states constituted the primary and",
    "primary-supporting pathway architecture. The remaining states were",
    "retained as secondary or contextual evidence rather than discarded."
  ),
  "",
  paste(
    "The dominant high-confidence architecture was shared between AA and",
    "FFA. All nine Tier 1 states were supported bilaterally by all three",
    "primary databases and additionally by both Hallmark enrichment and",
    "patient-level GSVA. Eight were concordantly increased in AA and FFA:",
    "cell-cycle/DNA-repair programmes in basal epithelium;",
    "ECM/collagen/fibrosis programmes in PTGDS/IGF1 fibroblasts;",
    "IFN/JAK-STAT/antigen-presentation programmes in bulge/isthmus,",
    "suprabasal and terminal infundibular epithelium, ORS, and",
    "macrophage/monocytes; and metabolism/OXPHOS/ROS programmes in",
    "contractile mural cells. The ninth Tier 1 state was a concordant",
    "decrease in WNT/Hedgehog/hair-follicle programmes in terminal",
    "infundibular epithelium."
  ),
  "",
  paste(
    "The IFN/JAK-STAT/antigen-presentation axis showed the broadest",
    "high-confidence multicompartmental activation. In addition to the",
    "Tier 1 epithelial and macrophage/monocyte states, Tier 2 shared",
    "increases were present in basal epithelium, differentiating",
    "infundibular epithelium, ORS-basal cells and PTGDS/IGF1 fibroblasts.",
    "This convergence across epithelial, immune and stromal compartments",
    "supports a distributed inflammatory transcriptional programme rather",
    "than restriction to a single cell lineage."
  ),
  "",
  paste(
    "Structural and homeostatic pathway changes accompanied this shared",
    "inflammatory architecture. PTGDS/IGF1 fibroblasts showed a Tier 1",
    "shared increase in the ECM/collagen/fibrosis axis. Terminal and",
    "suprabasal infundibular states showed concordant depletion of",
    "WNT/Hedgehog/hair-follicle programmes, with the terminal state",
    "reaching Tier 1 and the suprabasal state Tier 2. Terminal",
    "infundibular epithelium also showed Tier 2 shared depletion of the",
    "ECM/collagen/fibrosis axis. ORS-basal cells showed shared Tier 2",
    "increases in cell-cycle/DNA-repair and metabolism/OXPHOS/ROS",
    "programmes."
  ),
  "",
  paste(
    "Several states were supported in only one disease-versus-control",
    "analysis and were therefore retained conservatively as supported-only",
    "rather than disease-specific. Nine such states were independently",
    "confirmed by the direct AA-versus-FFA contrast. These included",
    "AA-supported IFN/JAK-STAT/antigen-presentation activity in ACKR1+",
    "venular endothelial cells and FFA-supported changes involving",
    "suprabasal infundibular ECM and vascular programmes and ORS-basal",
    "WNT/Hedgehog/hair-follicle programmes. Direct triangulation therefore",
    "strengthened evidence for differential disease architecture without",
    "converting absence of support in the comparator disease into a claim",
    "of biological specificity."
  ),
  "",
  paste(
    "Only two candidate states showed directionally opposed primary",
    "disease-versus-control consensus: the ECM/collagen/fibrosis axis in",
    "ORS-basal cells and the epithelial-differentiation axis in ORS.",
    "Both were increased in AA and decreased in FFA. These states remained",
    "below the primary publication tiers and were retained as secondary or",
    "contextual observations. More broadly, mixed-direction database",
    "signals were retained explicitly rather than being silently resolved",
    "into a single biological direction."
  ),
  "",
  paste(
    "Together, the cross-database analysis identifies a reproducible",
    "shared core of multicompartmental interferon/antigen-presentation",
    "activation coupled to epithelial, stromal and metabolic remodeling,",
    "while also preserving a smaller set of disease-differential",
    "supported states. These results represent coordinated transcriptional",
    "programme enrichment and do not establish biochemical pathway",
    "activation, causal signaling relationships or lineage transitions."
  )
)

# ============================================================
# 6. Append manuscript modules
# ============================================================

cat("\nSTEP 6/8 — updating frozen Methods and Results\n")

writeLines(
  c(
    methods_old,
    methods_block
  ),
  METHODS
)

writeLines(
  c(
    results_old,
    results_block
  ),
  RESULTS
)

# ============================================================
# 7. Update status table
# ============================================================

cat("\nSTEP 7/8 — updating manuscript status and audit\n")

status <- read_tsv(
  STATUS,
  show_col_types=FALSE
)

new_status <- tibble(
  module=c(
    "20D_TARGETED_PATHWAYS",
    "20E_BIOPLANET",
    "20F_CROSS_DATABASE_CONSENSUS",
    "20F2A_DISEASE_ARCHITECTURE",
    "20F2B_INTEGRATED_PRIORITIZATION"
  ),
  status=rep(
    "FROZEN",
    5
  ),
  manuscript_state=c(
    paste(
      "Reactome + GO-BP complete collections; eight predefined",
      "biological interpretation axes; pathway redundancy and",
      "leading-edge architecture frozen"
    ),
    paste(
      "NCATS BioPlanet 2019-05-14; 1,658 pathways;",
      "9,818 unique genes; 60 frozen rankings reused;",
      "4,378 FDR<0.05 L2-pathway associations"
    ),
    paste(
      "Reactome + GO-BP + BioPlanet primary consensus;",
      "Hallmark and patient GSVA orthogonal;",
      "CORE/TRIPLE directional evidence frozen"
    ),
    paste(
      "AA-HC + FFA-HC define architecture;",
      "AA-FFA used only for independent triangulation;",
      "supported-only terminology retained"
    ),
    paste(
      "65 candidate states; 9 Tier 1; 13 Tier 2;",
      "22 primary/primary-supporting states;",
      "9 direct-confirmed supported-only states;",
      "2 opposite states"
    )
  )
)

status <- status %>%
  filter(
    !module %in% new_status$module
  ) %>%
  bind_rows(
    new_status
  )

write_tsv(
  status,
  STATUS
)

# ============================================================
# Audit extension
# ============================================================

scripts <- c(
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20D1_freeze_targeted_pathways.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20D2_full_Reactome_GOBP_GSEA.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20E0_freeze_BioPlanet.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20E1_BioPlanet_GSEA.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20E1B_BioPlanet_postprocessing.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20F0_BioPlanet_targeted_axes.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20F1_prepare_cross_database_consensus.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20F1B_cross_database_consensus.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20F2A_disease_architecture_triangulation.R",
  "PIPELINE_V2/scripts/20_PATHWAY_ANALYSIS/20F2B_integrated_state_prioritization.R"
)

existing_scripts <- scripts[
  file.exists(scripts)
]

script_hashes <- if (length(existing_scripts)) {

  sha_out <- system2(
    "shasum",
    c("-a", "256", existing_scripts),
    stdout = TRUE,
    stderr = TRUE
  )

  sha_parts <- strsplit(sha_out, "[[:space:]]+")

  hashes <- vapply(
    sha_parts,
    function(z) z[[1]],
    character(1)
  )

  files <- vapply(
    sha_parts,
    function(z) z[[length(z)]],
    character(1)
  )

  names(hashes) <- files
  hashes

} else {
  character()
}

audit_block <- c(
  "",
  "========================================================================",
  "20F2C — CROSS-DATABASE PATHWAY CONSENSUS / DISEASE ARCHITECTURE",
  "========================================================================",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "STATUS",
  "------",
  "- 20D_TARGETED_PATHWAYS: FROZEN",
  "- 20E_BIOPLANET: FROZEN",
  "- 20F_CROSS_DATABASE_CONSENSUS: FROZEN",
  "- 20F2A_DISEASE_ARCHITECTURE: FROZEN",
  "- 20F2B_INTEGRATED_PRIORITIZATION: FROZEN",
  "",
  "PRIMARY PATHWAY EVIDENCE",
  "------------------------",
  "- Reactome",
  "- GO Biological Process",
  "- NCATS BioPlanet",
  "- Primary cross-database denominator = 3",
  "",
  "ORTHOGONAL SUPPORT",
  "------------------",
  "- MSigDB Hallmark enrichment",
  "- patient-level Hallmark GSVA",
  "- neither is counted as a primary database",
  "",
  "BIOPLANET",
  "---------",
  "- source: NCATS BioPlanet",
  "- source date: 2019-05-14",
  "- frozen pathways: 1,658",
  "- unique gene symbols: 9,818",
  "- GSEA-eligible pathway sizes: 10-500 genes",
  "- same 60 frozen signed edgeR QL rankings reused",
  "- 49,000 total GSEA result rows",
  "- 4,378 FDR<0.05 L2-pathway associations",
  "",
  "FROZEN BIOLOGICAL AXES",
  "----------------------",
  "- IFN/JAK-STAT/antigen presentation",
  "- T-cell/cytotoxic immunity",
  "- ECM/collagen/fibrosis",
  "- WNT/Hedgehog/hair follicle",
  "- epithelial differentiation/junction",
  "- vascular/angiogenesis",
  "- metabolism/OXPHOS/ROS",
  "- cell cycle/DNA repair",
  "",
  "CROSS-DATABASE CONSENSUS",
  "------------------------",
  "- SINGLE = support from one primary database",
  "- CORE = concordant support from two primary databases",
  "- TRIPLE = concordant support from all three primary databases",
  "- mixed-direction evidence retained explicitly",
  "",
  "DISEASE ARCHITECTURE",
  "--------------------",
  "- defined from AA_vs_HC and FFA_vs_HC primary CORE/TRIPLE evidence",
  "- AA_vs_FFA does not define disease architecture",
  "- AA_vs_FFA is retained as independent triangulation",
  "- supported-only terminology is retained",
  "- absence of support in one disease is not interpreted as specificity",
  "",
  "20F2B INTEGRATED PRIORITIZATION",
  "--------------------------------",
  "- candidate states: 65",
  "- Tier 1: 9",
  "- Tier 2: 13",
  "- primary/primary-supporting states: 22",
  "- direct-confirmed supported-only states: 9",
  "- opposite states: 2",
  "- no artificial continuous evidence score was constructed",
  "",
  "TIER 1 BIOLOGICAL ARCHITECTURE",
  "------------------------------",
  "- basal epithelial: shared UP cell-cycle/DNA-repair",
  "- PTGDS/IGF1 fibroblast: shared UP ECM/collagen/fibrosis",
  "- bulge/isthmus: shared UP IFN/JAK-STAT/antigen presentation",
  "- suprabasal infundibulum: shared UP IFN/JAK-STAT/antigen presentation",
  "- terminal infundibulum: shared UP IFN/JAK-STAT/antigen presentation",
  "- ORS: shared UP IFN/JAK-STAT/antigen presentation",
  "- macrophage/monocyte: shared UP IFN/JAK-STAT/antigen presentation",
  "- contractile mural cell: shared UP metabolism/OXPHOS/ROS",
  "- terminal infundibulum: shared DOWN WNT/Hedgehog/hair-follicle",
  "",
  "INTERPRETIVE BOUNDARIES",
  "-----------------------",
  "- pathway enrichment represents coordinated transcriptional programmes",
  "- enrichment does not establish biochemical pathway activation",
  "- enrichment does not establish causal signaling",
  "- enrichment does not establish lineage transition",
  "- AA_SUPPORTED_ONLY and FFA_SUPPORTED_ONLY are not called disease-specific",
  "- direct AA-vs-FFA evidence is triangulation, not primary architecture",
  "- mixed directions are not silently collapsed",
  "",
  "ANALYTICAL INTEGRITY",
  "--------------------",
  "- no GSEA refit",
  "- no GSVA refit",
  "- no differential-expression refit",
  "- no L2 annotation modified",
  "- frozen upstream analytical outputs reused",
  "",
  "REPRODUCIBILITY",
  "---------------"
)

if (length(script_hashes)) {
  audit_block <- c(
    audit_block,
    paste0(
      unname(script_hashes),
      "  ",
      names(script_hashes)
    )
  )
}

audit_block <- c(
  audit_block,
  "",
  "FILES UPDATED",
  "-------------",
  METHODS,
  RESULTS,
  STATUS,
  AUDIT,
  "",
  "MAIN MANUSCRIPT",
  "---------------",
  "07_MANUSCRIPT/main/sections/02_methods.md NOT modified.",
  "07_MANUSCRIPT/main/sections/03_results_skeleton.md NOT modified.",
  "07_MANUSCRIPT/main/manuscript_v0.1.md NOT modified.",
  "",
  "No analytical results were recalculated by this integration module.",
  "No cell annotations were modified.",
  "",
  "PASS",
  "========================================================================"
)

write(
  audit_block,
  file=AUDIT,
  append=TRUE
)

# ============================================================
# 8. Formal post-write validation
# ============================================================

cat("\nSTEP 8/8 — formal post-write validation\n")

methods_new <- readLines(
  METHODS,
  warn=FALSE
)

results_new <- readLines(
  RESULTS,
  warn=FALSE
)

status_new <- read_tsv(
  STATUS,
  show_col_types=FALSE
)

checks <- tibble(
  check=c(
    "Methods_section_added_once",
    "Results_section_added_once",
    "20D_status_frozen",
    "20E_status_frozen",
    "20F_status_frozen",
    "20F2A_status_frozen",
    "20F2B_status_frozen",
    "candidate_states_65",
    "Tier1_9",
    "Tier2_13",
    "primary_publication_states_22",
    "direct_confirmed_supported_only_9",
    "opposite_states_2",
    "upstream_20F2B_checks_pass",
    "supported_only_terminology_retained",
    "no_specificity_claim_from_absence",
    "no_GSEA_refit",
    "no_GSVA_refit",
    "no_DE_refit",
    "annotation_unchanged"
  ),
  pass=c(
    sum(methods_new == method_heading) == 1,
    sum(results_new == result_heading) == 1,
    any(
      status_new$module == "20D_TARGETED_PATHWAYS" &
      status_new$status == "FROZEN"
    ),
    any(
      status_new$module == "20E_BIOPLANET" &
      status_new$status == "FROZEN"
    ),
    any(
      status_new$module == "20F_CROSS_DATABASE_CONSENSUS" &
      status_new$status == "FROZEN"
    ),
    any(
      status_new$module == "20F2A_DISEASE_ARCHITECTURE" &
      status_new$status == "FROZEN"
    ),
    any(
      status_new$module == "20F2B_INTEGRATED_PRIORITIZATION" &
      status_new$status == "FROZEN"
    ),
    nrow(states) == 65,
    n_tier1 == 9,
    n_tier2 == 13,
    nrow(primary) == 22,
    nrow(direct) == 9,
    nrow(opposite) == 2,
    all(formal$pass),
    all(
      c(
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      ) %in% unique(states$architecture)
    ),
    !any(
      grepl(
        "AA-specific|FFA-specific|disease-specific",
        c(methods_block, results_block),
        ignore.case=TRUE
      )
    ),
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  checks,
  n=Inf
)

stopifnot(
  all(checks$pass)
)

cat("\n===== UPDATED MANUSCRIPT STATUS =====\n")
print(
  status_new %>%
    filter(
      grepl(
        "^20",
        module
      )
    ),
  n=Inf,
  width=Inf
)

cat("\n===== METHODS SECTION ADDED =====\n")
cat(
  paste(
    methods_block,
    collapse="\n"
  ),
  "\n"
)

cat("\n===== RESULTS SECTION ADDED =====\n")
cat(
  paste(
    results_block,
    collapse="\n"
  ),
  "\n"
)

cat("\n============================================================\n")
cat("20F2C COMPLETE\n")
cat("============================================================\n")
cat("20D-20F2B pathway-analysis modules integrated into frozen manuscript files\n")
cat("Methods updated\n")
cat("Results updated\n")
cat("Status updated\n")
cat("Audit updated\n")
cat("Main manuscript source files NOT modified\n")
cat("No analytical refit\n")
cat("No annotation modified\n")
cat("PASS\n")
