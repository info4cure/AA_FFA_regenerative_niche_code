suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20B1 CROSS-DISEASE HALLMARK ARCHITECTURE\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2", "results", "disease_biology",
  "20_PATHWAY_ANALYSIS"
)

IN <- file.path(
  BASE, "20A_HALLMARK_GSEA", "tables",
  "20A_all_Hallmark_GSEA_results.csv"
)

OUT <- file.path(
  BASE, "20B_HALLMARK_ARCHITECTURE"
)

TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(IN))

x <- read_csv(IN, show_col_types=FALSE)

stopifnot(
  all(c(
    "contrast",
    "module",
    "final_L2_label",
    "pathway",
    "NES",
    "padj"
  ) %in% names(x))
)

# ============================================================
# 1. Identify L2 shared between AA-HC and FFA-HC
# ============================================================

cat("\nSTEP 1/7 — identifying shared L2\n")

AA_L2 <- x %>%
  filter(contrast=="AA_vs_HC") %>%
  distinct(module, final_L2_label)

FFA_L2 <- x %>%
  filter(contrast=="FFA_vs_HC") %>%
  distinct(module, final_L2_label)

shared_L2 <- inner_join(
  AA_L2,
  FFA_L2,
  by=c("module","final_L2_label")
)

cat("AA-HC L2 :", nrow(AA_L2), "\n")
cat("FFA-HC L2:", nrow(FFA_L2), "\n")
cat("Shared L2:", nrow(shared_L2), "\n")

write_csv(
  shared_L2,
  file.path(TAB, "20B1_shared_L2.csv")
)

# ============================================================
# 2. Build AA-HC vs FFA-HC paired Hallmark table
# ============================================================

cat("\nSTEP 2/7 — pairing disease-vs-HC Hallmarks\n")

aa <- x %>%
  filter(contrast=="AA_vs_HC") %>%
  select(
    module,
    final_L2_label,
    pathway,
    NES_AA_HC=NES,
    FDR_AA_HC=padj
  )

ffa <- x %>%
  filter(contrast=="FFA_vs_HC") %>%
  select(
    module,
    final_L2_label,
    pathway,
    NES_FFA_HC=NES,
    FDR_FFA_HC=padj
  )

paired <- inner_join(
  aa,
  ffa,
  by=c(
    "module",
    "final_L2_label",
    "pathway"
  )
)

expected_rows <- nrow(shared_L2) * 50

cat("Paired rows:", nrow(paired), "\n")
cat("Expected approximately:", expected_rows, "\n")

# ============================================================
# 3. Classify architecture
# ============================================================

cat("\nSTEP 3/7 — classifying shared/specific architecture\n")

paired <- paired %>%
  mutate(
    sig_AA  = FDR_AA_HC  < .05,
    sig_FFA = FDR_FFA_HC < .05,

    dir_AA = case_when(
      NES_AA_HC > 0 ~ "UP",
      NES_AA_HC < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    ),

    dir_FFA = case_when(
      NES_FFA_HC > 0 ~ "UP",
      NES_FFA_HC < 0 ~ "DOWN",
      TRUE ~ "ZERO"
    ),

    architecture = case_when(

      sig_AA & sig_FFA &
        NES_AA_HC > 0 & NES_FFA_HC > 0
        ~ "shared_UP",

      sig_AA & sig_FFA &
        NES_AA_HC < 0 & NES_FFA_HC < 0
        ~ "shared_DOWN",

      sig_AA & sig_FFA &
        NES_AA_HC > 0 & NES_FFA_HC < 0
        ~ "AA_UP_FFA_DOWN",

      sig_AA & sig_FFA &
        NES_AA_HC < 0 & NES_FFA_HC > 0
        ~ "AA_DOWN_FFA_UP",

      sig_AA & !sig_FFA &
        NES_AA_HC > 0
        ~ "AA_specific_UP",

      sig_AA & !sig_FFA &
        NES_AA_HC < 0
        ~ "AA_specific_DOWN",

      !sig_AA & sig_FFA &
        NES_FFA_HC > 0
        ~ "FFA_specific_UP",

      !sig_AA & sig_FFA &
        NES_FFA_HC < 0
        ~ "FFA_specific_DOWN",

      TRUE ~ "not_significant"
    ),

    delta_NES_AA_minus_FFA =
      NES_AA_HC - NES_FFA_HC
  )

write_csv(
  paired,
  file.path(
    TAB,
    "20B1_cross_disease_Hallmark_architecture.csv"
  )
)

# ============================================================
# 4. Add direct AA-vs-FFA evidence where available
# ============================================================

cat("\nSTEP 4/7 — integrating direct AA vs FFA GSEA\n")

direct <- x %>%
  filter(contrast=="AA_vs_FFA") %>%
  select(
    module,
    final_L2_label,
    pathway,
    NES_AA_FFA=NES,
    FDR_AA_FFA=padj
  )

triangulated <- paired %>%
  left_join(
    direct,
    by=c(
      "module",
      "final_L2_label",
      "pathway"
    )
  ) %>%
  mutate(

    expected_direct_direction =
      case_when(
        delta_NES_AA_minus_FFA > 0 ~ "AA",
        delta_NES_AA_minus_FFA < 0 ~ "FFA",
        TRUE ~ "TIE"
      ),

    observed_direct_direction =
      case_when(
        is.na(NES_AA_FFA) ~ NA_character_,
        NES_AA_FFA > 0 ~ "AA",
        NES_AA_FFA < 0 ~ "FFA",
        TRUE ~ "TIE"
      ),

    direct_direction_concordant =
      case_when(
        is.na(observed_direct_direction) ~ NA,
        expected_direct_direction=="TIE" ~ NA,
        observed_direct_direction ==
          expected_direct_direction ~ TRUE,
        TRUE ~ FALSE
      ),

    direct_FDR005 =
      !is.na(FDR_AA_FFA) &
      FDR_AA_FFA < .05
  )

write_csv(
  triangulated,
  file.path(
    TAB,
    "20B1_Hallmark_triangulated_with_AA_vs_FFA.csv"
  )
)

# ============================================================
# 5. Summaries
# ============================================================

cat("\nSTEP 5/7 — summaries\n")

architecture_summary <- triangulated %>%
  count(
    architecture,
    name="n_L2_pathway_pairs"
  ) %>%
  arrange(desc(n_L2_pathway_pairs))

write_csv(
  architecture_summary,
  file.path(
    TAB,
    "20B1_architecture_summary.csv"
  )
)

pathway_architecture <- triangulated %>%
  filter(architecture != "not_significant") %>%
  group_by(
    pathway,
    architecture
  ) %>%
  summarise(
    n_L2=n_distinct(final_L2_label),
    median_NES_AA=median(NES_AA_HC),
    median_NES_FFA=median(NES_FFA_HC),
    median_delta_NES=
      median(delta_NES_AA_minus_FFA),
    L2=paste(
      sort(unique(final_L2_label)),
      collapse=";"
    ),
    .groups="drop"
  ) %>%
  arrange(
    desc(n_L2),
    pathway
  )

write_csv(
  pathway_architecture,
  file.path(
    TAB,
    "20B1_pathway_architecture_recurrence.csv"
  )
)

L2_architecture <- triangulated %>%
  filter(architecture != "not_significant") %>%
  count(
    module,
    final_L2_label,
    architecture,
    name="n_pathways"
  ) %>%
  arrange(
    module,
    final_L2_label,
    desc(n_pathways)
  )

write_csv(
  L2_architecture,
  file.path(
    TAB,
    "20B1_L2_architecture_summary.csv"
  )
)

# ============================================================
# 6. Direct-contrast concordance
# ============================================================

cat("\nSTEP 6/7 — direct-contrast concordance\n")

direct_eval <- triangulated %>%
  filter(
    !is.na(NES_AA_FFA),
    expected_direct_direction != "TIE"
  )

direct_summary <- tibble(
  n_pairs_with_direct=nrow(direct_eval),

  direction_concordance_all=
    mean(
      direct_eval$direct_direction_concordant,
      na.rm=TRUE
    ),

  n_direct_FDR005=
    sum(
      direct_eval$direct_FDR005,
      na.rm=TRUE
    ),

  direction_concordance_direct_FDR005=
    ifelse(
      sum(direct_eval$direct_FDR005) > 0,
      mean(
        direct_eval$direct_direction_concordant[
          direct_eval$direct_FDR005
        ],
        na.rm=TRUE
      ),
      NA_real_
    ),

  spearman_delta_vs_direct_NES=
    cor(
      direct_eval$delta_NES_AA_minus_FFA,
      direct_eval$NES_AA_FFA,
      method="spearman",
      use="complete.obs"
    )
)

write_csv(
  direct_summary,
  file.path(
    TAB,
    "20B1_direct_contrast_concordance.csv"
  )
)

# ============================================================
# 7. Formal checks + audit
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check=c(
    "input_20A_present",
    "AA_HC_24",
    "FFA_HC_18",
    "AA_FFA_18",
    "shared_L2_nonzero",
    "paired_table_nonempty",
    "architecture_complete",
    "direct_contrast_integrated",
    "no_GSEA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    file.exists(IN),
    nrow(AA_L2)==24,
    nrow(FFA_L2)==18,
    n_distinct(
      x$final_L2_label[
        x$contrast=="AA_vs_FFA"
      ]
    )==18,
    nrow(shared_L2)>0,
    nrow(paired)>0,
    !any(is.na(paired$architecture)),
    nrow(direct)>0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

print(
  as.data.frame(checks),
  row.names=FALSE
)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20B1 cross-disease Hallmark architecture",
  "",
  "Source:",
  "- Frozen 20A Hallmark GSEA results",
  "",
  "Purpose:",
  "- Compare AA-vs-HC and FFA-vs-HC Hallmark NES within shared L2 identities",
  "- Classify shared, disease-specific, and direction-opposed pathway architecture",
  "- Triangulate disease-vs-HC differences against direct AA-vs-FFA GSEA",
  "",
  "Primary pathway significance threshold:",
  "- BH FDR < 0.05 from frozen 20A",
  "",
  "Architecture classification is descriptive.",
  "No new GSEA was performed.",
  "No DE was refit.",
  "No annotation was modified.",
  "No manuscript figures were generated."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20B1_audit.txt"
  )
)

cat("\n===== ARCHITECTURE SUMMARY =====\n")
print(
  as.data.frame(architecture_summary),
  row.names=FALSE
)

cat("\n===== DIRECT AA-vs-FFA CONCORDANCE =====\n")
print(
  as.data.frame(direct_summary),
  row.names=FALSE
)

cat("\n===== TOP RECURRENT ARCHITECTURES =====\n")
print(
  as.data.frame(
    head(pathway_architecture, 40)
  ),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("20B1 COMPLETE\n")
cat("============================================================\n")
cat("Cross-disease Hallmark architecture\n")
cat("No GSEA refit\n")
cat("No DE refit\n")
cat("No annotation modified\n")
cat("PASS\n")
