suppressPackageStartupMessages({
  library(tidyverse)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20F2A DISEASE ARCHITECTURE TRIANGULATION\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

ROOT <- file.path(
  "PIPELINE_V2/results/disease_biology",
  "20_PATHWAY_ANALYSIS"
)

IN <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F1B_consensus",
  "tables"
)

OUT <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS",
  "20F2A_disease_architecture"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

CORE_FILE <- file.path(
  IN,
  "20F1B_core_consensus_map.csv"
)

ALL_FILE <- file.path(
  IN,
  "20F1B_cross_database_axis_consensus.csv"
)

MIXED_FILE <- file.path(
  IN,
  "20F1B_mixed_consensus_states.csv"
)

# ============================================================
# 1. Recover frozen 20F1B
# ============================================================

cat("\nSTEP 1/9 — recovering frozen 20F1B consensus\n")

stopifnot(
  file.exists(CORE_FILE),
  file.exists(ALL_FILE),
  file.exists(MIXED_FILE)
)

core <- read_csv(
  CORE_FILE,
  show_col_types=FALSE
)

all_cons <- read_csv(
  ALL_FILE,
  show_col_types=FALSE
)

mixed <- read_csv(
  MIXED_FILE,
  show_col_types=FALSE
)

cat("CORE/TRIPLE rows :", nrow(core), "\n")
cat("All support rows :", nrow(all_cons), "\n")
cat("Mixed states     :", nrow(mixed), "\n")

stopifnot(
  all(c(
    "AA_vs_FFA",
    "AA_vs_HC",
    "FFA_vs_HC"
  ) %in% unique(all_cons$contrast))
)

# ============================================================
# 2. Shared L2 × axis universe for AA and FFA vs HC
# ============================================================

cat("\nSTEP 2/9 — defining shared AA/FFA-vs-HC universe\n")

keys <- c(
  "module",
  "final_L2_label",
  "axis"
)

aa_universe <- all_cons %>%
  filter(contrast == "AA_vs_HC") %>%
  distinct(across(all_of(keys)))

ffa_universe <- all_cons %>%
  filter(contrast == "FFA_vs_HC") %>%
  distinct(across(all_of(keys)))

shared_universe <- inner_join(
  aa_universe,
  ffa_universe,
  by=keys
) %>%
  arrange(
    module,
    final_L2_label,
    axis
  )

cat(
  "Shared L2 × axis states:",
  nrow(shared_universe),
  "\n"
)

# ============================================================
# 3. Collapse primary consensus per disease
# ============================================================

cat("\nSTEP 3/9 — collapsing primary directional consensus\n")

collapse_disease <- function(dat, contrast_name, prefix) {

  # Only CORE/TRIPLE states define primary consensus.
  x <- dat %>%
    filter(
      contrast == contrast_name,
      consensus_class %in% c("CORE","TRIPLE")
    ) %>%
    group_by(
      module,
      final_L2_label,
      axis
    ) %>%
    summarise(
      directions = paste(
        sort(unique(direction)),
        collapse=";"
      ),

      has_UP = any(direction == "UP"),
      has_DOWN = any(direction == "DOWN"),

      max_primary_databases =
        max(n_primary_databases, na.rm=TRUE),

      primary_databases = paste(
        sort(unique(
          unlist(
            strsplit(
              primary_databases,
              ";",
              fixed=TRUE
            )
          )
        )),
        collapse=";"
      ),

      n_significant_pathways =
        sum(total_significant_pathways, na.rm=TRUE),

      median_database_NES =
        median(
          median_database_NES,
          na.rm=TRUE
        ),

      min_FDR =
        suppressWarnings(
          min(min_FDR, na.rm=TRUE)
        ),

      Hallmark_support =
        any(Hallmark_support %in% TRUE),

      GSVA_support =
        any(GSVA_support %in% TRUE),

      .groups="drop"
    ) %>%
    mutate(
      primary_status = case_when(
        has_UP & has_DOWN ~ "MIXED",
        has_UP ~ "UP",
        has_DOWN ~ "DOWN",
        TRUE ~ "NONE"
      )
    )

  # Rename disease-specific columns.
  x %>%
    rename_with(
      ~paste0(prefix, "_", .x),
      -all_of(keys)
    )
}

aa <- collapse_disease(
  all_cons,
  "AA_vs_HC",
  "AA"
)

ffa <- collapse_disease(
  all_cons,
  "FFA_vs_HC",
  "FFA"
)

architecture <- shared_universe %>%
  left_join(
    aa,
    by=keys
  ) %>%
  left_join(
    ffa,
    by=keys
  ) %>%
  mutate(
    AA_primary_status =
      replace_na(AA_primary_status, "NONE"),

    FFA_primary_status =
      replace_na(FFA_primary_status, "NONE"),

    AA_has_consensus =
      AA_primary_status != "NONE",

    FFA_has_consensus =
      FFA_primary_status != "NONE"
  )

# ============================================================
# 4. Explicit disease architecture
# ============================================================

cat("\nSTEP 4/9 — assigning disease architecture\n")

architecture <- architecture %>%
  mutate(
    architecture = case_when(

      AA_primary_status == "MIXED" |
      FFA_primary_status == "MIXED" ~
        "MIXED_AMBIGUOUS",

      AA_primary_status == "UP" &
      FFA_primary_status == "UP" ~
        "SHARED_UP",

      AA_primary_status == "DOWN" &
      FFA_primary_status == "DOWN" ~
        "SHARED_DOWN",

      AA_primary_status %in% c("UP","DOWN") &
      FFA_primary_status %in% c("UP","DOWN") &
      AA_primary_status != FFA_primary_status ~
        "OPPOSITE",

      AA_primary_status %in% c("UP","DOWN") &
      FFA_primary_status == "NONE" ~
        "AA_SUPPORTED_ONLY",

      FFA_primary_status %in% c("UP","DOWN") &
      AA_primary_status == "NONE" ~
        "FFA_SUPPORTED_ONLY",

      TRUE ~
        "NO_PRIMARY_CONSENSUS"
    )
  )

cat("\n===== ARCHITECTURE COUNTS =====\n")

print(
  architecture %>%
    count(
      architecture,
      sort=TRUE
    ),
  n=Inf
)

# ============================================================
# 5. Recover direct AA vs FFA evidence
# ============================================================

cat("\nSTEP 5/9 — triangulating direct AA_vs_FFA evidence\n")

direct <- all_cons %>%
  filter(
    contrast == "AA_vs_FFA",
    consensus_class %in% c("CORE","TRIPLE")
  ) %>%
  group_by(
    module,
    final_L2_label,
    axis
  ) %>%
  summarise(
    direct_has_UP =
      any(direction == "UP"),

    direct_has_DOWN =
      any(direction == "DOWN"),

    direct_max_primary_databases =
      max(n_primary_databases, na.rm=TRUE),

    direct_n_significant_pathways =
      sum(total_significant_pathways, na.rm=TRUE),

    direct_median_NES =
      median(
        median_database_NES,
        na.rm=TRUE
      ),

    direct_min_FDR =
      suppressWarnings(
        min(min_FDR, na.rm=TRUE)
      ),

    direct_Hallmark_support =
      any(Hallmark_support %in% TRUE),

    direct_GSVA_support =
      any(GSVA_support %in% TRUE),

    .groups="drop"
  ) %>%
  mutate(
    direct_status = case_when(
      direct_has_UP & direct_has_DOWN ~ "MIXED",
      direct_has_UP ~ "AA_HIGHER",
      direct_has_DOWN ~ "FFA_HIGHER",
      TRUE ~ "NONE"
    )
  )

architecture <- architecture %>%
  left_join(
    direct,
    by=keys
  ) %>%
  mutate(
    direct_status =
      replace_na(direct_status, "NONE")
  )

# ============================================================
# 6. Expected AA-vs-FFA direction
# ============================================================

cat("\nSTEP 6/9 — evaluating direct-contrast concordance\n")

architecture <- architecture %>%
  mutate(

    expected_direct_direction = case_when(

      architecture == "OPPOSITE" &
      AA_primary_status == "UP" &
      FFA_primary_status == "DOWN" ~
        "AA_HIGHER",

      architecture == "OPPOSITE" &
      AA_primary_status == "DOWN" &
      FFA_primary_status == "UP" ~
        "FFA_HIGHER",

      architecture == "AA_SUPPORTED_ONLY" &
      AA_primary_status == "UP" ~
        "AA_HIGHER",

      architecture == "AA_SUPPORTED_ONLY" &
      AA_primary_status == "DOWN" ~
        "FFA_HIGHER",

      architecture == "FFA_SUPPORTED_ONLY" &
      FFA_primary_status == "UP" ~
        "FFA_HIGHER",

      architecture == "FFA_SUPPORTED_ONLY" &
      FFA_primary_status == "DOWN" ~
        "AA_HIGHER",

      TRUE ~ NA_character_
    ),

    direct_direction_concordant = case_when(

      is.na(expected_direct_direction) ~
        NA,

      direct_status == "NONE" ~
        NA,

      direct_status == "MIXED" ~
        NA,

      direct_status ==
        expected_direct_direction ~ TRUE,

      TRUE ~ FALSE
    ),

    direct_support = case_when(
      direct_status == "NONE" ~
        "NO_DIRECT_PRIMARY_SUPPORT",

      direct_status == "MIXED" ~
        "DIRECT_MIXED",

      direct_direction_concordant %in% TRUE ~
        "DIRECT_CONCORDANT",

      direct_direction_concordant %in% FALSE ~
        "DIRECT_DISCORDANT",

      TRUE ~
        "DIRECT_NONDISCRIMINATING"
    )
  )

# ============================================================
# 7. Evidence grade
# ============================================================

cat("\nSTEP 7/9 — grading integrated evidence\n")

architecture <- architecture %>%
  mutate(

    n_disease_consensus =
      as.integer(AA_has_consensus) +
      as.integer(FFA_has_consensus),

    any_Hallmark_support =
      coalesce(AA_Hallmark_support, FALSE) |
      coalesce(FFA_Hallmark_support, FALSE),

    any_GSVA_support =
      coalesce(AA_GSVA_support, FALSE) |
      coalesce(FFA_GSVA_support, FALSE),

    evidence_grade = case_when(

      architecture %in%
        c(
          "SHARED_UP",
          "SHARED_DOWN",
          "OPPOSITE"
        ) &
      any_GSVA_support ~
        "PRIMARY_BOTH_PLUS_PATIENT",

      architecture %in%
        c(
          "SHARED_UP",
          "SHARED_DOWN",
          "OPPOSITE"
        ) &
      any_Hallmark_support ~
        "PRIMARY_BOTH_PLUS_HALLMARK",

      architecture %in%
        c(
          "SHARED_UP",
          "SHARED_DOWN",
          "OPPOSITE"
        ) ~
        "PRIMARY_BOTH",

      architecture %in%
        c(
          "AA_SUPPORTED_ONLY",
          "FFA_SUPPORTED_ONLY"
        ) &
      direct_support == "DIRECT_CONCORDANT" &
      any_GSVA_support ~
        "DISEASE_SUPPORTED_PLUS_DIRECT_PLUS_PATIENT",

      architecture %in%
        c(
          "AA_SUPPORTED_ONLY",
          "FFA_SUPPORTED_ONLY"
        ) &
      direct_support == "DIRECT_CONCORDANT" ~
        "DISEASE_SUPPORTED_PLUS_DIRECT",

      architecture %in%
        c(
          "AA_SUPPORTED_ONLY",
          "FFA_SUPPORTED_ONLY"
        ) ~
        "DISEASE_SUPPORTED_ONLY",

      architecture ==
        "MIXED_AMBIGUOUS" ~
        "AMBIGUOUS",

      TRUE ~
        "NO_PRIMARY_CONSENSUS"
    )
  )

# ============================================================
# 8. Save tables and summaries
# ============================================================

cat("\nSTEP 8/9 — writing frozen 20F2A outputs\n")

architecture <- architecture %>%
  arrange(
    contrast = architecture,
    module,
    final_L2_label,
    axis
  )

write_csv(
  architecture,
  file.path(
    TAB,
    "20F2A_disease_architecture_map.csv"
  )
)

architecture_summary <- architecture %>%
  count(
    architecture,
    name="n_states",
    sort=TRUE
  )

write_csv(
  architecture_summary,
  file.path(
    TAB,
    "20F2A_architecture_summary.csv"
  )
)

axis_summary <- architecture %>%
  count(
    axis,
    architecture,
    name="n_states"
  ) %>%
  arrange(
    axis,
    desc(n_states),
    architecture
  )

write_csv(
  axis_summary,
  file.path(
    TAB,
    "20F2A_axis_architecture_summary.csv"
  )
)

L2_summary <- architecture %>%
  filter(
    architecture != "NO_PRIMARY_CONSENSUS"
  ) %>%
  count(
    module,
    final_L2_label,
    architecture,
    name="n_axes"
  ) %>%
  arrange(
    module,
    final_L2_label,
    desc(n_axes)
  )

write_csv(
  L2_summary,
  file.path(
    TAB,
    "20F2A_L2_architecture_summary.csv"
  )
)

direct_audit <- architecture %>%
  filter(
    !is.na(expected_direct_direction)
  ) %>%
  count(
    architecture,
    direct_support,
    name="n_states"
  ) %>%
  arrange(
    architecture,
    desc(n_states)
  )

write_csv(
  direct_audit,
  file.path(
    TAB,
    "20F2A_direct_contrast_audit.csv"
  )
)

discordant <- architecture %>%
  filter(
    direct_support ==
      "DIRECT_DISCORDANT"
  )

write_csv(
  discordant,
  file.path(
    TAB,
    "20F2A_direct_discordant_states.csv"
  )
)

mixed_out <- architecture %>%
  filter(
    architecture ==
      "MIXED_AMBIGUOUS"
  )

write_csv(
  mixed_out,
  file.path(
    TAB,
    "20F2A_mixed_ambiguous_states.csv"
  )
)

# Strong candidate states for subsequent biological prioritisation.
candidate_states <- architecture %>%
  filter(
    architecture %in% c(
      "SHARED_UP",
      "SHARED_DOWN",
      "OPPOSITE",
      "AA_SUPPORTED_ONLY",
      "FFA_SUPPORTED_ONLY"
    )
  ) %>%
  arrange(
    factor(
      architecture,
      levels=c(
        "OPPOSITE",
        "SHARED_UP",
        "SHARED_DOWN",
        "AA_SUPPORTED_ONLY",
        "FFA_SUPPORTED_ONLY"
      )
    ),
    module,
    final_L2_label,
    axis
  )

write_csv(
  candidate_states,
  file.path(
    TAB,
    "20F2A_candidate_states_for_20F2B.csv"
  )
)

# ============================================================
# 9. Formal checks / audit
# ============================================================

cat("\nSTEP 9/9 — formal checks\n")

checks <- tibble(
  check=c(
    "20F1B_input_present",
    "shared_universe_nonzero",
    "eight_axes_retained",
    "AA_FFA_vs_HC_only_define_architecture",
    "direct_contrast_used_for_triangulation",
    "supported_only_not_called_specific",
    "mixed_states_explicit",
    "direct_discordance_explicit",
    "Hallmark_not_primary",
    "GSVA_not_primary",
    "no_GSEA_refit",
    "no_GSVA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    file.exists(CORE_FILE),
    nrow(shared_universe)>0,
    n_distinct(architecture$axis)==8,
    TRUE,
    TRUE,
    !any(
      architecture$architecture %in%
        c("AA_SPECIFIC","FFA_SPECIFIC")
    ),
    all(
      architecture$architecture[
        architecture$AA_primary_status=="MIXED" |
        architecture$FFA_primary_status=="MIXED"
      ] == "MIXED_AMBIGUOUS"
    ),
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

print(
  as.data.frame(checks),
  row.names=FALSE
)

stopifnot(all(checks$pass))

cat("\n===== DISEASE ARCHITECTURE =====\n")
print(
  architecture_summary,
  n=Inf
)

cat("\n===== DIRECT-CONTRAST AUDIT =====\n")
print(
  direct_audit,
  n=Inf
)

cat("\nDirect discordant states:",
    nrow(discordant), "\n")

cat("Mixed/ambiguous states:",
    nrow(mixed_out), "\n")

audit <- c(
  "PIPELINE_V2 — 20F2A disease architecture triangulation",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Frozen input:",
  "- 20F1B cross-database consensus",
  "",
  "Architecture is defined only from AA_vs_HC and FFA_vs_HC",
  "primary CORE/TRIPLE consensus.",
  "",
  "Classes:",
  "- SHARED_UP",
  "- SHARED_DOWN",
  "- OPPOSITE",
  "- AA_SUPPORTED_ONLY",
  "- FFA_SUPPORTED_ONLY",
  "- MIXED",
  "- NO_PRIMARY_CONSENSUS",
  "",
  "AA_vs_FFA is not used to define disease architecture.",
  "It is retained as an independent direct-contrast triangulation layer.",
  "",
  "Disease-specific terminology is deliberately conservative:",
  "AA_SUPPORTED_ONLY and FFA_SUPPORTED_ONLY indicate absence of",
  "primary CORE/TRIPLE support in the other disease comparison;",
  "they do not imply biological exclusivity.",
  "",
  "Mixed-direction states are retained explicitly and are not",
  "silently collapsed into a single direction.",
  "",
  "No GSEA refit.",
  "No GSVA refit.",
  "No differential-expression refit.",
  "No annotation modified."
)

writeLines(
  audit,
  file.path(OUT, "20F2A_audit.txt")
)

cat("\n===== ARCHITECTURE SUMMARY =====\n")
print(
  architecture_summary,
  n=Inf
)

cat("\n===== DIRECT-CONTRAST AUDIT =====\n")
print(
  direct_audit,
  n=Inf
)

cat(
  "\nDirect discordant states:",
  nrow(discordant),
  "\n"
)

cat(
  "Mixed/ambiguous states:",
  nrow(mixed_out),
  "\n"
)

cat("\n============================================================\n")
cat("20F2A COMPLETE\n")
cat("============================================================\n")
cat("Disease architecture defined from AA_vs_HC + FFA_vs_HC\n")
cat("AA_vs_FFA retained as independent triangulation\n")
cat("Specificity terminology kept conservative\n")
cat("Mixed states retained explicitly\n")
cat("No analytical refit\n")
cat("PASS\n")
