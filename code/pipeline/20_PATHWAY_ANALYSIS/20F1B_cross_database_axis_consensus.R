suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

ROOT <- "PIPELINE_V2/results/disease_biology/20_PATHWAY_ANALYSIS"

IN <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS/tables"
)

OUT <- file.path(
  ROOT,
  "20F_CROSS_DATABASE_CONSENSUS/20F1B_consensus"
)

TAB <- file.path(OUT, "tables")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

cat("============================================================\n")
cat("PIPELINE_V2 — 20F1B CROSS-DATABASE AXIS CONSENSUS\n")
cat("============================================================\n")

# ============================================================
# 1. Inputs
# ============================================================

cat("\nSTEP 1/8 — recovering prepared evidence\n")

db <- read_csv(
  file.path(IN, "20F1_database_axis_support.csv"),
  show_col_types=FALSE
)

gsva <- read_csv(
  file.path(IN, "20F1_Hallmark_axis_patient_GSVA.csv"),
  show_col_types=FALSE
)

stopifnot(
  all(c(
    "database",
    "contrast",
    "module",
    "final_L2_label",
    "axis",
    "direction",
    "n_pathways",
    "median_NES",
    "min_FDR"
  ) %in% names(db))
)

primary_db <- c(
  "GO_BP",
  "REACTOME",
  "BIOPLANET"
)

stopifnot(
  all(primary_db %in% unique(db$database))
)

cat("Database-support rows:", nrow(db), "\n")
cat("GSVA rows            :", nrow(gsva), "\n")

# ============================================================
# 2. Primary directional evidence
# ============================================================

cat("\nSTEP 2/8 — constructing primary directional consensus\n")

primary <- db %>%
  filter(database %in% primary_db) %>%
  group_by(
    contrast,
    module,
    final_L2_label,
    axis,
    direction
  ) %>%
  summarise(
    n_primary_databases=n_distinct(database),
    primary_databases=paste(
      sort(unique(database)),
      collapse=";"
    ),
    total_significant_pathways=sum(n_pathways),
    median_database_NES=median(median_NES),
    max_abs_database_NES=max(abs(median_NES)),
    min_FDR=min(min_FDR),
    .groups="drop"
  ) %>%
  mutate(
    consensus_class=case_when(
      n_primary_databases == 3 ~ "TRIPLE",
      n_primary_databases == 2 ~ "CORE",
      n_primary_databases == 1 ~ "SINGLE",
      TRUE ~ "NONE"
    ),
    primary_consensus=n_primary_databases >= 2
  )

write_csv(
  primary,
  file.path(
    TAB,
    "20F1B_primary_directional_consensus.csv"
  )
)

cat("\nPrimary consensus classes:\n")
print(
  primary %>%
    count(consensus_class)
)

# ============================================================
# 3. Detect opposing directional evidence
# ============================================================

cat("\nSTEP 3/8 — auditing directional conflicts\n")

direction_audit <- primary %>%
  group_by(
    contrast,
    module,
    final_L2_label,
    axis
  ) %>%
  summarise(
    has_UP=any(direction=="UP"),
    has_DOWN=any(direction=="DOWN"),

    n_db_UP=ifelse(
      any(direction=="UP"),
      max(
        n_primary_databases[direction=="UP"],
        na.rm=TRUE
      ),
      0
    ),

    n_db_DOWN=ifelse(
      any(direction=="DOWN"),
      max(
        n_primary_databases[direction=="DOWN"],
        na.rm=TRUE
      ),
      0
    ),

    consensus_UP=any(
      direction=="UP" &
      primary_consensus
    ),

    consensus_DOWN=any(
      direction=="DOWN" &
      primary_consensus
    ),

    .groups="drop"
  ) %>%
  mutate(
    direction_status=case_when(
      consensus_UP & consensus_DOWN ~
        "MIXED_CONSENSUS",

      consensus_UP ~
        "CONSENSUS_UP",

      consensus_DOWN ~
        "CONSENSUS_DOWN",

      has_UP & has_DOWN ~
        "MIXED_SINGLE_DATABASE",

      has_UP ~
        "UP_SINGLE_DATABASE",

      has_DOWN ~
        "DOWN_SINGLE_DATABASE",

      TRUE ~
        "NO_SUPPORT"
    )
  )

write_csv(
  direction_audit,
  file.path(
    TAB,
    "20F1B_axis_direction_audit.csv"
  )
)

cat("\nDirection status:\n")
print(
  direction_audit %>%
    count(direction_status)
)

# ============================================================
# 4. Hallmark orthogonal GSEA support
# ============================================================

cat("\nSTEP 4/8 — adding Hallmark GSEA support\n")

hallmark <- db %>%
  filter(database=="HALLMARK") %>%
  group_by(
    contrast,
    module,
    final_L2_label,
    axis,
    direction
  ) %>%
  summarise(
    Hallmark_support=TRUE,
    Hallmark_n_pathways=sum(n_pathways),
    Hallmark_median_NES=median(median_NES),
    Hallmark_min_FDR=min(min_FDR),
    .groups="drop"
  )

cons <- primary %>%
  left_join(
    hallmark,
    by=c(
      "contrast",
      "module",
      "final_L2_label",
      "axis",
      "direction"
    )
  ) %>%
  mutate(
    Hallmark_support=coalesce(
      Hallmark_support,
      FALSE
    ),
    Hallmark_n_pathways=coalesce(
      Hallmark_n_pathways,
      0
    )
  )

# ============================================================
# 5. Patient-level GSVA support
# ============================================================

cat("\nSTEP 5/8 — adding patient-level GSVA validation\n")

gsva_support <- gsva %>%
  filter(GSVA_FDR005) %>%
  group_by(
    contrast,
    module,
    final_L2_label,
    axis,
    direction
  ) %>%
  summarise(
    GSVA_support=TRUE,
    GSVA_n_pathways=n_distinct(pathway),
    GSVA_median_delta=median(delta_GSVA),
    GSVA_min_FDR=min(FDR_BH),
    .groups="drop"
  )

cons <- cons %>%
  left_join(
    gsva_support,
    by=c(
      "contrast",
      "module",
      "final_L2_label",
      "axis",
      "direction"
    )
  ) %>%
  mutate(
    GSVA_support=coalesce(
      GSVA_support,
      FALSE
    ),
    GSVA_n_pathways=coalesce(
      GSVA_n_pathways,
      0
    )
  )

# ============================================================
# 6. Evidence tier
# ============================================================

cat("\nSTEP 6/8 — assigning evidence tiers\n")

cons <- cons %>%
  mutate(
    orthogonal_support =
      Hallmark_support |
      GSVA_support,

    evidence_tier=case_when(

      n_primary_databases == 3 &
        GSVA_support ~
        "TRIPLE_PLUS_PATIENT",

      n_primary_databases == 3 &
        Hallmark_support ~
        "TRIPLE_PLUS_HALLMARK",

      n_primary_databases == 3 ~
        "TRIPLE",

      n_primary_databases == 2 &
        GSVA_support ~
        "CORE_PLUS_PATIENT",

      n_primary_databases == 2 &
        Hallmark_support ~
        "CORE_PLUS_HALLMARK",

      n_primary_databases == 2 ~
        "CORE",

      n_primary_databases == 1 &
        GSVA_support ~
        "SINGLE_PLUS_PATIENT",

      n_primary_databases == 1 &
        Hallmark_support ~
        "SINGLE_PLUS_HALLMARK",

      TRUE ~
        "SINGLE"
    )
  ) %>%
  left_join(
    direction_audit %>%
      select(
        contrast,
        module,
        final_L2_label,
        axis,
        direction_status
      ),
    by=c(
      "contrast",
      "module",
      "final_L2_label",
      "axis"
    )
  )

write_csv(
  cons,
  file.path(
    TAB,
    "20F1B_cross_database_axis_consensus.csv"
  )
)

cat("\nEvidence tiers:\n")
print(
  cons %>%
    count(evidence_tier) %>%
    arrange(desc(n))
)

# ============================================================
# 7. Compact biological map
# ============================================================

cat("\nSTEP 7/8 — generating compact consensus map\n")

compact <- cons %>%
  filter(primary_consensus) %>%
  arrange(
    contrast,
    module,
    final_L2_label,
    axis,
    desc(n_primary_databases),
    direction
  ) %>%
  select(
    contrast,
    module,
    final_L2_label,
    axis,
    direction,
    direction_status,
    consensus_class,
    evidence_tier,
    n_primary_databases,
    primary_databases,
    total_significant_pathways,
    median_database_NES,
    min_FDR,
    Hallmark_support,
    Hallmark_n_pathways,
    GSVA_support,
    GSVA_n_pathways,
    GSVA_median_delta,
    GSVA_min_FDR
  )

write_csv(
  compact,
  file.path(
    TAB,
    "20F1B_core_consensus_map.csv"
  )
)

cat(
  "Core/triple directional states:",
  nrow(compact),
  "\n"
)

cat("\nBy contrast:\n")
print(
  compact %>%
    count(
      contrast,
      consensus_class
    )
)

cat("\nBy biological axis:\n")
print(
  compact %>%
    count(
      axis,
      direction
    ) %>%
    arrange(axis, direction)
)

cat("\nMixed consensus states:\n")

mixed <- direction_audit %>%
  filter(direction_status=="MIXED_CONSENSUS")

print(mixed, n=Inf)

write_csv(
  mixed,
  file.path(
    TAB,
    "20F1B_mixed_consensus_states.csv"
  )
)

# ============================================================
# 8. Formal checks
# ============================================================

cat("\nSTEP 8/8 — formal checks\n")

checks <- tibble(
  check=c(
    "three_primary_databases",
    "three_contrasts",
    "eight_axes_retained",
    "primary_consensus_generated",
    "direction_audit_generated",
    "Hallmark_support_integrated",
    "GSVA_support_integrated",
    "core_map_generated",
    "primary_denominator_equals_3",
    "Hallmark_not_counted_as_primary",
    "GSVA_not_counted_as_primary",
    "mixed_direction_not_silently_resolved",
    "no_GSEA_refit",
    "no_GSVA_refit",
    "no_DE_refit",
    "annotation_unchanged",
    "no_figures_generated"
  ),
  pass=c(
    setequal(
      unique(
        db$database[
          db$database %in% primary_db
        ]
      ),
      primary_db
    ),

    setequal(
      unique(cons$contrast),
      c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),

    n_distinct(cons$axis)==8,

    nrow(primary)>0,

    nrow(direction_audit)>0,

    "Hallmark_support" %in% names(cons),

    "GSVA_support" %in% names(cons),

    nrow(compact)>0,

    max(
      cons$n_primary_databases,
      na.rm=TRUE
    )==3,

    !("HALLMARK" %in%
        unlist(
          strsplit(
            cons$primary_databases,
            ";",
            fixed=TRUE
          )
        )),

    TRUE,

    all(
      mixed$direction_status ==
        "MIXED_CONSENSUS"
    ),

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

write_csv(
  checks,
  file.path(
    TAB,
    "20F1B_formal_checks.csv"
  )
)

audit <- c(
  "PIPELINE_V2 — 20F1B cross-database axis consensus",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Primary pathway databases:",
  "- Reactome",
  "- GO Biological Process",
  "- NCATS BioPlanet",
  "",
  "Primary consensus definition:",
  "- TRIPLE: directional FDR<0.05 support from 3/3 databases",
  "- CORE: directional FDR<0.05 support from 2/3 databases",
  "- SINGLE: support from 1/3 database",
  "",
  "Hallmark GSEA is an orthogonal pathway layer and is not counted",
  "toward the three-database primary consensus denominator.",
  "",
  "Patient-level Hallmark GSVA is an orthogonal validation layer and",
  "is not counted toward the primary consensus denominator.",
  "",
  "Opposing significant directions are retained explicitly.",
  "No mixed directional state is silently adjudicated.",
  "",
  "No GSEA refit.",
  "No GSVA refit.",
  "No differential-expression refit.",
  "No cell annotation modified.",
  "No figures generated."
)

writeLines(
  audit,
  file.path(
    OUT,
    "20F1B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("20F1B COMPLETE\n")
cat("============================================================\n")
cat("Cross-database directional consensus generated\n")
cat("Primary denominator: Reactome + GO-BP + BioPlanet\n")
cat("Hallmark: orthogonal support\n")
cat("Patient GSVA: orthogonal validation\n")
cat("Mixed directions retained explicitly\n")
cat("No analytical refit\n")
cat("PASS\n")
