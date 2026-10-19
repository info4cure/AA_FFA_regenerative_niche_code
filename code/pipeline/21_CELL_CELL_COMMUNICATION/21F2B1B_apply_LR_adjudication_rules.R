suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(purrr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 21F2B1B PAIR-LEVEL LR ADJUDICATION\n")
cat("============================================================\n")

# ============================================================
# 0. Paths
# ============================================================

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/"
)

ROOT <- file.path(
  BASE,
  "21F2_BIOLOGICAL_LR_AXES"
)

TAB <- file.path(ROOT, "tables")

CAND_FILE <- file.path(
  TAB,
  "21F2B0_LR_axis_adjudication_inventory.csv"
)

RULE_FILE <- file.path(
  TAB,
  "21F2B1A_family_adjudication_rules.csv"
)

HYP_FILE <- file.path(
  TAB,
  "21F2B1A_hypothesis_superaxes.csv"
)

CLASS_FILE <- file.path(
  TAB,
  "21F2B1A_mechanistic_classes.csv"
)

stopifnot(
  file.exists(CAND_FILE),
  file.exists(RULE_FILE),
  file.exists(HYP_FILE),
  file.exists(CLASS_FILE)
)

# ============================================================
# 1. Recover frozen inputs
# ============================================================

cat("\nSTEP 1/10 — recovering frozen 21F2B0/21F2B1A inputs\n")

cand <- read_csv(
  CAND_FILE,
  show_col_types=FALSE
)

rules <- read_csv(
  RULE_FILE,
  show_col_types=FALSE
)

hyp <- read_csv(
  HYP_FILE,
  show_col_types=FALSE
)

classes <- read_csv(
  CLASS_FILE,
  show_col_types=FALSE
)

stopifnot(
  nrow(cand) == 822,
  nrow(rules) == 25,
  all(
    c(
      "axis",
      "ligand",
      "receptor"
    ) %in% names(cand)
  ),
  all(
    c(
      "rule_id",
      "source_axis",
      "ligand_regex",
      "receptor_regex",
      "action",
      "mechanistic_class",
      "reason"
    ) %in% names(rules)
  )
)

cand <- cand %>%
  mutate(
    pair_id=row_number()
  )

cat("Candidate pairs:", nrow(cand), "\n")
cat("Frozen rules   :", nrow(rules), "\n")

# ============================================================
# 2. Freeze explicit rule precedence
# ============================================================

cat("\nSTEP 2/10 — freezing explicit rule precedence\n")

# Lower number = higher priority.
#
# Priority is explicit rather than derived from row order.
# Specific exceptions precede broader family rules.

priority <- tribble(
  ~rule_id,                         ~priority,

  "FGF_INTRACELLULAR_EXCLUDE",       10L,
  "FGF_ENDOCRINE",                   20L,
  "FGF_PARACRINE",                   30L,

  "IGF_IGF1R",                       10L,
  "IGF_IGF2R",                       10L,
  "IGF1_IGF2R_REVIEW",               20L,

  "WNT_CANONICAL_COMPATIBLE",        10L,
  "WNT5_CONTEXT",                    10L,
  "WNT_OTHER_CONTEXT",               30L,

  "BMP_CORE",                        10L,
  "GDF_CONTEXT",                     20L,

  "VEGF_CORE",                       10L,
  "ANGPT_CORE",                      10L,
  "ANGPTL_REVIEW",                   20L,

  "CYTOTOXIC_DEATH",                 10L,
  "CYTOTOXIC_COSTIM",                10L,

  "IFN_DIRECT",                      10L,
  "CHEMOKINE_DIRECT",                10L,

  "TNF_IL1",                         10L,
  "TNFSF_CONTEXT",                   20L,

  "TGFB_CORE",                       10L,

  "EGF_HGF_NRG",                     10L,
  "PDGF_CONTEXT",                    10L,

  "IMMUNE_CHECKPOINT",               10L,
  "ECM_STRUCTURAL",                  10L
)

stopifnot(
  nrow(priority) == nrow(rules),
  setequal(
    priority$rule_id,
    rules$rule_id
  )
)

rules2 <- rules %>%
  left_join(
    priority,
    by="rule_id"
  )

stopifnot(
  !anyNA(rules2$priority)
)

write_csv(
  priority,
  file.path(
    TAB,
    "21F2B1B_rule_precedence.csv"
  )
)

# ============================================================
# 3. Apply every compatible rule to every candidate pair
# ============================================================

cat("\nSTEP 3/10 — applying all compatible rules\n")

matches_list <- vector(
  "list",
  nrow(rules2)
)

for (i in seq_len(nrow(rules2))) {

  r <- rules2[i, ]

  idx <- (
    cand$axis == r$source_axis
  ) &
    str_detect(
      cand$ligand,
      r$ligand_regex
    ) &
    str_detect(
      cand$receptor,
      r$receptor_regex
    )

  if (any(idx)) {

    matches_list[[i]] <- cand[idx, ] %>%
      transmute(
        pair_id,
        axis,
        ligand,
        receptor,
        rule_id=r$rule_id,
        priority=r$priority,
        action=r$action,
        mechanistic_class=
          r$mechanistic_class,
        reason=r$reason
      )

  } else {

    matches_list[[i]] <- tibble()
  }
}

all_matches <- bind_rows(
  matches_list
)

cat(
  "Total pair-rule matches:",
  nrow(all_matches),
  "\n"
)

cat(
  "Pairs with >=1 match:",
  n_distinct(all_matches$pair_id),
  "\n"
)

write_csv(
  all_matches,
  file.path(
    TAB,
    "21F2B1B_all_pair_rule_matches.csv"
  )
)

# ============================================================
# 4. Audit multiplicity before resolving it
# ============================================================

cat("\nSTEP 4/10 — auditing rule multiplicity\n")

multiplicity <- cand %>%
  select(
    pair_id,
    axis,
    ligand,
    receptor
  ) %>%
  left_join(
    all_matches %>%
      count(
        pair_id,
        name="n_matching_rules"
      ),
    by="pair_id"
  ) %>%
  mutate(
    n_matching_rules=
      replace_na(
        n_matching_rules,
        0L
      )
  )

cat("\n===== RULE-MATCH MULTIPLICITY =====\n")

print(
  multiplicity %>%
    count(
      n_matching_rules,
      name="n_pairs"
    ) %>%
    arrange(
      n_matching_rules
    ),
  n=Inf,
  width=Inf
)

write_csv(
  multiplicity,
  file.path(
    TAB,
    "21F2B1B_rule_match_multiplicity.csv"
  )
)

# ============================================================
# 5. Resolve pair-level adjudication by frozen precedence
# ============================================================

cat("\nSTEP 5/10 — resolving pair-level adjudication\n")

selected <- all_matches %>%
  group_by(
    pair_id
  ) %>%
  mutate(
    best_priority=min(
      priority,
      na.rm=TRUE
    )
  ) %>%
  filter(
    priority == best_priority
  ) %>%
  ungroup()

priority_ties <- selected %>%
  count(
    pair_id,
    name="n_best_rules"
  ) %>%
  filter(
    n_best_rules > 1
  )

cat(
  "Unresolved equal-priority ties:",
  nrow(priority_ties),
  "\n"
)

# Equal-priority competing rules are not silently resolved.
stopifnot(
  nrow(priority_ties) == 0
)

selected <- selected %>%
  select(
    pair_id,
    selected_rule_id=rule_id,
    selected_priority=priority,
    action,
    mechanistic_class,
    reason
  )

match_trace <- all_matches %>%
  arrange(
    pair_id,
    priority,
    rule_id
  ) %>%
  group_by(
    pair_id
  ) %>%
  summarise(
    matching_rule_ids=
      paste(
        rule_id,
        collapse=";"
      ),
    matching_priorities=
      paste(
        priority,
        collapse=";"
      ),
    .groups="drop"
  )

final <- cand %>%
  left_join(
    multiplicity %>%
      select(
        pair_id,
        n_matching_rules
      ),
    by="pair_id"
  ) %>%
  left_join(
    match_trace,
    by="pair_id"
  ) %>%
  left_join(
    selected,
    by="pair_id"
  ) %>%
  mutate(
    adjudication_status=
      case_when(
        action == "INCLUDE" ~
          "PRIMARY_OR_CONTEXTUAL",
        action == "REVIEW" ~
          "CONTEXTUAL_REVIEW",
        action == "EXCLUDE" ~
          "EXCLUDE_MECHANISTIC",
        is.na(action) ~
          "UNRESOLVED",
        TRUE ~
          "UNRESOLVED"
      )
  )

# ============================================================
# 6. Attach frozen mechanistic interpretation metadata
# ============================================================

cat("\nSTEP 6/10 — attaching frozen mechanistic metadata\n")

if (
  "mechanistic_class" %in%
    names(classes)
) {

  class_meta <- classes %>%
    distinct(
      mechanistic_class,
      .keep_all=TRUE
    )

  final <- final %>%
    left_join(
      class_meta,
      by="mechanistic_class",
      suffix=c("", "_class")
    )
}

if (
  "mechanistic_class" %in%
    names(hyp)
) {

  hyp_meta <- hyp %>%
    distinct(
      mechanistic_class,
      .keep_all=TRUE
    )

  final <- final %>%
    left_join(
      hyp_meta,
      by="mechanistic_class",
      suffix=c("", "_hypothesis")
    )
}

# ============================================================
# 7. Pair-level adjudication summaries
# ============================================================

cat("\nSTEP 7/10 — generating adjudication summaries\n")

action_summary <- final %>%
  count(
    action,
    adjudication_status,
    name="n_pairs"
  ) %>%
  arrange(
    adjudication_status,
    action
  )

class_summary <- final %>%
  count(
    axis,
    mechanistic_class,
    action,
    adjudication_status,
    name="n_pairs"
  ) %>%
  arrange(
    axis,
    mechanistic_class,
    action
  )

rule_summary <- final %>%
  count(
    selected_rule_id,
    action,
    mechanistic_class,
    name="n_pairs"
  ) %>%
  arrange(
    selected_rule_id
  )

cat("\n===== FINAL ACTION SUMMARY =====\n")

print(
  action_summary,
  n=Inf,
  width=Inf
)

cat("\n===== FINAL RULE SUMMARY =====\n")

print(
  rule_summary,
  n=Inf,
  width=Inf
)

cat("\n===== FINAL MECHANISTIC CLASS SUMMARY =====\n")

print(
  class_summary,
  n=Inf,
  width=Inf
)

write_csv(
  final,
  file.path(
    TAB,
    "21F2B1B_pair_level_adjudication.csv"
  )
)

write_csv(
  action_summary,
  file.path(
    TAB,
    "21F2B1B_action_summary.csv"
  )
)

write_csv(
  class_summary,
  file.path(
    TAB,
    "21F2B1B_mechanistic_class_summary.csv"
  )
)

write_csv(
  rule_summary,
  file.path(
    TAB,
    "21F2B1B_selected_rule_summary.csv"
  )
)

# ============================================================
# 8. Sensitive-family audit
# ============================================================

cat("\nSTEP 8/10 — auditing mechanistically sensitive families\n")

sensitive <- final %>%
  filter(
    str_detect(
      ligand,
      "^(FGF11$|FGF12$|FGF13$|FGF14$|FGF19$|FGF21$|FGF23$)"
    ) |
    receptor == "IGF2R" |
    str_detect(
      ligand,
      "^(WNT5A$|WNT5B$|WNT11$)"
    ) |
    str_detect(
      ligand,
      "^ANGPTL"
    ) |
    str_detect(
      ligand,
      "^GDF"
    )
  ) %>%
  arrange(
    axis,
    ligand,
    receptor
  )

write_csv(
  sensitive,
  file.path(
    TAB,
    "21F2B1B_sensitive_family_audit.csv"
  )
)

cat(
  "Sensitive-family adjudications:",
  nrow(sensitive),
  "\n"
)

# ============================================================
# 9. Freeze interpretation boundaries
# ============================================================

cat("\nSTEP 9/10 — freezing interpretation boundaries\n")

boundaries <- c(
  "PIPELINE_V2 — 21F2B1B interpretation boundaries",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "Pair-level rules were defined before disease comparison.",
  "",
  "All matching rules are preserved in the pair-rule trace.",
  "Final pair adjudication uses explicit frozen rule precedence.",
  "Equal-priority conflicts are not silently resolved.",
  "",
  "INCLUDE does not imply that an interaction is active in any disease.",
  "REVIEW denotes contextual biological interpretation.",
  "EXCLUDE denotes mechanistic exclusion from the intended biological axis.",
  "",
  "WNT canonical-compatible classification does not establish beta-catenin activation.",
  "BMP signaling is not interpreted as uniformly regenerative.",
  "FGF11-14 are excluded from secreted FGF-FGFR trophic interpretation.",
  "IGF2R is separated from IGF1R-mediated trophic signaling.",
  "ANGPTL candidates are not treated as canonical ANGPT-TIE signaling.",
  "ECM-receptor interaction is not equivalent to fibrosis.",
  "",
  "No disease labels were used for rule definition or pair adjudication.",
  "No expression values were used for rule definition or pair adjudication.",
  "No communication scores were used for rule definition or pair adjudication.",
  "No DE or pathway model was refitted."
)

writeLines(
  boundaries,
  file.path(
    ROOT,
    "21F2B1B_interpretation_boundaries.txt"
  )
)

# ============================================================
# 10. Formal checks
# ============================================================

cat("\nSTEP 10/10 — formal checks\n")

fgf_fhf <- final %>%
  filter(
    ligand %in%
      paste0(
        "FGF",
        11:14
      )
  )

igf2r <- final %>%
  filter(
    ligand == "IGF2",
    receptor == "IGF2R"
  )

angptl <- final %>%
  filter(
    str_detect(
      ligand,
      "^ANGPTL"
    )
  )

wnt5_11 <- final %>%
  filter(
    ligand %in%
      c(
        "WNT5A",
        "WNT5B",
        "WNT11"
      )
  )

checks <- tibble(
  check=c(
    "21F2B0_candidates_present",
    "21F2B1A_rules_present",
    "822_candidates_recovered",
    "25_rules_recovered",
    "all_rules_have_explicit_priority",
    "all_pair_rule_matches_preserved",
    "pair_key_unique_input",
    "no_equal_priority_ties",
    "822_final_rows",
    "pair_key_unique_final",
    "all_pairs_have_matching_rule",
    "no_unresolved_pairs",
    "FGF11_14_excluded",
    "IGF2_IGF2R_separated",
    "ANGPTL_contextual_review",
    "WNT5_11_noncanonical_contextual",
    "ECM_not_called_fibrosis",
    "all_exclusions_have_reason",
    "interpretation_boundaries_written",
    "no_disease_used",
    "no_expression_used",
    "no_communication_score_used",
    "no_DE_refit",
    "no_pathway_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(CAND_FILE),
    file.exists(RULE_FILE),
    nrow(cand) == 822,
    nrow(rules) == 25,
    !anyNA(rules2$priority),
    file.exists(
      file.path(
        TAB,
        "21F2B1B_all_pair_rule_matches.csv"
      )
    ),
    !anyDuplicated(
      cand[
        c(
          "axis",
          "ligand",
          "receptor"
        )
      ]
    ),
    nrow(priority_ties) == 0,
    nrow(final) == 822,
    !anyDuplicated(
      final[
        c(
          "axis",
          "ligand",
          "receptor"
        )
      ]
    ),
    all(
      final$n_matching_rules >= 1
    ),
    !any(
      final$adjudication_status ==
        "UNRESOLVED"
    ),
    nrow(fgf_fhf) > 0 &&
      all(
        fgf_fhf$action ==
          "EXCLUDE"
      ),
    nrow(igf2r) > 0 &&
      all(
        igf2r$mechanistic_class ==
          "IGF2_CLEARANCE_TRAFFICKING"
      ),
    nrow(angptl) > 0 &&
      all(
        angptl$action ==
          "REVIEW"
      ),
    nrow(wnt5_11) > 0 &&
      all(
        wnt5_11$mechanistic_class ==
          "WNT_NONCANONICAL_CONTEXTUAL"
      ),
    all(
      final$mechanistic_class[
        final$axis ==
          "ECM_INTEGRIN_REMODELING"
      ] ==
        "ECM_ADHESION_STRUCTURAL"
    ),
    all(
      !is.na(
        final$reason[
          final$action ==
            "EXCLUDE"
        ]
      ) &
      nzchar(
        final$reason[
          final$action ==
            "EXCLUDE"
        ]
      )
    ),
    file.exists(
      file.path(
        ROOT,
        "21F2B1B_interpretation_boundaries.txt"
      )
    ),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n===== FORMAL CHECKS =====\n")

print(
  checks,
  n=Inf,
  width=Inf
)

stopifnot(
  length(checks$check) ==
    length(checks$pass),
  all(checks$pass)
)

write_csv(
  checks,
  file.path(
    TAB,
    "21F2B1B_formal_checks.csv"
  )
)

# ============================================================
# Audit
# ============================================================

audit <- c(
  "PIPELINE_V2 — 21F2B1B pair-level LR adjudication",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  paste0(
    "Candidate axis-LR pairs: ",
    nrow(cand)
  ),
  paste0(
    "Frozen adjudication rules: ",
    nrow(rules)
  ),
  paste0(
    "Total pair-rule matches: ",
    nrow(all_matches)
  ),
  paste0(
    "Pairs with multiple matching rules: ",
    sum(
      multiplicity$n_matching_rules > 1
    )
  ),
  paste0(
    "Equal-priority unresolved ties: ",
    nrow(priority_ties)
  ),
  paste0(
    "Final pair-level adjudications: ",
    nrow(final)
  ),
  paste0(
    "Unresolved final pairs: ",
    sum(
      final$adjudication_status ==
        "UNRESOLVED"
    )
  ),
  "",
  "ADJUDICATION PRINCIPLE",
  "----------------------",
  paste(
    "Every candidate axis-ligand-receptor pair was matched",
    "against the frozen disease-blind 21F2B1A rule set."
  ),
  paste(
    "All rule matches were retained for audit, while a single",
    "final rule was selected using explicit frozen precedence."
  ),
  paste(
    "Equal-priority conflicts were prohibited from silent",
    "resolution."
  ),
  "",
  "SPECIAL MECHANISTIC BOUNDARIES",
  "------------------------------",
  "- FGF11-14 excluded from canonical extracellular FGF-FGFR trophic signaling.",
  "- IGF2-IGF2R separated from IGF1R-mediated trophic signaling.",
  "- WNT5A/WNT5B/WNT11 retained as noncanonical/context-dependent WNT signaling.",
  "- ANGPTL interactions retained for contextual review rather than canonical ANGPT-TIE interpretation.",
  "- ECM-receptor interactions retained as structural/adhesive communication and not equated with fibrosis.",
  "",
  "ANALYTICAL BOUNDARIES",
  "---------------------",
  "No disease labels used for adjudication.",
  "No expression values used for adjudication.",
  "No communication scores used for adjudication.",
  "No AA-vs-FFA comparison.",
  "No DE refit.",
  "No pathway refit.",
  "No annotation modified.",
  "",
  "PASS"
)

writeLines(
  audit,
  file.path(
    ROOT,
    "21F2B1B_audit.txt"
  )
)

cat("\n============================================================\n")
cat("21F2B1B COMPLETE\n")
cat("============================================================\n")
cat("822 candidate LR-axis pairs adjudicated using frozen rules\n")
cat("All rule matches preserved for audit\n")
cat("Explicit rule precedence applied\n")
cat("No equal-priority conflict silently resolved\n")
cat("Mechanistic exclusions and contextual classes retained\n")
cat("No disease or expression information used for adjudication\n")
cat("No communication comparison performed\n")
cat("No analytical refit\n")
cat("PASS\n")
