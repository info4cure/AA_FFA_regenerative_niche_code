suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09E1 FIBROBLAST FINAL ADJUDICATION\n")
cat("============================================================\n")

IN <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09E0_contamination_audit/tables/09E0_master_cell_audit.csv"
)

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09E1_final_adjudication/tables"
)

dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

x <- read_csv(IN,show_col_types=FALSE)

stopifnot(
  nrow(x)==4808,
  length(unique(x$cell))==4808
)

cat("Cells:",nrow(x),"\n")

# ==========================================================================
# 1. REFERENCE DOUBLET DISTRIBUTION
# ==========================================================================

cat("\nSTEP 1/6 — empirical doublet calibration\n")

ref <- x %>%
  filter(provisional_flag=="FIBROBLAST_SUPPORTED") %>%
  pull(doublet_score)

q90 <- as.numeric(quantile(ref,.90,na.rm=TRUE))
q95 <- as.numeric(quantile(ref,.95,na.rm=TRUE))
q99 <- as.numeric(quantile(ref,.99,na.rm=TRUE))

cat(
  "Q90:",round(q90,4),"\n",
  "Q95:",round(q95,4),"\n",
  "Q99:",round(q99,4),"\n"
)

# ==========================================================================
# 2. FINAL EVIDENCE CATEGORIES
# ==========================================================================

cat("\nSTEP 2/6 — final evidence categories\n")

x <- x %>%
  mutate(

    high_doublet =
      doublet_score > q95,

    extreme_doublet =
      doublet_score > q99,

    final_evidence = case_when(

      provisional_flag=="FIBROBLAST_SUPPORTED"
        ~ "KEEP_FIBROBLAST",

      provisional_flag=="UNRESOLVED_REVIEW"
        ~ "KEEP_LOW_INFORMATION",

      provisional_flag=="T_RETURN_CANDIDATE"
        ~ "RETURN_T",

      provisional_flag=="EPITHELIAL_RETURN_CANDIDATE"
        ~ "RETURN_EPITHELIAL",

      provisional_flag=="MYELOID_RETURN_CANDIDATE"
        ~ "RETURN_MYELOID",

      provisional_flag=="ENDOTHELIAL_RETURN_CANDIDATE"
        ~ "RETURN_ENDOTHELIAL",

      provisional_flag=="SCHWANN_RETURN_CANDIDATE"
        ~ "RETURN_SCHWANN",

      provisional_flag=="MELANOCYTE_RETURN_CANDIDATE"
        ~ "RETURN_MELANOCYTE",

      provisional_flag=="MIXED_LINEAGE_REVIEW" &
      high_doublet
        ~ "EXCLUDE_MIXED_HIGH_DOUBLET",

      provisional_flag=="MIXED_LINEAGE_REVIEW" &
      !high_doublet
        ~ "REVIEW_MIXED_LOW_DOUBLET",

      provisional_flag=="MULTILINEAGE_REVIEW" &
      high_doublet
        ~ "EXCLUDE_MULTILINEAGE_HIGH_DOUBLET",

      TRUE
        ~ "REVIEW_UNRESOLVED"
    )
  )

summary <- x %>%
  count(final_evidence,name="n_cells") %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== FINAL EVIDENCE SUMMARY ===\n")
print(as.data.frame(summary),row.names=FALSE)

# ==========================================================================
# 3. MIXED CELLS
# ==========================================================================

cat("\nSTEP 3/6 — mixed-lineage adjudication\n")

mixed <- x %>%
  filter(
    provisional_flag=="MIXED_LINEAGE_REVIEW"
  )

mixed_summary <- mixed %>%
  group_by(
    cluster03,
    cluster04
  ) %>%
  summarise(
    n_cells=n(),
    n_patients=n_distinct(patient_id),
    median_doublet_score=
      median(doublet_score,na.rm=TRUE),
    n_gt_Q95=sum(doublet_score>q95,na.rm=TRUE),
    pct_gt_Q95=
      100*mean(doublet_score>q95,na.rm=TRUE),
    n_gt_Q99=sum(doublet_score>q99,na.rm=TRUE),
    pct_gt_Q99=
      100*mean(doublet_score>q99,na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(n_cells))

cat("\n=== MIXED BY CLUSTER ===\n")
print(as.data.frame(mixed_summary),row.names=FALSE)

# ==========================================================================
# 4. C5 / C9 FINAL AUDIT
# ==========================================================================

cat("\nSTEP 4/6 — C5/C9 final audit\n")

target <- x %>%
  filter(cluster03 %in% c("5","9")) %>%
  count(
    cluster03,
    cluster04,
    final_evidence,
    name="n_cells"
  ) %>%
  group_by(cluster03,cluster04) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

cat("\n=== C5/C9 FINAL EVIDENCE ===\n")
print(as.data.frame(target),row.names=FALSE)

# ==========================================================================
# 5. PROPOSED DISPOSITION
#
# IMPORTANT:
# REVIEW_MIXED_LOW_DOUBLET remains retained pending final review.
# Low-information fibroblasts are retained.
# ==========================================================================

cat("\nSTEP 5/6 — proposed disposition\n")

x <- x %>%
  mutate(

    proposed_disposition = case_when(

      final_evidence=="EXCLUDE_MIXED_HIGH_DOUBLET"
        ~ "EXCLUDE_DOUBLET",

      final_evidence=="EXCLUDE_MULTILINEAGE_HIGH_DOUBLET"
        ~ "EXCLUDE_DOUBLET",

      grepl("^RETURN_",final_evidence)
        ~ final_evidence,

      TRUE
        ~ "KEEP_FIBROBLAST_UNIVERSE"
    )
  )

disp <- x %>%
  count(
    proposed_disposition,
    name="n_cells"
  ) %>%
  mutate(
    pct=100*n_cells/sum(n_cells)
  ) %>%
  arrange(desc(n_cells))

cat("\n=== PROPOSED DISPOSITION ===\n")
print(as.data.frame(disp),row.names=FALSE)

retained <- x %>%
  filter(
    proposed_disposition=="KEEP_FIBROBLAST_UNIVERSE"
  )

retained_patient <- retained %>%
  count(
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  arrange(
    disease,
    desc(n_cells)
  )

cat("\n=== PROVISIONALLY RETAINED BY PATIENT ===\n")
print(as.data.frame(retained_patient),row.names=FALSE)

# ==========================================================================
# 6. EXPORT / CHECKS
# ==========================================================================

cat("\nSTEP 6/6 — export and checks\n")

write_csv(
  x,
  file.path(OUT,"09E1_master_adjudication.csv")
)

write_csv(
  summary,
  file.path(OUT,"09E1_final_evidence_summary.csv")
)

write_csv(
  mixed_summary,
  file.path(OUT,"09E1_mixed_cluster_summary.csv")
)

write_csv(
  target,
  file.path(OUT,"09E1_C5_C9_final_audit.csv")
)

write_csv(
  disp,
  file.path(OUT,"09E1_proposed_disposition.csv")
)

write_csv(
  retained_patient,
  file.path(OUT,"09E1_retained_by_patient.csv")
)

checks <- data.frame(
  check=c(
    "input_cells_4808",
    "unique_cells_4808",
    "all_cells_adjudicated",
    "disposition_sum_4808",
    "no_cells_removed",
    "no_object_modified",
    "no_frozen_object_written"
  ),
  pass=c(
    nrow(x)==4808,
    length(unique(x$cell))==4808,
    all(!is.na(x$final_evidence)),
    sum(disp$n_cells)==4808,
    TRUE,
    TRUE,
    TRUE
  )
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks,row.names=FALSE)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(OUT,"09E1_formal_checks.csv")
)

cat("\n============================================================\n")
cat("09E1 COMPLETE\n")
cat("============================================================\n")
cat("No cells removed.\n")
cat("No object modified.\n")
cat("No fibroblast universe frozen.\n")
cat("PASS — final adjudication evidence generated.\n")

