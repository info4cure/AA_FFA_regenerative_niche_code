suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09C FIBROBLAST RESOLUTION / STABILITY AUDIT\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09B_reclustered.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09C_resolution_audit/tables"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==4808,
  nrow(x)==23412,
  length(unique(x$patient_id))==16
)

resolutions <- c(0.2,0.3,0.4,0.5,0.6)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")

# ==========================================================================
# 1. PATIENT / BATCH SUPPORT
# ==========================================================================

cat("\nSTEP 1/6 — patient and batch support\n")

support_list <- list()

for(r in resolutions){

  cn <- paste0("fibroblast_snn_res.",r)

  stopifnot(cn %in% colnames(x@meta.data))

  md <- x@meta.data %>%
    mutate(cluster=as.character(.data[[cn]]))

  by_patient <- md %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n_cells"
    )

  s <- by_patient %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n_cells),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      dominant_patient=patient_id[which.max(n_cells)],
      dominant_patient_cells=max(n_cells),
      dominant_patient_fraction=max(n_cells)/sum(n_cells),
      .groups="drop"
    ) %>%
    mutate(resolution=r)

  support_list[[as.character(r)]] <- s
}

support <- bind_rows(support_list)

write_csv(
  support,
  file.path(OUT,"09C_patient_support.csv")
)

summary_support <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_cluster_cells=min(cluster_cells),
    median_cluster_cells=median(cluster_cells),
    max_cluster_cells=max(cluster_cells),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt5_patients=sum(n_patients<5),
    clusters_single_batch=sum(n_batches<2),
    clusters_patient_gt50=sum(dominant_patient_fraction>0.5),
    max_patient_fraction=max(dominant_patient_fraction),
    .groups="drop"
  )

cat("\n=== SUPPORT SUMMARY ===\n")
print(as.data.frame(summary_support),row.names=FALSE)

write_csv(
  summary_support,
  file.path(OUT,"09C_support_summary.csv")
)

# ==========================================================================
# 2. PATIENT-DOMINATED / LOW-SUPPORT CLUSTERS
# ==========================================================================

cat("\nSTEP 2/6 — patient-dominated and low-support clusters\n")

flagged <- support %>%
  filter(
    dominant_patient_fraction>0.5 |
    n_patients<5
  ) %>%
  arrange(
    resolution,
    desc(dominant_patient_fraction)
  )

cat("\n=== FLAGGED CLUSTERS ===\n")

if(nrow(flagged)>0){
  print(as.data.frame(flagged),row.names=FALSE)
} else {
  cat("None.\n")
}

write_csv(
  flagged,
  file.path(OUT,"09C_flagged_clusters.csv")
)

flagged_patient_list <- list()

if(nrow(flagged)>0){

  for(i in seq_len(nrow(flagged))){

    r <- flagged$resolution[i]
    cl <- as.character(flagged$cluster[i])
    cn <- paste0("fibroblast_snn_res.",r)

    z <- x@meta.data %>%
      mutate(cluster=as.character(.data[[cn]])) %>%
      filter(cluster==cl) %>%
      count(
        patient_id,
        disease,
        batch,
        name="n_cells"
      ) %>%
      mutate(
        resolution=r,
        cluster=cl,
        fraction=n_cells/sum(n_cells)
      )

    flagged_patient_list[[i]] <- z
  }
}

if(length(flagged_patient_list)>0){

  flagged_patient <- bind_rows(flagged_patient_list)

  cat("\n=== FLAGGED CLUSTER PATIENT COMPOSITION ===\n")
  print(
    as.data.frame(
      flagged_patient %>%
        arrange(
          resolution,
          cluster,
          desc(n_cells)
        )
    ),
    row.names=FALSE
  )

  write_csv(
    flagged_patient,
    file.path(
      OUT,
      "09C_flagged_cluster_patient_composition.csv"
    )
  )

} else {

  flagged_patient <- tibble()

}

# ==========================================================================
# 3. DISEASE COMPOSITION
# Descriptive QC only.
# ==========================================================================

cat("\nSTEP 3/6 — disease composition\n")

disease_list <- list()

for(r in resolutions){

  cn <- paste0("fibroblast_snn_res.",r)

  d <- x@meta.data %>%
    mutate(cluster=as.character(.data[[cn]])) %>%
    count(
      cluster,
      disease,
      name="n_cells"
    ) %>%
    group_by(cluster) %>%
    mutate(
      cluster_total=sum(n_cells),
      fraction=n_cells/cluster_total
    ) %>%
    ungroup() %>%
    mutate(resolution=r)

  disease_list[[as.character(r)]] <- d
}

disease <- bind_rows(disease_list)

write_csv(
  disease,
  file.path(OUT,"09C_disease_composition.csv")
)

disease_dom <- disease %>%
  group_by(resolution,cluster) %>%
  slice_max(
    fraction,
    n=1,
    with_ties=FALSE
  ) %>%
  ungroup() %>%
  select(
    resolution,
    cluster,
    dominant_disease=disease,
    dominant_disease_fraction=fraction
  )

disease_summary <- disease_dom %>%
  group_by(resolution) %>%
  summarise(
    max_disease_fraction=max(dominant_disease_fraction),
    clusters_gt80pct=sum(dominant_disease_fraction>0.8),
    clusters_gt90pct=sum(dominant_disease_fraction>0.9),
    .groups="drop"
  )

cat("\n=== MAX DISEASE DOMINANCE ===\n")
print(as.data.frame(disease_summary),row.names=FALSE)

write_csv(
  disease_dom,
  file.path(OUT,"09C_disease_dominance.csv")
)

# ==========================================================================
# 4. RESOLUTION TRANSITIONS
# ==========================================================================

cat("\nSTEP 4/6 — resolution transitions\n")

transition_table <- function(from,to){

  col_from <- paste0("fibroblast_snn_res.",from)
  col_to   <- paste0("fibroblast_snn_res.",to)

  a <- as.character(x@meta.data[[col_from]])
  b <- as.character(x@meta.data[[col_to]])

  stopifnot(
    length(a)==ncol(x),
    length(b)==ncol(x)
  )

  as.data.frame(
    table(parent=a,child=b),
    stringsAsFactors=FALSE
  ) %>%
    filter(Freq>0) %>%
    rename(n_cells=Freq) %>%
    group_by(parent) %>%
    mutate(
      parent_total=sum(n_cells),
      fraction_parent=n_cells/parent_total
    ) %>%
    ungroup() %>%
    mutate(
      resolution_from=from,
      resolution_to=to
    )
}

pairs <- list(
  c(0.2,0.3),
  c(0.3,0.4),
  c(0.4,0.5),
  c(0.5,0.6)
)

transitions <- bind_rows(
  lapply(
    pairs,
    function(z){
      transition_table(z[1],z[2])
    }
  )
)

write_csv(
  transitions,
  file.path(OUT,"09C_resolution_transitions.csv")
)

cat("\n=== TRANSITIONS >=5% OF PARENT ===\n")

print(
  transitions %>%
    filter(fraction_parent>=0.05) %>%
    arrange(
      resolution_from,
      parent,
      desc(fraction_parent)
    ),
  n=Inf
)

# ==========================================================================
# 5. TRUE SPLITS
# ==========================================================================

cat("\nSTEP 5/6 — true split detection\n")

true_splits <- transitions %>%
  filter(fraction_parent>=0.05) %>%
  group_by(
    resolution_from,
    resolution_to,
    parent
  ) %>%
  summarise(
    n_children=n(),
    largest_child_fraction=max(fraction_parent),
    second_child_fraction=
      if(n()>=2)
        sort(fraction_parent,decreasing=TRUE)[2]
      else
        NA_real_,
    .groups="drop"
  ) %>%
  filter(n_children>=2)

cat("\n=== PARENTS WITH >=2 CHILDREN CONTRIBUTING >=5% ===\n")

if(nrow(true_splits)>0){
  print(true_splits,n=Inf)
} else {
  cat("None.\n")
}

write_csv(
  true_splits,
  file.path(OUT,"09C_true_splits.csv")
)

# Detailed child transitions for true splits
true_split_edges <- transitions %>%
  semi_join(
    true_splits,
    by=c(
      "resolution_from",
      "resolution_to",
      "parent"
    )
  ) %>%
  filter(fraction_parent>=0.05) %>%
  arrange(
    resolution_from,
    parent,
    desc(fraction_parent)
  )

write_csv(
  true_split_edges,
  file.path(OUT,"09C_true_split_edges.csv")
)

# ==========================================================================
# 6. FORMAL CHECKS
# ==========================================================================

cat("\nSTEP 6/6 — formal checks\n")

checks <- data.frame(
  check=c(
    "cells_4808",
    "genes_23412",
    "patients_16",
    "all_candidate_resolutions_present",
    "support_complete",
    "transition_pairs_4",
    "transitions_generated",
    "no_resolution_selected",
    "no_FibConsensus_used",
    "no_L2_labels_assigned",
    "no_cells_removed"
  ),

  pass=c(
    ncol(x)==4808,
    nrow(x)==23412,
    length(unique(x$patient_id))==16,

    all(
      paste0(
        "fibroblast_snn_res.",
        resolutions
      ) %in% colnames(x@meta.data)
    ),

    nrow(support)>0,

    length(
      unique(
        paste(
          transitions$resolution_from,
          transitions$resolution_to,
          sep="_"
        )
      )
    )==4,

    nrow(transitions)>0,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),

  stringsAsFactors=FALSE
)

cat("\n=== FORMAL CHECKS ===\n")
print(checks,row.names=FALSE)

stopifnot(all(checks$pass))

write_csv(
  checks,
  file.path(OUT,"09C_formal_checks.csv")
)

audit <- c(
  "PIPELINE_V2 — 09C fibroblast resolution/stability audit",
  paste0("Generated: ",Sys.time()),
  "Input: fibroblast_v2_09B_reclustered.rds",
  "Candidate resolutions: 0.2, 0.3, 0.4, 0.5, 0.6.",
  "Patient and batch support audited.",
  "Patient-dominated and low-support clusters explicitly flagged.",
  "Disease composition evaluated as descriptive QC only.",
  "Disease was not used to define cluster identity.",
  "Transitions 0.2->0.3->0.4->0.5->0.6 evaluated.",
  "True parent-child splits defined using >=5% parent contribution.",
  "FibConsensus was not used.",
  "No resolution selected.",
  "No L2 identities assigned.",
  "No cells removed.",
  "Global L1 atlas unchanged."
)

writeLines(
  audit,
  file.path(OUT,"09C_audit.txt")
)

cat("\n============================================================\n")
cat("09C COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(x),"\n")
cat("No fibroblast resolution selected.\n")
cat("No fibroblast identities assigned.\n")
cat("FibConsensus not yet used.\n")
cat("PASS — resolution/stability evidence generated.\n")

