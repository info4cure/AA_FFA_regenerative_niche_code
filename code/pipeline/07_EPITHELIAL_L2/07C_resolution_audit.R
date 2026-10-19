suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07C EPITHELIAL RESOLUTION AUDIT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07C_resolution_audit"
TDIR <- file.path(OUT,"tables")

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==63893,
  nrow(x)==23412,
  length(unique(x$patient_id))==16
)

resolutions <- c(0.2,0.3,0.4,0.5,0.6,0.8,1.0)

rcol <- function(r)
  paste0("epi_snn_res.",format(r,trim=TRUE,scientific=FALSE))

stopifnot(
  all(vapply(
    resolutions,
    function(r) rcol(r) %in% colnames(x@meta.data),
    logical(1)
  ))
)

# ============================================================
# 1. CORRECT PATIENT SUPPORT
# ============================================================

cat("STEP 1/5 — patient support\n")

support_all <- list()

for(r in resolutions){

  cc <- rcol(r)

  z <- x@meta.data %>%
    mutate(
      cluster=as.character(.data[[cc]])
    ) %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="patient_cells"
    )

  totals <- z %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(patient_cells),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      max_patient_cells=max(patient_cells),
      max_patient_fraction=
        max(patient_cells)/sum(patient_cells),
      .groups="drop"
    ) %>%
    mutate(resolution=r)

  support_all[[as.character(r)]] <- totals
}

support <- bind_rows(support_all)

write_csv(
  support,
  file.path(TDIR,"07C_patient_support.csv")
)

support_summary <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt5_patients=sum(n_patients<5),
    clusters_single_batch=sum(n_batches<2),
    clusters_patient_gt50=sum(max_patient_fraction>0.50),
    max_patient_fraction=
      max(max_patient_fraction),
    .groups="drop"
  )

cat("\n=== CORRECTED PATIENT SUPPORT ===\n")
print(
  as.data.frame(support_summary),
  row.names=FALSE
)

write_csv(
  support_summary,
  file.path(TDIR,"07C_patient_support_summary.csv")
)

# ============================================================
# 2. DISEASE COMPOSITION
# ============================================================

cat("\nSTEP 2/5 — disease composition\n")

disease_comp <- bind_rows(
  lapply(resolutions,function(r){

    cc <- rcol(r)

    x@meta.data %>%
      mutate(cluster=as.character(.data[[cc]])) %>%
      count(cluster,disease,name="n_cells") %>%
      group_by(cluster) %>%
      mutate(
        cluster_total=sum(n_cells),
        fraction=n_cells/cluster_total
      ) %>%
      ungroup() %>%
      mutate(resolution=r)
  })
)

write_csv(
  disease_comp,
  file.path(TDIR,"07C_disease_composition.csv")
)

# This is descriptive QC only, NOT differential abundance.
disease_flags <- disease_comp %>%
  group_by(resolution,cluster) %>%
  summarise(
    max_disease_fraction=max(fraction),
    dominant_disease=disease[which.max(fraction)],
    .groups="drop"
  )

write_csv(
  disease_flags,
  file.path(TDIR,"07C_disease_dominance.csv")
)

# ============================================================
# 3. BATCH MIXING IN IDENTIFIABLE SUBSETS
#
# AA cannot diagnose batch mixing because all AA are CABIMER.
# Evaluate HC and FFA, where both batches are represented.
# ============================================================

cat("STEP 3/5 — batch mixing in HC/FFA\n")

batch_mix <- bind_rows(
  lapply(resolutions,function(r){

    cc <- rcol(r)

    x@meta.data %>%
      filter(disease %in% c("HC","FFA")) %>%
      mutate(cluster=as.character(.data[[cc]])) %>%
      count(
        disease,
        cluster,
        batch,
        name="n_cells"
      ) %>%
      group_by(disease,cluster) %>%
      mutate(
        total=sum(n_cells),
        fraction=n_cells/total
      ) %>%
      ungroup() %>%
      mutate(resolution=r)
  })
)

write_csv(
  batch_mix,
  file.path(TDIR,"07C_HC_FFA_batch_mix.csv")
)

batch_support <- batch_mix %>%
  filter(n_cells>0) %>%
  group_by(resolution,disease,cluster) %>%
  summarise(
    n_batches=n_distinct(batch),
    total=max(total),
    max_batch_fraction=max(fraction),
    .groups="drop"
  )

write_csv(
  batch_support,
  file.path(TDIR,"07C_HC_FFA_batch_support.csv")
)

# ============================================================
# 4. RESOLUTION TRANSITIONS
#
# Quantify how each parent cluster splits at next resolution.
# ============================================================

cat("STEP 4/5 — resolution transitions\n")

pairs <- list(
  c(0.2,0.3),
  c(0.3,0.4),
  c(0.4,0.5),
  c(0.5,0.6),
  c(0.6,0.8),
  c(0.8,1.0)
)

transitions <- list()

for(pp in pairs){

  r1 <- pp[1]
  r2 <- pp[2]

  c1 <- rcol(r1)
  c2 <- rcol(r2)

  z <- x@meta.data %>%
    transmute(
      parent=as.character(.data[[c1]]),
      child=as.character(.data[[c2]])
    ) %>%
    count(parent,child,name="n_cells") %>%
    group_by(parent) %>%
    mutate(
      parent_total=sum(n_cells),
      fraction_parent=n_cells/parent_total
    ) %>%
    ungroup() %>%
    mutate(
      resolution_from=r1,
      resolution_to=r2
    )

  transitions[[paste(r1,r2,sep="_")]] <- z
}

trans <- bind_rows(transitions)

write_csv(
  trans,
  file.path(TDIR,"07C_resolution_transitions.csv")
)

transition_summary <- trans %>%
  group_by(
    resolution_from,
    resolution_to,
    parent
  ) %>%
  summarise(
    parent_cells=max(parent_total),
    n_children=sum(fraction_parent>=0.05),
    dominant_child_fraction=max(fraction_parent),
    entropy=-sum(
      ifelse(
        fraction_parent>0,
        fraction_parent*log(fraction_parent),
        0
      )
    ),
    .groups="drop"
  )

write_csv(
  transition_summary,
  file.path(TDIR,"07C_transition_summary.csv")
)

pair_summary <- transition_summary %>%
  group_by(resolution_from,resolution_to) %>%
  summarise(
    n_parent=n(),
    parents_split_ge2=sum(n_children>=2),
    median_dominant_child=
      median(dominant_child_fraction),
    median_entropy=median(entropy),
    .groups="drop"
  )

cat("\n=== RESOLUTION TRANSITIONS ===\n")
print(
  as.data.frame(pair_summary),
  row.names=FALSE
)

write_csv(
  pair_summary,
  file.path(TDIR,"07C_transition_pair_summary.csv")
)

# ============================================================
# 5. CANDIDATE-RESOLUTION REPORT
# ============================================================

cat("\nSTEP 5/5 — candidate resolution report\n")

candidate <- support_summary %>%
  left_join(
    disease_flags %>%
      group_by(resolution) %>%
      summarise(
        max_disease_dominance=
          max(max_disease_fraction),
        .groups="drop"
      ),
    by="resolution"
  )

write_csv(
  candidate,
  file.path(TDIR,"07C_candidate_resolution_metrics.csv")
)

cat("\n=== CANDIDATE RESOLUTION METRICS ===\n")
print(
  as.data.frame(candidate),
  row.names=FALSE
)

cat("\n============================================================\n")
cat("07C COMPLETE\n")
cat("============================================================\n")
cat("No resolution selected.\n")
cat("No markers calculated.\n")
cat("No L2 labels assigned.\n")
cat("PASS — corrected support/stability audit complete.\n")
