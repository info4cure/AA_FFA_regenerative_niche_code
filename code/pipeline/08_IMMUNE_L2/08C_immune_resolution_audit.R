suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08C IMMUNE RESOLUTION AUDIT\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/immune_v2_08B_reclustered.rds"
OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08C_resolution_audit"
TDIR <- file.path(OUT,"tables")
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==4762,
  nrow(x)==23412,
  length(unique(x$patient_id))==16
)

resolutions <- c(0.2,0.3,0.4,0.5,0.6,0.8,1.0)

res_cols <- paste0("immune_snn_res.",resolutions)

stopifnot(all(res_cols %in% colnames(x@meta.data)))

md <- x@meta.data
md$cell <- rownames(md)

cat("Cells:",nrow(md),"\n")
cat("Patients:",n_distinct(md$patient_id),"\n")

# ============================================================
# 1. PATIENT / BATCH SUPPORT
# ============================================================

cat("\nSTEP 1/5 — patient and batch support\n")

support_all <- list()

for(i in seq_along(resolutions)){

  res <- resolutions[i]
  cc <- res_cols[i]

  z <- md %>%
    mutate(cluster=as.character(.data[[cc]])) %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n_cells"
    ) %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n_cells),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      dominant_patient=patient_id[which.max(n_cells)],
      max_patient_cells=max(n_cells),
      max_patient_fraction=max(n_cells)/sum(n_cells),
      .groups="drop"
    ) %>%
    mutate(resolution=res)

  support_all[[i]] <- z
}

support <- bind_rows(support_all)

write_csv(
  support,
  file.path(TDIR,"08C_patient_support.csv")
)

support_summary <- support %>%
  group_by(resolution) %>%
  summarise(
    n_clusters=n(),
    min_cluster_cells=min(cluster_cells),
    min_patients=min(n_patients),
    median_patients=median(n_patients),
    clusters_lt5_patients=sum(n_patients<5),
    clusters_single_batch=sum(n_batches<2),
    clusters_patient_gt50=sum(max_patient_fraction>0.50),
    max_patient_fraction=max(max_patient_fraction),
    .groups="drop"
  )

cat("\n=== PATIENT SUPPORT SUMMARY ===\n")
print(as.data.frame(support_summary),row.names=FALSE)

write_csv(
  support_summary,
  file.path(TDIR,"08C_patient_support_summary.csv")
)

# ============================================================
# 2. DISEASE COMPOSITION
# ============================================================

cat("\nSTEP 2/5 — disease composition\n")

disease_all <- list()

for(i in seq_along(resolutions)){

  res <- resolutions[i]
  cc <- res_cols[i]

  z <- md %>%
    mutate(cluster=as.character(.data[[cc]])) %>%
    count(cluster,disease,name="n_cells") %>%
    group_by(cluster) %>%
    mutate(
      cluster_total=sum(n_cells),
      disease_fraction=n_cells/cluster_total
    ) %>%
    ungroup() %>%
    mutate(resolution=res)

  disease_all[[i]] <- z
}

disease_comp <- bind_rows(disease_all)

write_csv(
  disease_comp,
  file.path(TDIR,"08C_disease_composition.csv")
)

disease_dom <- disease_comp %>%
  group_by(resolution,cluster) %>%
  slice_max(
    disease_fraction,
    n=1,
    with_ties=FALSE
  ) %>%
  ungroup() %>%
  transmute(
    resolution,
    cluster,
    dominant_disease=disease,
    dominant_disease_fraction=disease_fraction
  )

write_csv(
  disease_dom,
  file.path(TDIR,"08C_disease_dominance.csv")
)

cat("\n=== MAX DISEASE DOMINANCE ===\n")

print(
  as.data.frame(
    disease_dom %>%
      group_by(resolution) %>%
      summarise(
        max_disease_fraction=max(dominant_disease_fraction),
        clusters_gt80pct=sum(dominant_disease_fraction>0.80),
        clusters_gt90pct=sum(dominant_disease_fraction>0.90),
        .groups="drop"
      )
  ),
  row.names=FALSE
)

# ============================================================
# 3. SOURCE L1 COMPOSITION
#
# Diagnostic only: does de novo immune clustering separate
# original T_NK and Myeloid_DC, and where do mixed clusters occur?
# ============================================================

cat("\nSTEP 3/5 — source L1 composition\n")

source_all <- list()

for(i in seq_along(resolutions)){

  res <- resolutions[i]
  cc <- res_cols[i]

  z <- md %>%
    mutate(cluster=as.character(.data[[cc]])) %>%
    count(
      cluster,
      immune_source_L1,
      name="n_cells"
    ) %>%
    group_by(cluster) %>%
    mutate(
      cluster_total=sum(n_cells),
      fraction=n_cells/cluster_total
    ) %>%
    ungroup() %>%
    mutate(resolution=res)

  source_all[[i]] <- z
}

source_comp <- bind_rows(source_all)

write_csv(
  source_comp,
  file.path(TDIR,"08C_source_L1_composition.csv")
)

# ============================================================
# 4. RESOLUTION TRANSITIONS
# ============================================================

cat("\nSTEP 4/5 — resolution transitions\n")

transition_rows <- list()
transition_summary <- list()

for(i in seq_len(length(resolutions)-1)){

  r1 <- resolutions[i]
  r2 <- resolutions[i+1]

  c1 <- res_cols[i]
  c2 <- res_cols[i+1]

  z <- data.frame(
    parent=as.character(md[[c1]]),
    child=as.character(md[[c2]])
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

  transition_rows[[i]] <- z

  ss <- z %>%
    group_by(parent) %>%
    summarise(
      dominant_child=max(fraction_parent),
      n_children_05=sum(fraction_parent>=0.05),
      entropy=-sum(
        ifelse(
          fraction_parent>0,
          fraction_parent*log(fraction_parent),
          0
        )
      ),
      .groups="drop"
    )

  transition_summary[[i]] <- data.frame(
    resolution_from=r1,
    resolution_to=r2,
    n_parent=nrow(ss),
    parents_split_ge2=sum(ss$n_children_05>=2),
    median_dominant_child=median(ss$dominant_child),
    median_entropy=median(ss$entropy)
  )
}

transitions <- bind_rows(transition_rows)
trans_summary <- bind_rows(transition_summary)

write_csv(
  transitions,
  file.path(TDIR,"08C_resolution_transitions.csv")
)

write_csv(
  trans_summary,
  file.path(TDIR,"08C_resolution_transition_summary.csv")
)

cat("\n=== RESOLUTION TRANSITIONS ===\n")
print(as.data.frame(trans_summary),row.names=FALSE)

# ============================================================
# 5. CANDIDATE RESOLUTION REPORT
# ============================================================

cat("\nSTEP 5/5 — candidate resolution report\n")

candidate <- support_summary %>%
  left_join(
    disease_dom %>%
      group_by(resolution) %>%
      summarise(
        max_disease_dominance=max(dominant_disease_fraction),
        clusters_disease_gt90=sum(dominant_disease_fraction>0.90),
        .groups="drop"
      ),
    by="resolution"
  )

write_csv(
  candidate,
  file.path(TDIR,"08C_candidate_resolution_metrics.csv")
)

cat("\n=== CANDIDATE RESOLUTION METRICS ===\n")
print(as.data.frame(candidate),row.names=FALSE)

# Detailed 0.4 / 0.5 / 0.6 reports for next marker step
for(res in c(0.4,0.5,0.6)){

  cat("\n------------------------------------------------------------\n")
  cat("RESOLUTION",res,"\n")

  zz <- support %>%
    filter(resolution==res) %>%
    arrange(as.numeric(cluster))

  print(as.data.frame(zz),row.names=FALSE)
}

audit <- c(
  "PIPELINE_V2 — 08C immune resolution audit",
  paste0("Generated: ",Sys.time()),
  "Input: immune_v2_08B_reclustered.rds",
  "Cells: 4762",
  "Resolutions audited: 0.2,0.3,0.4,0.5,0.6,0.8,1.0",
  "Metrics: patient support, batch support, disease composition, source-L1 composition, resolution transitions",
  "Disease composition is diagnostic only and is not used to define identity",
  "No resolution selected",
  "No markers calculated",
  "No immune identities assigned",
  "No object modified"
)

writeLines(
  audit,
  file.path(TDIR,"08C_audit.txt")
)

cat("\n============================================================\n")
cat("08C COMPLETE\n")
cat("============================================================\n")
cat("No resolution selected.\n")
cat("No immune identities assigned.\n")
cat("PASS — immune support/stability audit complete.\n")
