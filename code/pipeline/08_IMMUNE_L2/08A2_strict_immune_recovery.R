suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 08A2 STRICT IMMUNE RECOVERY\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/global_atlas_v2_L1_frozen.rds"
OUT <- "PIPELINE_V2/results/annotation/08_IMMUNE_L2/08A2_strict_recovery"
TDIR <- file.path(OUT,"tables")
dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==82893,
  nrow(x)==23412
)

GC <- "RNA_snn_res.0.6"
stopifnot(GC %in% colnames(x@meta.data))

md <- x@meta.data
md$cell <- rownames(md)
md$global_cluster <- as.character(md[[GC]])
md$current_immune <- as.character(md$L1) %in% c("T_NK","Myeloid_DC")

cat("Cells:",ncol(x),"\n")
cat("Current broad immune:",sum(md$current_immune),"\n")

# ============================================================
# 1. STRICT MARKER SETS
#
# Deliberately exclude promiscuous HLA/IFN markers from
# identity recovery.
# ============================================================

sets <- list(

  Leukocyte_core=c(
    "PTPRC","CORO1A","LCP1","CD37"
  ),

  T_specific=c(
    "CD3D","CD3E","CD3G","TRAC"
  ),

  NK_specific=c(
    "NKG7","GNLY","KLRD1","KLRF1"
  ),

  B_specific=c(
    "CD79A","CD79B","MS4A1","CD37"
  ),

  Plasma_specific=c(
    "MZB1","JCHAIN","TNFRSF17","DERL3","SDC1"
  ),

  Myeloid_specific=c(
    "LYZ","LST1","FCER1G","TYROBP","AIF1"
  ),

  cDC1_specific=c(
    "CLEC9A","XCR1","BATF3"
  ),

  cDC2_specific=c(
    "CD1C","CLEC10A","FCER1A","CD1E"
  ),

  pDC_specific=c(
    "CLEC4C","IL3RA","TCF4","GZMB"
  ),

  Mast_specific=c(
    "TPSAB1","TPSB2","CPA3","KIT","MS4A2","HDC"
  )
)

coverage <- bind_rows(
  lapply(names(sets),function(nm){
    g <- sets[[nm]]
    present <- intersect(g,rownames(x))
    data.frame(
      programme=nm,
      n_present=length(present),
      n_requested=length(g),
      coverage=length(present)/length(g),
      missing=paste(setdiff(g,present),collapse=";")
    )
  })
)

cat("\n=== STRICT MARKER COVERAGE ===\n")
print(coverage,row.names=FALSE)

write_csv(
  coverage,
  file.path(TDIR,"08A2_marker_coverage.csv")
)

genes <- intersect(unique(unlist(sets)),rownames(x))

counts <- LayerData(
  x,
  assay="RNA",
  layer="counts"
)[genes,,drop=FALSE]

# ============================================================
# 2. NUMBER OF DETECTED GENES PER PROGRAMME / CELL
# ============================================================

detected <- lapply(sets,function(g){

  gg <- intersect(g,rownames(counts))

  if(length(gg)==0)
    return(rep(0,ncol(x)))

  Matrix::colSums(
    counts[gg,,drop=FALSE] > 0
  )
})

detected <- as.data.frame(detected)
rownames(detected) <- colnames(x)

# ============================================================
# 3. STRICT CANDIDATE DEFINITIONS
#
# These are candidate flags, NOT annotations.
# ============================================================

flag <- data.frame(
  cell=colnames(x),

  leukocyte_core =
    detected$Leukocyte_core >= 2,

  T_candidate =
    detected$T_specific >= 2,

  NK_candidate =
    detected$NK_specific >= 2,

  B_candidate =
    detected$B_specific >= 2,

  Plasma_candidate =
    detected$Plasma_specific >= 2,

  Myeloid_candidate =
    detected$Myeloid_specific >= 2,

  cDC1_candidate =
    detected$cDC1_specific >= 2,

  cDC2_candidate =
    detected$cDC2_specific >= 2,

  pDC_candidate =
    detected$pDC_specific >= 2,

  Mast_candidate =
    detected$Mast_specific >= 2,

  stringsAsFactors=FALSE
)

flag <- flag %>%
  left_join(
    md %>%
      select(
        cell,L1,global_cluster,
        patient_id,disease,batch,
        current_immune
      ),
    by="cell"
  )

# ============================================================
# 4. GLOBAL SUMMARY
# ============================================================

candidate_cols <- c(
  "T_candidate",
  "NK_candidate",
  "B_candidate",
  "Plasma_candidate",
  "Myeloid_candidate",
  "cDC1_candidate",
  "cDC2_candidate",
  "pDC_candidate",
  "Mast_candidate"
)

summary_rows <- list()

for(v in candidate_cols){

  z <- flag[flag[[v]],,drop=FALSE]

  summary_rows[[length(summary_rows)+1]] <-
    data.frame(
      candidate=v,
      n_cells=nrow(z),
      n_current_immune=sum(z$current_immune),
      n_outside_current_immune=sum(!z$current_immune),
      n_patients=n_distinct(z$patient_id),
      n_patients_outside=n_distinct(
        z$patient_id[!z$current_immune]
      )
    )
}

summ <- bind_rows(summary_rows)

cat("\n=== STRICT CANDIDATE SUMMARY ===\n")
print(summ,row.names=FALSE)

write_csv(
  summ,
  file.path(TDIR,"08A2_candidate_summary.csv")
)

# ============================================================
# 5. OUTSIDE-CURRENT-IMMUNE CANDIDATES
# ============================================================

outside_rows <- list()

for(v in candidate_cols){

  z <- flag %>%
    filter(
      !current_immune,
      .data[[v]]
    ) %>%
    count(
      L1,
      global_cluster,
      patient_id,
      disease,
      batch,
      name="n_cells"
    ) %>%
    mutate(candidate=v)

  outside_rows[[length(outside_rows)+1]] <- z
}

outside <- bind_rows(outside_rows)

write_csv(
  outside,
  file.path(TDIR,"08A2_candidates_outside_current_immune.csv")
)

cat("\n=== OUTSIDE CURRENT IMMUNE — SUMMARY BY L1 ===\n")

outside_l1 <- outside %>%
  group_by(candidate,L1) %>%
  summarise(
    n_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_clusters=n_distinct(global_cluster),
    .groups="drop"
  ) %>%
  arrange(candidate,desc(n_cells))

print(outside_l1,row.names=FALSE)

write_csv(
  outside_l1,
  file.path(TDIR,"08A2_outside_summary_by_L1.csv")
)

cat("\n=== OUTSIDE CURRENT IMMUNE — SUMMARY BY CLUSTER ===\n")

outside_cluster <- outside %>%
  group_by(candidate,L1,global_cluster) %>%
  summarise(
    n_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    .groups="drop"
  ) %>%
  arrange(candidate,desc(n_cells))

print(outside_cluster,row.names=FALSE)

write_csv(
  outside_cluster,
  file.path(TDIR,"08A2_outside_summary_by_cluster.csv")
)

# ============================================================
# 6. LEUKOCYTE-CORE FILTERED OUTSIDE CANDIDATES
#
# Particularly important for B/plasma/mast/DC recovery.
# ============================================================

strict_outside <- flag %>%
  filter(
    !current_immune,
    leukocyte_core,
    B_candidate |
      Plasma_candidate |
      Myeloid_candidate |
      cDC1_candidate |
      cDC2_candidate |
      pDC_candidate |
      Mast_candidate
  )

cat("\n=== LEUKOCYTE-CORE POSITIVE OUTSIDE CANDIDATES ===\n")
cat("Cells:",nrow(strict_outside),"\n")
cat("Patients:",n_distinct(strict_outside$patient_id),"\n")

strict_summary <- strict_outside %>%
  count(
    L1,
    global_cluster,
    patient_id,
    disease,
    batch,
    name="n_cells"
  )

print(strict_summary,row.names=FALSE)

write_csv(
  strict_outside,
  file.path(TDIR,"08A2_strict_outside_candidate_cells.csv")
)

write_csv(
  strict_summary,
  file.path(TDIR,"08A2_strict_outside_candidate_support.csv")
)

# ============================================================
# 7. INSIDE CURRENT IMMUNE: WHAT LINEAGES ARE PRESENT?
# ============================================================

inside <- flag %>%
  filter(current_immune)

inside_summary <- bind_rows(
  lapply(candidate_cols,function(v){

    z <- inside[inside[[v]],,drop=FALSE]

    data.frame(
      candidate=v,
      n_cells=nrow(z),
      pct_current_immune=
        100*nrow(z)/nrow(inside),
      n_patients=n_distinct(z$patient_id),
      n_HC=n_distinct(z$patient_id[z$disease=="HC"]),
      n_AA=n_distinct(z$patient_id[z$disease=="AA"]),
      n_FFA=n_distinct(z$patient_id[z$disease=="FFA"])
    )
  })
)

cat("\n=== LINEAGE CANDIDATES WITHIN CURRENT 4,762 IMMUNE CELLS ===\n")
print(inside_summary,row.names=FALSE)

write_csv(
  inside_summary,
  file.path(TDIR,"08A2_inside_lineage_summary.csv")
)

# ============================================================
# 8. MULTI-LINEAGE FLAGS INSIDE IMMUNE
# ============================================================

flag$n_lineage_candidates <- rowSums(
  flag[,candidate_cols,drop=FALSE]
)

multi <- flag %>%
  filter(
    current_immune,
    n_lineage_candidates>=2
  ) %>%
  count(
    L1,
    global_cluster,
    n_lineage_candidates,
    name="n_cells"
  )

write_csv(
  multi,
  file.path(TDIR,"08A2_multilineage_flags_inside_immune.csv")
)

# ============================================================
# 9. AUDIT
# ============================================================

audit <- c(
  "PIPELINE_V2 — 08A2 strict immune recovery",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_L1_frozen.rds",
  "Global cells: 82893",
  "Current broad immune: T_NK + Myeloid_DC = 4762 cells",
  "Recovery based on co-detection of >=2 relatively lineage-informative genes",
  "Promiscuous HLA/IFN genes excluded from strict identity recovery panels",
  "Candidate flags are evidence only",
  "No cells reassigned",
  "No L1 labels changed",
  "No clustering performed",
  "No object written"
)

writeLines(
  audit,
  file.path(TDIR,"08A2_audit.txt")
)

cat("\n============================================================\n")
cat("08A2 COMPLETE\n")
cat("============================================================\n")
cat("No cells reassigned.\n")
cat("PASS — strict immune recovery evidence generated.\n")
