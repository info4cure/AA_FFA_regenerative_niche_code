suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 07E EPITHELIAL ANATOMICAL/STATE EVIDENCE\n")
cat("============================================================\n\n")

IN <- "PIPELINE_V2/objects/annotated/epithelial_v2_07B_reclustered.rds"
MARKERS <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07D_markers/07D_markers_res0.5_top15.csv"
OUT <- "PIPELINE_V2/results/annotation/07_EPITHELIAL_L2/07E_anatomical_state"
TDIR <- file.path(OUT,"tables")

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)

x <- readRDS(IN)
DefaultAssay(x) <- "RNA"

CC <- "epi_snn_res.0.5"

stopifnot(
  ncol(x)==63893,
  nrow(x)==23412,
  CC %in% colnames(x@meta.data)
)

clusters <- sort(unique(as.integer(x@meta.data[[CC]])))
stopifnot(identical(clusters,0:18))

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Resolution: 0.5\n")
cat("Clusters:",length(clusters),"\n\n")

# ============================================================
# 1. IDENTITY PROGRAMMES
# Same anatomical framework used in 06B3.
# ============================================================

identity_programmes <- list(

  Epithelial_core=c(
    "EPCAM","KRT8","KRT18","KRT19",
    "KRT5","KRT14","KRT17","DSP"
  ),

  Keratinocyte_core=c(
    "KRT5","KRT14","KRT1","KRT10",
    "KRT2","KRT15","KRT17","KRTDAP"
  ),

  Epidermal_basal=c(
    "KRT5","KRT14","KRT15","KRT17",
    "TP63","COL17A1","ITGA6","KRT19"
  ),

  Epidermal_suprabasal=c(
    "KRT1","KRT10","KRT2","IVL",
    "KRTDAP","DSC1","DSP","SFN"
  ),

  Epidermal_terminal=c(
    "FLG","FLG2","LOR","CDSN",
    "TGM1","CASP14","KLK7","KPRP"
  ),

  HF_ORS=c(
    "KRT5","KRT14","KRT15","KRT17",
    "KRT19","SOX9","KRT75","KRT79"
  ),

  HF_stem_bulge=c(
    "SOX9","LHX2","KRT15","KRT19",
    "CD200","NFATC1","KRT5","KRT14"
  ),

  HF_IRS_companion=c(
    "KRT25","KRT27","KRT28","KRT71",
    "KRT72","KRT73","KRT74","TCHH",
    "KRT75"
  ),

  HF_hair_shaft=c(
    "KRT31","KRT32","KRT33A","KRT33B",
    "KRT34","KRT35","KRT36",
    "KRT81","KRT83","KRT85"
  ),

  HF_matrix_differentiation=c(
    "LEF1","SHH","HOXC13","FOXN1",
    "MSX2","DLX3","WNT3","WNT10B"
  ),

  Sebaceous=c(
    "PPARG","FASN","MOGAT1","MOGAT2",
    "AWAT1","AWAT2","PLIN2","THRSP"
  )
)

# ============================================================
# 2. STATE PROGRAMMES
# Kept separate from anatomical identity.
# ============================================================

state_programmes <- list(

  Cycling=c(
    "MKI67","TOP2A","UBE2C","BIRC5",
    "CENPF","TYMS","RRM2","PCNA"
  ),

  Epithelial_stress=c(
    "KRT6A","KRT6B","KRT6C","KRT16",
    "KRT17","S100A7","S100A8","S100A9"
  ),

  IFN_response=c(
    "CXCL9","CXCL10","CXCL11",
    "ISG15","IFI6","IFIT1","IFIT2",
    "IFIT3","RSAD2","OASL","BST2"
  ),

  MHCII_AP=c(
    "HLA-DRA","HLA-DRB1","HLA-DPA1",
    "HLA-DPB1","HLA-DQA1","HLA-DQB1",
    "CD74","CIITA"
  )
)

programmes <- c(identity_programmes,state_programmes)

# ============================================================
# 3. COVERAGE
# ============================================================

coverage <- bind_rows(
  lapply(names(programmes),function(p){

    requested <- programmes[[p]]
    present <- intersect(requested,rownames(x))

    data.frame(
      programme=p,
      class=ifelse(
        p %in% names(identity_programmes),
        "identity",
        "state"
      ),
      n_present=length(present),
      n_requested=length(requested),
      coverage=length(present)/length(requested),
      missing=paste(
        setdiff(requested,present),
        collapse=";"
      )
    )
  })
)

write_csv(
  coverage,
  file.path(TDIR,"07E_programme_coverage.csv")
)

cat("=== PROGRAMME COVERAGE ===\n")
print(as.data.frame(coverage),row.names=FALSE)

# ============================================================
# 4. EXPRESSION + CELL SCORES
# ============================================================

genes <- intersect(
  unique(unlist(programmes)),
  rownames(x)
)

mat <- LayerData(
  x,
  assay="RNA",
  layer="data"
)[genes,,drop=FALSE]

scores <- lapply(
  programmes,
  function(g){

    gg <- intersect(g,rownames(mat))

    Matrix::colMeans(
      mat[gg,,drop=FALSE]
    )
  }
)

score_mat <- do.call(cbind,scores)
rownames(score_mat) <- colnames(x)

# ============================================================
# 5. CLUSTER EVIDENCE
# ============================================================

md <- x@meta.data
md$cell <- rownames(md)
md$cluster <- as.character(md[[CC]])

rows <- list()

for(cl in as.character(0:18)){

  cells <- md$cell[md$cluster==cl]

  for(p in colnames(score_mat)){

    v <- score_mat[cells,p]

    rows[[length(rows)+1]] <- data.frame(
      cluster=cl,
      programme=p,
      class=ifelse(
        p %in% names(identity_programmes),
        "identity",
        "state"
      ),
      n_cells=length(cells),
      mean_programme=mean(v),
      median_programme=median(v),
      pct_positive=100*mean(v>0)
    )
  }
}

ev <- bind_rows(rows)

# Correct z-score: within programme across clusters.
ev <- ev %>%
  group_by(programme) %>%
  mutate(
    programme_z={
      s <- sd(mean_programme,na.rm=TRUE)
      if(is.finite(s) && s>0){
        as.numeric(scale(mean_programme))
      } else {
        rep(0,n())
      }
    }
  ) %>%
  ungroup()

write_csv(
  ev,
  file.path(TDIR,"07E_cluster_programme_evidence.csv")
)

# ============================================================
# 6. PATIENT-LEVEL PROGRAMME SUPPORT
# ============================================================

patient_rows <- list()

for(cl in as.character(0:18)){

  m <- md[md$cluster==cl,,drop=FALSE]

  for(pt in unique(m$patient_id)){

    cells <- m$cell[m$patient_id==pt]

    for(p in colnames(score_mat)){

      patient_rows[[length(patient_rows)+1]] <-
        data.frame(
          cluster=cl,
          patient_id=pt,
          disease=unique(
            m$disease[m$patient_id==pt]
          )[1],
          programme=p,
          class=ifelse(
            p %in% names(identity_programmes),
            "identity",
            "state"
          ),
          n_cells=length(cells),
          mean_programme=
            mean(score_mat[cells,p])
        )
    }
  }
}

patient_ev <- bind_rows(patient_rows)

write_csv(
  patient_ev,
  file.path(TDIR,"07E_patient_programme_evidence.csv")
)

# ============================================================
# 7. TOP IDENTITY + STATE PROGRAMMES
# ============================================================

identity_rank <- ev %>%
  filter(class=="identity") %>%
  group_by(cluster) %>%
  arrange(desc(programme_z),.by_group=TRUE) %>%
  mutate(rank=row_number()) %>%
  filter(rank<=5) %>%
  ungroup()

state_rank <- ev %>%
  filter(class=="state") %>%
  group_by(cluster) %>%
  arrange(desc(programme_z),.by_group=TRUE) %>%
  mutate(rank=row_number()) %>%
  ungroup()

write_csv(
  identity_rank,
  file.path(TDIR,"07E_top_identity_programmes.csv")
)

write_csv(
  state_rank,
  file.path(TDIR,"07E_state_programmes.csv")
)

# ============================================================
# 8. TOP MARKERS FROM 07D
# ============================================================

mk <- read_csv(
  MARKERS,
  show_col_types=FALSE
)

mk_summary <- mk %>%
  mutate(cluster=as.character(cluster)) %>%
  group_by(cluster) %>%
  summarise(
    top_markers=paste(gene,collapse=", "),
    .groups="drop"
  )

# ============================================================
# 9. PATIENT / DISEASE SUPPORT
# ============================================================

support <- md %>%
  count(
    cluster,
    patient_id,
    disease,
    name="n_cells"
  ) %>%
  group_by(cluster) %>%
  summarise(
    n_cells=sum(n_cells),
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(patient_id[disease=="HC"]),
    n_AA=n_distinct(patient_id[disease=="AA"]),
    n_FFA=n_distinct(patient_id[disease=="FFA"]),
    dominant_patient_fraction=max(n_cells)/sum(n_cells),
    .groups="drop"
  )

# ============================================================
# 10. COMPACT EVIDENCE MATRIX
# ============================================================

id_compact <- identity_rank %>%
  group_by(cluster) %>%
  summarise(
    identity_top5=paste0(
      programme,
      "[z=",
      sprintf("%.2f",programme_z),
      ";mean=",
      sprintf("%.2f",mean_programme),
      "]",
      collapse=" | "
    ),
    .groups="drop"
  )

state_compact <- state_rank %>%
  group_by(cluster) %>%
  summarise(
    states=paste0(
      programme,
      "[z=",
      sprintf("%.2f",programme_z),
      ";mean=",
      sprintf("%.2f",mean_programme),
      "]",
      collapse=" | "
    ),
    .groups="drop"
  )

evidence_matrix <- support %>%
  left_join(mk_summary,by="cluster") %>%
  left_join(id_compact,by="cluster") %>%
  left_join(state_compact,by="cluster") %>%
  mutate(cluster_num=as.integer(cluster)) %>%
  arrange(cluster_num) %>%
  select(-cluster_num)

write_csv(
  evidence_matrix,
  file.path(TDIR,"07E_L2_evidence_matrix.csv")
)

# ============================================================
# 11. PRINT CLUSTER-BY-CLUSTER REPORT
# ============================================================

cat("\n============================================================\n")
cat("L2 EVIDENCE MATRIX\n")
cat("============================================================\n")

for(i in seq_len(nrow(evidence_matrix))){

  z <- evidence_matrix[i,]

  cat("\n------------------------------------------------------------\n")
  cat(
    "CLUSTER",z$cluster,
    "| n =",z$n_cells,
    "| patients =",z$n_patients,
    "\n"
  )

  cat("MARKERS:",z$top_markers,"\n")
  cat("IDENTITY:",z$identity_top5,"\n")
  cat("STATES:",z$states,"\n")
}

# ============================================================
# 12. AUDIT
# ============================================================

audit <- c(
  "PIPELINE_V2 — 07E epithelial anatomical/state evidence",
  paste0("Generated: ",Sys.time()),
  "Input: epithelial_v2_07B_reclustered.rds",
  "Working resolution: epi_snn_res.0.5",
  "Cells: 63893",
  "Clusters: 19",
  "Identity and state programmes scored separately",
  "Identity framework reuses 06B3 anatomical programmes",
  "States include Cycling, epithelial stress, IFN response and MHC-II/AP",
  "Disease not used to assign identities",
  "No L2 labels assigned",
  "No resolution frozen"
)

writeLines(
  audit,
  file.path(TDIR,"07E_audit.txt")
)

cat("\n============================================================\n")
cat("07E COMPLETE\n")
cat("============================================================\n")
cat("No L2 labels assigned.\n")
cat("No resolution frozen.\n")
cat("PASS — anatomical/state evidence matrix generated.\n")
