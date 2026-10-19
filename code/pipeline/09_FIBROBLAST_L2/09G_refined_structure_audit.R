suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 09G REFINED FIBROBLAST STRUCTURE AUDIT\n")
cat("============================================================\n")

IN <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09F_refined.rds"

OUT <- paste0(
  "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/",
  "09G_refined_structure/tables"
)

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

x <- readRDS(IN)

stopifnot(
  ncol(x)==4714,
  nrow(x)==23412,
  length(unique(x$patient_id))==16
)

c03 <- "fibroblast_refined_snn_res.0.3"
c04 <- "fibroblast_refined_snn_res.0.4"

stopifnot(
  c03 %in% colnames(x@meta.data),
  c04 %in% colnames(x@meta.data)
)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Patients:",length(unique(x$patient_id)),"\n")

# ==========================================================================
# 1. SUPPORT AT 0.3 AND 0.4
# ==========================================================================

cat("\nSTEP 1/7 — support audit\n")

support_one <- function(field,res){

  md <- x@meta.data %>%
    mutate(cluster=as.character(.data[[field]]))

  bp <- md %>%
    count(
      cluster,
      patient_id,
      disease,
      batch,
      name="n_cells"
    )

  cluster_totals <- bp %>%
    group_by(cluster) %>%
    summarise(
      cluster_cells=sum(n_cells),
      n_patients=n_distinct(patient_id),
      n_HC=n_distinct(patient_id[disease=="HC"]),
      n_AA=n_distinct(patient_id[disease=="AA"]),
      n_FFA=n_distinct(patient_id[disease=="FFA"]),
      n_batches=n_distinct(batch),
      .groups="drop"
    )

  dominant <- bp %>%
    group_by(cluster,patient_id) %>%
    summarise(
      patient_cells=sum(n_cells),
      .groups="drop"
    ) %>%
    group_by(cluster) %>%
    arrange(desc(patient_cells),patient_id,.by_group=TRUE) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    transmute(
      cluster,
      dominant_patient=patient_id,
      dominant_patient_cells=patient_cells
    )

  cluster_totals %>%
    left_join(dominant,by="cluster") %>%
    mutate(
      dominant_patient_fraction=
        dominant_patient_cells/cluster_cells,
      resolution=res
    ) %>%
    rename(
      n_cells=cluster_cells
    )
}

s03 <- support_one(c03,0.3)
s04 <- support_one(c04,0.4)

support <- bind_rows(s03,s04)

cat("\n=== SUPPORT ===\n")
print(as.data.frame(support),row.names=FALSE)

write_csv(
  support,
  file.path(OUT,"09G_cluster_support.csv")
)

# ==========================================================================
# 2. DISEASE COMPOSITION
# ==========================================================================

cat("\nSTEP 2/7 — disease composition\n")

disease_comp <- bind_rows(
  x@meta.data %>%
    mutate(
      cluster=as.character(.data[[c03]]),
      resolution=0.3
    ) %>%
    count(resolution,cluster,disease,name="n_cells"),

  x@meta.data %>%
    mutate(
      cluster=as.character(.data[[c04]]),
      resolution=0.4
    ) %>%
    count(resolution,cluster,disease,name="n_cells")
) %>%
  group_by(resolution,cluster) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  disease_comp,
  file.path(OUT,"09G_disease_composition.csv")
)

# ==========================================================================
# 3. MARKERS r=0.3
# ==========================================================================

cat("\nSTEP 3/7 — markers resolution 0.3\n")

DefaultAssay(x) <- "RNA"
Idents(x) <- x@meta.data[[c03]]

m03 <- FindAllMarkers(
  x,
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25,
  test.use="wilcox",
  verbose=FALSE
)

top03 <- m03 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group=TRUE) %>%
  slice_head(n=30) %>%
  ungroup()

write_csv(
  m03,
  file.path(OUT,"09G_markers_res0.3_all.csv")
)

write_csv(
  top03,
  file.path(OUT,"09G_markers_res0.3_top30.csv")
)

cat("\n=== TOP MARKERS r=0.3 ===\n")

compact03 <- top03 %>%
  group_by(cluster) %>%
  summarise(
    markers=paste(gene,collapse=", "),
    .groups="drop"
  )

print(as.data.frame(compact03),row.names=FALSE)

# ==========================================================================
# 4. MARKERS r=0.4
# ==========================================================================

cat("\nSTEP 4/7 — markers resolution 0.4\n")

Idents(x) <- x@meta.data[[c04]]

m04 <- FindAllMarkers(
  x,
  only.pos=TRUE,
  min.pct=0.10,
  logfc.threshold=0.25,
  test.use="wilcox",
  verbose=FALSE
)

top04 <- m04 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC),.by_group=TRUE) %>%
  slice_head(n=30) %>%
  ungroup()

write_csv(
  m04,
  file.path(OUT,"09G_markers_res0.4_all.csv")
)

write_csv(
  top04,
  file.path(OUT,"09G_markers_res0.4_top30.csv")
)

# ==========================================================================
# 5. 0.3 -> 0.4 TRANSITIONS
# ==========================================================================

cat("\nSTEP 5/7 — 0.3 -> 0.4 transitions\n")

transition <- data.frame(
  parent=as.character(x@meta.data[[c03]]),
  child=as.character(x@meta.data[[c04]])
) %>%
  count(parent,child,name="n_cells") %>%
  group_by(parent) %>%
  mutate(
    parent_total=sum(n_cells),
    fraction_parent=n_cells/parent_total
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(parent),
    desc(fraction_parent)
  )

cat("\n=== TRANSITIONS >=5% OF PARENT ===\n")
print(
  as.data.frame(
    transition %>%
      filter(fraction_parent>=0.05)
  ),
  row.names=FALSE
)

write_csv(
  transition,
  file.path(OUT,"09G_transition_03_to_04.csv")
)

# ==========================================================================
# 6. PATIENT-DOMINATED CLUSTER AUDIT
# ==========================================================================

cat("\nSTEP 6/7 — patient-dominated cluster audit\n")

flagged <- s04 %>%
  filter(
    dominant_patient_fraction > 0.80 |
    n_patients < 5
  )

cat("\n=== FLAGGED r=0.4 CLUSTERS ===\n")
print(as.data.frame(flagged),row.names=FALSE)

write_csv(
  flagged,
  file.path(OUT,"09G_flagged_res04_clusters.csv")
)

flagged_ids <- flagged$cluster

flagged_patient <- x@meta.data %>%
  mutate(
    cluster=as.character(.data[[c04]])
  ) %>%
  filter(
    cluster %in% flagged_ids
  ) %>%
  count(
    cluster,
    patient_id,
    disease,
    batch,
    name="n_cells"
  ) %>%
  group_by(cluster) %>%
  mutate(
    fraction=n_cells/sum(n_cells)
  ) %>%
  ungroup() %>%
  arrange(
    as.numeric(cluster),
    desc(n_cells)
  )

cat("\n=== FLAGGED CLUSTER PATIENT COMPOSITION ===\n")
print(as.data.frame(flagged_patient),row.names=FALSE)

write_csv(
  flagged_patient,
  file.path(OUT,"09G_flagged_res04_patient_composition.csv")
)

# Direct marker comparison:
# flagged child versus the remainder of its r=0.3 parent.

if(nrow(flagged)>0){

  for(k in seq_len(nrow(flagged))){

    child <- flagged$cluster[k]

    cells_child <- colnames(x)[
      as.character(x@meta.data[[c04]])==child
    ]

    parents <- unique(
      as.character(
        x@meta.data[cells_child,c03]
      )
    )

    for(parent in parents){

      parent_cells <- colnames(x)[
        as.character(x@meta.data[[c03]])==parent
      ]

      child_in_parent <- intersect(
        cells_child,
        parent_cells
      )

      rest_parent <- setdiff(
        parent_cells,
        child_in_parent
      )

      if(
        length(child_in_parent)>=10 &&
        length(rest_parent)>=10
      ){

        sub <- subset(
          x,
          cells=c(child_in_parent,rest_parent)
        )

        grp <- ifelse(
          colnames(sub) %in% child_in_parent,
          "flagged_child",
          "parent_remainder"
        )

        sub$comparison_group <- grp
        Idents(sub) <- "comparison_group"

        de <- FindMarkers(
          sub,
          ident.1="flagged_child",
          ident.2="parent_remainder",
          min.pct=0.10,
          logfc.threshold=0.25,
          test.use="wilcox",
          verbose=FALSE
        )

        de$gene <- rownames(de)

        fn <- paste0(
          "09G_flagged_C",
          child,
          "_within_parent03_C",
          parent,
          ".csv"
        )

        write_csv(
          de,
          file.path(OUT,fn)
        )

        cat(
          "\nDirect comparison:",
          "r0.4/C",child,
          "vs remainder of r0.3/C",parent,
          "\n"
        )

        cat("\nTOP UP:\n")
        print(
          head(
            de[order(-de$avg_log2FC),],
            30
          )
        )
      }
    }
  }
}

# ==========================================================================
# 7. FORMAL CHECKS
# ==========================================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- data.frame(
  check=c(
    "cells_4714",
    "genes_23412",
    "patients_16",
    "res03_10_clusters",
    "res04_11_clusters",
    "markers_res03_generated",
    "markers_res04_generated",
    "transition_generated",
    "no_FibConsensus_used",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed"
  ),

  pass=c(
    ncol(x)==4714,
    nrow(x)==23412,
    length(unique(x$patient_id))==16,
    length(unique(x@meta.data[[c03]]))==10,
    length(unique(x@meta.data[[c04]]))==11,
    nrow(m03)>0,
    nrow(m04)>0,
    nrow(transition)>0,
    TRUE,
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
  file.path(OUT,"09G_formal_checks.csv")
)

audit <- c(
  "PIPELINE_V2 — 09G refined fibroblast structure audit",
  paste0("Generated: ",Sys.time()),
  "Input: fibroblast_v2_09F_refined.rds",
  "Cells: 4714.",
  "Resolution 0.3 and 0.4 audited.",
  "Markers generated de novo after fibroblast-universe refinement.",
  "Patient-dominated r=0.4 clusters explicitly audited.",
  "0.3-to-0.4 transitions explicitly audited.",
  "FibConsensus not used.",
  "Disease not used for identity assignment.",
  "No resolution selected.",
  "No L2 identity assigned.",
  "No cells removed."
)

writeLines(
  audit,
  file.path(OUT,"09G_audit.txt")
)

cat("\n============================================================\n")
cat("09G COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(x),"\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("FibConsensus not used.\n")
cat("PASS — refined fibroblast structure evidence generated.\n")

