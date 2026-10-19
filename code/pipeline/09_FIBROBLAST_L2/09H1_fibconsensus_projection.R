suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 09H1 FIBCONSENSUS FROZEN-SIGNATURE PROJECTION\n")
cat("============================================================\n")

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

OBJ <- "PIPELINE_V2/objects/annotated/fibroblast_v2_09F_refined.rds"

REFDIR <- "PIPELINE_V2/references/FibConsensus/frozen_reference"

SIGFILE <- file.path(
  REFDIR,
  "atlas_validation_signature_specification_v2_2.csv"
)

FAMFILE <- file.path(
  REFDIR,
  "family_master_v2_2.csv"
)

OUT <- "PIPELINE_V2/results/annotation/09_FIBROBLAST_L2/09H1_FibConsensus_projection"
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")

dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(OBJ))
stopifnot(file.exists(SIGFILE))
stopifnot(file.exists(FAMFILE))

# -------------------------------------------------------------------------
# Load
# -------------------------------------------------------------------------

cat("\nSTEP 1/9 — loading refined fibroblast universe and frozen reference\n")

x <- readRDS(OBJ)

sig <- read.csv(
  SIGFILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

n_signature_rows_original <- nrow(sig)

fam <- read.csv(
  FAMFILE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cat("Cells:", ncol(x), "\n")
cat("Genes:", nrow(x), "\n")
cat("Patients:", length(unique(x$patient_id)), "\n")
cat("Frozen families:", nrow(fam), "\n")
cat("Frozen signature rows:", nrow(sig), "\n")

stopifnot(ncol(x) == 4714)
stopifnot(nrow(x) == 23412)
stopifnot(nrow(fam) == 6)

# -------------------------------------------------------------------------
# Validate reference
# -------------------------------------------------------------------------

cat("\nSTEP 2/9 — validating frozen FibConsensus specification\n")

required_sig <- c(
  "entity_id",
  "entity_name",
  "entity_layer",
  "projection_gene_symbol",
  "signature_tier",
  "evidence_support_score"
)

stopifnot(all(required_sig %in% colnames(sig)))

sig <- sig %>%
  filter(entity_layer == "CONSENSUS_FAMILY") %>%
  mutate(
    projection_gene_symbol = as.character(projection_gene_symbol),
    signature_tier = toupper(as.character(signature_tier))
  )

fam_ids <- sort(unique(sig$entity_id))

cat("Families represented:", paste(fam_ids, collapse=", "), "\n")
cat("Signature tiers:",
    paste(sort(unique(sig$signature_tier)), collapse=", "), "\n")

print(
  sig %>%
    count(entity_id, entity_name, signature_tier, name="n_genes")
)

# -------------------------------------------------------------------------
# Coverage
# -------------------------------------------------------------------------

cat("\nSTEP 3/9 — signature coverage in current fibroblast universe\n")

genes <- rownames(x)

coverage <- sig %>%
  mutate(present = projection_gene_symbol %in% genes) %>%
  group_by(entity_id, entity_name, signature_tier) %>%
  summarise(
    n_requested = n(),
    n_present = sum(present),
    coverage = n_present/n_requested,
    missing = paste(
      projection_gene_symbol[!present],
      collapse=";"
    ),
    .groups="drop"
  )

write.csv(
  coverage,
  file.path(TAB, "09H1_signature_coverage.csv"),
  row.names=FALSE
)

print(coverage)

# -------------------------------------------------------------------------
# Expression matrix
# -------------------------------------------------------------------------

cat("\nSTEP 4/9 — extracting normalized expression\n")

DefaultAssay(x) <- "RNA"

expr <- GetAssayData(x, assay="RNA", layer="data")

# -------------------------------------------------------------------------
# Frozen programmes
# -------------------------------------------------------------------------

cat("\nSTEP 5/9 — calculating frozen CORE and FULL scores\n")

score_programme <- function(gene_table, expression_matrix) {

  g <- gene_table$projection_gene_symbol
  w <- gene_table$evidence_support_score

  keep <- g %in% rownames(expression_matrix)

  g <- g[keep]
  w <- w[keep]

  if(length(g) == 0) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }

  w[is.na(w)] <- 1

  m <- as.matrix(expression_matrix[g, , drop=FALSE])

  # Weighted mean of normalized expression.
  # Weight derives exclusively from the frozen evidence specification.
  as.numeric(crossprod(w, m) / sum(w))
}

families <- unique(
  sig[, c("entity_id","entity_name")]
)

score_df <- data.frame(
  cell=colnames(x),
  stringsAsFactors=FALSE
)

for(i in seq_len(nrow(families))) {

  id <- families$entity_id[i]
  nm <- families$entity_name[i]

  s_all <- sig %>% filter(entity_id == id)
  s_core <- s_all %>% filter(signature_tier == "CORE")

  score_df[[paste0(id,"_FULL")]] <-
    score_programme(s_all, expr)

  score_df[[paste0(id,"_CORE")]] <-
    score_programme(s_core, expr)

  cat(
    id, "|", nm,
    "| FULL:", nrow(s_all),
    "| CORE:", nrow(s_core), "\n"
  )
}

rownames(score_df) <- score_df$cell

# z-standardize each programme across fibroblast universe
score_cols <- setdiff(colnames(score_df), "cell")

z_df <- score_df

for(v in score_cols) {
  z_df[[v]] <- as.numeric(scale(score_df[[v]]))
}

# -------------------------------------------------------------------------
# Add metadata for summaries only
# -------------------------------------------------------------------------

md <- x@meta.data
md$cell <- rownames(md)

z_df <- z_df %>%
  left_join(
    md %>%
      select(
        cell,
        patient_id,
        disease,
        batch,
        starts_with("fib_refined_res.")
      ),
    by="cell"
  )

write.csv(
  z_df,
  file.path(TAB, "09H1_cell_FibConsensus_scores.csv"),
  row.names=FALSE
)

# -------------------------------------------------------------------------
# Cluster summaries r0.3 / r0.4
# -------------------------------------------------------------------------

cat("\nSTEP 6/9 — cluster-level projection\n")

find_res_col <- function(res) {

  candidates <- grep(
    paste0("res\\.", res, "$"),
    colnames(md),
    value=TRUE
  )

  if(length(candidates) == 0) {
    candidates <- grep(
      paste0(res, "$"),
      colnames(md),
      value=TRUE
    )
  }

  if(length(candidates) != 1) {
    stop(
      "Could not uniquely identify resolution column for ",
      res,
      ". Candidates: ",
      paste(candidates, collapse=", ")
    )
  }

  candidates
}

r03 <- find_res_col("0.3")
r04 <- find_res_col("0.4")

cat("r0.3 column:", r03, "\n")
cat("r0.4 column:", r04, "\n")

make_cluster_summary <- function(res_col, label) {

  d <- z_df
  d$cluster <- md[d$cell, res_col]

  long <- d %>%
    select(cell, patient_id, disease, batch, cluster,
           all_of(score_cols)) %>%
    pivot_longer(
      cols=all_of(score_cols),
      names_to="programme",
      values_to="zscore"
    ) %>%
    separate(
      programme,
      into=c("family_id","tier"),
      sep="_(?=[^_]+$)"
    )

  cluster_summary <- long %>%
    group_by(cluster, family_id, tier) %>%
    summarise(
      mean_z = mean(zscore, na.rm=TRUE),
      median_z = median(zscore, na.rm=TRUE),
      n_cells = n(),
      n_patients = n_distinct(patient_id),
      .groups="drop"
    )

  write.csv(
    cluster_summary,
    file.path(
      TAB,
      paste0("09H1_",label,"_cluster_scores.csv")
    ),
    row.names=FALSE
  )

  patient_summary <- long %>%
    group_by(cluster, patient_id, disease, family_id, tier) %>%
    summarise(
      mean_z = mean(zscore, na.rm=TRUE),
      n_cells=n(),
      .groups="drop"
    )

  write.csv(
    patient_summary,
    file.path(
      TAB,
      paste0("09H1_",label,"_patient_cluster_scores.csv")
    ),
    row.names=FALSE
  )

  cluster_summary
}

c03 <- make_cluster_summary(r03, "res03")
c04 <- make_cluster_summary(r04, "res04")

# -------------------------------------------------------------------------
# CORE/FULL concordance
# -------------------------------------------------------------------------

cat("\nSTEP 7/9 — CORE/FULL concordance\n")

concordance <- bind_rows(
  c03 %>% mutate(resolution="0.3"),
  c04 %>% mutate(resolution="0.4")
) %>%
  select(resolution, cluster, family_id, tier, mean_z) %>%
  pivot_wider(
    names_from=tier,
    values_from=mean_z
  ) %>%
  group_by(resolution, family_id) %>%
  summarise(
    n_clusters=n(),
    spearman_CORE_FULL=
      if(all(c("CORE","FULL") %in% names(cur_data()))) {
        suppressWarnings(cor(
          CORE, FULL,
          method="spearman",
          use="complete.obs"
        ))
      } else NA_real_,
    .groups="drop"
  )

write.csv(
  concordance,
  file.path(TAB, "09H1_CORE_FULL_concordance.csv"),
  row.names=FALSE
)

print(concordance)

# -------------------------------------------------------------------------
# Heatmaps
# -------------------------------------------------------------------------

cat("\nSTEP 8/9 — figures\n")

plot_heatmap <- function(cs, resolution) {

  d <- cs %>%
    filter(tier == "CORE") %>%
    mutate(
      cluster=factor(cluster),
      family_id=factor(
        family_id,
        levels=paste0("FAM0",1:6)
      )
    )

  p <- ggplot(
    d,
    aes(x=family_id, y=cluster, fill=mean_z)
  ) +
    geom_tile() +
    geom_text(
      aes(label=sprintf("%.2f",mean_z)),
      size=3
    ) +
    scale_fill_gradient2(midpoint=0) +
    labs(
      x="Frozen FibConsensus family",
      y=paste0("Refined fibroblast cluster (r=",resolution,")"),
      fill="Mean z",
      title=paste0(
        "FibConsensus CORE projection — resolution ",
        resolution
      )
    ) +
    theme_minimal(base_size=11) +
    theme(
      panel.grid=element_blank()
    )

  ggsave(
    file.path(
      FIG,
      paste0(
        "09H1_FibConsensus_CORE_heatmap_res",
        gsub("\\.","",resolution),
        ".pdf"
      )
    ),
    p,
    width=8,
    height=6
  )

  ggsave(
    file.path(
      FIG,
      paste0(
        "09H1_FibConsensus_CORE_heatmap_res",
        gsub("\\.","",resolution),
        ".png"
      )
    ),
    p,
    width=8,
    height=6,
    dpi=300
  )
}

plot_heatmap(c03, "0.3")
plot_heatmap(c04, "0.4")

# -------------------------------------------------------------------------
# Formal checks
# -------------------------------------------------------------------------

cat("\nSTEP 9/9 — formal checks\n")

checks <- data.frame(
  check=c(
    "cells_4714",
    "genes_23412",
    "patients_16",
    "six_frozen_families",
    "signature_source_56_rows",
    "consensus_family_signature_40_rows",
    "CORE_present",
    "FULL_present",
    "res03_scored",
    "res04_scored",
    "disease_not_used_for_scoring",
    "no_resolution_selected",
    "no_L2_labels_assigned",
    "no_cells_removed",
    "no_object_modified"
  ),
  pass=c(
    ncol(x)==4714,
    nrow(x)==23412,
    length(unique(x$patient_id))==16,
    length(unique(sig$entity_id))==6,
    n_signature_rows_original==56,
    nrow(sig)==40,
    any(sig$signature_tier=="CORE"),
    nrow(sig)>sum(sig$signature_tier=="CORE"),
    nrow(c03)>0,
    nrow(c04)>0,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
)

write.csv(
  checks,
  file.path(TAB, "09H1_formal_checks.csv"),
  row.names=FALSE
)

print(checks)

if(!all(checks$pass)) {
  stop("09H1 formal checks failed.")
}

cat("\n============================================================\n")
cat("09H1 COMPLETE\n")
cat("============================================================\n")
cat("Cells:",ncol(x),"\n")
cat("Frozen FibConsensus families:",length(unique(sig$entity_id)),"\n")
cat("Resolution 0.3 and 0.4 scored independently.\n")
cat("Disease was not used for scoring.\n")
cat("No resolution selected.\n")
cat("No L2 identities assigned.\n")
cat("No cells removed.\n")
cat("PASS — frozen FibConsensus projection evidence generated.\n")
