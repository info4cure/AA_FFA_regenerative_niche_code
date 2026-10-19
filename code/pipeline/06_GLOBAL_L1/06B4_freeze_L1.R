suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

cat("\n============================================================\n")
cat("PIPELINE_V2 — 06B4 FREEZE L1\n")
cat("============================================================\n\n")

V2 <- "PIPELINE_V2"
IN <- file.path(V2,"objects","ambient_corrected",
                "global_atlas_v2_SoupX_expression.rds")
REF <- file.path(V2,"objects","frozen",
                 "global_atlas_v2_frozen.rds")
OUTDIR <- file.path(V2,"results","annotation",
                    "06_GLOBAL_L1","06B4_L1_frozen")
TDIR <- file.path(OUTDIR,"tables")
OBJDIR <- file.path(V2,"objects","annotated")
OUT <- file.path(OBJDIR,"global_atlas_v2_L1_frozen.rds")
CC <- "RNA_snn_res.0.6"

dir.create(TDIR,recursive=TRUE,showWarnings=FALSE)
dir.create(OBJDIR,recursive=TRUE,showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Load + canonical invariance
# ------------------------------------------------------------

cat("STEP 1/6 — canonical object verification\n")

x <- readRDS(IN)
ref <- readRDS(REF)
DefaultAssay(x) <- "RNA"

stopifnot(
  ncol(x)==82893,
  nrow(x)==23412,
  ncol(ref)==82893,
  CC %in% colnames(x@meta.data),
  identical(colnames(x),colnames(ref))
)

cl <- as.character(x@meta.data[[CC]])

stopifnot(
  identical(sort(unique(as.integer(cl))),0:26)
)

cat("Cells:",ncol(x),"\n")
cat("Genes:",nrow(x),"\n")
cat("Cell identity/order: EXACT\n\n")

# ------------------------------------------------------------
# 2. Frozen broad L1 map
# ------------------------------------------------------------

cat("STEP 2/6 — L1 adjudication map\n")

map <- data.frame(
  cluster=as.character(0:26),
  L1=c(
    "Epithelial",      # 0
    "Epithelial",      # 1
    "Epithelial",      # 2
    "Epithelial",      # 3
    "Fibroblast",      # 4
    "Epithelial",      # 5
    "Epithelial",      # 6
    "Epithelial",      # 7
    "Epithelial",      # 8
    "T_NK",            # 9
    "Epithelial",      # 10
    "Vascular_EC",     # 11
    "Epithelial",      # 12
    "Epithelial",      # 13
    "Epithelial",      # 14
    "Epithelial",      # 15
    "Epithelial",      # 16
    "Epithelial",      # 17
    "Melanocyte",      # 18
    "Myeloid_DC",      # 19
    "Sebocyte",        # 20
    "Pericyte_SMC",    # 21
    "Epithelial",      # 22
    "Sweat_gland",     # 23
    "Sweat_gland",     # 24
    "Smooth_muscle",   # 25
    "Lymphatic_EC"     # 26
  ),
  status=c(
    rep("CONFIDENT",25),
    "CONFIDENT_FLAGGED",
    "CONFIDENT"
  ),
  stringsAsFactors=FALSE
)

stopifnot(
  nrow(map)==27,
  identical(map$cluster,as.character(0:26)),
  !anyDuplicated(map$cluster)
)

write_csv(map,file.path(TDIR,"06B4_cluster_to_L1_map.csv"))

# ------------------------------------------------------------
# 3. Assign metadata
# ------------------------------------------------------------

cat("STEP 3/6 — assigning L1\n")

ii <- match(cl,map$cluster)
stopifnot(!anyNA(ii))

x$L1 <- map$L1[ii]
x$L1_status <- map$status[ii]

x$L1_flag <- NA_character_
x$L1_flag[cl=="25"] <-
  "cluster25_KRT_high_contractile_subpopulation_retained"

stopifnot(
  !anyNA(x$L1),
  length(x$L1)==82893
)

# ------------------------------------------------------------
# 4. Cell-count audit
# ------------------------------------------------------------

cat("STEP 4/6 — cell-count audit\n")

cluster_tab <- data.frame(
  cluster=as.character(0:26),
  n_cells=as.integer(table(factor(cl,levels=as.character(0:26))))
) %>%
  left_join(map,by="cluster")

l1_tab <- as.data.frame(table(L1=x$L1))
names(l1_tab)[2] <- "n_cells"
l1_tab <- l1_tab %>% arrange(desc(n_cells))

stopifnot(
  sum(cluster_tab$n_cells)==82893,
  sum(l1_tab$n_cells)==82893
)

write_csv(
  cluster_tab,
  file.path(TDIR,"06B4_cluster_L1_cell_counts.csv")
)
write_csv(
  l1_tab,
  file.path(TDIR,"06B4_L1_cell_counts.csv")
)

cat("\n=== L1 CELL COUNTS ===\n")
print(l1_tab,row.names=FALSE)
cat("\nTOTAL:",sum(l1_tab$n_cells),"\n\n")

# ------------------------------------------------------------
# 5. Freeze metadata + save NEW object
# ------------------------------------------------------------

cat("STEP 5/6 — saving new annotated object\n")

x$L1_frozen <- TRUE
x$annotation_version <- "PIPELINE_V2_06B4_L1_frozen"
x$annotation_expression_source <- "SoupX_corrected"

# Final invariance before save
stopifnot(
  identical(colnames(x),colnames(ref)),
  identical(as.character(x@meta.data[[CC]]),cl),
  ncol(x)==82893,
  nrow(x)==23412
)

saveRDS(x,OUT,compress=FALSE)

stopifnot(file.exists(OUT))

# Re-open saved object and verify
y <- readRDS(OUT)

stopifnot(
  ncol(y)==82893,
  nrow(y)==23412,
  identical(colnames(y),colnames(x)),
  identical(as.character(y$L1),as.character(x$L1)),
  identical(as.character(y@meta.data[[CC]]),cl),
  all(y$L1_frozen)
)

rm(y)

cat("Saved:",OUT,"\n\n")

# ------------------------------------------------------------
# 6. Audit
# ------------------------------------------------------------

cat("STEP 6/6 — audit trail\n")

audit <- c(
  "PIPELINE_V2 — 06B4 L1 FREEZE",
  paste0("Generated: ",Sys.time()),
  "Input: global_atlas_v2_SoupX_expression.rds",
  "Canonical reference: global_atlas_v2_frozen.rds",
  "Output: global_atlas_v2_L1_frozen.rds",
  "Cells: 82893",
  "Genes: 23412",
  "Canonical clusters: RNA_snn_res.0.6, 0-26",
  "Expression: SoupX-corrected",
  "L1 evidence: 06B1 + 06B2 + corrected 06B3",
  "Cycling and follicular/anatomical states deferred to L2",
  "Cluster 25 broad L1: Smooth_muscle",
  "Cluster 25 KRT-high subset retained and flagged",
  "No cells removed",
  "No expression recalculated",
  "No PCA/Harmony/neighbors/UMAP/clusters recalculated",
  "Original SoupX object not overwritten",
  "L1_frozen: TRUE"
)

writeLines(audit,file.path(TDIR,"06B4_audit.txt"))

cat("\n============================================================\n")
cat("06B4 COMPLETE\n")
cat("============================================================\n")
cat("Output:",OUT,"\n")
cat("Cells:",ncol(x),"\n")
cat("L1 classes:",nrow(l1_tab),"\n")
cat("PASS — L1 annotation frozen; canonical atlas preserved.\n")
