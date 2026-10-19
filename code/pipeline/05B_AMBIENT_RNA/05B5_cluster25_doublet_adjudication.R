suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})
options(stringsAsFactors=FALSE)
set.seed(20260924)
source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 05B5 CLUSTER 25 DOUBLET ADJUDICATION\n")
cat("============================================================\n\n")

PRE_FILE <- file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
POST_FILE <- file.path(V2_ROOT,"objects","ambient_corrected","global_atlas_v2_SoupX_expression.rds")
OUT <- file.path(V2_ROOT,"results","QC","ambient_RNA","05B5_cluster25_doublet_adjudication")
TD <- file.path(OUT,"tables")
FD <- file.path(OUT,"figures")
dir.create(TD,recursive=TRUE,showWarnings=FALSE)
dir.create(FD,recursive=TRUE,showWarnings=FALSE)

pre <- readRDS(PRE_FILE)
post <- readRDS(POST_FILE)
DefaultAssay(pre) <- "RNA"
DefaultAssay(post) <- "RNA"

stopifnot(
  identical(colnames(pre),colnames(post)),
  identical(rownames(pre),rownames(post)),
  identical(as.character(pre$RNA_snn_res.0.6),as.character(post$RNA_snn_res.0.6))
)

cl <- as.character(post$RNA_snn_res.0.6)
cells25 <- colnames(post)[cl=="25"]
stopifnot(length(cells25)==330)

KRT_PANEL <- intersect(c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10"),rownames(post))
counts <- LayerData(post,assay="RNA",layer="counts")
dat <- LayerData(post,assay="RNA",layer="data")
krt_counts <- counts[KRT_PANEL,cells25,drop=FALSE]
krt_data <- dat[KRT_PANEL,cells25,drop=FALSE]
nKRT <- Matrix::colSums(krt_counts>0)
KRT_score <- Matrix::colMeans(krt_data)
KRT_high <- nKRT>=3
stopifnot(sum(KRT_high)==147)

meta <- data.frame(
  cell_id=cells25,KRT_high=KRT_high,nKRT=nKRT,KRT_score=KRT_score,
  patient_id=as.character(post$patient_id[cells25]),
  disease=as.character(post$disease[cells25]),
  batch=as.character(post$batch[cells25]),
  doublet_score=as.numeric(pre$doublet_score[cells25]),
  nFeature_RNA=as.numeric(pre$nFeature_RNA[cells25]),
  nCount_RNA=as.numeric(pre$nCount_RNA[cells25]),
  stringsAsFactors=FALSE
)
write.csv(meta,file.path(TD,"05B5_cluster25_cell_metadata.csv"),row.names=FALSE)

cat("STEP 1/6 — full-transcriptome exploratory DE\n")
sub <- subset(post,cells=cells25)
Idents(sub) <- factor(ifelse(KRT_high,"KRT_high","KRT_low"),levels=c("KRT_low","KRT_high"))
de <- FindMarkers(sub,ident.1="KRT_high",ident.2="KRT_low",assay="RNA",slot="data",
                  test.use="wilcox",logfc.threshold=0,min.pct=0,only.pos=FALSE,verbose=FALSE)
de$gene <- rownames(de); rownames(de) <- NULL
de$is_defining_KRT <- de$gene %in% KRT_PANEL
de <- de[,c("gene","is_defining_KRT",setdiff(colnames(de),c("gene","is_defining_KRT")))]
write.csv(de,file.path(TD,"05B5_celllevel_DE_all_genes.csv"),row.names=FALSE)

de_nondef <- de[!de$is_defining_KRT,]
top_up <- head(de_nondef[order(de_nondef$avg_log2FC,decreasing=TRUE),],100)
top_dn <- head(de_nondef[order(de_nondef$avg_log2FC,decreasing=FALSE),],100)
write.csv(top_up,file.path(TD,"05B5_top100_up_KRThigh_excluding_definingKRT.csv"),row.names=FALSE)
write.csv(top_dn,file.path(TD,"05B5_top100_down_KRThigh_excluding_definingKRT.csv"),row.names=FALSE)

cat("\nTop 20 UP excluding defining KRT panel:\n")
print(top_up[seq_len(min(20,nrow(top_up))),c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")],row.names=FALSE)

cat("\nSTEP 2/6 — doublet and complexity associations\n")
scor <- function(x,y) suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE))
ct_ds <- scor(meta$KRT_score,meta$doublet_score)
ct_nf <- scor(meta$KRT_score,meta$nFeature_RNA)
ct_nc <- scor(meta$KRT_score,meta$nCount_RNA)
assoc <- data.frame(
  comparison=c("KRT_score_vs_doublet_score","KRT_score_vs_nFeature_RNA","KRT_score_vs_nCount_RNA"),
  spearman_rho=c(unname(ct_ds$estimate),unname(ct_nf$estimate),unname(ct_nc$estimate)),
  p_value=c(ct_ds$p.value,ct_nf$p.value,ct_nc$p.value)
)
assoc$p_adj_BH <- p.adjust(assoc$p_value,"BH")
write.csv(assoc,file.path(TD,"05B5_continuous_associations.csv"),row.names=FALSE)

grp <- do.call(rbind,lapply(split(meta,meta$KRT_high),function(z)data.frame(
  KRT_high=z$KRT_high[1],n_cells=nrow(z),
  median_doublet=median(z$doublet_score,na.rm=TRUE),
  median_nFeature=median(z$nFeature_RNA,na.rm=TRUE),
  median_nCount=median(z$nCount_RNA,na.rm=TRUE),
  median_KRT_score=median(z$KRT_score,na.rm=TRUE)
)))
write.csv(grp,file.path(TD,"05B5_KRThigh_complexity_summary.csv"),row.names=FALSE)

cat("\nSTEP 3/6 — patient-level reproducibility\n")
patients <- sort(unique(meta$patient_id))
patient_effects <- do.call(rbind,lapply(patients,function(pt){
  ii <- meta$patient_id==pt; hi <- ii & meta$KRT_high; lo <- ii & !meta$KRT_high
  data.frame(
    patient_id=pt,disease=unique(meta$disease[ii]),batch=unique(meta$batch[ii]),
    n_cluster25=sum(ii),n_KRThigh=sum(hi),n_KRTlow=sum(lo),pct_KRThigh=100*mean(meta$KRT_high[ii]),
    median_doublet_high=if(sum(hi)>0) median(meta$doublet_score[hi],na.rm=TRUE) else NA,
    median_doublet_low=if(sum(lo)>0) median(meta$doublet_score[lo],na.rm=TRUE) else NA,
    delta_median_doublet=if(sum(hi)>0&&sum(lo)>0) median(meta$doublet_score[hi],na.rm=TRUE)-median(meta$doublet_score[lo],na.rm=TRUE) else NA,
    median_nFeature_high=if(sum(hi)>0) median(meta$nFeature_RNA[hi],na.rm=TRUE) else NA,
    median_nFeature_low=if(sum(lo)>0) median(meta$nFeature_RNA[lo],na.rm=TRUE) else NA
  )
}))
write.csv(patient_effects,file.path(TD,"05B5_patient_effects.csv"),row.names=FALSE)
eligible <- patient_effects[is.finite(patient_effects$delta_median_doublet),]
direction_summary <- data.frame(
  patients_with_both_states=nrow(eligible),
  patients_doublet_higher_in_KRThigh=sum(eligible$delta_median_doublet>0),
  patients_doublet_lower_in_KRThigh=sum(eligible$delta_median_doublet<0),
  patients_equal=sum(eligible$delta_median_doublet==0),
  median_within_patient_delta=median(eligible$delta_median_doublet,na.rm=TRUE)
)
write.csv(direction_summary,file.path(TD,"05B5_patient_direction_summary.csv"),row.names=FALSE)

cat("\nSTEP 4/6 — individual marker audit\n")
MARKERS <- intersect(c(
  "KRT5","KRT14","KRT15","KRT17","KRT1","KRT10","EPCAM","KRT8","KRT18","KRT19",
  "ACTG2","TAGLN","CNN1","MYH11","MYL9","ACTA2","RGS5","CSPG4","PDGFRB","NOTCH3","ABCC9",
  "DES","ACTA1","ACTC1","PCP4"
),rownames(post))
marker_tab <- do.call(rbind,lapply(MARKERS,function(g){
  v <- dat[g,cells25]; cts <- counts[g,cells25]
  data.frame(gene=g,pct_KRTlow=100*mean(cts[!KRT_high]>0),pct_KRThigh=100*mean(cts[KRT_high]>0),
             mean_KRTlow=mean(v[!KRT_high]),mean_KRThigh=mean(v[KRT_high]),
             delta_mean=mean(v[KRT_high])-mean(v[!KRT_high]))
}))
write.csv(marker_tab,file.path(TD,"05B5_marker_audit.csv"),row.names=FALSE)

cat("\nSTEP 5/6 — canonical UMAP neighborhood audit\n")
emb <- Embeddings(post,"umap_harmony")[cells25,,drop=FALSE]
stopifnot(identical(rownames(emb),cells25))
k <- 15L
dm <- as.matrix(dist(emb)); diag(dm) <- Inf
nn <- t(apply(dm,1,function(v) order(v)[seq_len(k)]))
local_frac_high <- apply(nn,1,function(ix) mean(KRT_high[ix]))
umap_tab <- cbind(meta,UMAP1=emb[,1],UMAP2=emb[,2],local_KRThigh_fraction=local_frac_high)
write.csv(umap_tab,file.path(TD,"05B5_cluster25_UMAP_localization.csv"),row.names=FALSE)
local_summary <- do.call(rbind,lapply(split(umap_tab,umap_tab$KRT_high),function(z)data.frame(
  KRT_high=z$KRT_high[1],n_cells=nrow(z),
  mean_local_KRThigh_fraction=mean(z$local_KRThigh_fraction),
  median_local_KRThigh_fraction=median(z$local_KRThigh_fraction)
)))
write.csv(local_summary,file.path(TD,"05B5_local_neighborhood_summary.csv"),row.names=FALSE)

cat("\nSTEP 6/6 — diagnostic figures\n")
p1 <- ggplot(meta,aes(KRT_score,doublet_score))+geom_point(size=1.5,alpha=.75)+
  geom_smooth(method="lm",se=FALSE,linewidth=.6)+theme_classic(base_size=11)+
  labs(x="SoupX-corrected KRT programme",y="scDblFinder score",title="Cluster 25: KRT programme vs doublet score")
ggsave(file.path(FD,"05B5_KRT_vs_doublet.pdf"),p1,width=6,height=5,bg="white")
ggsave(file.path(FD,"05B5_KRT_vs_doublet.png"),p1,width=6,height=5,dpi=600,bg="white")

p2 <- ggplot(umap_tab,aes(UMAP1,UMAP2,color=KRT_high))+geom_point(size=2,alpha=.9)+
  coord_equal()+theme_void()+scale_color_manual(values=c("FALSE"="grey75","TRUE"="black"))+
  labs(color="KRT-high",title="Cluster 25: residual KRT-high cells")
ggsave(file.path(FD,"05B5_cluster25_KRThigh_UMAP.pdf"),p2,width=6,height=5,bg="white")
ggsave(file.path(FD,"05B5_cluster25_KRThigh_UMAP.png"),p2,width=6,height=5,dpi=600,bg="white")

p3 <- ggplot(umap_tab,aes(UMAP1,UMAP2,color=doublet_score))+geom_point(size=2,alpha=.9)+
  coord_equal()+theme_void()+scale_color_viridis_c()+
  labs(color="Doublet score",title="Cluster 25: scDblFinder score")
ggsave(file.path(FD,"05B5_cluster25_doublet_UMAP.pdf"),p3,width=6,height=5,bg="white")
ggsave(file.path(FD,"05B5_cluster25_doublet_UMAP.png"),p3,width=6,height=5,dpi=600,bg="white")

topgenes <- unique(head(top_up$gene,30))
pp <- list()
if(length(topgenes)>0){
  for(pt in patients) for(state in c(FALSE,TRUE)){
    cc <- meta$cell_id[meta$patient_id==pt & meta$KRT_high==state]
    if(length(cc)>0){
      mm <- Matrix::rowMeans(dat[topgenes,cc,drop=FALSE])
      pp[[length(pp)+1]] <- data.frame(patient_id=pt,KRT_high=state,gene=names(mm),mean_expr=as.numeric(mm))
    }
  }
}
if(length(pp)>0) write.csv(do.call(rbind,pp),file.path(TD,"05B5_topDE_patient_state_means.csv"),row.names=FALSE)

evidence <- data.frame(
  metric=c("n_cluster25","n_KRThigh","pct_KRThigh","median_doublet_KRThigh","median_doublet_KRTlow",
           "spearman_KRT_doublet","patients_with_both_states","patients_doublet_higher_in_KRThigh",
           "median_local_KRThigh_fraction_KRThigh","median_local_KRThigh_fraction_KRTlow"),
  value=c(nrow(meta),sum(meta$KRT_high),100*mean(meta$KRT_high),
          median(meta$doublet_score[meta$KRT_high],na.rm=TRUE),
          median(meta$doublet_score[!meta$KRT_high],na.rm=TRUE),
          unname(ct_ds$estimate),
          direction_summary$patients_with_both_states,
          direction_summary$patients_doublet_higher_in_KRThigh,
          local_summary$median_local_KRThigh_fraction[local_summary$KRT_high],
          local_summary$median_local_KRThigh_fraction[!local_summary$KRT_high])
)
write.csv(evidence,file.path(TD,"05B5_evidence_summary.csv"),row.names=FALSE)

writeLines(c(
  "PIPELINE_V2 — 05B5 cluster 25 adjudication audit",
  paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("Cluster 25 cells: ",nrow(meta)),
  paste0("KRT-high cells: ",sum(meta$KRT_high)),
  "Expression: SoupX-corrected RNA counts/data.",
  "Cluster membership and UMAP: canonical frozen atlas.",
  "Cell-level DE is exploratory/diagnostic, not patient-level inference.",
  "Defining KRT genes are flagged and excluded from unbiased top-DE ranking.",
  "No cells removed and no atlas identities changed."
),file.path(TD,"05B5_audit.txt"))

cat("\n=== CONTINUOUS ASSOCIATIONS ===\n");print(assoc,row.names=FALSE)
cat("\n=== KRT-HIGH COMPLEXITY ===\n");print(grp,row.names=FALSE)
cat("\n=== PATIENT DIRECTION SUMMARY ===\n");print(direction_summary,row.names=FALSE)
cat("\n=== LOCAL NEIGHBORHOOD SUMMARY ===\n");print(local_summary,row.names=FALSE)
cat("\nPASS — adjudication audit only; canonical atlas unchanged.\n")
