# =============================================================================
# PIPELINE_V2 — 05B4_cluster25_dissection.R
# Focused dissection of global cluster 25 after SoupX correction.
# Tests SMC/mural vs epithelial/glandular coexpression and residual doublets.
# READ-ONLY: no cells are removed and no atlas identities are changed.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2);library(patchwork)})
options(stringsAsFactors=FALSE);set.seed(20260924);source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\nPIPELINE_V2 — 05B4 CLUSTER 25 DISSECTION\n============================================================\n\n")

PRE_FILE<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
POST_FILE<-file.path(V2_ROOT,"objects","ambient_corrected","global_atlas_v2_SoupX_expression.rds")
OUT<-file.path(V2_ROOT,"results","QC","ambient_RNA","05B4_cluster25")
TD<-file.path(OUT,"tables");FD<-file.path(OUT,"figures");dir.create(TD,recursive=TRUE,showWarnings=FALSE);dir.create(FD,recursive=TRUE,showWarnings=FALSE)

pre<-readRDS(PRE_FILE);post<-readRDS(POST_FILE);DefaultAssay(pre)<-"RNA";DefaultAssay(post)<-"RNA"
stopifnot(identical(colnames(pre),colnames(post)),identical(as.character(pre$RNA_snn_res.0.6),as.character(post$RNA_snn_res.0.6)))
cells25<-colnames(post)[as.character(post$RNA_snn_res.0.6)=="25"];stopifnot(length(cells25)==330)

programs<-list(
 KRT=c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10"),
 Epithelial=c("EPCAM","KRT8","KRT18","KRT19"),
 Glandular=c("AQP5","DCD","PIP","FOXA1","MUC7","FOXI1"),
 SMC=c("ACTG2","TAGLN","CNN1","MYH11","MYL9","ACTA2"),
 Pericyte=c("RGS5","CSPG4","PDGFRB","NOTCH3","ABCC9")
)
programs<-lapply(programs,function(g)intersect(g,rownames(post)))
write.csv(do.call(rbind,lapply(names(programs),function(n)data.frame(program=n,gene=programs[[n]]))),file.path(TD,"05B4_program_genes.csv"),row.names=FALSE)

C<-LayerData(post,assay="RNA",layer="counts")[unique(unlist(programs)),cells25,drop=FALSE]
D<-LayerData(post,assay="RNA",layer="data")[unique(unlist(programs)),cells25,drop=FALSE]
score<-sapply(programs,function(g)Matrix::colMeans(D[g,,drop=FALSE]))
ndet<-sapply(programs,function(g)Matrix::colSums(C[g,,drop=FALSE]>0))

tab<-data.frame(cell_id=cells25,patient_id=post$patient_id[cells25],disease=post$disease[cells25],batch=post$batch[cells25],doublet_score=pre$doublet_score[cells25],nFeature_RNA=pre$nFeature_RNA[cells25],nCount_RNA=pre$nCount_RNA[cells25],stringsAsFactors=FALSE)
for(n in colnames(score)){tab[[paste0(n,"_score")]]<-score[,n];tab[[paste0(n,"_n_detected")]]<-ndet[,n]}
tab$KRT_high<-tab$KRT_n_detected>=3
tab$SMC_high<-tab$SMC_n_detected>=3
tab$epi_contractile_coexpression<-tab$KRT_high & tab$SMC_high
write.csv(tab,file.path(TD,"05B4_cluster25_cells.csv"),row.names=FALSE)

summary_group<-do.call(rbind,lapply(split(tab,tab$KRT_high),function(z)data.frame(KRT_high=z$KRT_high[1],n_cells=nrow(z),median_doublet_score=median(z$doublet_score,na.rm=TRUE),median_nFeature=median(z$nFeature_RNA),median_nCount=median(z$nCount_RNA),mean_KRT=mean(z$KRT_score),mean_Epithelial=mean(z$Epithelial_score),mean_Glandular=mean(z$Glandular_score),mean_SMC=mean(z$SMC_score),mean_Pericyte=mean(z$Pericyte_score),pct_SMC_high=100*mean(z$SMC_high),pct_epi_contractile=100*mean(z$epi_contractile_coexpression))))
write.csv(summary_group,file.path(TD,"05B4_KRThigh_vs_low_summary.csv"),row.names=FALSE)

patient<-do.call(rbind,lapply(split(tab,tab$patient_id),function(z)data.frame(patient_id=z$patient_id[1],batch=z$batch[1],n_cells=nrow(z),pct_KRT_high=100*mean(z$KRT_high),pct_SMC_high=100*mean(z$SMC_high),pct_epi_contractile=100*mean(z$epi_contractile_coexpression),median_doublet=median(z$doublet_score,na.rm=TRUE))))
write.csv(patient,file.path(TD,"05B4_patient_distribution.csv"),row.names=FALSE)

# Marker expression table, useful for inspection without plotting.
genes<-unique(unlist(programs))
means<-do.call(rbind,lapply(genes,function(g)data.frame(gene=g,KRT_high_FALSE=mean(D[g,!tab$KRT_high]),KRT_high_TRUE=mean(D[g,tab$KRT_high]))))
write.csv(means,file.path(TD,"05B4_marker_means_KRThigh_vs_low.csv"),row.names=FALSE)

# Global UMAP location of cluster 25, with corrected expression / canonical geometry.
emb<-Embeddings(post,"umap_harmony");u<-data.frame(cell_id=rownames(emb),UMAP1=emb[,1],UMAP2=emb[,2],cluster25=rownames(emb)%in%cells25)
p0<-ggplot(u,aes(UMAP1,UMAP2))+geom_point(data=u[!u$cluster25,],size=.08,color="grey85")+geom_point(data=u[u$cluster25,],size=.8,color="black")+coord_equal()+theme_void()+ggtitle("Global location of cluster 25")
ggsave(file.path(FD,"05B4_cluster25_global_location.pdf"),p0,width=6,height=5,bg="white")

# Feature plots on canonical UMAP for cluster 25 only.
focus<-c("KRT14","KRT17","KRT19","EPCAM","AQP5","PIP","DCD","ACTG2","TAGLN","CNN1","MYH11","RGS5")
focus<-focus[focus%in%rownames(post)]
pfeat<-FeaturePlot(post,cells=cells25,features=focus,reduction="umap_harmony",ncol=4,order=TRUE,keep.scale="all")
ggsave(file.path(FD,"05B4_cluster25_featureplots_corrected.pdf"),pfeat,width=12,height=3*ceiling(length(focus)/4),bg="white",limitsize=FALSE)
ggsave(file.path(FD,"05B4_cluster25_featureplots_corrected.png"),pfeat,width=12,height=3*ceiling(length(focus)/4),dpi=400,bg="white",limitsize=FALSE)

# Key paired relationships.
p1<-ggplot(tab,aes(KRT_score,SMC_score,color=doublet_score))+geom_point(size=1.3,alpha=.8)+scale_color_viridis_c()+labs(x="Corrected KRT programme",y="Corrected SMC programme",color="Doublet score",title="Cluster 25: KRT vs contractile programme")+theme_classic()
ggsave(file.path(FD,"05B4_KRT_vs_SMC_doublet.pdf"),p1,width=6,height=5,bg="white");ggsave(file.path(FD,"05B4_KRT_vs_SMC_doublet.png"),p1,width=6,height=5,dpi=600,bg="white")

p2<-ggplot(tab,aes(factor(KRT_high),doublet_score))+geom_violin(trim=TRUE,fill="grey85")+geom_boxplot(width=.15,outlier.shape=NA)+labs(x="≥3 corrected KRT detected",y="scDblFinder score",title="Residual KRT burden and doublet score")+theme_classic()
ggsave(file.path(FD,"05B4_KRThigh_doublet_violin.pdf"),p2,width=5,height=4.5,bg="white")

writeLines(c(
 "PIPELINE_V2 — 05B4 cluster 25 dissection",
 paste0("Cluster 25 cells: ",length(cells25)),
 "Expression: SoupX-corrected RNA/data and counts.",
 "Geometry and cluster membership: canonical frozen atlas.",
 "KRT-high definition: >=3 detected corrected KRT genes.",
 "SMC-high definition: >=3 detected SMC genes.",
 "No cells removed. No identities changed.",
 "This analysis distinguishes ambient-resistant KRT signal, residual doublets, and coherent contractile/glandular phenotypes."
),file.path(TD,"05B4_audit.txt"))

cat("=== KRT-HIGH vs KRT-LOW ===\n");print(summary_group,row.names=FALSE)
cat("\n=== PATIENT DISTRIBUTION ===\n");print(patient,row.names=FALSE)
cat("\nCells KRT-high:",sum(tab$KRT_high),"/",nrow(tab),"\n")
cat("Cells SMC-high:",sum(tab$SMC_high),"/",nrow(tab),"\n")
cat("Cells KRT-high + SMC-high:",sum(tab$epi_contractile_coexpression),"\n")
cat("\nPASS — focused audit only; atlas unchanged.\n")
