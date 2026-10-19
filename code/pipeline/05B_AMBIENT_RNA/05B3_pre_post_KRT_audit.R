# =============================================================================
# PIPELINE_V2 — 05B3_pre_post_KRT_audit.R
# Paired original-vs-SoupX KRT audit on the exact same 82,893 cells/clusters.
# =============================================================================
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2)})
options(stringsAsFactors=FALSE);source("PIPELINE_V2/config/paths.R")

cat("\n============================================================\n")
cat("PIPELINE_V2 — 05B3 PRE/POST KRT AUDIT\n")
cat("============================================================\n\n")

PRE_FILE<-file.path(V2_FROZEN,"global_atlas_v2_frozen.rds")
POST_FILE<-file.path(V2_ROOT,"objects","ambient_corrected","global_atlas_v2_SoupX_expression.rds")
OUT<-file.path(V2_ROOT,"results","QC","ambient_RNA","05B3_pre_post_KRT")
TD<-file.path(OUT,"tables");FD<-file.path(OUT,"figures")
dir.create(TD,recursive=TRUE,showWarnings=FALSE);dir.create(FD,recursive=TRUE,showWarnings=FALSE)

pre<-readRDS(PRE_FILE);post<-readRDS(POST_FILE)
DefaultAssay(pre)<-"RNA";DefaultAssay(post)<-"RNA"
stopifnot(ncol(pre)==82893,ncol(post)==82893,nrow(pre)==23412,nrow(post)==23412)
stopifnot(identical(colnames(pre),colnames(post)),identical(rownames(pre),rownames(post)))
stopifnot(identical(as.character(pre$RNA_snn_res.0.6),as.character(post$RNA_snn_res.0.6)))

KRT<-c("KRT5","KRT14","KRT15","KRT17","KRT1","KRT10")
stopifnot(all(KRT%in%rownames(pre)))
cl<-as.character(pre$RNA_snn_res.0.6)
cells<-colnames(pre)

preC<-LayerData(pre,assay="RNA",layer="counts")[KRT,,drop=FALSE]
postC<-LayerData(post,assay="RNA",layer="counts")[KRT,,drop=FALSE]
preD<-LayerData(pre,assay="RNA",layer="data")[KRT,,drop=FALSE]
postD<-LayerData(post,assay="RNA",layer="data")[KRT,,drop=FALSE]

preN<-Matrix::colSums(preC>0);postN<-Matrix::colSums(postC>0)
preMean<-Matrix::colMeans(preD);postMean<-Matrix::colMeans(postD)

# Cell-level paired transitions.
cell<-data.frame(
 cell_id=cells,cluster=cl,patient_id=pre$patient_id,disease=pre$disease,batch=pre$batch,
 pre_nKRT=preN,post_nKRT=postN,pre_KRT_mean=preMean,post_KRT_mean=postMean,
 delta_KRT_mean=postMean-preMean,doublet_score=pre$doublet_score,nFeature_RNA=pre$nFeature_RNA,
 stringsAsFactors=FALSE
)
for(g in KRT){
 cell[[paste0("pre_",g,"_det")]]<-as.numeric(preC[g,]>0)
 cell[[paste0("post_",g,"_det")]]<-as.numeric(postC[g,]>0)
}
con<-gzfile(file.path(TD,"05B3_cell_paired_KRT.csv.gz"),"wt");write.csv(cell,con,row.names=FALSE);close(con)

# Cluster summary.
summ<-do.call(rbind,lapply(0:26,function(k){
 ii<-which(cl==as.character(k))
 r<-data.frame(
  cluster=k,n_cells=length(ii),
  pre_pct_ge1=100*mean(preN[ii]>=1),post_pct_ge1=100*mean(postN[ii]>=1),
  pre_pct_ge2=100*mean(preN[ii]>=2),post_pct_ge2=100*mean(postN[ii]>=2),
  pre_pct_ge3=100*mean(preN[ii]>=3),post_pct_ge3=100*mean(postN[ii]>=3),
  pre_pct_ge4=100*mean(preN[ii]>=4),post_pct_ge4=100*mean(postN[ii]>=4),
  pre_mean_KRT=mean(preMean[ii]),post_mean_KRT=mean(postMean[ii]),
  stringsAsFactors=FALSE
 )
 for(g in KRT){
  r[[paste0("pre_pct_",g)]]<-100*mean(preC[g,ii]>0)
  r[[paste0("post_pct_",g)]]<-100*mean(postC[g,ii]>0)
 }
 r
}))
for(metric in c("pct_ge1","pct_ge2","pct_ge3","pct_ge4","mean_KRT"))
 summ[[paste0("delta_",metric)]]<-summ[[paste0("post_",metric)]]-summ[[paste0("pre_",metric)]]
write.csv(summ,file.path(TD,"05B3_cluster_pre_post_KRT.csv"),row.names=FALSE)

non_epi<-c(4,9,11,18,19,21,25,26)
write.csv(summ[summ$cluster%in%non_epi,],file.path(TD,"05B3_non_epithelial_pre_post_KRT.csv"),row.names=FALSE)

# Patient × cluster summary.
pc<-do.call(rbind,lapply(split(seq_len(ncol(pre)),list(pre$patient_id,cl),drop=TRUE),function(ii){
 data.frame(patient_id=pre$patient_id[ii[1]],disease=pre$disease[ii[1]],batch=pre$batch[ii[1]],cluster=as.integer(cl[ii[1]]),n_cells=length(ii),
 pre_pct_ge3=100*mean(preN[ii]>=3),post_pct_ge3=100*mean(postN[ii]>=3),
 pre_mean_KRT=mean(preMean[ii]),post_mean_KRT=mean(postMean[ii]))
}))
write.csv(pc,file.path(TD,"05B3_patient_cluster_pre_post_KRT.csv"),row.names=FALSE)

# Paired transition table for stringent burden.
trans<-as.data.frame(table(cluster=factor(cl,levels=as.character(0:26)),pre_ge3=preN>=3,post_ge3=postN>=3),stringsAsFactors=FALSE)
names(trans)[4]<-"n_cells";write.csv(trans,file.path(TD,"05B3_ge3KRT_paired_transitions.csv"),row.names=FALSE)

# Heatmaps with identical scale.
hm<-rbind(
 do.call(rbind,lapply(KRT,function(g)data.frame(cluster=0:26,gene=g,state="PRE",pct=sapply(0:26,function(k)100*mean(preC[g,cl==as.character(k)]>0))))),
 do.call(rbind,lapply(KRT,function(g)data.frame(cluster=0:26,gene=g,state="POST",pct=sapply(0:26,function(k)100*mean(postC[g,cl==as.character(k)]>0)))))
)
mx<-max(hm$pct)
for(st in c("PRE","POST")){
 z<-hm[hm$state==st,];z$cluster<-factor(z$cluster,levels=rev(0:26))
 p<-ggplot(z,aes(gene,cluster,fill=pct))+geom_tile()+scale_fill_gradient(low="white",high="black",limits=c(0,mx))+labs(x=NULL,y="Global cluster",fill="% positive",title=paste("KRT detection",st))+theme_classic(base_size=10)
 ggsave(file.path(FD,paste0("05B3_KRT_heatmap_",st,".pdf")),p,width=7,height=8,bg="white")
 ggsave(file.path(FD,paste0("05B3_KRT_heatmap_",st,".png")),p,width=7,height=8,dpi=600,bg="white")
}

# Delta heatmap.
wide<-merge(hm[hm$state=="PRE",c("cluster","gene","pct")],hm[hm$state=="POST",c("cluster","gene","pct")],by=c("cluster","gene"),suffixes=c("_pre","_post"))
wide$delta<-wide$pct_post-wide$pct_pre;wide$cluster<-factor(wide$cluster,levels=rev(0:26));lim<-max(abs(wide$delta))
pd<-ggplot(wide,aes(gene,cluster,fill=delta))+geom_tile()+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-lim,lim))+labs(x=NULL,y="Global cluster",fill="POST-PRE\npercentage points",title="Change in KRT detection after SoupX")+theme_classic(base_size=10)
ggsave(file.path(FD,"05B3_KRT_delta_heatmap.pdf"),pd,width=7,height=8,bg="white");ggsave(file.path(FD,"05B3_KRT_delta_heatmap.png"),pd,width=7,height=8,dpi=600,bg="white")

# Cluster-level PRE vs POST mean KRT scatter.
ps<-ggplot(summ,aes(pre_mean_KRT,post_mean_KRT,label=cluster))+geom_abline(slope=1,intercept=0,linetype=2)+geom_point(size=2)+geom_text(nudge_y=.01,size=3)+labs(x="PRE mean KRT programme",y="POST mean KRT programme",title="Cluster-level KRT programme before vs after SoupX")+theme_classic(base_size=10)
ggsave(file.path(FD,"05B3_cluster_KRT_pre_post_scatter.pdf"),ps,width=6,height=5,bg="white");ggsave(file.path(FD,"05B3_cluster_KRT_pre_post_scatter.png"),ps,width=6,height=5,dpi=600,bg="white")

# Cluster 25 detailed paired metrics.
c25<-cell[cell$cluster=="25",]
c25sum<-data.frame(
 n_cells=nrow(c25),
 pre_ge1=sum(c25$pre_nKRT>=1),post_ge1=sum(c25$post_nKRT>=1),
 pre_ge3=sum(c25$pre_nKRT>=3),post_ge3=sum(c25$post_nKRT>=3),
 pre_ge4=sum(c25$pre_nKRT>=4),post_ge4=sum(c25$post_nKRT>=4),
 cells_losing_ge3=sum(c25$pre_nKRT>=3 & c25$post_nKRT<3),
 cells_remaining_ge3=sum(c25$pre_nKRT>=3 & c25$post_nKRT>=3),
 mean_pre_KRT=mean(c25$pre_KRT_mean),mean_post_KRT=mean(c25$post_KRT_mean),
 median_doublet_ge3_pre=median(c25$doublet_score[c25$pre_nKRT>=3],na.rm=TRUE),
 median_doublet_lt3_pre=median(c25$doublet_score[c25$pre_nKRT<3],na.rm=TRUE)
)
write.csv(c25sum,file.path(TD,"05B3_cluster25_detailed.csv"),row.names=FALSE)

writeLines(c(
 "PIPELINE_V2 — 05B3 paired PRE/POST KRT audit",
 paste0("Cells: ",ncol(pre)),
 paste0("Genes: ",nrow(pre)),
 paste0("KRT panel: ",paste(KRT,collapse=", ")),
 "PRE and POST cell IDs/order: EXACT.",
 "Clusters are canonical RNA_snn_res.0.6 and were not recalculated.",
 "Detection is based on counts > 0; programme magnitude on LogNormalized data.",
 "No filtering or annotation performed."
),file.path(TD,"05B3_audit.txt"))

cat("\n=== NON-EPITHELIAL PRE vs POST ===\n")
print(summ[summ$cluster%in%non_epi,c("cluster","n_cells","pre_pct_ge1","post_pct_ge1","pre_pct_ge3","post_pct_ge3","pre_pct_ge4","post_pct_ge4","pre_mean_KRT","post_mean_KRT")],row.names=FALSE)
cat("\n=== CLUSTER 25 ===\n");print(c25sum,row.names=FALSE)
cat("\nPASS — paired sensitivity analysis; canonical atlas unchanged.\n")
