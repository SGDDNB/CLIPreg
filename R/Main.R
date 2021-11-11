# rm(list=ls(all=TRUE));
# options(stringsAsFactors=FALSE);

library(pheatmap)
library(data.table)
library(fastmatch)
library(doParallel)
library(foreach)
library(ggnet)
library(ALL)
library(topGO)
library(ggplot2)
library(grid)
library(stringr)

# input
load("data/clusters.Rdata")

# head(clusters)
folder=getwd()

load("data/RBP_POSTAR.RData")
load("data/RBP_ENCODE.RData")


res_Postar=CLIPReg_V4(RBP_data=RBP_POSTAR,cluster=clusters) # run takes several minutes
res_Encode=CLIPReg_V4(RBP_data=RBP_ENCODE,cluster=clusters) # run takes several minutes
#save(res_Encode,file=paste0(folder,"/Fibroblasts/Res_Encode.RData"))
#save(res_Postar,file=paste0(folder,"/Fibroblasts/Res_Postar.RData"))
load("data/Res_Postar.RData") # load the results of the example set instead of running CLIP
load("data/Res_Encode.RData") #

# Combine RBP results from Postar and Encode
all_RBPs=unique(c(res_Encode[[1]]$RBP,res_Postar[[1]]$RBP))

res_both=res_Postar
for (i in 1:nrow(res_Encode[[1]])) {
  r=res_Encode[[1]]$RBP[i]
  if (r%in%res_both[[1]]$RBP) {
    for (n in names(res_both)) {
      for (j in 2:6) {
        res_both[[n]][which(res_both[[1]]$RBP==r),j]=(res_both[[n]][which(res_both[[1]]$RBP==r),j]+res_Encode[[n]][i,j])/2
      }
    }
  } else {
    for (n in names(res_both)) {
      res_both[[n]]=rbind(res_both[[n]],res_Encode[[n]][i,])
    }
  }
}


# keep only significant RBPs

to_keep=c()
for (n in names(res_both)) {
  to_keep=c(to_keep,res_both[[n]]$RBP[res_both[[n]]$pval<0.05])
}
to_keep=unique(to_keep)

for (n in names(res_both)) {
  res_both[[n]]=res_both[[n]][res_both[[n]]$RBP%in%to_keep,]
}
res=res_both



# load fold change and tpm

load("data/tpm_ribo.Rdata")
load("data/ensembl.Rdata")

# IDENTIFIER and geneID loading
load("data/ribo_lfc.Rdata")
ribo_lfc=ribo_lfc[rownames(ribo_lfc)%in%ensembl$geneID,]
index_geneID=match(rownames(ribo_lfc),ensembl$geneID)
ribo_lfc$IDENTIFIER=ensembl$IDENTIFIER[index_geneID]


# Change of RBPs

rbp_lfc=ribo_lfc[ribo_lfc$IDENTIFIER%in%res[[1]]$RBP,2:6]
rbp_lfc=rbp_lfc[!duplicated(rbp_lfc$IDENTIFIER),]
rownames(rbp_lfc)=rbp_lfc$IDENTIFIER
rbp_lfc=rbp_lfc[,-5]

for (n in names(res_both)) {
  res_both[[n]]=res_both[[n]][res_both[[n]]$RBP%in%ribo_lfc$IDENTIFIER,]
}


# Heatmap of RBP scores

e=HeatmapRBP(res=res,RBP_change=rbp_lfc,grid=T)
e # Plot the heatmap

# Save the heatmap
location=paste0(folder,"heatmap.pdf")
n=length(e$tree_row$order)
pdf(location,8,3+n*0.15)
e
dev.off()


# Gather targets from ENCODE and POSTAR that are in the clusters
TargetsEncode=getTarget(RBP_data = RBP_ENCODE,background = clusters$geneID)
TargetsPOSTAR=getTarget(RBP_data = RBP_POSTAR,background = clusters$geneID)
Targets=c(TargetsEncode,TargetsPOSTAR)
Targets=split(unlist(Targets, use.names = FALSE), rep(names(Targets), lengths(Targets)))
for (t in names(Targets)) {
  Targets[[t]]=unique(Targets[[t]])
}


# Bubble plot clusters if clusters are

BubbleRBPs(RBP_res = res,clusters = clusters)

# Draw network

Draw_network(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],clusters=clusters,n=5)
Draw_network_by_group(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],Targets=Targets,clusters=clusters,n=5)

# plot GO of nodes of size >th.
Plot_GO(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],Targets=Targets,clusters=clusters,n=5,
  all_genes=rownames(tpm_ribo),th=200,GO_to_show=3)
Plot_GO_node_name(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],Targets=Targets,clusters=clusters,n=5,
                  all_genes=rownames(tpm_ribo),Nodes_to_keep=c(19,15),GO_to_show=3)
Plot_GO(rbp_of_interest="QKI",Targets=Targets,clusters=clusters,GO_to_show=3)



























