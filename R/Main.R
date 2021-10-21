# rm(list=ls(all=TRUE));
# options(stringsAsFactors=FALSE);

library(pheatmap)
library(data.table)
library(fastmatch)
library(doParallel)
library(ggnet)
library(ggplot2)
library(grid)
library(stringr)

# input

folder="C:/Users/e0545037/Desktop/Baptiste/PhD/CLIPreg"
# RBP_data="rbp_gene_sel.txt"
# RBP_data="rbp_gene_postar.txt"

clusters_file="Fibroblasts/Clusters.csv"

iterations=100000
#res=CLIPReg_V3(folder=folder,RBP_data=RBP_data,cluster=clusters_file,iterations=100000)
#save(res_Encode,file=paste0(folder,"/Fibroblasts/Res_RBP_Encode_V2.RData"))
#save(res_Postar,file=paste0(folder,"/Fibroblasts/Res_RBP_Postar_V2.RData"))
load(paste0(folder,"/Fibroblasts/Res_RBP_Postar.RData"))
load(paste0(folder,"/Fibroblasts/Res_RBP_Encode_V2.RData"))



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

setwd("C:/Users/e0545037/Desktop/Baptiste/PhD/CLIPreg/Fibroblasts")

rna_lfc = as.data.frame(read.delim(paste(sep="","Sig_genes_all3_RNA_time_series_pw_LRT.lfc")))
ribo_lfc = as.data.frame(read.delim(paste(sep="","Sig_genes_all3_Ribo_time_series_pw_LRT.lfc")))
ribOnly_lfc = as.data.frame(read.delim(paste(sep="","Sig_genes_all3_RibOnly_time_series_pw_LRT.lfc")))

tpm_rna = read.delim("tpm_rna.tab")
tpm_ribo = read.delim("tpm_ribo.tab ")

# IDENTIFIER and geneID loading
ensembl=read.csv("C:/Users/e0545037/Desktop/Baptiste/Bioinfo/Source/Homo_sapiens.GRCh38.95_ID_L_Name.csv")
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

e=HeatmapRBP(res=res,RBP_change=rbp_lfc,grid=F)
e # Plot the heatmap

# Save the heatmap
location="C:/Users/e0545037/Desktop/Baptiste/PhD/CLIPreg/Fibroblasts/Heatmap_fibroblasts_V4.pdf"
n=length(e$tree_row$order)
pdf(location,8,3+n*0.15)
e
dev.off()


# Get targets from POSTAR and ENCODE and combine them

TargetsEncode=getTarget(folder = folder,RBP_data = "rbp_gene_sel.txt",background = clusters$geneID)
TargetsPOSTAR=getTarget(folder = folder,RBP_data = "rbp_gene_postar.txt",background = clusters$geneID)
Targets=c(TargetsEncode,TargetsPOSTAR)
Targets=split(unlist(Targets, use.names = FALSE), rep(names(Targets), lengths(Targets)))
for (t in names(Targets)) {
  Targets[[t]]=unique(Targets[[t]])
}


# Heatmap of targets of 1 RBP

clusters=fread(paste0(folder,"/",clusters_file),header = T)
index_geneID_clusters=match(clusters$geneID,ensembl$geneID)
clusters$IDENTIFIER=ensembl$IDENTIFIER[index_geneID_clusters]

HeatmapTargets(folder="Folder",Clusters=clusters,RBP_name="YTHDF3",targets=Targets,
               counts_rna=counts_rna, counts_ribo=counts_ribo[,-16],scale=T,Not_DE=F)

HeatmapTargets(folder="Folder",Clusters=clusters,RBP_name="YTHDF3",targets=Targets,
               counts_rna=LFC_RNA[,c(2,8,14)], counts_ribo=LFC_RIBO[,c(3,8,13)],scale=F)


# Bubble plot clusters if clusters are

BubbleRBPs(RBP_res = res,clusters = clusters)

# Draw network

Draw_network(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],clusters=clusters,n=5)
Draw_network_by_group(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],Targets=Targets,clusters=clusters,n=5)

# plot GO of nodes of size >th.
Plot_GO(rbp_lfc=rbp_lfc,res=res_both[-c(1,2)],Targets=Targets,clusters=clusters,n=5,
  all_genes=rownames(tpm_ribo),th=200,GO_to_show=3)




























