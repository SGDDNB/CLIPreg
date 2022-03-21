# Processing of CLIPseq summary files


# Step 1 : Download all RBP clip-seq bed files from encode and combine them into 1 be
# Step 2 : Filter to keep high score peaks only. >8fold enrichment and P value <10e-5
# Step 3 : Intersect bed file with Ensembl genes gtf file to find target genes
# Step 4 : Shortlist columns to keep RBP and targets

# For Postar data, download from POSTAR3 website CLIPdb
# Start from step3 of ENCODE


# Processing of miRNA file

library(data.table)
library(miRBaseConverter)
library(stringr)

# miR file was downloaded from TargetScan context++
miR=fread("C:/Users/e0545037/Downloads/Predicted_Targets_Context_Scores.default_predictions.txt/Predicted_Targets_Context_Scores.default_predictions.txt")

index=str_sub(miR$miRNA,-1,-1)=="p"
miR$miRNA[index]=gsub('.{3}$', '',miR$miRNA[index] )

miRBase=getAllMiRNAs()
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "mirbase_id", attributes= c("ensembl_gene_id","mirbase_id"),values=miRBase$Name,mart= mart)

index=str_sub(G_list$mirbase_id,-1,-1)=="p"
G_list$mirbase_id[index]=gsub('.{3}$', '',G_list$mirbase_id[index] )
index=str_sub(G_list$mirbase_id,-2,-2)=="-"
G_list$mirbase_id[index]=gsub('.{2}$', '',G_list$mirbase_id[index] )

index=str_sub(miRBase$Name,-1,-1)=="p"
miRBase$Name[index]=gsub('.{3}$', '',miRBase$Name[index] )
index=str_sub(miRBase$Name,-2,-2)=="-"
miRBase$Name[index]=gsub('.{2}$', '',miRBase$Name[index] )

# miR$`Predicted relative KD`[miR$`Predicted relative KD`=="NULL"]=0
# miR$`Predicted relative KD`=as.numeric(miR$`Predicted relative KD`)
# miR=miR[abs(miR$`Predicted relative KD`)>0]

Thres=-0.5
miR=miR[tolower(miR$miRNA)%in%tolower(miRBase$Name),]
miR=miR[miR$`weighted context++ score`<Thres,]
miR_names=unique(miR$miRNA)


miR_data=list()

for (i in miR_names) {
  miR_i=miR[miR$miRNA==i,]
  miR_data[[i]]=unique(miR_i$`Gene ID`)
  miR_data[[i]]=gsub("\\..*","",miR_data[[i]])
}
miR_data=miR_data[sort(names(miR_data))]
names(miR_data)
miR_data=subset(miR_data,tolower(names(miR_data))%in%tolower(G_list$mirbase_id))

save(miR_data,file = "data/miR_data0.1.RData")
save(miR_data,file = "data/miR_data0.5.RData")

G_list=G_list[G_list$mirbase_id%in%c(miR_names,tolower(miR_names)),]
colnames(G_list)=c("geneID","miRBase_ID")
miR_info=G_list

save(miR_info,file = "data/miR_info.RData")

### Example of processing data from HeLa cells stimulated by EGF

# Get gene groups from existing data

file="HeLa_EGF_TE.csv"
df=read.csv(file)
sig_30=df$IDENTIFIER[df$pval_TE_30<0.05]
sig_60=df$IDENTIFIER[df$pval_TE_60<0.05]
sig_90=df$IDENTIFIER[df$pval_TE_90<0.05]

TE_genes=unique(c(sig_30,sig_60,sig_90))

df=df[df$IDENTIFIER%in%TE_genes,]
df=df[!duplicated(df$IDENTIFIER),]

gene_groups=data.frame(IDENTIFIER=df$IDENTIFIER,Gene_group=0)

for (i in 1:nrow(df)) {
  ID=df$IDENTIFIER[i]
  lowest_pval=which(df[i,4:6]==min(df[i,4:6]))
  TE_sign=sign(df[i,lowest_pval+1])
  if (TE_sign>0) {gene_groups$Gene_group[i]="TE_up"} else {gene_groups$Gene_group[i]="TE_down"}
}


library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- df$IDENTIFIER
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)

colnames(G_list)=c("geneID","IDENTIFIER")

gene_groups=merge(G_list,gene_groups,by.y="IDENTIFIER")
save(gene_groups,file="data/gene_groups_HeLa_EGF.RData")

# prepare ribo_lfc data

file="HeLa_EGF_ribo_lfc.csv"
df=read.csv(file)
df[which(is.na(df),arr.ind = T)]=1

geneNames=unique(df$IDENTIFIER)
ribo_lfc=data.frame(IDENTIFIER=geneNames,FoldChange=0)

for (i in 1:nrow(df)) {
  lowest_pval=min(df[i,6:8])
  if (lowest_pval<0.05) {
    lowest_pval_index=which(df[i,6:8]==min(df[i,6:8]))[1]
    ribo_lfc[ribo_lfc$IDENTIFIER==df$IDENTIFIER[i],"FoldChange"]=df[i,lowest_pval_index+2]
  }
}


library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- ribo_lfc$IDENTIFIER
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)

colnames(G_list)=c("geneID","IDENTIFIER")
ribo_lfc=merge(G_list,ribo_lfc,by.y="IDENTIFIER")
ribo_lfc=ribo_lfc[!duplicated(ribo_lfc$IDENTIFIER),]
rownames(ribo_lfc)=ribo_lfc$IDENTIFIER
save(ribo_lfc,file="data/ribo_lfc_HeLa_EGF.RData")

# TPM

file="HeLa_EGF_Ribo_TPM.csv"
df=read.csv(file)

df=df[rowSums(df[,2:7] > 1) >= 1,]
df=df[!duplicated(df$IDENTIFIER),]

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- df$IDENTIFIER
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(G_list)

colnames(G_list)=c("geneID","IDENTIFIER")
tpm_ribo=merge(G_list,df,by.y="IDENTIFIER")
tpm_ribo=tpm_ribo[!duplicated(tpm_ribo$IDENTIFIER),]
rownames(tpm_ribo)=tpm_ribo$IDENTIFIER
save(tpm_ribo,file="data/tpm_ribo_HeLa_EGF.RData")


###################### to delete

df=res_miR$intensified_up
df=df[order(df$padj),]
miR_sig=df$RBP[1:3]
geneID_sig=miR_info$ensembl_gene_id[miR_info$mirbase_id%in%miR_sig]

ggplot(tpm_all_RNA[geneID_sig,])+geom_point()

df1=as.numeric(tpm_all_RNA[geneID_sig[1],])
df1=data.frame(tpm=df1,
               timepoint=factor(rep(c("24h","2h","45min","6h","BSL"),4),
                                levels = c("BSL","45min","2h","6h","24h")))
ggplot(df1,aes(x=timepoint,y=tpm))+geom_point()+theme_bw()

#

df2=as.numeric(tpm_all_RNA[geneID_sig[2],])
df2=data.frame(tpm=df2,
               timepoint=factor(rep(c("24h","2h","45min","6h","BSL"),4),
                                levels = c("BSL","45min","2h","6h","24h")))
ggplot(df2,aes(x=timepoint,y=tpm))+geom_point()+theme_bw()


#

df3=as.numeric(tpm_all_RNA[geneID_sig[3],])
df3=data.frame(tpm=df3,
               timepoint=factor(rep(c("24h","2h","45min","6h","BSL"),4),
                                levels = c("BSL","45min","2h","6h","24h")))
ggplot(df3,aes(x=timepoint,y=tpm))+geom_point()+theme_bw()

#

df4=as.numeric(tpm_all_RNA[geneID_sig[4],])
df4=data.frame(tpm=df4,
               timepoint=factor(rep(c("24h","2h","45min","6h","BSL"),4),
                                levels = c("BSL","45min","2h","6h","24h")))
ggplot(df4,aes(x=timepoint,y=tpm))+geom_point()+theme_bw()






