# Processing of miRNA file

library(data.table)

# miR file was downloaded from TargetScan context++
miR=fread("C:/Users/e0545037/Downloads/Predicted_Targets_Context_Scores.default_predictions.txt/Predicted_Targets_Context_Scores.default_predictions.txt")
miR$`Predicted relative KD`[miR$`Predicted relative KD`=="NULL"]=0
miR$`Predicted relative KD`=as.numeric(miR$`Predicted relative KD`)
miR=miR[abs(miR$`Predicted relative KD`)>0]

miR_names=unique(miR$miRNA)

miR_data=list()

for (i in miR_names) {
  miR_i=miR[miR$miRNA==i,]
  miR_data[[i]]=unique(miR_i$`Gene ID`)
  miR_data[[i]]=gsub("\\..*","",miR_data[[i]])
}

save(miR_data,file = "data/miR_data.RData")


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






