rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);


df=read.csv("HeLa_EGF_raw_counts.csv",header = T) # raw counts of RIBO and RNA
coldata=data.frame(SampleID=colnames(df[,-1]),Condition=rep("1",24),SeqType=c(rep("Ribo-seq",12),rep("RNA-seq",12)))
coldata$Condition[grepl("EGF",coldata$SampleID)]=2

coldata_30=coldata[grepl("30_min",coldata$SampleID),]
counts_30=df[,c("geneID",coldata_30$SampleID)]


coldata_60=coldata[grepl("60_min",coldata$SampleID),]
counts_60=df[,c("geneID",coldata_60$SampleID)]


coldata_90=coldata[grepl("90_min",coldata$SampleID),]
counts_90=df[,c("geneID",coldata_90$SampleID)]


gene_groups_30=deltaTE_gene_groups(counts=counts_30,coldata = coldata_30,batch=0)
gene_groups_60=deltaTE_gene_groups(counts=counts_60,coldata = coldata_60,batch=0)
gene_groups_90=deltaTE_gene_groups(counts=counts_90,coldata = coldata_90,batch=0)

gene_groups=rbind(gene_groups_30,gene_groups_60,gene_groups_90)
gene_groups=unique(gene_groups)
