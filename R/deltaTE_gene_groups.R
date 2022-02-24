#' @title deltaTE_gene_groups
#'
#' @description Computes deltaTE from raw RNA and RIBO counts
#'
#' @param symbol counts, coldata, batch
#'
#' @return
#'
#' @examples deltaTE_gene_groups(counts=counts,coldata=coldata,batch=0)
#'
#' @export
#'
#'
#'
deltaTE_gene_groups <-function(counts=counts,coldata=coldata,batch=0)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  # Read and merge count matrices
  merged=counts
  rownames(merged)=merged[,1]
  merged=merged[,-1]

  # Sample information file
  coldata <- as.data.frame(coldata)
  coldata <- as.data.frame(apply(coldata,2,as.factor))
  rownames(coldata)=coldata$SampleID

  #head(coldata)

  ribo_sample=coldata$SampleID[which(coldata$SeqType=="Ribo-seq")]
  ribo=merged[,ribo_sample]

  rna_sample=coldata$SampleID[which(coldata$SeqType=="RNA-seq")]
  rna=merged[,rna_sample]

  ## Detecting differential translation regulation
  ### DESeq2 object with batch and interaction term in the design

  if(batch == 1){
    ddsMat_merged <- DESeqDataSetFromMatrix(countData = merged,
                                     colData = coldata, design = ~ Batch + Condition + SeqType + Condition:SeqType)
  }else if(batch == 0){
    ddsMat_merged <- DESeqDataSetFromMatrix(countData = merged,
                                     colData = coldata, design =~ Condition + SeqType + Condition:SeqType)
  }else{
    stop("Batch presence should be indicated by 0 or 1 only", call.=FALSE)
  }

  ddsMat_merged$SeqType = relevel(ddsMat_merged$SeqType,"RNA-seq")
  ddsMat_merged <- DESeq(ddsMat_merged)
  resultsNames(ddsMat_merged)

  system("mkdir Results")
  setwd("Results")
  system("mkdir fold_changes")
  system("mkdir gene_lists")

  res_merged <- results(ddsMat_merged, contrast=list("Condition2.SeqTypeRibo.seq"))
  summary(res_merged)
  length(which(res_merged$padj < 0.05))
  write.table(rownames(res_merged)[which(res_merged$padj < 0.05)],"gene_lists/DTEGs.txt",quote=F,sep="\t",col.names = F,row.names = F)
  write.table(res_merged,"fold_changes/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)

  ribOnly_lfc=as.data.frame(res_merged@listData)
  rownames(ribOnly_lfc)=res_merged@rownames


  ### DESeq2 object with batch for Ribo-seq
  ind = which(coldata$SeqType == "Ribo-seq")
  coldata_ribo = coldata[ind,]

  # PCA
  if(batch == 1){
    ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                          colData = coldata_ribo, design =~ Batch + Condition)
    }else if(batch ==0){
    ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                          colData = coldata_ribo, design =~ Condition)
    }

  ddsMat_ribo <- DESeq(ddsMat_ribo)
  res_ribo <- results(ddsMat_ribo, contrast=c("Condition","2","1"))
  res_ribo <- lfcShrink(ddsMat_ribo, contrast=c("Condition","2","1"),res=res_ribo,type = "ashr")
  write.table(res_ribo,"fold_changes/deltaRibo.txt",quote=F,sep="\t",col.names = T,row.names = T)

  ribo_lfc=as.data.frame(res_ribo@listData)
  rownames(ribo_lfc)=res_ribo@rownames

  ### DESeq2 object with batch for RNA-seq
  ind = which(coldata$SeqType == "RNA-seq")
  coldata_rna = coldata[ind,]
  # PCA
  if(batch == 1){
    ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                         colData = coldata_rna, design =~ Batch + Condition)
    }else if(batch ==0){
    ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                         colData = coldata_rna, design =~ Condition)
    vsd <- vst(ddsMat_rna)
  }

  ddsMat_rna <- DESeq(ddsMat_rna)


  res_rna <- results(ddsMat_rna, contrast=c("Condition","2","1"))
  res_rna <- lfcShrink(ddsMat_rna, contrast=c("Condition","2","1"),res=res_rna,type = "ashr")
  write.table(res_rna,"fold_changes/deltaRNA.txt",quote=F,sep="\t",col.names = T,row.names = T)
  write.table(rownames(res_rna)[which(res_rna$padj < 0.05)],"gene_lists/DTG.txt",quote=F,sep="\t",col.names = F,row.names = F)

  rna_lfc=as.data.frame(res_rna@listData)
  rownames(rna_lfc)=res_rna@rownames

  ## Classes of genes
  rna_lfc$padj[which(is.na(rna_lfc$padj))]=1
  ribo_lfc$padj[which(is.na(ribo_lfc$padj))]=1
  ribOnly_lfc$padj[which(is.na(ribOnly_lfc$padj))]=1

  fw_up=c()
  fw_dn=c()
  it_up=c()
  it_dn=c()
  bd_up=c()
  bd_dn=c()
  ee_up=c()
  ee_dn=c()

  for (i in 1:nrow(rna_lfc)) {
    if (rna_lfc$padj[i]<0.05 & ribo_lfc$padj[i] < 0.05 & ribOnly_lfc$padj[i] > 0.05) {
      if (rna_lfc$log2FoldChange[i]>0 & ribo_lfc >0) {
        fw_up=c(fw_up,rownames(rna_lfc)[i])
      } else {fw_dn=c(fw_dn,rownames(rna_lfc)[i])}
    }
    if (rna_lfc$padj[i]<0.05 & ribo_lfc$padj[i] < 0.05 & ribOnly_lfc$padj[i] < 0.05) {
      if (rna_lfc$log2FoldChange[i]>0 & ribOnly_lfc$log2FoldChange[i] > 0) {it_up=c(it_up,rownames(rna_lfc)[i])}
      if (rna_lfc$log2FoldChange[i]<0 & ribOnly_lfc$log2FoldChange[i] > 0) {it_dn=c(it_dn,rownames(rna_lfc)[i])}
      if (rna_lfc$log2FoldChange[i]>0 & ribOnly_lfc$log2FoldChange[i] < 0) {bd_dn=c(bd_dn,rownames(rna_lfc)[i])}
      if (rna_lfc$log2FoldChange[i]<0 & ribOnly_lfc$log2FoldChange[i] < 0) {bd_up=c(bd_up,rownames(rna_lfc)[i])}
    }
    if (rna_lfc$padj[i]>0.05 & ribo_lfc$padj[i] < 0.05 & ribOnly_lfc$padj[i] < 0.05) {
      if (ribo_lfc$log2FoldChange[i]>0) {ee_up=c(ee_up,rownames(rna_lfc)[i])}
      else {ee_dn=c(ee_dn,rownames(rna_lfc)[i])}
    }
    if (rna_lfc$padj[i]<0.05 & ribo_lfc$padj[i] > 0.05 & ribOnly_lfc$padj[i] < 0.05) {
      if (rna_lfc$log2FoldChange[i]>0) {bd_dn=c(bd_dn,rownames(rna_lfc)[i])}
      else {bd_up=c(bd_up,rownames(rna_lfc)[i])}
    }


  }

  gene_groups=data.frame(geneID=c(fw_up,fw_dn,ee_up,ee_dn,bd_up,bd_dn,it_up,it_dn),
                              Gene_group=c(rep("forwarded_up",length(fw_up)),rep("forwarded_down",length(fw_dn)),
                                        rep("exclusive_up",length(ee_up)),rep("exclusive_down",length(ee_dn)),
                                        rep("buffered_up",length(bd_up)),rep("buffered_down",length(bd_dn)),
                                        rep("intensified_up",length(it_up)),rep("intensified_down",length(it_dn))))

  print(table(gene_groups$Gene_group))


  write.csv(gene_groups,"gene_lists/gene_groups.csv")
  return(gene_groups)

}
