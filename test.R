library(data.table)

Folder="C:/Users/e0545037/Desktop/Baptiste/PhD/CLIPreg/Fibroblasts/"

Ribo1=fread(paste0(Folder,"lfc_Ribo_H00M45.tab"))
Ribo2=fread(paste0(Folder,"lfc_Ribo_H02M00.tab"))
Ribo3=fread(paste0(Folder,"lfc_Ribo_H06M00.tab"))
Ribo4=fread(paste0(Folder,"lfc_Ribo_H24M00.tab"))

RNA1=fread(paste0(Folder,"lfc_RNA_H00M45.tab"))
RNA2=fread(paste0(Folder,"lfc_RNA_H02M00.tab"))
RNA3=fread(paste0(Folder,"lfc_RNA_H06M00.tab"))
RNA4=fread(paste0(Folder,"lfc_RNA_H24M00.tab"))

RNA=cbind(RNA1,RNA2[,-1],RNA3[,-1],RNA4[,-1])
colnames(RNA)=c("geneID",
                paste0(colnames(RNA1)[-1],".0.45"),
                paste0(colnames(RNA1)[-1],".2"),
                paste0(colnames(RNA1)[-1],".6"),
                paste0(colnames(RNA1)[-1],".24"))

Ribo=cbind(Ribo1,Ribo2[,-1],Ribo3[,-1],Ribo4[,-1])
colnames(Ribo)=c("geneID",
                paste0(colnames(Ribo1)[-1],".0.45"),
                paste0(colnames(Ribo1)[-1],".2"),
                paste0(colnames(Ribo1)[-1],".6"),
                paste0(colnames(Ribo1)[-1],".24"))

miR_geneIDs=data.frame(geneID=c("ENSG00000199165","ENSG00000198974","ENSG00000207721",
                                "ENSG00000208014","ENSG00000207708","ENSG00000207730",
                                "ENSG00000211514","ENSG00000212027","ENSG00000207996",
                                "ENSG00000207808","ENSG00000199121","ENSG00000207980",
                                "ENSG00000207568","ENSG00000207585","ENSG00000283819"),
                       miR_name=c("hsa-let-7a-5p","hsa−miR−30e−5p","hsa−miR−186−5p",
                                  "hsa−miR−653−5p","hsa−miR−141−3p","hsa−miR−200a−3p",
                                  "hsa−miR−454−3p","hsa−miR−374b−5p","hsa−miR−301b−3p",
                                  "hsa−miR−27a−3p","hsa−miR−26b−5p","hsa−miR−23a−3p",
                                  "hsa−miR−203a−3p.1","hsa−miR−181d−5p","hsa−miR−144−3p"))

mir_RNA_lfc=RNA[RNA$geneID%in%miR_geneIDs$geneID,]








