#' @title Visualise
#' @description Generate visual output
#' @param symbol
#' @return
#' @examples
#' @export
#'
Visualise=function(results=results,folder=getwd()){
  gene_groups=results$gene_groups
  tpm_ribo=results$tpm_ribo
  rbp_lfc=results$rbp_lfc
  ribo_lfc=results$ribo_lfc
  Targets=results$Targets
  res=results$res

  save(results,file=paste0(folder,"/Results_CLIPreg.RData"))

  A=HeatmapRBP(res=res,rbp_lfc = rbp_lfc)
  n=length(A@ht_list$`RBP direction`@matrix)
  pdf(paste0(folder,"/Heatmap.pdf"),7,3+n*0.15,useDingbats = F)
  print(A)
  dev.off()

  B=BubbleRBPs(res=res,gene_groups=gene_groups,FDR=0.1)
  pdf(paste0(folder,"/Bubble_plot.pdf"),6,4.5,useDingbats = F)
  print(B)
  dev.off()

  C=Draw_network_by_group(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,forwarded = F)
  pdf(paste0(folder,"/Network_5_most_DE_RBP.pdf"),6,4.5,useDingbats = F)
  print(C)
  dev.off()

  D=Plot_GO(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,
          tpm_ribo = tpm_ribo,th=200,GO_to_show=3,forwarded=F)
  n=nrow(D$data)
  pdf(paste0(folder,"/GO_from_nodes_200+.pdf"),12,3+0.3*n,useDingbats = F)
  print(D)
  dev.off()

  rbp_of_interest=as.character(rbp_lfc$IDENTIFIER)[1]
  E=Plot_GO_RBP(rbp_of_interest=rbp_of_interest,tpm_ribo = tpm_ribo,Targets=Targets,gene_groups=gene_groups,GO_to_show=3)
  pdf(paste0(folder,"/GO_RBP_",rbp_of_interest,".pdf"),12,3+n*0.3,useDingbats = F)
  print(E)
  dev.off()


}
