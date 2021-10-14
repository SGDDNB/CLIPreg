#' @title HeatmapTargets
#'
#' @description Plot the counts of the targets of a specified RBP
#'
#' @param symbol folder, Clusters, RBP_name, targets, counts_rna, counts_ribo, scale, Not_DE)
#'
#' @return Plot a heatmap
#'
#' @examples HeatmapTargets(folder="Folder",Clusters=clusters,RBP_name="RBM47",targets=Targets,counts_rna=counts_rna, counts_ribo=counts_ribo,scale=T,Not_DE=F)
#'
#' @export
#'
#'
#'
HeatmapTargets <-function(folder="Folder",Clusters=clusters,RBP_name="RBM47",targets=Targets,
                      counts_rna=counts_rna, counts_ribo=counts_ribo,scale=T,Not_DE=F)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  if (scale==T) {
    counts_ribo=scale(counts_ribo)
    counts_rna=scale(counts_rna)
    counts_ribo=t(scale(t(counts_ribo)))
    counts_rna=t(scale(t(counts_rna)))
  }

  annotation=data.frame(Cluster=rep("Not DE",length(targets[[RBP_name]])))
  j=0
  for (i in targets[[RBP_name]]) {
    j=j+1
    if (i%in%clusters$geneID) {
      annotation[j,1]=clusters$Cluster[clusters$geneID==i]
      }
  }

  annotation=as.data.frame(annotation[order(annotation$Cluster),])
  rownames(annotation)=targets[[RBP_name]]
  colnames(annotation)="Cluster"

  if (Not_DE==F) {
    annotation=subset(annotation,annotation$Cluster!="Not DE")
  }

  df=cbind(counts_rna,counts_ribo)
  df=df[rownames(annotation),]

  pheatmap(df,cluster_cols = F,cluster_rows = F,annotation_row = annotation )
}
