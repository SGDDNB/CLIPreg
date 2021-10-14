#' @title Draw_network
#'
#' @description This package draw a network from CLIPreg output
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples
#'
#' @export
#'
#'
#'
Draw_network <-function(rbp_lfc=rbp_lfc,res=res_both,clusters=clusters,n=5)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);


  RBPs=c()
  for (c in names(res)) {
    res[[c]]=res[[c]][res[[c]]$pval<0.05,]
    if (length(res[[c]][,1])==0) {
      res=res[-which(names(res)==c)]
    }
    RBPs=c(RBPs,res[[c]]$RBP)
  }
  RBPs=unique(RBPs)


  rbp_lfc=rbp_lfc[RBPs,]
  rbp_lfc$max=apply(abs(rbp_lfc), 1, max)
  rbp_lfc=rbp_lfc[order(-rbp_lfc$max),]
  rbp_lfc=rbp_lfc[1:n,]

  RBP_kept=rownames(rbp_lfc)

  clusters=clusters[clusters$Cluster%in%names(res),]
  adjacency_matrix=as.data.frame(matrix(0,nrow = n+nrow(clusters),ncol = n+nrow(clusters)))
  colnames(adjacency_matrix)=c(RBP_kept,clusters$geneID)
  rownames(adjacency_matrix)=c(RBP_kept,clusters$geneID)
  for (r in RBP_kept) {
    for (j in clusters$geneID) {
      if (j%in%Targets[[r]]) {
        adjacency_matrix[r,j]=1
        adjacency_matrix[j,r]=1
      }
    }
  }

  adjacency_matrix=adjacency_matrix[rowSums(adjacency_matrix)>0,colSums(adjacency_matrix)>0]


  node_color=rep(0,nrow(adjacency_matrix))
  node_color=clusters$Cluster[match(rownames(adjacency_matrix),clusters$geneID)]
  node_color[1:n]="RBP"

  node_size=rep(0.5,nrow(adjacency_matrix))
  node_size[1:n]=10
  ggnet2(adjacency_matrix,size = node_size,color = node_color,label = RBP_kept,
         palette="Set1")+ guides(size = FALSE,color=guide_legend(title="Clusters"))


}
