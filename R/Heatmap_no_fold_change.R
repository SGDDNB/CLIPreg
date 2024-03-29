#' @title Heatmap_no_fold_change
#'
#' @description Plot a heatmap of the pvalues of enrichment of each miRNA in each gene group
#'
#' @param symbol res
#'
#' @return Plot a heatmap
#'
#' @examples Heatmap_no_fold_change(res=res)
#'
#' @export
#'
#'
#'
Heatmap_no_fold_change <-function(res=res)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  pvalues = data.frame(value1 = res[[1]]$padj)

  for (i in names(res)[-1]) {
    pvalues = cbind(pvalues, res[[i]]$padj)
  }

  pvalues[which(pvalues>0.05,arr.ind = T)]=1

  pvalues = -log10(pvalues)
  pvalues[which(pvalues == Inf, arr.ind = T)] = 5
  rownames(pvalues) = res[[1]]$RBP
  colnames(pvalues) = names(res)
  pvalues = pvalues[rowSums(pvalues > -log10(0.05)) > 0, ]
  pvalues=round(pvalues)
  colnames(pvalues)=sub("_"," ",colnames(pvalues))


  HT1=Heatmap(pvalues,column_title = "Regulators enrichment per gene group",
              col = colorRampPalette(c("oldlace", "darkred"))(max(pvalues)),name = "-log(FDR)",
              column_names_rot=45,show_row_names = T,width = unit(ncol(pvalues),"cm"),
              heatmap_legend_param =list(at = 0:max(pvalues)))

  return(HT1)

}


