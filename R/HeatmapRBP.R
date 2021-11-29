#' @title HeatmapRBP
#'
#' @description Plot a heatmap of the pvalues of enrichment of each RBP in each gene group
#'
#' @param symbol res, rbp_lfc
#'
#' @return Plot a heatmap
#'
#' @examples HeatmapRBP(res=res,RBP_change=rbp_lfc)
#'
#' @export
#'
#'
#'
HeatmapRBP <-function(res=res,rbp_lfc=rbp_lfc)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  pvalues = data.frame(value1 = res[[1]]$padj)

  for (i in names(res)[-1]) {
    pvalues = cbind(pvalues, res[[i]]$padj)
  }

  pvalues = -log10(pvalues)
  pvalues[which(pvalues == Inf, arr.ind = T)] = 5
  rownames(pvalues) = res[[1]]$RBP
  colnames(pvalues) = names(res)
  pvalues = pvalues[rowSums(pvalues > -log10(0.05)) > 0, ]
  pvalues=trunc(pvalues)
  pvalues$status = NA
  pvalues$status[match(rbp_lfc$IDENTIFIER, rownames(pvalues))] = sign(rbp_lfc$FoldChange)
  pvalues$colors = ifelse(pvalues$status == 1, "dodgerblue3",
                          "darkorange1")
  pvalues = pvalues[!is.na(pvalues$colors), ]
  colnames(pvalues)=sub("_"," ",colnames(pvalues))

  colorRBP=as.character(pvalues$status)
  names(colorRBP)=rownames(pvalues)

  HT1=Heatmap(pvalues[,-c(9,10)],column_title = "RBP enrichment per gene group",
              col = colorRampPalette(c("oldlace", "darkred"))(5),name = "-log(P)",
              column_names_rot=45,show_row_names = T,width = unit(10,"cm"),
              heatmap_legend_param =list(at = 0:5))

  HT2=Heatmap(colorRBP,show_column_names=F,
              show_row_names = T,name = "RBP direction",
              col=c("dodgerblue3","darkorange1"),width = unit(0.5,"cm"),
              heatmap_legend_param = list(at=c("-1","1"),labels=c("Down","Up")))

  A=HT1+HT2
  return(A)

}


