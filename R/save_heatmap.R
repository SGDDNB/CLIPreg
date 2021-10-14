#' @title save_heatmap
#'
#' @description Generates a pdf from the heatmap function
#'
#' @param symbol res, RBP_change,location
#'
#' @return pdf is created at the location
#'
#' @examples save_heatmap(res=res,RBP_change=rbp_lfc,location=location)
#'
#' @export
#'
#'
#'
save_heatmap <-function(res=res,RBP_change=rbp_lfc,location=location)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  RBP_change=rowMeans(RBP_change)

  pvalues=data.frame(value1=res[[1]]$pval)
  for (i in names(res)[-1]) {
    pvalues=cbind(pvalues,res[[i]]$pval)
  }
  pvalues=-log10(pvalues)
  pvalues[which(pvalues==Inf,arr.ind = T)]=5
  rownames(pvalues)=res[[1]]$RBP
  colnames(pvalues)=names(res)
  pvalues=pvalues[rowSums(pvalues>-log10(0.05))>0,]

  pvalues$status=sign(RBP_change[rownames(pvalues)])
  pvalues$colors=ifelse(pvalues$status==1,"dodgerblue3","darkorange1")
  pvalues=pvalues[!is.na(pvalues$colors),]
  e=pheatmap(pvalues[,1:(ncol(pvalues)-2)],breaks =0:5,color = colorRampPalette(c("oldlace","darkred"))(5),
             angle_col = 45,legend = T,border_color = FALSE)#,fontsize_row = 6)
  cols=pvalues[order(match(rownames(pvalues), e$gtable$grobs[[5]]$label)), ]$colors
  e$gtable$grobs[[5]]$gp=gpar(col=cols)

  n=nrow(pvalues)

  pdf(location,8,3+n*0.15)
  e
  dev.off()

}


