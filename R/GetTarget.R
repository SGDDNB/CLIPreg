#' @title GetTarget
#'
#' @description Summarize RBPs in a list with their targets
#'
#' @param symbol folder, RBP_data, background
#'
#' @return List of RBPs each containing 1 vector of targets in form of characters
#'
#' @examples GetTarget(folder="Folder",RBP_data="rbp_gene_postar.txt",background="bg.txt")
#'
#' @export
#'
#'
#'
getTarget <-function(folder="Folder",RBP_data="rbp_gene_postar.txt",background=clusters$geneID)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  rbp=fread(paste0(folder,"/",RBP_data))
  rbp=rbp[rbp$V3%in%background,]

  RBP=list()
  for (i in unique(rbp$V1)) {
    RBP[[i]]=rbp$V3[rbp$V1==i]
    RBP[[i]]=unique(RBP[[i]])
  }
  return(RBP)
}
