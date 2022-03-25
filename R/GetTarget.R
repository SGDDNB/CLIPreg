#' @title GetTarget
#'
#' @description Summarize RBPs in a list with their targets
#'
#' @param symbol folder, RBP_data, background
#'
#' @return List of RBPs each containing 1 vector of targets in form of characters
#'
#' @examples GetTarget(RBP_data=RBP_POSTAR, background=gene_groups$geneID)
#'
#' @export
#'
#'
#'
GetTarget <-function(Regulator_data=RBP_POSTAR, background=gene_groups$geneID)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  for (i in names(Regulator_data)) {
    Regulator_data[[i]]=subset(Regulator_data[[i]],Regulator_data[[i]]%in%background)
    }
  return(Regulator_data)
}
