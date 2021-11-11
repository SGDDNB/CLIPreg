#' @title GetTarget
#'
#' @description Summarize RBPs in a list with their targets
#'
#' @param symbol folder, RBP_data, background
#'
#' @return List of RBPs each containing 1 vector of targets in form of characters
#'
#' @examples GetTarget(RBP_data=RBP_POSTAR, background=clusters$geneID)
#'
#' @export
#'
#'
#'
getTarget <-function(RBP_data=RBP_POSTAR, background=clusters$geneID)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  for (i in names(RBP_data)) {
    RBP_data[[i]]=subset(RBP_data[[i]],RBP_data[[i]]%in%background)
    }
  return(RBP_data)
}
