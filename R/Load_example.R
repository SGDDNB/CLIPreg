#' @title Load_example
#' @description
#' @param symbol
#' @return NULL
#' @examples
#' @export
#'
Load_example <-function(){
  clusters=data("clusters")
  tpm_ribo=data("tpm_ribo")
  ribo_lfc=data("ribo_lfc")
  RBP_ENCODE=data("RBP_ENCODE")
  RBP_POSTAR=data("RBP_POSTAR")
  L=list(clusters,tpm_ribo,ribo_lfc,RBP_ENCODE,RBP_POSTAR)
}
