#' @title Load_example
#' @description
#' @param symbol
#' @return NULL
#' @examples
#' @export
#'
Load_example <-function(){
  data("clusters")
  data("tpm_ribo")
  data("ribo_lfc")
  data("RBP_ENCODE")
  data("RBP_POSTAR")

  L=list(clusters=clusters,
         tpm_ribo=tpm_ribo,
         ribo_lfc=ribo_lfc,
         RBP_ENCODE=RBP_ENCODE,
         RBP_POSTAR=RBP_POSTAR)
  return(L)
}

#' @title Load_input_files
#' @description
#' @param symbol
#' @return NULL
#' @examples
#' @export
#'
Load_input_files <-function(folder="Path/to/folder"){
  clusters=fread(paste0(folder,"/clusters.txt"))
  tpm_ribo=fread(paste0(folder,"/tpm_ribo.txt"))
  ribo_lfc=fread(paste0(folder,"/ribo_lfc.txt"))
  data("RBP_ENCODE")
  data("RBP_POSTAR")
  L=list(clusters=clusters,
         tpm_ribo=tpm_ribo,
         ribo_lfc=ribo_lfc,
         RBP_ENCODE=RBP_ENCODE,
         RBP_POSTAR=RBP_POSTAR)
  return(L)
}
