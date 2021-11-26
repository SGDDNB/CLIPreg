#' @title Load_example
#' @description
#' @param symbol
#' @return NULL
#' @examples
#' @export
#'
Load_example <-function(){
  L=list(clusters=CLIPreg::clusters,
         tpm_ribo=CLIPreg::tpm_ribo,
         ribo_lfc=CLIPreg::ribo_lfc,
         RBP_ENCODE=CLIPreg::RBP_ENCODE,
         RBP_POSTAR=CLIPreg::RBP_POSTAR)
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
  L=list(clusters=clusters,
         tpm_ribo=tpm_ribo,
         ribo_lfc=ribo_lfc,
         RBP_ENCODE=CLIPreg::RBP_ENCODE,
         RBP_POSTAR=CLIPreg::RBP_POSTAR)
  return(L)
}
