#' @title run_CLIPreg
#' @description
#' @param symbol
#' @return NULL
#' @examples
#' @export
#'
run_CLIPreg=function(input_data=Example, is.example=T)
{
  gene_groups=input_data$gene_groups
  tpm_ribo=input_data$tpm_ribo
  ribo_lfc=input_data$ribo_lfc
  RBP_ENCODE=input_data$RBP_ENCODE
  RBP_POSTAR=input_data$RBP_POSTAR

  if (is.example==T) {
    res_Encode=CLIPreg::res_Encode
    res_Postar=CLIPreg::res_Postar
  } else {
    res_Encode=CLIPReg_V4(RBP_data=RBP_ENCODE,gene_groups=gene_groups)
    res_Postar=CLIPReg_V4(RBP_data=RBP_POSTAR,gene_groups=gene_groups)
  }

  Targets=combine_targets(RBP_list1=RBP_ENCODE,RBP_list2=RBP_POSTAR,background=gene_groups$geneID)
  res=CLIPreg::combine(res1=res_Encode,res2=res_Postar)
  rbp_lfc=rbp_change(res=res,ribo_lfc=ribo_lfc)
  res=cure_res(res=res,rbp_lfc=rbp_lfc)

  L=list(gene_groups=gene_groups,
         tpm_ribo=tpm_ribo,
         rbp_lfc=rbp_lfc,
         ribo_lfc=ribo_lfc,
         Targets=Targets,
         res=res)
  return(L)
}
