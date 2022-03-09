#' @title load_gene_groups
#' @description Load the gene group file from user
#' @param symbol gene_groups_file
#' @return NULL
#' @examples
#' @export
#'
load_gene_groups <-function(gene_groups_file=gene_groups_file) {
  gene_groups=as.data.frame(fread(gene_groups_file))
  return(gene_groups)
}

#' @title load_ribo_lfc
#' @description Load the ribo_lfc file from user
#' @param symbol ribo_lfc
#' @return NULL
#' @examples
#' @export
#'
load_ribo_lfc <-function(ribo_lfc_file=ribo_lfc_file) {
  ribo_lfc=as.data.frame(fread(ribo_lfc_file))
  return(ribo_lfc)
}

#' @title load_ribo_tpm
#' @description Load the ribo_tpm file from user
#' @param symbol ribo_tpm
#' @return NULL
#' @examples
#' @export
#'
load_ribo_tpm <-function(ribo_tpm_file=ribo_tpm_file) {
  ribo_tpm=as.data.frame(fread(ribo_tpm_file))
  return(ribo_tpm)
}


#' @title combine
#' @description Load the gene_groups file from user
#' @param symbol gene_groups_file
#' @return NULL
#' @examples
#' @export
#'
combine <- function(res1=res_Encode,res2=res_Postar,FDR=0.1){
  if (!exists("res2")) {
    res2=res1
  }
  res_both=res2
  for (i in 1:nrow(res1[[1]])) {
    r=res1[[1]]$RBP[i]
    if (r%in%res_both[[1]]$RBP) {
      for (n in names(res_both)) {
        if (res_both[[n]][which(res_both[[1]]$RBP==r),7]>res1[[n]][i,7]) {
          for (j in 2:7) {
            res_both[[n]][which(res_both[[1]]$RBP==r),j]=res1[[n]][i,j]
          }
        }
      }
    } else {
      for (n in names(res_both)) {
        res_both[[n]]=rbind(res_both[[n]],res1[[n]][i,])
      }
    }
  }

  to_keep=c()
  for (n in names(res_both)) {
    to_keep=c(to_keep,res_both[[n]]$RBP[res_both[[n]]$padj<FDR])
  }
  to_keep=unique(to_keep)

  for (n in names(res_both)) {
    res_both[[n]]=res_both[[n]][res_both[[n]]$RBP%in%to_keep,]
  }
  return(res_both)

}

#' @title rbp_change
#' @description Subset ribo_lfc to keep only rbp_lfc
#' @param symbol res, ribo_lfc
#' @return NULL
#' @examples
#' @export
#'
rbp_change=function(res=res,ribo_lfc=ribo_lfc){
  options(stringsAsFactors=FALSE)
  rbp_lfc=ribo_lfc[ribo_lfc$IDENTIFIER%in%res[[1]]$RBP,1:3]
  rbp_lfc=rbp_lfc[!duplicated(rbp_lfc$IDENTIFIER),]
  rownames(rbp_lfc)=rbp_lfc$IDENTIFIER
  return(rbp_lfc)
}

#' @title cure_res
#' @description Remove RBPs that are not in rbp_lfc from res
#' @param symbol res, ribo_lfc
#' @return NULL
#' @examples
#' @export
#'
cure_res=function(res=res,rbp_lfc=rbp_lfc){
  for (n in 1:length(names(res))) {
    res[[n]]=res[[n]][res[[n]]$RBP%in%as.character(rbp_lfc$IDENTIFIER),]
  }
  return(res)
}

#' @title combine_targets
#' @description Combine targets from POSTAR and ENCODE and keep only genes from background
#' @param symbol RBP_list1, RBP_list2, background
#' @return NULL
#' @examples
#' @export
#'
combine_targets=function(RBP_list1=RBP_ENCODE,RBP_list2=RBP_POSTAR,background=gene_groups$geneID){
  TargetsEncode=GetTarget(RBP_data = RBP_list1,background = background)
  if (exists("RBP_list2")) {
    TargetsPOSTAR=GetTarget(RBP_data = RBP_list2,background = background)
  } else {TargetsPOSTAR=c()}
  Targets=c(TargetsEncode,TargetsPOSTAR)
  Targets=split(unlist(Targets, use.names = FALSE), rep(names(Targets), lengths(Targets)))
  for (t in names(Targets)) {
    Targets[[t]]=unique(Targets[[t]])
  }
  return(Targets)
}






