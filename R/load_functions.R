#' @title load_clusters
#' @description Load the cluster file from user
#' @param symbol cluster_file
#' @return NULL
#' @examples
#' @export
#'
load_clusters <-function(cluster_file=cluster_file) {
  clusters=as.data.frame(fread(cluster_file))
  return(clusters)
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
#' @description Load the cluster file from user
#' @param symbol cluster_file
#' @return NULL
#' @examples
#' @export
#'
combine <- function(res1=res_Encode,res2=res_Postar){
  res_both=res2
  for (i in 1:nrow(res1[[1]])) {
    r=res1[[1]]$RBP[i]
    if (r%in%res_both[[1]]$RBP) {
      for (n in names(res_both)) {
        for (j in 2:6) {
          res_both[[n]][which(res_both[[1]]$RBP==r),j]=(res_both[[n]][which(res_both[[1]]$RBP==r),j]+res1[[n]][i,j])/2
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
    to_keep=c(to_keep,res_both[[n]]$RBP[res_both[[n]]$pval<0.05])
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
}


