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


#' @title load_clusters
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
















