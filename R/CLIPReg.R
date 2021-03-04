#' @title CLIPReg
#'
#' @description This package identifies...
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples RbpReg('AAPL')
#'
#' @export RbpReg
#'
#'
#'
CLIPReg <-function(folder="Folder",RBP_data="rbp_gene_postar.txt",background="bg.txt",cluster="Tx_down.txt",iterations=10)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  library(data.table)
  library(fastmatch)

  # Function fastmatch
  `%fin%` <- function(x, table) {
    stopifnot(require(fastmatch))
    fmatch(x, table, nomatch = 0L) > 0L
  }

  # open files
  bg=fread(paste0(folder,"/",background),header = F)
  rbp=fread(paste0(folder,"/",RBP_data))
  cluster=fread(paste0(folder,"/",cluster),header = F)

  rbp_names=unique(rbp$V1)
  RBP=list()
  nb_genes=nrow(cluster)


  # Build the list of RBP and their targets
  for (r in rbp_names) {
    RBP[[r]]=rbp$V3[rbp$V1==r]
  }

  # Prepare the output table
  overlap=data.frame(RBP=rbp_names,real_overlap=0,simulated_overlap_mean=0,simulated_overlap_sd=0,z=0,pval=0)

  # Get the real overlap
  for (r in 1:nrow(overlap)) {
    rbp_r=overlap$RBP[r]
    overlap$real_overlap[r]=sum(cluster$V1%fin%RBP[[rbp_r]])
  }


  # Get simulated data
  simulations=list()
  for (i in 1:iterations) {
    bg_shuffled=sample(bg$V1,nb_genes)
    for (r in rbp_names) {
      simulations[[r]]=c(simulations[[r]],sum(bg_shuffled%fin%RBP[[r]]))
    }
  }

  # summarize the simulated data
  overlap$simulated_overlap_mean=sapply(simulations,mean)
  overlap$simulated_overlap_sd=sapply(simulations,sd)

  # Calculate z-score and pvalue
  z=(overlap$real_overlap-overlap$simulated_overlap_mean)/overlap$simulated_overlap_sd
  overlap$z=z

  for (i in 1:nrow(overlap)) {
    rbp_i=overlap$RBP[i]
    count_lower=sum(simulations[[i]]>overlap$real_overlap[i])
    overlap$pval[i]=count_lower/iterations
  }
  return(overlap)
}
