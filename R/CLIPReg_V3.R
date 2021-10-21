#' @title CLIPReg_V3
#'
#' @description Calculates the enrichment of each RBP for each cluster
#'
#' @param symbol folder, RBP_data, cluster
#'
#' @return List of clusters containing dataframe with all the enrichment score for each RBP
#'
#' @examples CLIPreg_V3(folder="Folder",RBP_data="rbp_gene_postar.txt",cluster="Tx_down.txt")
#'
#' @export
#'
#'
#'
CLIPReg_V3 <-function(folder="Folder",RBP_data="rbp_gene_postar.txt",cluster="Tx_down.txt")
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  numCores <- detectCores()
  registerDoParallel(numCores-1)  # use multicore, set to the number of our cores


  # Function fastmatch
  `%fin%` <- function(x, table) {
    stopifnot(require(fastmatch))
    fmatch(x, table, nomatch = 0L) > 0L
  }

  # open files
  rbp=fread(paste0(folder,"/",RBP_data))
  clusters=fread(paste0(folder,"/",cluster),header = T)
  #bg=fread(paste0(folder,"/",background),header = F)

  rbp_names=unique(rbp$V1)
  categories=unique(clusters$Cluster)
  out=list()
  for (cl in categories) {
    print(cl)
    RBP=list()
    cluster=clusters$geneID[clusters$Cluster==cl]
    nb_genes=length(cluster)

    # Build the list of RBP and their targets
    for (r in rbp_names) {
      RBP[[r]]=rbp$V3[rbp$V1==r]
    }

    # Prepare the output table
    overlap=data.frame(RBP=rbp_names,real_overlap=0,simulated_overlap_mean=0,simulated_overlap_sd=0,z=0,pval=0)

    # Get the real overlap
    for (r in 1:nrow(overlap)) {
      rbp_r=overlap$RBP[r]
      overlap$real_overlap[r]=sum(cluster%fin%RBP[[rbp_r]])
    }


    # Function for the loop in parallel
    Get_simulation=function(clusters,nb_genes,simulations,rbp_names,RBP,iterations=1000){
      for(j in 1:iterations) {
        bg_shuffled=sample(clusters$geneID,nb_genes)
        for (r in rbp_names) {
          simulations[[r]]=c(simulations[[r]],sum(bg_shuffled%fin%RBP[[r]]))
        }
      }
      return(simulations)
    }

    simulations=list()
    # Get simulated data
    finalMatrix <- foreach(i=1:100, .combine='c') %dopar% {
      tempMatrix =  Get_simulation(clusters,nb_genes,simulations,rbp_names,RBP,iterations = 1000)#calling a function
      tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }

    simulations=split(unlist(finalMatrix, use.names = FALSE), rep(names(finalMatrix), lengths(finalMatrix)))


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

    out[[cl]]=overlap


  }
  stopImplicitCluster()
  return(out)
}
