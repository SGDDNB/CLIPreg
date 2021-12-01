#' @title CLIPreg
#'
#' @description Calculates the enrichment of each RBP for each gene group
#'
#' @param symbol folder, RBP_data, gene_groups
#'
#' @return List of gene_groups containing dataframe with all the enrichment score for each RBP
#'
#' @examples CLIPreg(RBP_data=RBP_POSTAR,gene_groups=gene_groups)
#' @export
#'
#'
#'
CLIPreg <-function(RBP_data=RBP_POSTAR,gene_groups=gene_groups)
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

  rbp_names=names(RBP_data)
  categories=unique(gene_groups$Gene_group)
  out=list()
  for (cl in categories) {
    print(cl)
    RBP=RBP_data
    gene_group=gene_groups$geneID[gene_groups$Gene_group==cl]
    nb_genes=length(gene_group)

    # Prepare the output table
    overlap=data.frame(RBP=rbp_names,real_overlap=0,simulated_overlap_mean=0,simulated_overlap_sd=0,z=0,pval=0,padj=0)

    # Get the real overlap
    for (r in 1:nrow(overlap)) {
      rbp_r=overlap$RBP[r]
      overlap$real_overlap[r]=sum(gene_group%fin%RBP[[rbp_r]])
    }


    # Function for the loop in parallel
    Get_simulation=function(gene_groups,nb_genes,simulations,rbp_names,RBP,iterations=1000){
      for(j in 1:iterations) {
        bg_shuffled=sample(gene_groups$geneID,nb_genes)
        for (r in rbp_names) {
          simulations[[r]]=c(simulations[[r]],sum(bg_shuffled%fin%RBP[[r]]))
        }
      }
      return(simulations)
    }

    simulations=list()
    # Get simulated data
    finalMatrix <- foreach(i=1:100, .combine='c') %dopar% {
      tempMatrix =  Get_simulation(gene_groups,nb_genes,simulations,rbp_names,RBP,iterations = 1000)#calling a function
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
      count_lower=sum(simulations[[i]]>overlap$real_overlap[i])
      if (count_lower=0) {
        count_lower=1
      }
      overlap$pval[i]=count_lower/100000
    }
    overlap$padj=p.adjust(overlap$pval,method="fdr")
    out[[cl]]=overlap


  }
  stopImplicitCluster()
  return(out)
}
