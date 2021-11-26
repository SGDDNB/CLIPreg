#' @title BubbleRBPs
#'
#' @description Use gene groups from DeltaTE output to generate enrichement plot of RBPs
#'
#' @param symbol RBP_res, gene_groups, rbp_lfc
#'
#' @return NULL. Generates a bubble plot
#'
#' @examples BubbleRBPs(RBP_res=res,gene_groups=gene_groups,rbp_lfc=rbp_lfc)
#'
#' @export
#'
#'
#'
BubbleRBPs <-function(res=res,gene_groups=gene_groups,rbp_lfc=rbp_lfc)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  perc_in <- function(a, b)
    {
    A=length(a)
    return(as.numeric(table(a%in%b)[2]*100/A))
  }

  sig_res=res
  for (cl in names(res)) {
    sig_res[[cl]]=sig_res[[cl]][sig_res[[cl]]$pval<0.05,]
  }


  # Set up the vectors
  RNA <- c("Trancription\nup","Transcription\ndown","Transcription\nunchanged")
  RIBO <- c("TE unchanged","TE down","TE up")

  # Create the data frame
  df <- expand.grid(RNA, RIBO)
  df$value=c(mean(sig_res$forwarded_up$z),
             mean(sig_res$forwarded_down$z),
             NA,
             mean(sig_res$buffered_down$z),
             mean(sig_res$intensified_down$z),
             mean(sig_res$exclusive_down$z),
             mean(sig_res$intensified_up$z),
             mean(sig_res$buffered_up$z),
             mean(sig_res$exclusive_up$z))
  df$value=round(df$value,1)

  df$Regulators=c(perc_in(gene_groups$geneID[gene_groups$Gene_group=="forwarded_up"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="forwarded_down"],rbp_lfc$geneID),
                  NA,
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="buffered_down"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="intensified_down"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="exclusive_down"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="intensified_up"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="buffered_up"],rbp_lfc$geneID),
                  perc_in(gene_groups$geneID[gene_groups$Gene_group=="exclusive_up"],rbp_lfc$geneID))

  df$Regulators[is.na(df$Regulators)]=0
  df$label=c("forwarded_up","forwarded_down",NA,"buffered_down","intensified_down",
             "exclusive_down","intensified_up","buffered_up","exclusive_up")

  #Plot the Data
  ggplot(df, aes(Var1, Var2,label=value)) + geom_point(aes(size = value, colour = Regulators)) +
    scale_size_continuous(range = c(1,10))+
    scale_colour_gradient(low = "white", high = "darkblue")+
    theme_bw() + xlab("") + ylab("")+  geom_text(hjust=0.5, vjust=-2)+
    ggtitle("Mean z-score per gene group")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    guides(col=guide_legend(title="Percentage of RBP"),
           size=guide_legend(title="mean z-score"))
}
