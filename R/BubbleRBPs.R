#' @title BubbleRBPs
#'
#' @description Use clusters from DeltaTE output to generate enrichement plot of RBPs
#'
#' @param symbol RBP_res, clusters
#'
#' @return NULL. Generates a bubble plot
#'
#' @examples BubbleRBPs(RBP_res=res,clusters=clusters)
#'
#' @export
#'
#'
#'
BubbleRBPs <-function(RBP_res=res,clusters=clusters)
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

  for (c in names(res)) {
    res[[c]]=res[[c]][res[[c]]$pval<0.05,]
  }


  # Set up the vectors
  RNA <- c("Trancription up","Transcription down","Transcription unchanged")
  RIBO <- c("TE unchanged","TE down","TE up")

  # Create the data frame
  df <- expand.grid(RNA, RIBO)
  df$value=c(mean(res$forwarded_down$z),
             mean(res$forwarded_up$z),
             NA,
             mean(res$intensified_down$z),
             mean(res$buffered_down$z),
             mean(res$exclusive_down$z),
             mean(res$buffered_up$z),
             mean(res$intensified_up$z),
             mean(res$exclusive_up$z))

  df$Regulators=c(perc_in(clusters$IDENTIFIER[clusters$Cluster=="forwarded_down"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="forwarded_up"],RBP_res[[1]][,1]),
                    NA,
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="buffered_down"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="intensified_down"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="exclusive_down"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="intensified_up"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="buffered_up"],RBP_res[[1]][,1]),
                    perc_in(clusters$IDENTIFIER[clusters$Cluster=="exclusive_up"],RBP_res[[1]][,1]))

  df$Regulators[is.na(df$Regulators)]=0
  df$label=c("forwarded_down","forwarded_up",NA,"buffered_down","intensified_down",
             "exclusive_down","intensified_up","buffered_up","exclusive_up")

  #Plot the Data
  ggplot(df, aes(Var1, Var2,label=label)) + geom_point(aes(size = value, colour = Regulators)) +
    scale_size_continuous(range = c(4, 8))+
    scale_colour_gradient(low = "white", high = "darkblue")+
    theme_bw() + xlab("") + ylab("")+  geom_text(hjust=0.5, vjust=-3)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    guides(col=guide_legend(title="Percentage of RBP"),
           size=guide_legend(title="mean z-score"))
}
