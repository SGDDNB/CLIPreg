#' @title Plot_GO
#'
#' @description This function draws a network from CLIPreg output. Precise the number of RBP to consider with n
#' the thershold for size of node you want the GO done with th. The number of GO term to show per node with GO_to_show
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples
#'
#' @export
#'
#'
#'
Plot_GO_RBP <-function(rbp_of_interest="QKI",tpm_ribo = tpm_ribo,Targets=Targets,gene_groups=gene_groups,GO_to_show=3)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  Targets_RBP=Targets[[rbp_of_interest]]

  UpOrDown=rep("Down",nrow(gene_groups))
  UpOrDown[which(grepl("up",gene_groups$Gene_group,fixed = T))]="Up"
  gene_groups$Direction=UpOrDown

  RBP_t_up=subset(gene_groups$geneID,gene_groups$Direction=="Up")
  RBP_t_up=subset(RBP_t_up,RBP_t_up%in%Targets_RBP)

  RBP_t_down=subset(gene_groups$geneID,gene_groups$Direction=="Down")
  RBP_t_down=subset(RBP_t_up,RBP_t_down%in%Targets_RBP)


  selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")


  tpm_ribo=tpm_ribo[!duplicated(tpm_ribo$IDENTIFIER),]
  genes=tpm_ribo$IDENTIFIER


  GOterms=c()
  GOscore=c()
  NodeNames=c()
  ToSelect=rep(1,length(genes))
  names(ToSelect)=genes
  ToSelect[as.character(genes[tpm_ribo$geneID%in%RBP_t_down])]=0
  GOdata <- new("topGOdata",ontology="BP",allGenes=ToSelect,annot=annFUN.GO2genes,GO2genes=allGO2genes,geneSel=selection,
                  nodeSize=10)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  goEnrichment <- GenTable(GOdata, classicFisher = resultFisher, topNodes=GO_to_show)
  goEnrichment$classicFisher=as.numeric(goEnrichment$classicFisher)
  GOterms=c(GOterms,goEnrichment$Term)
  GOscore=c(GOscore,goEnrichment$classicFisher)
  NodeNames=c(NodeNames,rep(paste0(rbp_of_interest," down"),length(goEnrichment$GO.ID)))

  ToSelect=rep(1,length(genes))
  names(ToSelect)=genes
  ToSelect[as.character(genes[tpm_ribo$geneID%in%RBP_t_up])]=0
  GOdata <- new("topGOdata",ontology="BP",allGenes=ToSelect,annot=annFUN.GO2genes,GO2genes=allGO2genes,geneSel=selection,
                nodeSize=10)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  goEnrichment <- GenTable(GOdata, classicFisher = resultFisher, topNodes=GO_to_show)
  goEnrichment$classicFisher=as.numeric(goEnrichment$classicFisher)
  GOterms=c(GOterms,goEnrichment$Term)
  GOscore=c(GOscore,goEnrichment$classicFisher)
  NodeNames=c(NodeNames,rep(paste0(rbp_of_interest," up"),length(goEnrichment$GO.ID)))

  GOtermsNamed=paste0(NodeNames," : ",GOterms)
  df=data.frame(Term=GOtermsNamed,Pval=GOscore)

  ggplot(df, aes(x=Term, y=-log10(Pval),fill="1",color="1")) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("") +
    ylab("Enrichment") +
    ggtitle(paste0("GO of ",rbp_of_interest,"'s targets")) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(df$Pval)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
      axis.title=element_text(size=24, face="bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=18),  #Text size
      title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()+
    scale_color_manual(values ="red")+
    scale_fill_manual(values = "red")
}
