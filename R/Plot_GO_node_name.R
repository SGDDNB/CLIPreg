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
Plot_GO_node_name <-function(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,
                             tpm_ribo = tpm_ribo,Nodes_to_keep=c(19,15),GO_to_show=3,forwarded=F)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  options(stringsAsFactors=FALSE);

  if (forwarded==F) {
    fw=which(grepl("forwarded",names(res)))
    if (length(fw)>0) {
      res=res[-fw]
    }
  }

  RBPs=c()
  for (c in names(res)) {
    res[[c]]=res[[c]][res[[c]]$pval<0.05,]
    if (length(res[[c]][,1])==0) {
      res=res[-which(names(res)==c)]
    }
    RBPs=c(RBPs,res[[c]]$RBP)
  }
  RBPs=unique(RBPs)
  rm(c)

  rbp_lfc=rbp_lfc[RBPs,]
  rbp_lfc=rbp_lfc[order(-abs(rbp_lfc$FoldChange)),]
  rbp_lfc=rbp_lfc[1:n,]

  RBP_kept=rownames(rbp_lfc)

  gene_groups=gene_groups[gene_groups$Gene_group%in%names(res),]
  adjacency_matrix=as.data.frame(matrix(0,nrow = n+nrow(gene_groups),ncol = n+nrow(gene_groups)))
  colnames(adjacency_matrix)=c(RBP_kept,gene_groups$geneID)
  rownames(adjacency_matrix)=c(RBP_kept,gene_groups$geneID)
  for (r in RBP_kept) {
    for (j in gene_groups$geneID) {
      if (j%in%Targets[[r]]) {
        adjacency_matrix[r,j]=1
        adjacency_matrix[j,r]=1
      }
    }
  }

  adjacency_matrix=adjacency_matrix[rowSums(adjacency_matrix)>0,colSums(adjacency_matrix)>0]

  type=apply( adjacency_matrix[,1:n] , 1 , paste , collapse = "" )

  all_types=unique(type)[-1]

  UpOrDown=rep("Down",nrow(gene_groups))
  UpOrDown[which(grepl("up",gene_groups$Gene_group,fixed = T))]="Up"
  gene_groups$Direction=UpOrDown

  #adj_Down
  size_down=rep(1,n)
  size_up=c()
  TypesList=list()
  TypesList_up=list()
  for (i in all_types) {
    genes_in_group=names(which(type==i))
    genes_Down=genes_in_group[gene_groups$Direction[gene_groups$geneID==genes_in_group]=="Down"]
    genes_Up=genes_in_group[!(genes_in_group%in%genes_Down)]
    TypesList[[i]]=genes_Down
    TypesList_up[[i]]=genes_Up
    RBP_concerned=str_locate_all(pattern="1",i)[[1]][,1]
    if (length(RBP_concerned)>0) {
      size_down=c(size_down,length(genes_Down))
      size_up=c(size_up,length(genes_Up))
    }
  }
  names(size_down)=c(RBP_kept,all_types)
  names(size_up)=paste0(names(size_down)[-c(1:n)],"_up")
  sizes=c(size_down,size_up)

  sizes=sizes[sizes>0]

  names(TypesList_up)=paste0(names(TypesList_up),"_up")
  TypesListAll=c(TypesList,TypesList_up)
  color_nodes=rep("blue",length(TypesList))
  color_nodes=c(color_nodes,rep("orange",length(TypesList_up)))
  names(color_nodes)=names(TypesListAll)
  TypesListAll=subset(TypesListAll,names(TypesListAll)%in%names(sizes))
  color_nodes=color_nodes[names(color_nodes)%in%names(sizes)]

  sizes[1:n]=1.5*max(sizes)


  selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")


  tpm_ribo=tpm_ribo[!duplicated(tpm_ribo$IDENTIFIER),]
  genes=tpm_ribo$IDENTIFIER

  names(TypesListAll)=paste0("ID",1:length(TypesListAll))
  names(color_nodes)=names(TypesListAll)
  TypesToKeep=TypesListAll[paste0("ID",Nodes_to_keep)]

  # Nodes_to_keep=c(19,16)
  GOterms=c()
  GOscore=c()
  NodeNames=c()
  GOcolor=c()
  for (j in 1:length(TypesToKeep)) {
    i=names(TypesToKeep)[j]
    ToSelect=rep(1,length(genes))
    names(ToSelect)=genes
    ToSelect[as.character(genes[tpm_ribo$geneID%in%TypesListAll[[i]]])]=0
    GOdata <- new("topGOdata",ontology="BP",allGenes=ToSelect,annot=annFUN.GO2genes,GO2genes=allGO2genes,geneSel=selection,
                  nodeSize=10)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    goEnrichment <- GenTable(GOdata, classicFisher = resultFisher, topNodes=GO_to_show)
    goEnrichment$classicFisher=as.numeric(goEnrichment$classicFisher)
    GOterms=c(GOterms,goEnrichment$Term)
    GOscore=c(GOscore,goEnrichment$classicFisher)
    NodeNames=c(NodeNames,rep(Nodes_to_keep[j],length(goEnrichment$GO.ID)))
    GOcolor=rep(color_nodes[names(TypesToKeep[j])],length(goEnrichment$GO.ID))
  }




  GOtermsNamed=paste0("ID",NodeNames," : ",GOterms)
  df=data.frame(Term=GOtermsNamed,Pval=GOscore,Color=GOcolor)

  ggplot(df, aes(x=Term, y=-log10(Pval),fill=Color,colour=Color)) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("") +
    ylab("Enrichment") +
    ggtitle("GO of selected nodes") +
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
    #guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()+
    scale_color_manual(values = unique(GOcolor))+
    scale_fill_manual(values = unique(GOcolor))
}
