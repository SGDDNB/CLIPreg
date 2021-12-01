#' @title Draw_network_by_group
#'
#' @description This function draw a network from CLIPreg output. It plots the n most changing RBP and their targets group by TE category
#'
#' @param symbol rbp_lfc, res, Targets, gene_groups, n
#'
#' @return Plot network
#'
#' @examples Draw_network(rbp_lfc=rbp_lfc,res=res_both,gene_groups=gene_groups,n=5)
#'
#' @export
#'
#'
#'
Draw_network_by_group <-function(rbp_lfc=rbp_lfc,res=res,Targets=Targets,gene_groups=gene_groups,n=5,forwarded=F)
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
    res[[c]]=res[[c]][res[[c]]$padj<0.1,]
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
  new_adj=as.data.frame(matrix(0,nrow = n+length(all_types),ncol = n+length(all_types)))
  colnames(new_adj)=c(RBP_kept,all_types)
  rownames(new_adj)=c(RBP_kept,all_types)

  UpOrDown=rep("Down",nrow(gene_groups))
  UpOrDown[which(grepl("up",gene_groups$Gene_group,fixed = T))]="Up"
  gene_groups$Direction=UpOrDown

  adj_up=new_adj
  adj_down=new_adj

  #adj_Down
  size_down=rep(1,n)
  size_up=c()
  for (i in all_types) {
    genes_in_group=names(which(type==i))
    genes_Down=genes_in_group[gene_groups$Direction[gene_groups$geneID==genes_in_group]=="Down"]
    genes_Up=genes_in_group[!(genes_in_group%in%genes_Down)]
    RBP_concerned=str_locate_all(pattern="1",i)[[1]][,1]
    if (length(RBP_concerned)>0) {
        adj_down[i,RBP_concerned]=1
        adj_down[RBP_concerned,i]=1
        adj_up[i,RBP_concerned]=1
        adj_up[RBP_concerned,i]=1
        size_down=c(size_down,length(genes_Down))
        size_up=c(size_up,length(genes_Up))
    }
  }
  rownames(adj_up)=paste0(rownames(adj_up),"_up")
  colnames(adj_up)=paste0(colnames(adj_up),"_up")
  sizes=c(size_down,size_up)

  new_adj=cbind(adj_down,matrix(0,nrow=nrow(adj_down),ncol=ncol(adj_up)-n))
  new_adj2=matrix(0,nrow=ncol(adj_up)-n,ncol=ncol(adj_down)+ncol(adj_up)-n)
  colnames(new_adj2)=colnames(new_adj)
  new_adj=rbind(new_adj,new_adj2)
  new_adj[1:n,(ncol(adj_down)+1):ncol(new_adj)]=adj_up[1:n,-c(1:n)]
  new_adj[(nrow(adj_down)+1):nrow(new_adj),1:n]=adj_up[-c(1:n),1:n]
  rownames(new_adj)[(nrow(adj_down)+1):nrow(new_adj)]=rownames(adj_up)[-(1:n)]
  names(sizes)=rownames(new_adj)
  sizes=sizes[sizes>0]
  new_adj=new_adj[rownames(new_adj)%in%names(sizes),rownames(new_adj)%in%names(sizes)]


  node_color=rep("TE down",nrow(new_adj))
  node_color[!is.na(str_locate(pattern="up",rownames(new_adj))[,1])]="TE up"
  node_color[1:n]="RBP"


  forscale=c(min(sizes),
             sort(sizes)[round(length(sizes)/2)],
             sort(sizes)[round(length(sizes)/4)*3],
             max(sizes))
  sizes[1:n]=1.5*max(sizes)

  ggnet2(new_adj,size = sizes,color = node_color,label = c(RBP_kept,paste0("ID",1:(length(sizes)-n))),
         palette="Set1",vjust=-0.6,edge.color = "gray80")+ guides(color=guide_legend(title="Regulation"))+
    scale_size_discrete(breaks=forscale)+labs(size="Node size")+
    scale_color_manual(values=c("red","blue","orange"))


}
