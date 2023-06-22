
## homemade GO enrichment function

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(RColorBrewer)
library(org.Mm.eg.db,verbose=F,quietly=T)
#library(HGNChelper) #https://rdrr.io/cran/HGNChelper/man/mouse.table.html
#data("mouse.table", package="HGNChelper")
#head(mouse.table)
#length(genes)
#sum(genes %in% mouse.table$Symbol)
#tmp=mouse.table[mouse.table$Symbol %in% genes,]

GOenrich<-function(test.genes,category='BP',topn=5,cutoff=0.7){
  gene.df <- clusterProfiler::bitr(test.genes, fromType = "SYMBOL",
                  toType = c("ENTREZID","GENENAME"),
                  OrgDb = org.Mm.eg.db)
  #genes2=genes[!genes %in% gene.df$SYMBOL]
  #mouse.table[mouse.table$Symbol %in% genes2,]
  ego <- enrichGO(gene          = gene.df$ENTREZID,
                  OrgDb         = org.Mm.eg.db,
                  #keyType  = 'SYMBOL',
                  ont           = category,
                  #ont           = "BP",
                  #ont           = "MF", 
                  #ont           = "CC", 
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if(is.null(ego)){return(NULL)}

  # https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
  # https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html
  # larger cutoff, smaller returned GO terms
  x<-simplify(ego, cutoff=cutoff, by="p.adjust", select_fun=min)
  result<-x@result[order(x@result$p.adjust),]
  
  return(result)
}

dotplot_GOenrich<-function(x1,count.cut=5){  
  x1=x1[order(x1$p.adjust),]
  #x1=x1[x1$Count>=5,]
  x1$desp=factor(x1$Description,levels=rev(x1$Description))
  
  #count.cut=count.cut; #minimal GO term gene count
  p=ggplot(subset(x1,Count>=count.cut),aes(x=GeneRatio,y=desp,size=Count,col=p.adjust))+
    geom_point()+theme_bw(base_size=14)+
    scale_color_gradient(low="blue", high="red")+
    scale_size(range = c(2,8))+
    #scale_y_discrete(labels=y.lab.text)+
    #scale_size(breaks = seq(10,120,30))+
    #ylab('desp, Molecular Function')+
    #ylab(basename(file))+
    #ylab('Biological Process GO enrichment analysis of tissue specific genes')+
    ylab('')+
    theme(
      plot.title =element_text(size=20, face='bold'),
      panel.grid = element_blank(),
      axis.text=element_text(size=14),
      axis.title=element_text(size=20),
      axis.text.x=element_text(size=20,angle=45, hjust=1),
      axis.text.y=element_text(size=12,angle=0, hjust=1),
      axis.ticks.y = element_blank())
  return(p)
}

barplot_GOenrich<-function(x1, only.plot.top.n.GO=6){
  x1=x1[order(x1$p.adjust),]
  #x1=x1[x1$Count>=5,]
  x1$desp=factor(x1$Description,levels=rev(x1$Description))
  
  x1$log.p.adjust=-1*log(x1$p.adjust,base=10)
  x2=x1;
  x2$desp=as.character(x2$desp)
  x2$x.axis=x2$desp
  #if(nrow(x2)>=6){
  if(nrow(x2)>=only.plot.top.n.GO){
    x2=x2[1:only.plot.top.n.GO,]; #x2=x2[1:6,];
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }else{
    x2$x.axis=factor(x2$x.axis,levels=rev(x2$x.axis))
  }
  #plots2[[i]]<- ggplot(subset(x2,Count>=count.cut),aes(x=desp,y=log.p.adjust))+
  p<-ggplot(x2,aes(x=x.axis,y=log.p.adjust))+
    geom_bar(fill=adjustcolor('grey10',alpha.f =0.2),stat='identity',width=0.9)+
    coord_flip()+
    geom_text(aes(label=desp),y=max(x2$log.p.adjust*0.01),vjust=0,hjust = 0,size=8)+
    theme_bw(base_size=24)+ylab('-log10(p.adjust)')+xlab('')+    
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),limits = c(0, NA))+
    theme(legend.position = 'none',
          axis.title = element_text(size=22),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=20),
          plot.title = element_text(size = 20, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    return(p)
}


