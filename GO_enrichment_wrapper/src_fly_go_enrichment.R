
## homemade GO enrichment function

library(ggplot2)
library(gridExtra)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(RColorBrewer)
library(org.Dm.eg.db,verbose=F,quietly=T)


GOenrich<-function(test.genes,category='BP',topn=5,cutoff=0.7){
  gene.df <- clusterProfiler::bitr(test.genes, fromType = "SYMBOL",
                  toType = c("ENTREZID","FLYBASE","GENENAME"),
                  OrgDb = org.Dm.eg.db)
  ego <- enrichGO(gene          = gene.df$ENTREZID,
                  OrgDb         = org.Dm.eg.db,
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

