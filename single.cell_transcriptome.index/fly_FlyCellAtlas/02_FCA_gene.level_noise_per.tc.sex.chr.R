library(tidyverse)
library(Seurat)
library(ggplot2);library(gridExtra)
##########################################################
## read in gene meta info for FCA
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt')
dim(gene.meta) #16276     11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
dim(gene.meta) #16276    11

gene.meta$LOCATION_MAX=as.numeric(gene.meta$LOCATION_MAX)
gene.meta$LOCATION_MIN=as.numeric(gene.meta$LOCATION_MIN)
gene.meta[is.na(gene.meta$LOCATION_MAX),] #one gene, no information for chr arm

#############################################################################################
## read gene.expr per tc per sex with chr info from `FCA_gene.noise_per.tc.per.sex.per.chr.R`
files=Sys.glob('chr_gene.var/*.rds')

for(file in files){
 #(file=files[7])
  tissue=gsub('_gene.var.rds','',basename(file))
  x=readRDS(file)
  x=as.data.frame(Reduce(`rbind`,x))
  x=x[x$LOCATION_ARM %in% c('X','2L','2R','3L','3R'),]
  x$tc=unlist(lapply(strsplit(x$tc_sex,';'),'[',1))
  x$sex=unlist(lapply(strsplit(x$tc_sex,';'),'[',2))
  #table(x$sex);table(x$tc);
  
  pdf(paste0('FCA_',tissue,'_mean-cv.pdf'),useDingbats = T,width=12,height = 8)
  
  for(tc in unique(x$tc)){
    tc.gene=x[x$tc==tc,]
    summary(tc.gene$mean);
    tc.gene=tc.gene[tc.gene$mean!=0,]; #only look at expr genes
    
    if(length(unique(tc.gene$sex))<2){next}
    if(sum(tc.gene$n.cell<20)>0){next} #per tc_sex, >=20 cells in computing gene.mean.expr
    
    p1= ggplot(tc.gene,aes(x=LOCATION_ARM,y=mean,col=sex,group=paste(LOCATION_ARM,sex)))+
      geom_violin()+scale_y_log10()+
      geom_point(size=0.2,pch=21, position = position_jitterdodge())+
      stat_summary(geom = "point", 
                   fun.y = "median", 
                   position=position_dodge(width=1),
                   size = 3, 
                   col = "black", 
                   shape = "X")+
      theme_classic()+ggtitle(paste0(tissue,'\n',tc))
    
    p2= ggplot(tc.gene,aes(x=LOCATION_ARM,y=var^0.5/mean,col=sex,group=paste(LOCATION_ARM,sex)))+
      geom_violin()+
      geom_point(size=0.2,pch=21, position = position_jitterdodge())+
      stat_summary(geom = "point", 
                   fun.y = "median", 
                   position=position_dodge(width=1),
                   size = 3, 
                   col = "black", 
                   shape = "X")+
      theme_classic()+ggtitle(paste0(tissue,'\n',tc))
    print( grid.arrange(p1,p2,nrow=2) )
  }
  dev.off()
}




