library(tidyverse)
library(Seurat)

##########################################################
## read in gene meta info for FCA
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt')
dim(gene.meta) #16276     11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
dim(gene.meta) #16276    11

#####################################################################
## read in FCA data, average gene expr per tissue
folders=folders=Sys.glob('~/Documents/Data_fly_FCA/*FCA*')
out.path='chr_gene.var';
dir.create(out.path)

tc.expr.var=list()
for(folder in folders){
  tissue=strsplit(basename(folder),'\\_')[[1]][2]
  input.file=Sys.glob(paste0(folder,'/*valid.rds'))
  cat('begin',tissue,input.file,'\n')
  
  ## read in data
  #tissue='gut'
  #input.file = "../../single.cell_datasets/FCA_gut/whole_gut_filtered_valid.rds"
  sc = readRDS(input.file)
  dim(sc)
  
  # antenna and body has female, male, mix, three class labels
  # only keep cell types with both female and male cells profiled
  sc=sc[,sc$sex!='mix' & sc$annotation!='unannotated']
  sc$tc_sex=paste(sc$annotation,sc$sex,sep=';')
  
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 1e6,verbose = F) #ln(CPM+1)
  cat(tissue,'tc_sex #',length(unique(sc$tc_sex)),'\n');
  out=lapply(unique(sc$tc_sex),function(i){
    x=sc[,sc$tc_sex==i];
    x.mean=Matrix::rowMeans(x@assays$RNA@data)
    x.median=apply(x@assays$RNA@data,1,median)
    x.var=apply(x@assays$RNA@data,1,var)
    x=data.frame(SYMBOL=names(x.mean),n.cell=ncol(x),
                 mean=x.mean,median=x.median,var=x.var,tc_sex=i)
    merge(x,gene.meta)
  })
  length(out)
  
  #ggplot(out[[1]],aes(x=LOCATION_ARM,y=var^0.5/mean))+geom_violin()+geom_jitter(size=0.2)+theme_classic()
  out.file=paste0(out.path,'/',tissue,'_gene.var.rds')
  saveRDS(out,out.file)
  cat(tissue,'is done\n')
}

x=out[[2]]
head(x);unique(x$LOCATION_ARM)
x=x[x$LOCATION_ARM %in% c('X','2L','2R','3L','3R'),]
x$tc=unlist(lapply(strsplit(x$tc_sex,';'),'[',1))
x$sex=unlist(lapply(strsplit(x$tc_sex,';'),'[',2))
ggplot(x,aes(x=LOCATION_ARM,y=var^0.5/mean,col=sex))+geom_violin()+geom_jitter(size=0.2)+theme_classic()

## plot
files=Sys.glob('chr_gene.var/*.rds')
pdf('chr_gene.var.pdf',useDingbats = T)
for(file in files){
  tissue=gsub('_gene.var.rds','',basename(file))
  x=readRDS(file)
  x=as.data.frame(Reduce(`rbind`,x))
  x=x[x$LOCATION_ARM %in% c('X','2L','2R','3L','3R'),]
  x$tc=unlist(lapply(strsplit(x$tc_sex,';'),'[',1))
  x$sex=unlist(lapply(strsplit(x$tc_sex,';'),'[',2))
  table(x$tc);
  for(tc in unique(x$tc)){
    tmp=x[x$tc==tc,]
    print( ggplot(tmp,aes(x=LOCATION_ARM,y=var^0.5/mean,col=sex,group=paste(LOCATION_ARM,sex)))+
      geom_violin()+
      geom_point(size=0.2,pch=21, position = position_jitterdodge())+
      stat_summary(geom = "point", 
                   fun.y = "median", 
                   position=position_dodge(width=1),
                   size = 3, 
                   col = "black", 
                   shape = "X")+
      theme_classic()+ggtitle(paste0(tissue,'\n',tc))
    )
  }
}
dev.off()



