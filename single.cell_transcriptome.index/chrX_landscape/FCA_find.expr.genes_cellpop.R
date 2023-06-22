library(tidyverse)
library(Seurat)

##########################################################
## read in gene meta info for FCA
gene.meta=data.table::fread('~/Documents/Data_fly_FCA/gene.meta.txt')
dim(gene.meta) #16276     11
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$current_symbol
dim(gene.meta) #16276    11

#################################################################################
## read in FCA data, gene.expr>max(20%cells,20cell), average gene expr per tissue
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
    x=data.frame(SYMBOL=names(x.mean),
                 n.cell=ncol(x),n.expr.cell=Matrix::rowSums(x@assays$RNA@counts>0),
                 mean=x.mean,median=x.median,var=x.var,tc_sex=i)
    merge(x,gene.meta)
  })
  length(out)
  
  #ggplot(out[[1]],aes(x=LOCATION_ARM,y=var^0.5/mean))+geom_violin()+geom_jitter(size=0.2)+theme_classic()
  out.file=paste0(out.path,'/',tissue,'_gene.var.rds')
  saveRDS(out,out.file)
  cat(tissue,'is done\n')
}

# n.expr.gene per tissue across chrs
# define: gene.expr>max(20%cells,20cell)
files=Sys.glob('chr_gene.var/*.rds')

out=list()
for(file in files){
  tissue=gsub('_gene.var.rds','',basename(file))
  x=readRDS(file)
  x=as.data.frame(Reduce(`rbind`,x))
  x=x[x$LOCATION_ARM %in% c('X','2L','2R','3L','3R'),]
  x$tc=unlist(lapply(strsplit(x$tc_sex,';'),'[',1))
  x$sex=unlist(lapply(strsplit(x$tc_sex,';'),'[',2))
  table(x$tc);
  for(tc_sex in unique(x$tc_sex)){
    tc.gene=x[x$tc_sex==tc_sex,] 
    tc.gene$tissue_tc_sex=paste(tissue,tc_sex,sep=';')
   
     tc.gene=tc.gene[which(tc.gene$n.expr.cell>max(tc.gene$ncell * 0.2 ,20)),] #filter genes
    if(nrow(tc.gene)==0){out[[tc.gene$tissue_tc_sex[1]]]<-NULL;next }
   
    out[[tc.gene$tissue_tc_sex[1]]]<-tc.gene
  }
  cat('tissue',tissue,'is done\n')
}
saveRDS(out,'expr.genes_cellpop.rds')

# plot tissue, ngene on 2L, 2R, 3L, 3R, X
length(out); #per tc per sex
names(out)
tissues=sapply(strsplit(names(out),';'),'[',1)
genes=lapply(out,'[[',1)
length(tissues);length(genes)
sapply(genes,length)
tissue.genes=lapply(unique(tissues),function(i){
  x=genes[tissues==i]
  unique(unlist(x))
})
length(tissue.genes) #17 tissues
sapply(tissue.genes,length)
names(tissue.genes)=unique(tissues)
df.out=lapply(1:length(tissue.genes),function(i){
  x=data.frame(SYMBOL=tissue.genes[[i]])
  df=merge(x,gene.meta)
  df$tissue=names(tissue.genes)[i]
  df
})
df.out=as.data.frame(Reduce(`rbind`,df.out))
head(df.out)

x=data.frame(SYMBOL=unique(df.out$SYMBOL))
df=merge(x,gene.meta)
df$tissue='All'

df.out=rbind(df.out,df)
pdf('expr.genes_cellpop.pdf',useDingbats = T)
ggplot(df.out,aes(x=LOCATION_ARM))+geom_bar()+
  facet_wrap(.~tissue)+theme_classic()
dev.off()
