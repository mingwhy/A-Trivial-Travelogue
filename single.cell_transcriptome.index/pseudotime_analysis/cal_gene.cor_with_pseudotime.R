
setwd("~/Documents/sc_transcriptome.index/pseudotime_analysis/")

#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

library(BiocGenerics)
library(DelayedArray)
library(DelayedMatrixStats)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(batchelor)
library(Matrix.utils)
#devtools::install_github('cole-trapnell-lab/leidenbase')
library(leidenbase)
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

library(ggplot2)
library(dplyr)
library(ggpubr)
#install.packages('Seurat')
library(Seurat)

##########################################################################
## read in tc MA pattern 
pat=data.table::fread('~/Documents/sc_transcriptome.index/dnds_1_mouse_aging_atlas/FACS_3-18-24_up.down.txt')
head(pat)
unique(pat$sign)
#pat=pat[pat$trend=='-',]
#dim(pat) #16

##################################################################
## reading in gene id.mapping and dn.ds infomation
id.mapping=data.table::fread('~/Documents/Data_mouse_aging_atlas/fac_20449genes_id.mapping.txt')  #readin_h5ad.R
filter_mouse_rat=data.table::fread('~/Documents/sc_transcriptome.index/gene.properties/dNdS/mouse_rat.dnds.txt')
head(filter_mouse_rat)
dim(filter_mouse_rat) #20697 genes
dim(id.mapping) #20449
sum(filter_mouse_rat$ensembl_gene_id %in% id.mapping$ensembl_gene_id) #20697

mouse_rat_dnds=filter_mouse_rat[!duplicated(filter_mouse_rat$ensembl_gene_id),]
head(mouse_rat_dnds)
head(id.mapping)
gene.meta=merge(mouse_rat_dnds,id.mapping,by.x='ensembl_gene_id',by.y='ensembl_gene_id')
dim(gene.meta) #17029
gene.meta$omega=gene.meta$rnorvegicus_homolog_dn/gene.meta$rnorvegicus_homolog_ds
gene.meta[is.infinite(gene.meta$omega),]$omega=NA
gene.meta[is.nan(gene.meta$omega),]$omega=NA
summary(gene.meta$omega) 

#remove genes with dn/ds>1 (more likely to be under positive selection)
length(which(gene.meta$omega<1)) #16192 genes
gene.meta=gene.meta[which(gene.meta$omega<1),]
dim(gene.meta) #16192    12

########################################################
## read in dataset
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)

## read in data
sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
class(assay(sce))
sce


dim(assay(sce))
df.expr=assay(sce)
df.expr[1:10,1:10] #umi count
dim(df.expr) #22966 110824

## sample meta information
cell.meta=colData(sce)
head(cell.meta)
dim(cell.meta) #110824     13
colnames(cell.meta)


table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month. fac: 3,18,21,24 month
table(cell.meta$tissue) #23 tissue
length(table(cell.meta$cell_ontology_class)) #120 cell types
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #207 unique #consistent with `elife-62293-supp1-v2.xlsx`
x=as.data.frame(cell.meta) %>% group_by(tissue) %>% summarise(n.cell.type=length(unique(cell_ontology_class)))
x$n.cell.type

tissue_cell.type=paste(cell.meta$tissue,cell.meta$cell_ontology_class,sep=':')
cell.meta$tissue_cell.type=tissue_cell.type

table(cell.meta[cell.meta$age=='21m',]$tissue_cell.type) #remove these cells
i=cell.meta$age!='21m'
cell.meta=cell.meta[i,]
df.expr=df.expr[,i]
cell.meta$age=factor(cell.meta$age,levels = c('3m','18m','24m'))
dim(cell.meta) # 110096     14
dim(df.expr) # 22966 110096

## create seurat obj and normalize data
obj=CreateSeuratObject(df.expr)
obj=NormalizeData(obj)
dim(obj);dim(cell.meta)
df.expr=obj@assays$RNA@data; #log1p normalized
#######################################################################
## read in elife-62293-supp1-v2.xlsx, filter to 115 tc
library(readxl)
df=read_xlsx('~/Documents/Data_mouse_aging_atlas/SuppTables1-3/elife-62293-supp1-v2.xlsx',
             sheet='all FACS tissue-cell')
colnames(df)
dim(df) #207 tc
head(df)
unique(df$tissue) #23
unique(df$cell_ontology_class) #120
tc=paste(df$tissue,df$cell_ontology_class,sep=':')
length(unique(tc)) #207

## filter tc, contain>=20 cells in all three age groups
min.cell=20
x=df[df$`3m`>=min.cell & df$`18m`>=min.cell & df$`24m`>=min.cell,]
dim(x) #115 tc
table(x$tissue)
tc=paste(x$tissue,x$cell_ontology_class,sep=':')

df.expr=df.expr[,cell.meta$tissue_cell.type %in% tc]
cell.meta=cell.meta[cell.meta$tissue_cell.type %in% tc,]
dim(cell.meta) # 100817     14
dim(df.expr) # 22966 100817

###################################################################
## plot percentage of young, middle, old cells per ps time bin
ps.ages=readRDS('./100pc_result/pseudotime_115tc.rds')
length(ps.ages);
head(names(ps.ages));

gene.cor=list();
is=which(names(ps.ages) %in% pat[pat$trend=='-',]$tissue_cell.type)
#is=which(names(ps.ages) %in% pat[pat$trend=='+',]$tissue_cell.type)
#i=5;i=14
length(is)
pdf('nonMA_16tc.pdf',useDingbats = T)
#pdf('MA_59tc.pdf',useDingbats = T)
for(i in is){
  #plots=lapply(1:length(ps.ages),function(i){
  tc=as.data.frame(ps.ages[[i]])
  tc.name=names(ps.ages)[i]
  head(tc)
  summary(tc$pseudotime)
  tc$pseudotime[is.infinite(tc$pseudotime)]=NA
  br=seq(0,max(tc$pseudotime,na.rm = T)*1.01,len=20)
  gr=cut(tc$pseudotime,breaks=br,include.lowest = T)
  tc$gr=gr;
  x=tc %>% group_by(gr,age) %>% summarise(n=n())
  x=x %>% group_by(gr) %>% mutate(ncell=sum(n))
  x$proportion=x$n/x$ncell
  head(x)
  print( ggplot(x,aes(x=gr,y=proportion,group=age,col=age))+geom_line()+theme_classic()+
    xlab('pseudotime')+theme(axis.text.x = element_blank())+ggtitle(tc.name) )
  #})
  
  ###################################################################
  ## calculate gene spearman.cor with cell pseudotime  
  
  tc=as.data.frame(ps.ages[[i]])
  tc.name=names(ps.ages)[i]
  
  df.sub=df.expr[,cell.meta$tissue_cell.type==tc.name]
  dim(df.sub)
  dim(tc)
  
  tc$pseudotime[is.infinite(tc$pseudotime)]=NA
  sp.cors=cor(t(as.matrix(df.sub)),tc$pseudotime,method='spearman',use = "pairwise.complete.obs")
  length(sp.cors) #number of genes
  summary(sp.cors) #there may be NA due to 0 sd for some genes
  sp.cors=as.numeric(sp.cors)
  names(sp.cors)=rownames(df.sub)
  sp.cors=sort(sp.cors)
  tail(sp.cors)
  
  
  ################################################################
  ## for each tc, plot gene cor.value ~ dn/ds
  x=sp.cors[gene.meta$mgi_symbol]
  y=gene.meta$omega
  summary(x);summary(y)
  summary(lm(y~x))
  #plot(x,y,log='y',pch=16,cex=0.3,xlab='spearman.cor',ylab='dn/ds')
  tmp=data.frame(gene=gene.meta$mgi_symbol,sp.cors=x,omega=y)
  tmp=tmp[!is.na(tmp$sp.cors),]
  tmp=tmp[order(tmp$sp.cors),]
  print(ggplot(tmp,aes(x=sp.cors,y=omega))+geom_jitter(size=0.2,col=rgb(0,0,1,0.2))+
    scale_y_log10()+theme_classic()+
    xlab('spearman.cor')+ylab('dn/ds')+ggtitle(tc.name)+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x,lwd=0.5) +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x=0,label.y = 0.8)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
    stat_regline_equation(label.x=0,label.y = 1)) #this means at 30th unit regresion line equation will be shown
  #https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
  gene.cor[[tc.name]]=tmp
  cat('tc',tc.name,'is done\n')
}
saveRDS(gene.cor,'gene.cor_nonMA_16tc.rds')
#saveRDS(gene.cor,'gene.cor_MA_59tc.rds')
dev.off()

# shuffle x and y
tmp$omega.shuffle=sample(tmp$omega,nrow(tmp),replace = F)
ggplot(tmp,aes(x=sp.cors,y=omega.shuffle))+geom_jitter(size=0.2)+
  scale_y_log10()+theme_classic()+
  xlab('spearman.cor')+ylab('dn/ds')+ggtitle(tc.name)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  stat_cor(label.x=0,label.y = 0.8)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  stat_regline_equation(label.x=0,label.y = 1) #this means at 30th unit regresion line equation will be shown

#gene.cor=readRDS('gene.cor_MA_59tc.rds')
gene.cor=readRDS('gene.cor_nonMA_16tc.rds')
x=lapply(gene.cor,function(i){head(i,200)})
y=sort(table(unlist(sapply(x,'[[',1))))
tail(y)
table(y)

####################################################
## GSEA functional analysis
library(clusterProfiler) #https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
library(org.Mm.eg.db) #https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
keytypes(org.Mm.eg.db)
organism='org.Mm.eg.db'

gsea.out=list();
#pdf('gsea_nonMA_16tc.pdf',useDingbats = T,width = 16)
pdf('gsea_MA_59tc.pdf',useDingbats = T,width = 16)
for(i in 1:length(gene.cor)){
  head(gene.cor[[i]])
  tc.name=names(gene.cor)[i]
  
  gene_list=gene.cor[[1]]$sp.cors
  names(gene_list)=gene.cor[[1]]$gene
  gene_list<-na.omit(gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  
  print(dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)+ggtitle(tc.name)+theme_classic(base_size = 10))
  #ridgeplot(gse) + labs(x = "enrichment distribution")+ggtitle(tc.name)
  gsea.out[[tc.name]]=gse
  cat('tc',tc.name,'is done\n')
}
#saveRDS(gsea.out,'gsea_nonMA_16tc.rds')
saveRDS(gsea.out,'gsea_MA_59tc.rds')
dev.off()


