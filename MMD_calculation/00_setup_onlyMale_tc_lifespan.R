
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(zellkonverter)
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation

######################################################################
## read in mouse turnover rate data
#cell.lifespan=data.table::fread('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cellLifespan.txt')
#cell.lifespan=readxl::read_xlsx('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cellLifespan_20220829.xlsx');
cell.lifespan=data.table::fread('~/Documents/aging_cell.turnover/2021_RonMilo_data/mouse_tc_cell.lifespan_20221026.csv');
dim(cell.lifespan) #139
colnames(cell.lifespan);
head(cell.lifespan)

cell.lifespan=cell.lifespan[cell.lifespan$include==1,] #filter based on Kim's annotation
dim(cell.lifespan) #39 tc remain 
cell.lifespan$lifespan=cell.lifespan$lifespan.used.in.this.study

cell.lifespan[grep('-',cell.lifespan$lifespan),]
#https://onlinelibrary.wiley.com/doi/10.1111/imr.12693
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature T cells',]$lifespan='70'; #T cell 
cell.lifespan[cell.lifespan$lifespan=='40-70' & cell.lifespan$human_tc=='Blood:mature B cells',]$lifespan='40'; #B cell 
#https://jamanetwork.com/journals/jamainternalmedicine/article-abstract/565579
cell.lifespan[cell.lifespan$lifespan=='200-400',]$lifespan='400'; #hepatocyte
cell.lifespan$lifespan=as.numeric(cell.lifespan$lifespan)
dim(cell.lifespan)  #39 tc

tmp=cell.lifespan[!duplicated(cell.lifespan$human_tc),]
dim(tmp); #unique 21 human_tc
table(tmp$species) #1 human and 20 rodents

tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders #39 tc

######################################################################
if(!file.exists('select.tc.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,4)]) #3m and 24m
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & x[,3]>=50, ] #both age groups contain >=50 cells
  dim(y) #72
  data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  y[grep('Brain',y$tc),]
  y$tc=as.character(y$tc)
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',1))) #22 unique tissues
  unique(unlist(lapply(strsplit(y$tc,'\\:'),'[',2))) #52 unique tissues
  
  sum(y$tc %in% tc.orders) #39, only keep those with lifespan data 
  pick.cell.types=y[y$tc %in% tc.orders,]$tc
  
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  unique(sce$age) #'3m','18m','21m','24m'
  sce$age=droplevels(sce$age)
  unique(sce$age)
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  writeH5AD(sce, 'select.tc.h5ad')
}

y=data.table::fread('ncell_per.age_per.tc.txt')
y$tissue=sapply(strsplit(y$tc,':'),'[',1)
y$cell.type=sapply(strsplit(y$tc,':'),'[',2)
dim(y) #72
table(y$tissue) #22 tissues
length(table(y$cell.type)) #52 cell.type

######################################################################
sce=readH5AD('select.tc.h5ad')
cell.meta=colData(sce)
df.out=cell.meta %>% as.data.frame() %>% group_by(tissue_cell.type,age) %>% summarise(n=n())
df.out=df.out %>% spread(age,n)
df.out$tissue_cell.type #39

sum(df.out$tissue_cell.type %in% tc.orders) #39

