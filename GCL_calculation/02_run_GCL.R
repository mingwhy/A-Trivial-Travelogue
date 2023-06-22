#38 tc
library(zellkonverter)
library(SingleCellExperiment)
library(tidyverse)
library(ggExtra) #
library(scRNAseq)
library(ggplot2);theme_set(theme_classic())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
plotcol <- brewer.pal(8,'Dark2')
library(glmpca)
library(scry)
# for grid.table (https://stackoverflow.com/questions/31776557/how-to-adjust-the-font-size-of-tablegrob)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.6)),
  colhead = list(fg_params=list(cex = 0.6)),
  rowhead = list(fg_params=list(cex = 0.6)))

out=readRDS('sce_minCellCounts.rds')
out.folder='sce.minCellCounts_gcl'; dir.create(out.folder)
out.plot='sce.minCellCounts_gcl.pdf'
out.file='sce.minCellCounts_gcl.txt'
######################################################################
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
py_config()
sys=import('sys')
sys$path
# make sure you have 'gcl_library.py' in the working env
sys$path=c(sys$path,'./') #sys.path.append('./decibel/module/')
sys$path 
gcl_lib=import('gcl_library') 
gcl_lib$jackknife

py_run_string('from pathlib import Path')
py_run_string('import sys')
py_run_string('from pathlib import Path')
py_run_string('from os import path')
# for threading:
py_run_string('import threading')
# math tools:
np=import('numpy')
py_run_string('from scipy.spatial.distance import cdist')
py_run_string('import math')
py_run_string('import random')
# gcl library import:
py_run_string('from matplotlib import pyplot as plt')

######################################################################
## run
start=proc.time();

tc.names=names(out)
#tc.names=tc.names[c(1,2,5)]
#tc.names=tc.names[c(18,36,19,44,45)]; # "Thymus:DN4 thymocyte"
jack_knifes=as.integer(100)
num_divisions=as.integer(10)
jack_knife_percentage=0.8 

for(tc in tc.names){
  jack_knife_percentage=ifelse(tc=="Brain_Myeloid:microglial cell",0.5,0.8);
  # 2k~4k cells for microglial cell
  
  sce.tc=out[[tc]]
  #ages=names(sce.tc)
  ages=unique(sce.tc$age)
  
  (out.file=paste0(out.folder,'/',gsub(':','__',tc),'.txt'))
  if(file.exists(out.file)){cat('tc',tc,'out.file exists\n');next}
  
  res.df=data.frame();
  for(age in ages){
    tmp=sce.tc[,sce.tc$age==age]
    #tmp=sce.tc[[age]]
    
    if(T){
      csv_mat=assay(tmp,'counts')
      csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)>=5,] #filter gene
      csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100] #filter cell
      #csv_mat=csv_mat^0.5;
      dim(csv_mat) #(3000, 202), 3000 gene by 202 cell
    }
    
    if(F){
      #http://bioconductor.org/books/3.15/OSCA.basic/feature-selection.html
      tmp=logNormCounts(tmp)
      dec.tmp <- modelGeneVar(tmp)
      #pick.genes=rownames(dec.tmp[order(dec.tmp$bio, decreasing=TRUE),])[1:3000]
      pick.genes=rownames(dec.tmp[order(dec.tmp$bio, decreasing=TRUE),])[1:ngene]
      csv_mat=assay(tmp,'logcounts')
      csv_mat=csv_mat[pick.genes,]
    }
    
    #csv_mat=as.matrix(csv_mat)
    csv_mat=t(as.matrix(csv_mat)) #transpose here, the result makes more sense
    res.values=gcl_lib$jackknife(csv_mat, jack_knifes, jack_knife_percentage, num_divisions)
    
    res=data.frame(tc=tc,gcl=unlist(res.values),age=rep(age,jack_knifes))
    res.df=rbind(res.df,res)
    cat('tc',tc,'at age',age,'is done\n')
  }
  data.table::fwrite(res.df,out.file)
}
print(proc.time()-start) # neuron, 9166geneX265cell,0.5min


######################################################################
(files=Sys.glob(paste0(out.folder,'/*txt')))
x=lapply(files,function(file){x=data.table::fread(file)})
res.df=as.data.frame(Reduce(`rbind`,x))

#res.df$age=factor(res.df$age,levels=c('3m','18m','24m'))
res.df$age=factor(res.df$age,levels=c('3m','24m'))
ggplot(res.df,aes(x=age,y=1-gcl))+
  facet_wrap(.~tc ,scale='free')+
  geom_violin()+geom_jitter(size=0.1)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()

data.table::fwrite(res.df,out.file)


######################################################################
## plot 
## read in all gene result
#library(reticulate)
#pd<-import('pandas')
#pickle_data <- pd$read_pickle("scallop_out.pkl") 

## or read in shared gene result
#pickle_data=readRDS('sce.minCellCounts_gcl_top3kgene.rds')
pickle_data=data.table::fread(out.file)
pickle_data=pickle_data[pickle_data$tc!='Heart:fibroblast of cardiac tissue',]

colnames(pickle_data)[1]='tissue_cell.type'
pickle_data$tissue_cell.type=factor(pickle_data$tissue_cell.type,levels=sort(unique(pickle_data$tissue_cell.type)))
pickle_data$gcl
pickle_data$age=factor(pickle_data$age,levels=c('3m','18m','24m'))

unique(pickle_data$tissue_cell.type) #38 tc

p=ggplot(pickle_data,aes(x=age,y=1-gcl))+
  facet_wrap(.~tissue_cell.type ,scale='free')+
  geom_violin()+geom_jitter(size=0.1)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()

#pdf('sce.minCellCounts_gcl_top3kgene.pdf',useDingbats = T,width = 16,height = 12)
pdf(out.plot,useDingbats = T,width = 16,height = 12)
print(p)
dev.off()


######################################################################
######################################################################
######################################################################
if(F){
  ## get lm per cell type
  df=readRDS('sce.minCellCounts_gcl_top2kgene.rds')
  df$age=factor(df$age,levels=c('3m','18m','24m'))
  df$numeric.age=as.numeric(gsub('m','',df$age))
  df$noise=1-df$gcl
  
  #lmer_fit1<-lmer(scallop_noise~numeric.age+(1|mouse.id),data=tmp) #only random intercept
  #lm(scallop_noise~numeric.age,data=tmp)
  library(dplyr)
  library(lmerTest)
  library(broom)
  library(gridExtra)
  
  fitted_models=df %>% group_by(tc) %>%
    #do(lmer(noise~numeric.age+(1|mouse.id),data=.))
    do(model=lm(noise~numeric.age,data=.))
  #names(fitted_models); fitted_models$model
  fit.coef=as.data.frame(t(sapply(fitted_models$model,function(model){
    coef(model)
  })))
  
  fit.coef$tissue_cell.type=fitted_models$tc
  colnames(fit.coef)=c('intercept','beta','tissue_cell.type')
  saveRDS(fit.coef,'sce.minCellCounts_gcl_top2kgene_fitCoef.rds')
}

