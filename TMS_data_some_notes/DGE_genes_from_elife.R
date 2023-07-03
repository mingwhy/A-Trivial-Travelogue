library(tidyverse)
library(ggplot2)

inp.folder='~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/DGE_result_release_sig/'
inp.files=Sys.glob(paste0(inp.folder,'/facs*'))
length(inp.files) #76 tissue-cell types as presented in elife paper

tc=sapply(inp.files,function(i){
  gsub('facs\\.|\\.gz','',basename(i))
})
head(tc) #tissue.celltype names

all.dat=lapply(inp.files,function(file){
  data.table::fread(file)});
names(all.dat)=tc
sapply(all.dat,dim)
tc[1]
head(all.dat[[1]])
summary(abs(all.dat[[1]]$`coef (age.logFC)`)) #all > 0.005 as elife method section stated
summary(all.dat[[1]]$`fdr (based on age.H_p)`) #all <0.01 as elife method section stated

sum(all.dat[[1]]$`coef (age.logFC)`>0) #up genes with age
sum(all.dat[[1]]$`coef (age.logFC)`<0) #down genes with age

##