
library(Seurat)
library(tidyverse)
library(Matrix);
library(ggplot2);library(gridExtra)
library(RColorBrewer)
source("./gcl_r_version.R")

out=readRDS('sce_minCellCounts.rds')

######################################################################
## run
start=proc.time();
res.df=data.frame()

jack_knifes=as.integer(70) # 70 GCL value
jack_knife_percentage=0.7 #each time, sample 70% cell
num_divisions=as.integer(10) #for each sampled 70% cell, report Mean(10 division) 

#outfile=paste0(sex,'_GCL_',num_divisions,'_sub',sub.sample,'.rds');
(outfile=paste0('GCL_rep',jack_knifes,'.rds'));
(out.plot=gsub('.rds','pdf',outfile))

tc.names=names(out)
if(!file.exists(outfile)){
  
  for(tc in tc.names){ 
    sce.tc=out[[tc]]
    #ages=names(sce.tc)
    ages=unique(sce.tc$age)
    
    for(age in ages){
      tmp=sce.tc[,sce.tc$age==age]
      assayNames(tmp) #"counts"
      csv_mat=assay(tmp,'counts') #gene by cell matrix
      csv_mat=csv_mat[Matrix::rowSums(csv_mat>0)>=5,] #filter gene
      csv_mat=csv_mat[,Matrix::colSums(csv_mat>0)>=100] #filter cell
      #csv_mat=csv_mat^0.5;
      dim(csv_mat) #(3000, 202), 3000 gene by 202 cell
      
      nsample.cell=floor(jack_knife_percentage*ncol(csv_mat))
      
      gcl.values<-lapply(1:jack_knifes,function(j){
        input.mat=csv_mat[,sample(1:ncol(csv_mat),nsample.cell,replace = F)]
        res.values=gcl(input.mat,num_divisions)
        mean(res.values)
      })
      res=data.frame(tc=tc,gcl=unlist(gcl.values),age=rep(age,jack_knifes))
        
      res.df=rbind(res.df,res)
      cat('tc',tc,'at age',age,'is done\n')
    }
  }
  print(proc.time()-start) #38 tc, 3.5hr. neuron, 9166geneX265cell,3min
  saveRDS(res.df,outfile)
}  

res.df=readRDS(outfile)
res.df=res.df[res.df$tc!='Heart:fibroblast of cardiac tissue',]

colnames(res.df)[1]='tissue_cell.type'
res.df$tissue_cell.type=factor(res.df$tissue_cell.type,levels=sort(unique(res.df$tissue_cell.type)))
res.df$gcl
#res.df$age=factor(res.df$age,levels=c('3m','18m','24m'))
res.df$age=factor(res.df$age,levels=c('3m','24m'))

unique(res.df$tissue_cell.type) #38 tc

p=ggplot(res.df,aes(x=age,y=1-gcl))+
  facet_wrap(.~tissue_cell.type ,scale='free')+
  geom_violin()+geom_jitter(size=0.1)+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_classic()
p

pdf(out.plot,useDingbats = T,width = 16,height = 12)
print(p)
dev.off()


