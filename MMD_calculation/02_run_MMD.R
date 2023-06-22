if(F){
# in python
import os
import scmmd
import anndata
import numpy as np

adata=anndata.read_h5ad('sce_minCellCounts_PCA50.h5ad')
adata

type(adata.obsm['PCA']) #pandas.core.frame.DataFrame
adata.obsm['PCA'].to_numpy()
adata.obsm['PCA']=adata.obsm['PCA'].to_numpy() #need array here for using scmmd
adata.obsm['PCA'].shape #ncell X nPC

distances, p_values = scmmd.compute_mmd_contrast(
    adata=adata, # [Cells, Genes] object
    representation='PCA', # representation to use, "X" or key in `adata.obsm`.
    groupby='tissue_cell_type', # a categorical grouping variable in `adata.obs`
    contrast='age', # a binary contrast in `adata.obs`
    n_iters=100, # number of random sampling iterations
    sample_size=50, # sample size for random samples
    n_permutations=100, # permutations for p-val calculations
)

np.save('MMD_dist.npy', distances)
np.save('MMD_pval.npy', p_values)
}
###################################################
library(tidyverse);
library(ggplot2);library(gridExtra)
library(viridis)
library(ggplot2) 
library(grid)
library(gridExtra) 
library(reticulate)
np <- import("numpy") #https://cran.r-project.org/web/packages/RcppCNPy/vignettes/UsingReticulate.pdf

############################################################
## scallop_out.pkl
#pickle_data <- pd$read_pickle("MMD_out.pkl") #$ python python_run_scallop_batch.py
pickle_data <- np$load("MMD_dist.npy") #$ python python_run_scallop_batch.py
dim(pickle_data) #38tc, 10 reps,  3 (between,within,within)
dim(pickle_data[1,,])
out<-lapply(1:dim(pickle_data)[1],function(i) pickle_data[i,,])

# get tc names
sce.minCellCounts=readRDS('sce_minCellCounts.rds')
tc.names=sort(names(sce.minCellCounts))
names(out)<-tc.names

df.out=as.data.frame(Reduce(`rbind`,out))
df.out$tc=rep(tc.names,sapply(out,nrow))
colnames(df.out)[1:3]=c('3m-24m','3m-3m','24m-24m')


## plot
plots<-lapply(tc.names,function(tc){
  tmp=df.out[df.out$tc==tc,]
  tmp1=reshape2::melt(tmp,id='tc')
  tmp1$value=as.numeric(tmp1$value)
  tmp1$variable=factor(tmp1$variable)
  
  tmp1 %>% #filter(variable!='3m-24m') %>%
  ggplot(aes(x=variable,y=value,col=variable))+
    geom_violin()+
    #geom_point(aes(x=age_pair,y=median,group=variable), shape=18, size=3)+
    theme_classic()+ggtitle(tc)+
    stat_summary(fun.y=median, geom="point", size=2, color="red")+
    theme(axis.text.x = element_text(size=6))
})

length(plots)
pdf('MMD_out.pdf',useDingbats = T,width = 16,height = 12)
gridExtra::grid.arrange(grobs=plots,ncol=6)
dev.off()

##########################################################
## aging magnitude
tmp=reshape2::melt(df.out,id='tc')
df.res=tmp %>%
  group_by(tc,variable) %>%
  summarise_at(vars(matches("value")), funs(mean=mean,median=median))
head(df.res)

dim(df.res) #38 x 3=114
df.res=df.res[order(df.res$`median`,decreasing = T),]
df.res2=df.res[df.res$variable=='3m-24m',]

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

tc.orders.in.mouse=cell.lifespan[order(cell.lifespan$lifespan),]$`tissue: cell.type in mouse`
tc.orders.in.mouse #38 tc
tc.orders.in.human=cell.lifespan[order(cell.lifespan$lifespan),]$human_tc
tc.orders.in.human=unique(tc.orders.in.human)
######################################################################

df2=merge(df.res2,cell.lifespan,by.x='tc',by.y='tissue: cell.type in mouse')
head(df2)
colnames(df2) 

ggplot(df2,aes(x=lifespan,y=`mean`,col=human_tc))+
  geom_point()+scale_x_log10()+theme_bw()

test.out=cor.test(df2$lifespan,df2$mean,method='spearman') 
spearman.cor.value=test.out$estimate;
pval=test.out$p.value

test.out=cor.test(df2$lifespan,df2$mean,method='pearson') 
pearson.cor.value=test.out$estimate;
pearson.pval=test.out$p.value

df2$cell.type=df2$human_tc  
df2$cell.type=factor(df2$cell.type,levels=tc.orders.in.human)
p1<-ggplot(df2,aes(x=lifespan,y=mean,col=cell.type))+
  geom_point(size=4)+scale_x_log10()+theme_classic(base_size = 20)+
  xlab('cellular lifepan')+ylab('Aging magnitude')+#ylab( 'log2(ratio)')+
  #ggtitle(paste0('spearman.cor.coeff=',round(spearman.cor.value,3)))
  ggtitle(paste0('MMD',
                 ', Spearman.cor.coeff=',round(spearman.cor.value,3),
                 ', P value=',round(pval,6),
                 '\nPearson.cor.coeff=',round(pearson.cor.value,3),
                 ', P value=',round(pearson.pval,6)))+
  theme(legend.position = 'bottom')+
  geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),
              se=TRUE,level=0.90) + #https://r-graph-gallery.com/50-51-52-scatter-plot-with-ggplot2.html#:~:text=ggplot2%20provides%20the%20geom_smooth(),like%20glm%2C%20loess%20and%20more.
  scale_color_viridis(option='turbo',discrete=T)
print(p1);

pdf('MMD_out_cellular.lifespan.pdf',useDingbats = T,width = 16,height = 12)
print(p1);
dev.off()
  