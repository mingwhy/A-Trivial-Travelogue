
library(Matrix)
library(ggplot2);library(gridExtra)
library(tidyverse)

# from mz to gene, or from gene to mz from `00_Fruitfly-GEM_metabolite_neighbors.R`
mz_gene_list=readRDS('mz_gene_list.rds')
names(mz_gene_list)
########################################################
## read in gene symbol mapping
gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.id.rds')
#gene.meta=readRDS('~/Documents/Data_AgingFlyCellAtlas/AFCA_gene.meta.rds')
gene.meta=gene.meta[gene.meta$FLYBASE!='',]
head(gene.meta)
dup.names=names(which(table(gene.meta$original.id)>1))

##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]
metabolome.dat=dat;
sample.info=metabolome.dat[,1:6]

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86

#########################################################
## calculate P value
obs.result=readRDS('obs_metabolome_sc_trajectory.rds')
length(obs.result) # n.mz
head(obs.result[[1]])

permu.result=readRDS('shuffle_metabolome_out.rds')
length(permu.result) # n.rep
head(permu.result[[1]])

mz.out=list()
for(mz in names(obs.result)){
  obs.one.mz=obs.result[[mz]]
  permu.mz<-permu.result[[mz]]
  cat(mz,'\n')
  
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  
  #calculate P value across cell types
  df.cor.p<-lapply( unique(as.character(obs.one.mz$cell.type)), function(cell.type){
    cat(cell.type,'\n')
    x=obs.one.mz[obs.one.mz$cell.type==cell.type,]
    
    permu.cor.values<-permu.mz[permu.mz$cell.type==cell.type & permu.mz$gene==x$gene,]
      
    permu.cor.values=permu.cor.values$cor.coeff
    
    pval=min((sum(permu.cor.values>x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1),
             (sum(permu.cor.values<x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1));
    c(mz,cell.type,x$cor.coeff,pval)
  })
  df.cor.p=as.data.frame(Reduce(`rbind`,df.cor.p))
  colnames(df.cor.p)=c('mz','cell.type','cor.coeff','Pval')
  mz.out[[mz]]<-df.cor.p
}

df.mz.out=as.data.frame(Reduce(`rbind`,mz.out))
colnames(df.mz.out)
df.mz.out$cor.coeff=as.numeric(as.character(df.mz.out$cor.coeff))
df.mz.out$Pval=as.numeric(as.character(df.mz.out$Pval))
df.mz.out$FDR=p.adjust(df.mz.out$Pval,method='BH')
colnames(df.mz.out)
summary(df.mz.out$FDR) #minimal,0.16
summary(df.mz.out$Pval) #minimal,0.009901
saveRDS(df.mz.out,'metabolite_cell.type_Pvalue.rds')

#https://stackoverflow.com/questions/5890584/how-to-reshape-data-from-long-to-wide-format
# co.coeff mat
tmp=reshape(df.mz.out[,1:3], idvar = "cell.type", timevar = "mz", direction = "wide")
mat=tmp[,-1]
rownames(mat)=tmp[,1]
colnames(mat)=gsub('cor.coeff.','',colnames(mat))
mat=as.matrix(mat)

# p value mat
tmp2=reshape(df.mz.out[,c(1,2,4)], idvar = "cell.type", timevar = "mz", direction = "wide")
mat2=tmp2[,-1]
rownames(mat2)=tmp2[,1]
colnames(mat2)=gsub('cor.coeff.','',colnames(mat2))
mat2=as.matrix(mat2)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#cell-fun
#https://github.com/jokergoo/ComplexHeatmap/issues/675

pdf("heatmap_shuffle_metabolome_Pval0.05.pdf",width=10,height=5)
#png("heatmap_shuffle_metabolome_Pval0.05.png",width=10,height=5,units="in",res=1200)
Heatmap(mat,  name='Pearson\'s r',
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(!is.na(mat2[i,j]) & mat2[i, j] < 0.05)
      # grid.text(sprintf("%.1f", mat2[i, j]), x, y, gp = gpar(fontsize = 10))
      grid.text('*',x,y,gp = gpar(fontsize = 15))
  },
  col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 8)
)
dev.off()


