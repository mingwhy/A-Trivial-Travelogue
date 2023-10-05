
library(Matrix)
#https://stackoverflow.com/questions/51467276/how-to-find-the-column-means-for-a-sparse-matrix-excluding-0-values
# Note: rows or columsn with all 0 are dropped
colMeans_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), mean)
}
rowMeans_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  sapply(split(dgCMat@x, RowInd), mean)
}
# number of non-zero: https://stackoverflow.com/questions/51560456/r-package-matrix-get-number-of-non-zero-entries-per-rows-columns-of-a-sparse

##########
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

######################################################################
## calculate mean gene expr score per cell type per age
## only run once!!
library(SingleCellExperiment)
library(zellkonverter)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(SummarizedExperiment)

sce=readH5AD('sce_filteredBy_ncell100.h5ad')

assayNames(sce)='counts'
if(!file.exists('sizeFactors_sce.rds')){
  #Pooling normalization
  #https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/norm.html#3_Pooling_normalization
  clusters <- quickCluster(sce)
  sce <- computePooledFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  saveRDS(sizeFactors(sce),'sizeFactors_sce.rds')
}
#size.factors=readRDS('sizeFactors_sce.rds')
#length(size.factors);dim(sce)
#sce=logNormCounts(sce,size.factors=size.factors,pseudo.count =1)#https://rdrr.io/github/LTLA/scuttle/man/normalizeCounts.html

sce=logNormCounts(sce);
assayNames(sce) # "counts"    "logcounts"

out.file='bootstrap_log1p_expr_remove0.rds'
unique(sce$tc)
bootstrap=50;
prop.cell=0.8; #sample 80%

all.mz.out=list();

for(mz in names(mz_gene_list)){
  
  # obs metabolome data
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  if(F){
    obs.values=metabolome.dat[, c(1:6,which(colnames(metabolome.dat)==mz.name))]
    p.obs<-ggplot(obs.values,aes(x=AgeNum,y=obs.values[,7]))+geom_jitter()+
      theme_classic()+ggtitle(paste0('mz:',mz.name))+ylab('abundance')+
      scale_x_continuous(breaks=obs.values$AgeNum)+
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+theme(legend.position = 'none')# +scale_y_log10()
    #p.obs
  }
  # extract this mz associated gene expr from scRNA-seq
  df.g=mz_gene_list[[mz]]
  genes=c(df.g$from,df.g$to)
  genes=gsub("^C[0-9]{5}", "", genes)
  genes=genes[genes!='']  
  genes=unlist(lapply(genes,function(i){unlist(strsplit(i,'or'))}))
  genes=gsub('^\\s+|\\s+$','',genes)
  genes=unique(genes)
  
  #sce.sub=sce[intersect(genes,rownames(sce)),,drop=F]
  sce.sub=sce[intersect(genes,rownames(sce)),]
  
  one.mz.out=list()
  
  for(gene.name in rownames(sce.sub)){
    out<-lapply(as.character(unique(sce.sub$tc)),function(tc){
      sce.tc=sce.sub[gene.name,sce.sub$tc==tc]
      expr.m=assay(sce.tc,'logcounts')
      #expr.m[expr.m==0]=NA #remove zero expr gene, transform NA here too slow
      #class(expr.m) # "dgCMatrix"
      #table(sce.tc$age)
      
      x=lapply(unique(sce.tc$age),function(age){
        sce.tc.age=expr.m[,sce.tc$age==age,drop=F]
        expr.values=c()
        try=0;
        while(length(expr.values)<bootstrap){
          try=try+1;
          if(try>500){expr.values=rep(0,bootstrap);break}
          tmp=sce.tc.age[,sample(1:ncol(sce.tc.age),ceiling(ncol(sce.tc.age)*prop.cell),replace = F),drop=F]
          tmp1=tmp[rowSums(tmp)!=0,] # 1xncell matrix, rowsum==0
          if(nrow(tmp1)==0){
            #expr.values=c(expr.values,0)
            next
          } #sampled values are all 0
          #gene.na=tabulate(tmp1@i + 1L, nrow(tmp1)) ### nnz per row,number of non-zeros
          exprs=rowMeans_drop0(tmp)
          names(exprs)=rownames(tmp1) #only non-zero genes expr values would be returned
          expr.values=c(expr.values,exprs)
        }
        data.frame(age=age,expr.values=expr.values)
      })
      
      df=as.data.frame(Reduce(`rbind`,x))
      df$cell.type=tc
      cat('cell.type',tc,'is done\n')
      return(df)
    })
      
    df.out=as.data.frame(Reduce(`rbind`,out))
    df.out$gene=gene.name;
    one.mz.out[[gene.name]]=df.out
  }
  df.one.mz=as.data.frame(Reduce(`rbind`,one.mz.out))
  all.mz.out[[mz]]=df.one.mz
}

saveRDS(all.mz.out,out.file);

names(all.mz.out)
###################################################################  
## per mz, per cell.type, select the connected genes with the max cor with metabolome
## record the chosen gene, and record observed P value
all.mz.out=readRDS('bootstrap_log1p_expr_remove0.rds')
names(all.mz.out) #53 mz

pdf('obs_metabolome_sc_trajectory_remove0.pdf',height = 12,width = 15)
obs.result=list();
for(mz in names(all.mz.out)){
#for(mz in names(all.mz.out)[1:2]){
  cat(mz,'\n')
  
  ## obs metabolome data
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  df.met=data.frame(cell.type='metabolome',age=metabolome.dat$AgeNum,
                 mz.name=mz.name,value=metabolome.dat[,mz.name])
  fit.metabolome=loess(value ~ age, data=df.met, span=1)
  newd <- data.frame(age=0:90)
  newd$pred <- predict(fit.metabolome, newd)
  newd$dat='mz'
  
  p.obs<-ggplot(df.met,aes(x=age,y=value))+
    geom_jitter(size=1)+geom_violin(aes(group=age),fill=NA)+
    theme_classic()+ggtitle(paste0(mz,', ',mz.name))+ylab('abundance')+
    scale_x_continuous(breaks=df.met$age)+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+theme(legend.position = 'none',
          plot.title = element_text(size=10))+ # +scale_y_log10()
    geom_line(data=newd,aes(x=age,y=pred),color='blue')
  #p.obs
  
  ## obs per cell type per associated gene expr value
  ## for each mz, record cell.type, gene, obs.r
  df.one.mz=all.mz.out[[mz]]
  tmp=expand.grid(unique(df.one.mz$cell.type),unique(df.one.mz$gene))
  cor.values=apply(tmp,1,function(row){
    df1=subset(df.one.mz,cell.type==row[[1]] & gene==row[[2]])
    df1$age=as.numeric(as.character(df1$age))
    fit.sc=loess(expr.values ~ age, data=df1, span=1)
    if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NA)}
    
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    #newdf=rbind(newd,newd2)
    #ggplot(newdf,aes(x=age,y=pred,color=dat))+geom_point()+facet_wrap(.~dat,scale='free_y')+theme_classic(base_size = 15)
    if(sd(newd2$pred,na.rm=T)<1e-22 | sd(newd$pred,na.rm=T)<1e-22){return(NA)}
    cor.value=cor(newd2$pred,newd$pred,use='complete')
    return(cor.value)
  })
  tmp$cor.coeff=cor.values
  colnames(tmp)=c('cell.type','gene','cor.coeff')
  
  # for each cell type, select the gene with larger absolute cor.coeff value
  pick.gene<-tmp %>% group_by(cell.type) %>% 
    dplyr::summarise(gene[which.max(abs(cor.coeff))],cor.coeff[which.max(abs(cor.coeff))])
  
  colnames(pick.gene)=c('cell.type','gene','cor.coeff')
  pick.gene=pick.gene %>% arrange(desc(abs(cor.coeff)))
  obs.result[[mz]]<-pick.gene
  
  # plot
  plots=lapply(1:nrow(pick.gene),function(i){
    row=pick.gene[i,]
    df1=subset(df.one.mz,cell.type==row[[1]] & gene==row[[2]])
    df1$age=as.numeric(as.character(df1$age))
    fit.sc=loess(expr.values ~ age, data=df1, span=1)
    #if(is.null(fit.sc) || is.nan(fit.sc$residuals[1])){return(NA)}
    newd2 <- data.frame(age=0:90)
    newd2$pred <- predict(fit.sc, newd2)
    newd2$dat='sc'
    
    df1$age=as.numeric(as.character(df1$age))
    p=ggplot(df1,aes(x=age,y=expr.values))+
      facet_wrap(.~cell.type, scale='free')+#geom_violin(aes(group=age),fill=NA)+
      #facet_wrap(.~title,scale='free')+
      ylab(paste0('mz:',mz.name))+
      geom_jitter(size=1)+theme_classic()+#+scale_y_log10()
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+theme(legend.position = 'none') +
      ggtitle(paste0('Pearson r = ',round(row[[3]],3),', gene ',row[[2]]))
    p+geom_vline(xintercept = as.numeric(as.character(unique(df.met$age))),linetype = "dashed")+
      theme(axis.title.y = element_text(size=5),
            plot.title = element_text(size=10))+
      geom_line(data=newd2,aes(x=age,y=pred),color='blue')+
      scale_x_continuous(breaks=sort(unique(df1$age,df.met$age)))
  }) 
  plots[[length(plots)+1]]<-p.obs
  plots2=plots[c(length(plots),1:(length(plots)-1))]
  grid.arrange(grobs=plots2,ncol=5)
  
}

dev.off()

saveRDS(obs.result,'obs_metabolome_sc_trajectory_remove0.rds')



