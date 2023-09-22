
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

###################################################################  
## per mz, per cell.type, select the connected genes with the max cor with metabolome
## record the chosen gene, and record observed P value
all.mz.out=readRDS('bootstrap_log1p_expr.rds')
names(all.mz.out) #53 mz

obs.result=list();
mz='C00147'
#for(mz in names(all.mz.out)){
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
    theme_classic(base_size = 15)+
    ggtitle(paste0(mz,', ',mz.name))+ylab('Metabolite abundance')+
    scale_x_continuous(breaks=df.met$age)+
    stat_summary(
      geom = "point",fun.y = "median",
      col = "black",size = 2,
      shape = 24,fill = "red"
    )+ theme(legend.position = 'none',plot.title = element_text(size=10))+ # +scale_y_log10()
    geom_line(data=newd,aes(x=age,y=pred),color='blue')
  p.obs
  
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
  #plots=lapply(1:nrow(pick.gene),function(i){
  i=1
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
      #facet_wrap(.~cell.type, scale='free')+#geom_violin(aes(group=age),fill=NA)+
      #facet_wrap(.~title,scale='free')+
      #ylab(paste0('mz:',mz.name))+
      ylab('Gene expression level')+
      geom_jitter(size=1)+theme_classic(base_size = 15)+#+scale_y_log10()
      stat_summary(
        geom = "point",fun.y = "median",
        col = "black",size = 2,
        shape = 24,fill = "red"
      )+theme(legend.position = 'none') +
      ggtitle(paste0('Pearson r = ',round(row[[3]],3),', gene ',row[[2]]))
    p.one.cell.type<-p+geom_vline(xintercept = as.numeric(as.character(unique(df.met$age))),linetype = "dashed")+
      #theme(axis.title.y = element_text(size=5),
      #      plot.title = element_text(size=10))+
      geom_line(data=newd2,aes(x=age,y=pred),color='blue')+
      scale_x_continuous(breaks=sort(unique(df1$age,df.met$age)))
    
    
    pdf('one_mz_example_trajectory.pdf',height = 3.5,width = 8)
    grid.arrange(grobs=list(p.obs+xlab('Age (day)'),
                            p.one.cell.type+xlab('Age (day)')),
                 ncol=2)
#}
dev.off()

permu=readRDS('shuffle_metabolome_out.rds')
df.permu.out=permu[[mz]]
x=row
x$cor.coeff
permu.cor.values=df.permu.out[df.permu.out$cell.type==cell.type,]$cor.coeff

pval=min((sum(permu.cor.values>x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1),
         (sum(permu.cor.values<x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1));
pval

hist(permu.cor.values,
     #main=paste0(cell.type,'\nobs.cor=',round(x$cor.coeff,3),
    #             #'\none tail P=',formatC(pval, format = "e", digits = 2)),
    #             '\none tail P=',round(pval,3)),
     #https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals
     cex.main=1,cex.axis=1.2,
     breaks=13,#xlab='Pearson\'s r with metabolome data')
     xlab=paste0(mz,', ',mz.name),
     xlim=c(min(permu.cor.values,x$cor.coeff,na.rm=T),max(permu.cor.values,x$cor.coeff,na.rm=T)))
abline(v=x$cor.coeff,col='darkred',lwd=6)

