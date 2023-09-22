
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
#########################################################
## per mz, per cell.type, select the connected genes with the max cor with metabolome
## record the chosen gene, and record observed P value
## only run once
obs=readRDS('obs_metabolome_sc_trajectory.rds')
length(obs) # n.mz
head(obs[[1]])

(files=Sys.glob('shuffle_out/rep*rds'))

for(file in files){
  #all.mz.out=readRDS('shuffle_out/rep_1.rds')
  all.mz.out=readRDS(file)
  names(all.mz.out) #53 mz
  
  out.file.name=gsub('rep','cor_rep',file)
  if(file.exists(out.file.name)){next}
    
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
    
    if(F){
      p.obs<-ggplot(df.met,aes(x=age,y=value))+
        geom_jitter(size=1)+geom_violin(aes(group=age),fill=NA)+
        theme_classic()+ggtitle(paste0('mz:',mz.name))+ylab('abundance')+
        scale_x_continuous(breaks=df.met$age)+
        stat_summary(
          geom = "point",fun.y = "median",
          col = "black",size = 2,
          shape = 24,fill = "red"
        )+theme(legend.position = 'none',
                plot.title = element_text(size=10))+ # +scale_y_log10()
        geom_line(data=newd,aes(x=age,y=pred),color='blue')
      #p.obs
    }
    
    ## obs per cell type per associated gene expr value
    ## for each mz, record cell.type, gene, obs.r
    df.one.mz=all.mz.out[[mz]]
    df.one.mz$expr.values=as.numeric(df.one.mz$expr.values)
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
    #pick.gene<-tmp %>% group_by(cell.type) %>% 
    #  dplyr::summarise(gene[which.max(abs(cor.coeff))],cor.coeff[which.max(abs(cor.coeff))])
    
    # look at the same gene as in observed data
    pick.gene=lapply(unique(as.character(tmp$cell.type)),function(cell.type){
      #cat(cell.type,'\n')
      if(sum(obs[[mz]]$cell.type==cell.type)==0){return(NULL)}
      gene=as.character(obs[[mz]][obs[[mz]]$cell.type==cell.type,]$gene)
      r=tryCatch(
        same.gene.out<-tmp[tmp$cell.type==cell.type & tmp$gene==gene,],
        error=function(c) 'error')
      if(r=='error'){same.gene.out=NA
      }else{same.gene.out=r$cor.coeff}
      c(cell.type,gene,same.gene.out)
    })
    pick.gene=Filter(Negate(is.null), pick.gene)
    pick.gene=as.data.frame(Reduce(`rbind`,pick.gene))
    colnames(pick.gene)=c('cell.type','gene','cor.coeff')
    pick.gene$cor.coeff=as.numeric(as.character(pick.gene$cor.coeff))
    #pick.gene=pick.gene %>% arrange(desc(abs(cor.coeff)))
    obs.result[[mz]]<-pick.gene
    
    if(F){
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
  }
  
  saveRDS(obs.result,out.file.name)
  
}

#########################################################
## calculate P value
obs.result=readRDS('obs_metabolome_sc_trajectory.rds')
length(obs.result) # n.mz
head(obs.result[[1]])

(files=Sys.glob('shuffle_out/cor*rds'))
permu.result<-lapply(1:length(files),function(i){
  permu=readRDS(files[i])
  #length(permu) # n.mz
  #sapply(permu,dim)
  permu
})
length(permu.result) # n.rep


pdf('simu_head_body_out.pdf',useDingbats = T,height = 9,width = 14)
for(mz in names(obs.result)){
  obs.one.mz=obs.result[[mz]]
  permu.mz<-lapply(permu.result,function(i) i[[mz]])
  cat(mz,'\n')
  
  mz.name=mz.age.betas[mz.age.betas$KEGGid==mz,]$mz
  
  #calculate P value across cell types
  par(mfrow=c(5,6))
  for(cell.type in unique(obs.one.mz$cell.type)){
    cat(cell.type,'\n')
    x=obs.one.mz[obs.one.mz$cell.type==cell.type,]
    
    permu.cor.values<-lapply(permu.mz,function(i){
      i[i$cell.type==cell.type & i$gene==x$gene,]$cor.coeff
    })
    permu.cor.values=unlist(permu.cor.values)
    
    pval=min((sum(permu.cor.values>x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1),
             (sum(permu.cor.values<x$cor.coeff,na.rm=T)+1)/(length(permu.cor.values)+1));
    #pval
    
    hist(permu.cor.values,
         main=paste0(cell.type,'\nobs.cor=',round(x$cor.coeff,3),
                     #'\none tail P=',formatC(pval, format = "e", digits = 2)),
                     '\none tail P=',round(pval,3)),
         #https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals
         cex.main=1,
         breaks=12,#xlab='Pearson\'s r with metabolome data')
         xlab=paste0(mz,', ',mz.name),
         xlim=c(min(permu.cor.values,x$cor.coeff,na.rm=T),max(permu.cor.values,x$cor.coeff,na.rm=T)))
    abline(v=x$cor.coeff,col='darkred',lwd=2)
  }
}
dev.off()

