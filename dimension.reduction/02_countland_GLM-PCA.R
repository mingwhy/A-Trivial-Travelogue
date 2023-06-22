
library(tidyverse)
library(ggExtra) #
library(scRNAseq)
library(ggplot2);library(gridExtra); theme_set(theme_bw())
library(scran)
library(ggpubr)
library(Seurat)
library(grDevices);library(RColorBrewer)
library(SingleCellExperiment)
plotcol <- brewer.pal(8,'Dark2')
library(glmpca)
library(scry)

########################################################
## read in dataset
sce.minCellCounts=readRDS('sce.minCellCounts.rds')
tcs=names(sce.minCellCounts)

out=list();
for(tc in tcs){
  start.time=Sys.time()
  sce.ages=sce.minCellCounts[[tc]]
  assayNames(sce.ages)
  
  cell.meta=colData(sce.ages)
  counts=assay(sce.ages,'counts')
  
  # clean up rows
  i=rowSums(counts) > 0
  cat(nrow(counts)-sum(i),' with all zero umi\n'); #0
  Y <- as.matrix(counts[i,])
  
  dim(Y); # 9877  675
  sz<-colSums(Y) #cell lib size
  Ycpm<-1e6*t(t(Y)/sz)
  Yl2<-log2(1+Ycpm)
  z<-log10(sz)
  pz<-1-colMeans(Y>0)
  cm<-data.frame(total_counts=sz,zero_frac=pz,
                 age=cell.meta$age,batch=cell.meta$mouse.id,sex=cell.meta$sex)
  
  rank=20;
  #Negative binomial likelihood
  set.seed(202)
  res2<-glmpca(Y,L=rank,X=as.matrix(cell.meta[,'age',drop=F]), fam="nb") 
  res1<-glmpca(Y,L=rank, fam="nb")
  res3<-glmpca(Y,L=rank,X=as.matrix(cell.meta[,c('age','sex'),drop=F]), fam="nb")
  res4<-glmpca(Y,L=rank,X=as.matrix(cell.meta[,c('age','sex','mouse.id'),drop=F]), fam="nb")
  
  #check optimizer decreased deviance
  #plot(res2$dev,type="l",xlab="iterations",ylab="negative binomial deviance")
  
  #visualize results
  #pd2<-cbind(cm,res2$factors,dimreduce="glmpca-nb")
  #print( ggplot(pd2,aes(x=dim1,y=dim2,colour=age))+geom_point(size=4)+facet_wrap(~sex,scales="free",nrow=3)+ggtitle(tc) )
  #out[[tc]]<-list(poi.res=res1,nb.res=res2,pca.res=res3,cell.meta=cell.meta)
  out[[tc]]=list(res1=res1,res2=res2,res3=res3,res4=res4,cell.meta=cell.meta)
  
  end.time=Sys.time()
  print(start.time);print(end.time)
  cat(tc,'is done\n') 
  #rank=20, 2min for neuron. 5min for intestine
  #rank=50, 5min for neuron. 8min for intestine
}
saveRDS(out,'glmPCA_out.rds')


#########################################################
## mutual information between each PC and cell age vector
glm.out=readRDS('glmPCA_out.rds')
names(glm.out)

par(mfrow=c(2,2))
sapply(1:4,function(i) plot(glm.out[[1]][[i]]$dev,type="l",xlab="iterations",ylab="negative binomial deviance",main=i) )
#?glmpca 
#The objective is to minimize the deviance between Y and M. 
#The deviance quantifies the goodness-of-fit of the GLM-PCA model to the data (smaller=better).

#The deviance goodness of fit test
#https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
devs<-lapply(1:4,function(i){
  x=glm.out[[1]][[i]]$dev
  x[length(x)]
  #pchisq(x[length(x)],df=484-1,lower.tail = F)
})
#https://stats.stackexchange.com/questions/237702/comparing-models-using-the-deviance-and-log-likelihood-ratio-tests
pchisq(devs[[1]]-devs[[2]],df=1,lower.tail=F)
pchisq(devs[[2]]-devs[[3]],df=1,lower.tail=F)

library(mpmi) #mutual information for mixed distributions
sig.pc=list()
for(tc in names(glm.out)){
  all.res<-glm.out[[tc]]
  cell.meta=all.res$cell.meta
  
  out<-lapply(1:4,function(i){
    nbres=all.res[[i]]
    cat(sum(rownames(nbres$factors)==rownames(cell.meta)))
    df=cbind(nbres$factors,cell.meta)
    #ggplot(df,aes(x=age,y=dim1))+geom_violin()+geom_jitter(size=0.05)
    
    #mmi(nbres$factors[,1,drop=F],cell.meta[,'age',drop=F]) #mi: raw MI estimates, bcmi: Jackknife bias corrected MI estimates (BCMI). 
    pc.mi<-mmi(nbres$factors,cell.meta[,'age',drop=F])
    length(pc.mi$bcmi) #15 pc
    
    permu.mi<-lapply(1:100,function(i){
      mmi(nbres$factors,cell.meta[sample(1:nrow(cell.meta),nrow(cell.meta),replace = F),'age',drop=F])$bcmi
    })
    df.permu.mi<-as.data.frame(Reduce(`cbind`,permu.mi))
    pvals<-lapply(1:nrow(df.permu.mi),function(i){
      (sum(df.permu.mi[i,]>pc.mi$bcmi[i])+1)/(ncol(df.permu.mi)+1)
    })
    cat(tc,'model',i,'\n')
    print(summary(unlist(pvals)))
    qvals<-p.adjust(pvals,method='BH')
    print(summary(qvals))
    cat(tc,sum(unlist(pvals)<0.05),sum(unlist(qvals)<0.05),'\n')
    df=data.frame('model'=paste0('model',i),pc=paste0('pc',length(pc.mi$bcmi)),
                  pc.mi=pc.mi$bcmi,pvals=unlist(pvals),qvals=unlist(qvals))
    return(df)
  })
  df.out=as.data.frame(Reduce(`rbind`,out))
  sig.pc[[tc]]=df.out
}
saveRDS(sig.pc,'sig_pc.rds')


#######################################################################################
## use res3 ( age and sex) to select age-related PC and calcualte Mahalonobis distance

glm.out=readRDS('glmPCA_out.rds')
sig.pc=readRDS('sig_pc.rds')

all.distances=list();
for(tc in names(glm.out)){
  all.res=glm.out[[tc]]
  res=all.res[[3]]
  cell_metadata=all.res$cell.meta;
  
  pc.out=sig.pc[[tc]]
  select.i=which(pc.out[pc.out$model=='model3',]$pval<0.05)
  cat(tc,'select.pc',select.i,'\n')
  cell.embeddings<-res$factors[,select.i]
  
  out.cells=lapply(unique(cell_metadata$age),function(i){
    cell.embeddings[cell_metadata$age==i,]
  })
  dist.out.cells=lapply(out.cells,function(x){
    if(nrow(x)<20){return(NA)}
    mahalanobis(x,colMeans(x),cov(x)) #https://www.geeksforgeeks.org/how-to-calculate-mahalanobis-distance-in-r/
  })
  names(dist.out.cells)=unique(cell_metadata$age)
  df=data.frame(age=rep( names(dist.out.cells),sapply(dist.out.cells,length)),
                dist=unlist(dist.out.cells))
  all.distances[[tc]]=df
  cat('tc',tc,'is done\n')
}
length(all.distances) #
saveRDS(all.distances,'MH.dist.rds') 


all.OT.distances=readRDS('MH.dist.rds')
df=Reduce(`rbind`,all.OT.distances)
df$tc=rep(names(all.OT.distances),sapply(all.OT.distances,nrow))
head(df)
df$age=factor(df$age,levels=c('3m','18m','24m'))
df$tc2=gsub(':','\n',df$tc)
df$tissue=sapply(strsplit(df$tc,':'),'[',1)
df[is.na(df$dist),]
df=df[!is.na(df$dist),]

plotcol <- brewer.pal(8,'Dark2')
pdf('MH.dist.pdf',useDingbats = T,width = 24,height=14)
ggplot(df,aes(x=age,y=dist,col=age))+
  scale_color_manual(values=plotcol)+scale_y_log10()+
  facet_wrap(.~tc2,scale='free')+geom_violin()+geom_jitter(size=0.5)+theme_classic()
dev.off()

plotcol <- brewer.pal(8,'Dark2')
pdf('MH.dist.pdf2.pdf',useDingbats = T,width = 12,height=8)
for(tissue in unique(df$tissue)){
  x=df[df$tissue==tissue,]
  print( ggplot(x,aes(x=age,y=dist,col=age))+
           scale_color_manual(values=plotcol)+scale_y_log10()+
           ylab('Mahalanobis distance')+
           facet_wrap(.~tc2,scale='free')+geom_violin()+geom_jitter(size=0.5)+theme_classic()
  )
}
dev.off()



