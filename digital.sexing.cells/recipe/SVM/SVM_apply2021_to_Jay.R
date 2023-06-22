
library(ggplot2)
library(gridExtra)
library(Seurat)
library(e1071)
library(caTools)
###############################
classifier=readRDS('svm_2021sc_train80samples.rds')
classifier.gene.names=colnames(classifier$SV)
length(classifier.gene.names) #6073 genes


###############################
# read in gene.meta info which contains chr info
gene.meta=data.table::fread('~/Documents/Data_Jay_fly_development/validated_15392genes.txt')
y.genes=gene.meta[gene.meta$chromosome_name=='Y',]$SYMBOL
length(y.genes) #34 y genes
gene.meta=as.data.frame(gene.meta)
rownames(gene.meta)=gene.meta$SYMBOL

# read in embryo data by pred time window
(files=Sys.glob('~/Documents/Data_Jay_fly_development/RNA_seurat_object/pred_windows/*.rds'))

start.time=Sys.time()

out=c()
for(file in files){
  #x=readRDS('./RNA_seurat_object/pred_windows/18-20hrs_finished_processing.rds')
  x=readRDS(file)
  expr.mat=x@assays$RNA@data
  
  tmp.gene.names=rownames(expr.mat)
  rownames(expr.mat)=gene.meta[tmp.gene.names,]$FLYBASE
  tmp.gene.names=rownames(expr.mat)
  
  overlap.genes=tmp.gene.names[tmp.gene.names %in% classifier.gene.names]
  length(overlap.genes) #6071
  
  inp=(matrix(0,nrow=ncol(expr.mat),ncol=length(classifier.gene.names)))
  rownames(inp)=colnames(expr.mat)
  colnames(inp)=classifier.gene.names
  
  inp[,overlap.genes]<-t(as.matrix(expr.mat[overlap.genes,]))
  
  test=scale(inp); 
  dim(test) #sample by feature
  sum(is.nan(test))
  test[is.nan(test)]=0
  
  preds = predict(classifier, test,probability = TRUE)
  table(preds)
  
  time=gsub('_finished_processing.rds','',basename(file))
  out[[time]]=preds;
  cat('time',time,'finished\n')
}
saveRDS(out,'SVM_apply2021_to_Jay.out.rds')

end.time=Sys.time()
print(end.time-start.time) #3hr

out=readRDS('SVM_apply2021_to_Jay.out.rds')
names(out)
out=out[order(as.numeric(lapply(strsplit(names(out),'\\-'),'[[',1)))]

pdf('test.pdf',useDingbats = T)
par(mfrow=c(3,4))
for(i in 1:length(out)){
  x=out[[i]]
  preds=attr(x,'probabilities')
  larger.probs=apply(preds,1,max)
  hist(larger.probs,main=names(out)[i],xlim=c(0.5,1))
  #p=hist(larger.probs,main=names(out)[i],xlim=c(0.5,1),type='N')
  #p$freq=p$counts/sum(p$counts)
  #p$counts=p$freq;
  #plot(p,main=names(out)[i],xlim=c(0.5,1))  
}
dev.off();


lapply(out,function(i){
  preds=attr(i,'probabilities')
  #predict.labels=ifelse(preds[,1]>0.5,'embryoFemale','embryoMale')
  predict.labels=ifelse(preds[,1]>0.9,'embryoFemale','embryoMale')
  table(predict.labels)
})
