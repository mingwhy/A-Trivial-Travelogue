#https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md
files= Sys.glob('*.beta')
files=files[-grep('hg38',files)]
res<-lapply(files,function(fname){
 #fname <- PATH
 sample.id= strsplit(fname,'_')[[1]][[1]]
 print(sample.id)
 N <- file.info(fname)$size
 content <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
 i=content[,2]>10 
 return(c(sample.id,mean(content[i,1]/content[i,2])))
})

df=as.data.frame(Reduce(`rbind`,res))
colnames(df)=c('sample','mean');
df$mean=as.numeric(df$mean)
data.table::fwrite(df,'../mean_methylation_per_tc.txt')

##
library("GEOquery")
gse<-getGEO(filename='../GSE186458_series_matrix.txt.gz',GSEMatrix = TRUE,getGPL = FALSE) 
df2=pData(gse) #253 samples x 64 cols
colnames(df2)

df3=df2[,c('age:ch1','cell type:ch1','tissue:ch1')]
colnames(df3)=c('age','cell.type','tissue');
df3$sample=rownames(df3)
dfc=merge(df3,df)

dfc$tissue_tc=paste(dfc$tissue,dfc$cell.type,sep=':')

library(tidyverse)
#x=dfc %>% group_by(cell.type) %>% summarise(mean=mean(mean))
x=dfc %>% group_by(tissue_tc) %>% summarise(mean=mean(mean))
x=x[order(x$mean),]
print(x,n=39) #39 cell.type
print(x,n=79) #79 tc
