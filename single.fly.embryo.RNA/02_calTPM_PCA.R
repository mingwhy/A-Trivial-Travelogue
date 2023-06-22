#################################################
## read in meta information
df=read.table('GSE127176_series_matrix.txt',skip=30,sep='\t',fill=TRUE)
head(df)
df=t(df[c(1,10,39),])
df=df[-1,]
colnames(df)=c('sample.id','stage','GSM.id')
df.meta=df;
dim(df.meta) #54 samples
###################################################
## read in STAR output tab-limited read count data
files=Sys.glob('GSE127176_RAW/*tab')
files
raw.expr=list();
for(file in files){
  gsm=strsplit(basename(file),'_')[[1]][1]
  x=read.table(file,skip=4)
  raw.expr[[gsm]]=x
}
length(raw.expr) #54

# extract 'expressed genes', defined by >=21 counts in >=5 embryos
per.embryo.gene=lapply(raw.expr,function(x){
  x[x[,2]>=21,1]
})
sapply(per.embryo.gene,length)
expr.genes=names(which(table(unlist(per.embryo.gene))>=5))
length(expr.genes) #8983 genes
filtered.expr=lapply(raw.expr,function(x){
  x[x[,1] %in% expr.genes,]
})
sapply(filtered.expr,nrow) #8983 genes

#################################################
## TPM: read.count/gene.length
## then divided/ sum(length.corrented.count) x 10^6
df.info=data.table::fread('dme7_geneLength_chr.txt')
#query.genes=scan('all.gene.names.txt',what='')
query.genes=expr.genes
head(query.genes)
length(query.genes) #8983
sum(df.info$GENEID %in% query.genes) #8934
tmp=merge(filtered.expr[[1]],df.info,by.x='V1',by.y='GENEID')
length.corrected.counts=lapply(filtered.expr,function(x){
  tmp=merge(x,df.info,by.x='V1',by.y='GENEID')
  tmp$V2/tmp$exonic.gene.sizes/1000 #in kilobase
})
sapply(length.corrected.counts,length)
TPM=lapply(length.corrected.counts,function(x){
  x/sum(x)*10^6
})
log2.TPM=lapply(TPM,function(i) log(i+1,base=2))
#################################################
## PCA on log2(TPM)
log2TPM.df=Reduce(`cbind`,log2.TPM)
dim(log2TPM.df) #8934 x 54
rownames(log2TPM.df)=tmp[,1]
colnames(log2TPM.df)=names(log2.TPM)
log2TPM.df[1:3,1:3]

all.genes.var=apply(log2TPM.df,1,var)
quantile(all.genes.var,0.90)
log2TPM.df.genes=log2TPM.df[all.genes.var<quantile(all.genes.var,0.90),]
dim(log2TPM.df.genes) #8040   54

x=colnames(log2TPM.df.genes)
rownames(df.meta)=df.meta[,3]
df.meta=as.data.frame(df.meta[x,])
sum(df.meta[,3]==colnames(log2TPM.df.genes)) #54
df.meta.keep=df.meta;

library(factoextra)
res.pca=prcomp(t(log2TPM.df.genes),scale=T)

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
basic_plot <- fviz_pca_ind(res.pca, label="none")

ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=stage))+geom_point(shape=19)+theme_bw()

##############################################
## sex samples using k-means clustering
# sxl, msl-2, roX1, roX2
markers=c('FBgn0264270','FBgn0005616','FBgn0019661','FBgn0019660')
sum( markers %in% expr.genes) #4
sub.df=log2TPM.df[rownames(log2TPM.df) %in% markers,]
dim(sub.df) #4 x 54 samples
df=scale(t(sub.df)) #sample by feature matrix
head(df)
out=kmeans(df, centers=3, iter.max = 10, nstart = 1)
table(out$cluster)
# 1  2  3 
# 20 29  5

sum(df.meta$GSM.id==names(out$cluster))
df.meta$cluster=out$cluster
table(df.meta$stage,df.meta$cluster)
# compare with embr201948138-sup-0006-tableev1.xlsx
x=factor(df.meta$cluster);
levels(x)=c('male','female','PB')
df.meta$cluster=x
table(df.meta$cluster)

head((res.pca$sdev^2)/sum(res.pca$sdev^2))

basic_plot <- fviz_pca_ind(res.pca, label="none")
ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=stage,shape=cluster))+
  xlab('PC1 (37.5%)')+ylab('PC2 (23.3%)')+
  geom_point()+theme_bw()

out$centers
#https://towardsdatascience.com/how-to-use-and-visualize-k-means-clustering-in-r-19264374a53c
fviz_cluster(out,data=df,geom=c('point'))

# sxl, msl-2, roX1, roX2
markers=c('FBgn0264270','FBgn0005616','FBgn0019661','FBgn0019660')
rownames(sub.df)=c('msl-2','lncRNA:roX2','lncRNA:roX1','Sxl')
saveRDS(list(data=sub.df,meta=df.meta),'2019_bulk.rds')
###########################################################
# remove PB samples and perform PCA
df.meta=df.meta.keep
test=log2TPM.df.genes[,df.meta$cluster!='PB']
df.meta=df.meta[df.meta$cluster!='PB',]
res.pca=prcomp(t(test),scale=T)

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
basic_plot <- fviz_pca_ind(res.pca,axes = c(1,2), label="none")
ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=stage,shape=cluster))+geom_point()+theme_bw()

