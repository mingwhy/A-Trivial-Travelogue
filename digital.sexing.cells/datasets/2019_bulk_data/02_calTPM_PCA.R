
library(ggplot2)
library(gridExtra)
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
## check for chr Y genes
df.gene.table=data.table::fread('validate.id_2019_paper_data.txt')
head(df.gene.table)
table(df.gene.table$chromosome_name)
# no gene on chrY 

#################################################
## TPM: read.count/gene.length
## then divided/ sum(length.corrented.count) x 10^6
df.info=data.table::fread('dmel_geneLength_chr.txt')
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

data.table::fwrite(as.data.frame(log2TPM.df),'gene_by_sample_log2TPM.txt',row.names =TRUE,col.names=TRUE)

# remove genes with var larger than 90% quantile
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
res.pca=prcomp(t(log2TPM.df.genes),scale=T) #PCA using all filtered genes 8040 

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
var.explained=res.pca$sdev^2/sum(res.pca$sdev^2)
var.explained=round(var.explained*100,2)

df.meta$stage=as.numeric(gsub('Stage: Stage ','',df.meta$stage))
basic_plot <- fviz_pca_ind(res.pca, label="none")
ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=factor(stage)))+geom_point(shape=19)+theme_bw()+
  xlab(paste0('PC1 (',var.explained[1],'%)'))+
  ylab(paste0('PC2 (',var.explained[2],'%)'))

##############################################
## sex samples using k-means clustering
df.gene.table=data.table::fread('validate.id_2019_paper_data.txt')
pick.genes1=df.gene.table[grep('sxl|msl|mof|mle|roX',ignore.case = T,df.gene.table$SYMBOL),]
pick.genes1; #8 genes
#"mle" "msl-3"  "msl-2"   "msl-1" "mof" "lncRNA:roX2" "lncRNA:roX1" 'Sxl'
#markers=pick.genes1$FLYBASE

# only use Sxl, msl-2, roX1, rox2
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
levels(x)=c('PB','male','female')
df.meta$cluster=x
table(df.meta$cluster)

data.table::fwrite(df.meta,'sample.meta_sex.label.txt')

head((res.pca$sdev^2)/sum(res.pca$sdev^2))

pdf('all.stages.pdf',useDingbats = T,height=4)
basic_plot <- fviz_pca_ind(res.pca, label="none")
ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=factor(stage),shape=cluster))+
  geom_point(size=3)+theme_bw()+
  xlab(paste0('PC1 (',var.explained[1],'%)'))+
  ylab(paste0('PC2 (',var.explained[2],'%)'))
dev.off()

out$centers
#https://towardsdatascience.com/how-to-use-and-visualize-k-means-clustering-in-r-19264374a53c
fviz_cluster(out,data=df,geom=c('point'))

sum(colnames(sub.mat)==df.meta$GSM.id)
sub.mat=as.matrix(sub.df)
plot(sub.mat[1,],sub.mat[4,])
rownames(sub.mat)
x=pick.genes1[match(rownames(sub.mat),pick.genes1$FLYBASE),]
sum(x$FLYBASE==rownames(sub.mat))
rownames(sub.mat)=x$SYMBOL

# check roX2 expr
names(out$cluster)==colnames(sub.mat)
table(out$cluster) #3 female, 1 male, 2 PB
sub.mat[,out$cluster==3] #female cells
sub.mat[,out$cluster==1] #male cells

library(RColorBrewer)
my.col=brewer.pal(8,'Spectral')
my.col[1]='grey'
plots=list();iplot=0;
for(i in 1:(nrow(sub.mat)-1)){
  for(j in (i+1):nrow(sub.mat)){
    df.plot=data.frame(sub.mat[i,],sub.mat[j,],
                       stage=df.meta$stage,
                       sex=df.meta$cluster)
    gene1=rownames(sub.mat)[i];
    gene2=rownames(sub.mat)[j];
    colnames(df.plot)=c('gene1','gene2','stage','sex')
    iplot=iplot+1
    plots[[iplot]]<-ggplot(df.plot,aes(x=gene1,y=gene2,col=factor(stage),shape=sex))+
      geom_point(size=2)+
      xlab(paste0(gene1,', log2(TPM+1)'))+ylab(gene2)+
      #scale_colour_gradientn(colours = rev(terrain.colors(10)))+
      scale_color_manual(values=my.col)+
      #geom_text(aes(label=stage),nudge_x =0,nudge_y = 0.5,size=2)+
      theme_classic()
  }
}
length(plots)
grid.arrange(grobs=plots,ncol=3)

pdf('4gene.pair.pdf',height=8,width=12)
grid.arrange(grobs=plots,ncol=3)
dev.off()

###########################################################
###########################################################
# sxl, msl-2, roX1, roX2
markers=c('FBgn0264270','FBgn0005616','FBgn0019661','FBgn0019660')
rownames(sub.df)=c('msl-2','lncRNA:roX2','lncRNA:roX1','Sxl')
#saveRDS(list(data=sub.df,meta=df.meta),'2019_bulk.rds')

# whole transcriptome is informative for staging.
# 4 gene markers are informative for sexing.
# remove PB samples and perform PCA
df.meta=df.meta.keep
test=log2TPM.df.genes[,df.meta$cluster!='PB']
df.meta=df.meta[df.meta$cluster!='PB',]
res.pca=prcomp(t(test),scale=T)

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
basic_plot <- fviz_pca_ind(res.pca,axes = c(1,2), label="none")
ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,col=stage,shape=cluster))+geom_point()+theme_bw()

