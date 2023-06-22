
library(ggplot2)
library(gridExtra)
#################################################
## read in meta information
df=read.table('GSE25180_series_matrix.txt',skip=30,sep='\t',fill=TRUE)
head(df)
df=t(df[c(1,40,44),])
df=df[-1,]
df
tmp=basename(df[,2])
tmp=gsub('.bedgraph.gz','',tmp)
id.check=unlist(lapply(strsplit(tmp,'\\_'),'[',1))
sum(id.check==df[,3]) #24 checked
sample.id=gsub('GSM\\d.+?_','',tmp)
sample.id

df=cbind(df[,c(1,3)],sample.id)
colnames(df)=c('sample.info','GSM.id','sample.id')
df.meta=as.data.frame(df);
dim(df.meta) #24 samples

###################################################
## read in STAR output tab-limited read count data
dat=read.table('2011-embryo_normalized.read.count.data.txt',
               header=T,fill=T)
dim(dat) #12353    79
dat[1:3,1:9]

colnames(dat)
mat=dat[,-(1:7)]

#We were thus able to partition the overall expression of any mRNA containing w1-CaS differences into its maternal and zygotic component
#7meta+24embryo.sample+24 w1 mother sample+24 Canton-S father sample
x=mat[,c(1,1+24,1+24+24)]
tmp=x[,2]+x[,3]
cbind(x[,1],tmp)

expr.mat=mat[,1:24];
colnames(expr.mat)
sum(colnames(expr.mat) %in% df.meta$sample.id) #24 matched

rownames(expr.mat)=dat[,1]
apply(expr.mat,2,function(i){sum(i,na.rm=T)})

#################################################
## check for chr Y genes
df.gene.table=data.table::fread('validate.id_2011_paper_data.txt')
dim(df.gene.table) #12003     7
head(df.gene.table)
table(df.gene.table$chromosome_name)
# 26 gene on chrY 
df.gene.table=df.gene.table[!is.na(df.gene.table$chromosome_name),]
dim(df.gene.table) #12003 
y.genes=df.gene.table[df.gene.table$chromosome_name=='Y',]$query
expr.mat[y.genes,]
# both female and male expr Y genes, may due to mapping problem
colSums(expr.mat[y.genes,],na.rm=T)

gene.names=rownames(expr.mat);
length(gene.names) #12353
sum(gene.names %in% df.gene.table$query) #12003
rownames(df.gene.table)=df.gene.table$query
expr.mat=expr.mat[df.gene.table$query,]
rownames(expr.mat)=df.gene.table$FLYBASE

#################################################
## RPKM
RPKM=expr.mat
log2.RPKM=log(RPKM+1,base=2)

log2.RPKM[1:3,1:3]
head(df.meta)
tmp=log2.RPKM
tmp=cbind(rownames(log2.RPKM),tmp)
colnames(tmp)[1]='gene.name'
data.table::fwrite(tmp,'2011_log2.RPKM.txt')
data.table::fwrite(df.meta,'2011_sample.meta.txt')

#################################################
## PCA on log2(RPKM)
# remove NA
x=apply(log2.RPKM,1,function(i) sum(is.na(i)))
sum(x==0) #3401 genes
table(x)
log2.RPKM.df=log2.RPKM[x==0,]
dim(log2.RPKM.df) #3401 x 54
log2.RPKM.df[1:3,1:3]

all.genes.var=apply(log2.RPKM.df,1,var)
summary(all.genes.var)
quantile(all.genes.var,0.90)
log2.RPKM.df.genes=log2.RPKM.df[all.genes.var<quantile(all.genes.var,0.90),]
dim(log2.RPKM.df.genes) #3060   54

x=colnames(log2.RPKM.df.genes)
rownames(df.meta)=df.meta$sample.id
df.meta=as.data.frame(df.meta[x,])
dim(df.meta) #24 
sum(df.meta[,3]==colnames(log2.RPKM.df.genes)) #24
df.meta.keep=df.meta;

library(factoextra)
res.pca=prcomp(t(log2.RPKM.df.genes),scale=T) #PCA using all filtered genes 8040 

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
var.explained=res.pca$sdev^2/sum(res.pca$sdev^2)
var.explained=round(var.explained*100,2)

#df.meta$stage=as.numeric(gsub('Stage: Stage ','',df.meta$stage))
df.meta$stage=gsub('F|M|U|_r2','',df.meta$sample.id)
df.meta$sex=substr(df.meta$sample.id,0,1)
table(df.meta$sex)
df.meta[df.meta$sex=='U',]$sex='F'
table(df.meta$sex)
table(df.meta$stage)
basic_plot <- fviz_pca_ind(res.pca, label="none")
pca.stage=ggplot(cbind(basic_plot$data,df.meta),
       aes(x=x,y=y,shape=sex,col=factor(stage)))+
  geom_point(size=3)+theme_bw()+
  xlab(paste0('PC1 (',var.explained[1],'%)'))+
  ylab(paste0('PC2 (',var.explained[2],'%)'))+
  geom_text(aes(label=stage),nudge_x = 2,nudge_y = 2)
pca.stage

expr.genes=rownames(log2.RPKM.df.genes)

##############################################
## sex samples using k-means clustering
df.gene.table=data.table::fread('validate.id_2011_paper_data.txt')
pick.genes1=df.gene.table[grep('Sxl|msl|mof|mle|roX',ignore.case = F,df.gene.table$SYMBOL),]
pick.genes1$query; #8 genes
#"mle" "msl-3"  "msl-2"   "msl-1" "mof" "lncRNA:roX2" "lncRNA:roX1" 'Sxl'
RPKM[pick.genes1$SYMBOL,]
#markers=pick.genes1$FLYBASE
markers=pick.genes1$query[c(3,6,7,8)]
markers
sum( markers %in% rownames(log2.RPKM)) #4
#sub.df=log2.RPKM[rownames(log2.RPKM) %in% markers,]
sub.df=RPKM[rownames(RPKM) %in% markers,]
dim(sub.df) #4 x 24 samples

df=scale(t(sub.df[,-c(23,24)])) #2sample has NA,sample by feature matrix
head(df)
out=kmeans(df, centers=3, iter.max = 10, nstart = 1)
table(out$cluster)
fviz_cluster(out,data=df,geom=c('point'))

cluster.out=data.frame(sample.id=names(out$cluster),cluster=out$cluster)
df.tmp=merge(cluster.out,df.meta)
table(df.tmp$cluster,df.tmp$stage)
table(df.tmp$cluster,df.tmp$sex,df.tmp$stage)
dim(df.tmp) #22

out$centers
#https://towardsdatascience.com/how-to-use-and-visualize-k-means-clustering-in-r-19264374a53c
fviz_cluster(out,data=df,geom=c('point'))

sub.mat=as.matrix(sub.df)
plot(sub.mat[1,],sub.mat[4,])

library(RColorBrewer)
length(unique(df.meta$stage))
my.col=brewer.pal(8,'Spectral')
#my.col=brewer.pal(8,'YlGnBu')
#my.col=brewer.pal(8,'Dark2')

colnames(sub.df)=gsub('U','F',colnames(sub.df))
plots=list();iplot=0;
for(i in 1:(nrow(sub.mat)-1)){
  for(j in (i+1):nrow(sub.mat)){
    df.plot=data.frame(sub.mat[i,],sub.mat[j,],
                     stage=gsub('F|M|U|_r2','',colnames(sub.df)),
                     sex=substr(colnames(sub.df),0,1))
    gene1=rownames(sub.mat)[i];
    gene2=rownames(sub.mat)[j];
    colnames(df.plot)=c('gene1','gene2','stage','sex')
    iplot=iplot+1
    plots[[iplot]]<-ggplot(df.plot,aes(x=gene1,y=gene2,col=stage,shape=sex))+
      geom_point(size=2)+
      xlab(paste0(gene1,', RPKM'))+ylab(gene2)+
      #scale_colour_gradientn(colours = rev(terrain.colors(10)))+
      scale_color_manual(values=my.col)+
      #geom_text(aes(label=stage),nudge_x = 0,nudge_y = 0.1,size=2)+
      theme_classic()
  }
}
length(plots)
grid.arrange(grobs=plots,ncol=3)


#pdf('4gene.pair.pdf',height=6,width=12)
#grid.arrange(grobs=plots,ncol=3) 
pdf('3gene.pair.pdf',height=4,width=12)
grid.arrange(grobs=plots[c(1,3,5)],ncol=3) #remove roX2
print(pca.stage)
dev.off()

#https://plotly.com/r/3d-scatter-plots/
library(plotly)
rownames(sub.df)
df.plot=data.frame(sub.mat[1,],sub.mat[2,],sub.mat[4,],
                   stage=gsub('F|M|U','',colnames(sub.df)),
                   sex=substr(colnames(sub.df),0,1))
colnames(df.plot)=c('gene1','gene2','gene3','stage','sex')
fig <- plot_ly(df.plot, x = ~gene1, y = ~gene2, z = ~gene3, color = ~sex, 
               colors = c('#BF382A', '#0C4B8E','grey'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'msl-2'),
                                   yaxis = list(title = 'roX1'),
                                   zaxis = list(title = 'Sxl')))

fig
