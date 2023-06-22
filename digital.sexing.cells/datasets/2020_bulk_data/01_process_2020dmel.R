
setwd("~/Documents/Jay_data/external_data/2020_paper_data")
#file='expression-noise-across-fly-embryogenesis-master/150samples_retained.txt'
file='expression-noise-across-fly-embryogenesis-master/155samples_retained.txt'

x=data.table::fread(file)
dim(x)
head(x)
colnames(x)
genes=x[,1]$Ensembl.Gene.ID
dat=x[,-1]
rownames(dat)=genes
dim(dat) # 9898 gene x 154 embryo

library(org.Dm.eg.db);
x=AnnotationDbi::select(org.Dm.eg.db,keys=genes,keytype='FLYBASE',columns=c('SYMBOL'))
x[grep('Sxl|msl|roX',ignore.case = F,x$SYMBOL),]
sum(is.na(x$SYMBOL)) #45
############################
gene.meta=as.data.frame(x)
## get chr location
t2g=readRDS('../t2g_chr.coord.rds')
dim(t2g) #23932 
head(t2g)
sum(gene.meta$FLYBASE %in% t2g$ensembl_gene_id) #9854 overlap
df.gene.meta=merge(gene.meta,t2g,by.x='FLYBASE',by.y='ensembl_gene_id',all.x=T)
dim(df.gene.meta) #9898
head(df.gene.meta)

table(df.gene.meta$chromosome_name) #no Y chr genes

############################
stage=unlist(lapply(strsplit(colnames(dat),'\\.'),'[',1))
table(stage)


all.genes.var=apply(dat,1,var)
quantile(all.genes.var,0.90)
min(all.genes.var) #remove 0 variance genes
dat2=dat[all.genes.var!=0,]
dim(dat2) # 9898  154

library(factoextra)
res.pca=prcomp(t(dat2),scale=T) #PCA using all filtered genes 8040 

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
var.explained=res.pca$sdev^2/sum(res.pca$sdev^2)
var.explained=round(var.explained*100,2)

basic_plot <- fviz_pca_ind(res.pca, label="none")
ggplot(cbind(basic_plot$data,stage),
       aes(x=x,y=y,col=factor(stage)))+geom_point(shape=19)+theme_bw()+
  xlab(paste0('PC1 (',var.explained[1],'%)'))+
  ylab(paste0('PC2 (',var.explained[2],'%)'))

#########################################
# only use Sxl, msl-2, roX1, rox2
length(genes) #9898
markers=c('FBgn0264270','FBgn0005616','FBgn0019661','FBgn0019660')
genes[genes %in% markers] # 3 exist, rox2 didn't exist

sub.df=dat[genes %in% markers,]
colnames(sub.df)
rownames(sub.df)=genes[genes %in% markers]

dim(sub.df) #3 x 150
df=scale(t(sub.df)) #sample by feature matrix
dim(df); #sample by feature matrix
head(df)
out=kmeans(df, centers=3, iter.max = 10, nstart = 1)
table(out$cluster)
table(out$cluster,stage)

library(org.Dm.eg.db)
x=AnnotationDbi::select(org.Dm.eg.db,keys=rownames(sub.df),
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL #"msl-2"       "lncRNA:roX1" "Sxl"
test_set=as.matrix(sub.df)
rownames(test_set)=x$SYMBOL
df.plot=data.frame(t(test_set))
df.plot$stage=stage
ggplot(df.plot,aes(x=`msl.2`,y=`Sxl`,col=stage))+geom_point()+theme_bw()
sum(df.plot$msl.2==0) #38
dim(df.plot) #154  

#########################################
## try 2019 SVM classifier (05_SVM_on_embryo.R in 2019 folder)
library(e1071)
classifier=readRDS('../2019_paper_reproduce.result/svm_classifier_3genes_all49samples.rds')
classifier

# Predicting the Test set results
library(org.Dm.eg.db)
x=AnnotationDbi::select(org.Dm.eg.db,keys=rownames(sub.df),
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL #"msl-2"       "lncRNA:roX1" "Sxl"
test_set=as.matrix(sub.df)
rownames(test_set)=x$SYMBOL

y_pred = predict(classifier, newdata = scale(t(test_set)))
table(y_pred) #83 female, 67 male
table(y_pred, stage) ##1=female, 2=male

## plot msl-2 vs sxl
rownames(test_set)
df=as.data.frame(t(test_set))
head(df)
sum(rownames(df)==names(y_pred))
df$pred_cluster=y_pred

stage=unlist(lapply(strsplit(names(y_pred),'\\.'),'[',1))
table(stage)
df$stage=stage

plots=list();iplot=0;
for(i in 1:3){
  for(j in (i+1):4){
    iplot=iplot+1;
    #g<-ggplot(df,aes(x=df[,1],y=df[,3],col=stage,shape=pred_cluster))+
    g<-ggplot(df,aes(x=df[,i],y=df[,j],col=pred_cluster))+
      geom_point(size=2)+
      #geom_text(label=stage,nudge_x = 0.5)+
      xlab(paste0(colnames(df)[i],' log2(CPM+0.5)'))+
      ylab(colnames(df)[j])+theme_bw()
    g
    plots[[iplot]] <- ggplotGrob(g)
  }
}

length(plots)
library(gridExtra)

pdf('log2CPM_Sxl_msl-2_embryos.pdf',useDingbats =T,width = 16)
grid.arrange(grobs=plots,ncol=3)
dev.off()
