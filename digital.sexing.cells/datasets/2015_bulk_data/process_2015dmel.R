
#file='GSE68062/GSE68062_Gene_abundances_after_Q_normalization.txt.gz'
file='GSE68062/GSE68062_Gene_abundances_after_FPKM_normalization.txt.gz';
x=data.table::fread(file,header = T)
dim(x)
head(x)
colnames(x)
dat=x[,grep('mel',colnames(x))]
colnames(dat)
length(unique(x$mel)) #5834
dim(x) # 6003   64


grep('lnc',x$NAME) #no long non-coding  genes
x$NAME[grep('Sxl|roX|msl',x$NAME)] # Sxl and msl-2 exist
dat1=as.matrix(dat[,5:10]) #6 dmel samples
colnames(dat1)
rownames(dat1)=x$mel #use FBgn
dim(dat1) #6003    6
tmp=as.data.frame(dat1);
tmp=cbind(rownames(dat1),tmp)
colnames(tmp)[1]='gene.name'
data.table::fwrite(tmp,file='2015_FPKM_normalization_6mel.txt')

# remove genes with var larger than 90% quantile
all.genes.var=apply(dat1,1,var)
quantile(all.genes.var,0.90)
dat2=dat1[all.genes.var!=0,]
dim(dat2) #5981

library(factoextra)
res.pca=prcomp(t(dat2),scale=T) #PCA using all filtered genes 8040 

head((res.pca$sdev^2)/sum(res.pca$sdev^2))
var.explained=res.pca$sdev^2/sum(res.pca$sdev^2)
var.explained=round(var.explained*100,2)

sex=unlist(lapply(strsplit(colnames(dat2),'\\.'),'[',2))
basic_plot <- fviz_pca_ind(res.pca, label="none")
ggplot(cbind(basic_plot$data,sex),
       aes(x=x,y=y,col=factor(sex)))+geom_point(shape=19)+theme_bw()+
  xlab(paste0('PC1 (',var.explained[1],'%)'))+
  ylab(paste0('PC2 (',var.explained[2],'%)'))

sub.df=dat2[grep('Sxl|roX|msl-2',rownames(dat2)),]
rownames(sub.df)
df=as.data.frame(t(sub.df))
head(df)
df$sex=sex

pdf('Sxl_msl-2_6dmel.embryos.pdf',useDingbats =T)
ggplot(df,aes(x=df[,2],y=df[,1],col=sex))+geom_point(size=4)+
  ylab('Sxl (FPKM)')+xlab('msl-2(FPKM)')+theme_bw()
dev.off()

