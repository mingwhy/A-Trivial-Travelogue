options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)

# chech mz missing status 

## read in RDS data
input=readRDS('mz.dat.rds')
names(input)
attach(input)
table(genotype.id,age,batch.id)
# 60 samples: 5 genotype, 2 age, 3 batch, 2 reps each condition

## combine data by running order 
run.order=protein.abundance$`Original Sample ID`
run.order
sum(sample.id %in% run.order) #all sample ID match
sum(colnames(control) %in% run.order) #14 control samples' ID match
sum(rownames(control)==rownames(dat))
sum(control$`Current MS Compounds NA` == mz$COMPOUND) #393 mz row names match

tmp=cbind(control,dat)
df.com=tmp[,run.order]
dim(df.com)
x=lapply(df.com,as.numeric)
df=as.data.frame(Reduce(`cbind`,x))
colnames(df)=colnames(df.com)
rownames(df)=control$`Current MS Compounds NA`

## remove SSILIS
df=df[-which(is.na(control$`KEGG ID NA`)),]
dim(df) #361  74

## remove mz with >5% missing measures
tmp=df[,-grep('QC',colnames(df))]
dim(tmp)
x=apply(tmp,1,function(i) sum(is.na(i)))
table(x)
# 0   56  59  60 #NA
#170   1  1  189  
table(x/length(x))
# 170 mz with none missing values, 191 with >5% missing values
df=df[x/ncol(tmp)==0,]
dim(df) #170 x 74

## plot individual mz distribution in samples by running order
colnames(df)
group=rep('sample',ncol(df))
group[grep('QC\\(I\\)',colnames(df))]='QC1';
group[grep('QC\\(S\\)',colnames(df))]='QC2';
table(group)

if(F){
  plots=lapply(1:nrow(df),function(i){
    mz.name=rownames(df)[i]
    tmp=data.frame('sample.id'=colnames(df),'mz'=as.numeric(df[i,]))
    tmp$group=group;
    tmp$sample.id=factor(tmp$sample.id,levels=tmp$sample.id)
    ggplot(tmp,aes(x=sample.id,y=mz,col=group))+geom_point(size=2)+theme_bw(base_size = 6)+
      theme(axis.text.x = element_text(angle = 45,hjust=1))+ggtitle(mz.name)
  })
  length(plots)
  plots[[2]]
  
  pdf('check.batch.effect.pdf',width=9,useDingbats = T)
  (n=ceiling(length(plots) / 4)) #4 plots one page
  for(i in 1:n){
    start=i*4-3;
    end=i*4;
    if(i*4>=length(plots)){end=length(plots)}
    do.call(grid.arrange,c(plots[start:end],ncol=1))
  }
  dev.off()
}

## separate control and sample data
colnames(df)
df.control=df[,grep('QC',colnames(df))]
df.mat=df[,-grep('QC',colnames(df))]
dim(df.control)
dim(df.mat)
sum(is.na(df.mat)) #no missing values

if(F){
  ## no need! impute values
  library(impute) #https://bioconductor.org/packages/release/bioc/html/impute.html
  df.mat.raw=df.mat
  if(exists('.Random.seed')){rm(.Random.seed)}
  df.mat.imputed <- impute.knn(as.matrix(df.mat.raw))
  names(df.mat.imputed)
  df.mat.imputed$rng.seed #should be 362436069, according to manual
  df.mat.imputed$rng.state #should be NULL
  df.mat=df.mat.imputed$data
}

pheno=data.frame(sample=sample.id,batch=batch.id,age=age,genotype=genotype.id);
edata=df.mat;
saveRDS(list(pheno=pheno,edata=edata),file='mz.filter.dat.rds')

