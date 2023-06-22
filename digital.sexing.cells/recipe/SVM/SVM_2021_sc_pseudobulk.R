
library(ggplot2)
library(Seurat)
library(gridExtra)
########################################################
## read in gene chro info
df.gene.table=data.table::fread('../../single.cell_datasets/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol

########################################################
## read in dataset
if(F){
  library(Seurat)
  #tissue='embryoUnsex'
  #dat=readRDS('../../single.cell_datasets/embryo_germline/whole_embryo.germline_unsex_filtered_valid.rds')
  tissue='embryoSexed'
  dat1=readRDS('../../single.cell_datasets/embryo_germline/whole_embryo.germline_female_filtered_valid.rds')
  dat2=readRDS('../../single.cell_datasets/embryo_germline/whole_embryo.germline_male_filtered_valid.rds')
  
  # gene numbers are not the same
  sum(colnames(dat1) %in% colnames(dat2)) #17 same ID
  mat1=dat1@assays$RNA@counts
  colnames(mat1)=paste0('female_',colnames(mat1))
  mat2=dat2@assays$RNA@counts
  colnames(mat2)=paste0('female_',colnames(mat2))
  obj1=CreateSeuratObject(mat1)
  obj2=CreateSeuratObject(mat2)
  obj1$sex='embryoFemale';
  obj2$sex='embryoMale';
  dat<- merge(obj1, obj2, project = "embryo.sexed")
  table(dat$sex)
  
  ##
  mat <- dat@assays$RNA@counts
  gene.names <- rownames(mat)
  cell.names <- colnames(mat)
  #cell.types<-dat$annotation
  #sex.labels<-dat$sex
  sex=dat$sex
  cell.type='unknown';
  
  ##
  dim(mat) #12224gene, 17142cell
  length(gene.names)
  length(cell.names)
  #table(cell.types)
  #table(sex.labels)
  rownames(mat)=gene.names
  colnames(mat)=cell.names
  mat[1:3,1:3]
  
  umi.mat=as.matrix(mat)
  # filter out some genes
  #filter <- rowSums(umi.mat>2)>=5
  #table(filter) #FASE: 4646, TRUE:7578   
  #umi.mat<- umi.mat[filter,]
  
  ########################################################
  # use gene chr info to filter genes
  # (some genes annotated to non A or X chr, discard these genes)
  i=rownames(umi.mat) %in% df.gene.table$current_symbol
  sum(i);#7574
  dim(umi.mat) #7578 17142
  umi.mat=umi.mat[i,]
  
  
  #########################################
  ## use Seurat to integrate data 
  dat.both=CreateSeuratObject(umi.mat)
  dim(dat.both) # 12185 17142
  sum(colnames(dat.both)==colnames(dat)) #17142
  dat.both$sex=dat$sex
  
  dat.both<- NormalizeData(dat.both, verbose = FALSE)
  dat.both <- FindVariableFeatures(dat.both, selection.method = "vst", 
                                   nfeatures = 10000, 
                                   verbose = TRUE)
  #output of nfeatures=5000 or 10000 is consistent with monocle3 UMAP result
  dat.both <- ScaleData(dat.both, features =rownames(dat.both))
  dat.both <- RunPCA(dat.both,features=VariableFeatures(object = dat.both), npcs = 100)
  dat.both <- RunUMAP(dat.both, dims = 1:100)
  
  saveRDS(dat.both,'integrated.sexed.samples_seurat.obj.rds')
}

dat.both=readRDS('../../single.cell_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
grep("Sxl|msl-2",rownames(dat.both))
df=dat.both@assays$RNA@counts; #gene by cell matrix
dim(df)

########################################################
## pseudobulk: combine 20 cells into one pseudocell
class(df)
df.list=list()
for(sex in unique(dat.both$sex)){
  tmp=df[,dat.both$sex==sex]
  dim(tmp)
  
  ncell=floor(ncol(tmp)/20)
  df.pse=c()
  for(i in 1:(ncell+1)){
    start=20*(i-1)+1;
    end=start+20;
    if(end>ncol(tmp)){end=ncol(tmp)}
    j=Matrix::rowSums(tmp[,start:end])
    df.pse=cbind(df.pse,j)
  }
  df.pse=as.data.frame(df.pse)
  df.pse$sex=sex;
  df.list[[sex]]=df.pse
}

df.pse=as.data.frame(cbind(df.list[[1]],df.list[[2]])) #gene by cell data.frame
rownames(df.pse)
dim(df.pse) #12185   860
sex.label=rep(unique(dat.both$sex),c(ncol(df.list[[1]]),ncol(df.list[[2]])))
table(sex.label);sapply(df.list,ncol)


# construct SVM
gene.names=rownames(df.pse)
## change gene symbol to FBgn
df.gene.table=data.table::fread('../../single.cell_datasets/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
table(df.gene.table$LOCATION_ARM) #15 Y-chromosome genes

df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol
sum(gene.names %in% df.gene.table$submitted_item) #11740
length(gene.names) #12185

overlap.genes=gene.names[gene.names %in% df.gene.table$submitted_item]
df.pse=df.pse[overlap.genes,]
dim(df.pse) #11740   860

rownames(df.pse)<-df.gene.table[overlap.genes,]$FBID_KEY
df.pse=data.matrix(df.pse)
typeof(df.pse) #make sure it's numeric

## filter out zero variance genes or genes with >90% zero
i=apply(df.pse,1,sd)
sum(i==0) #0 genes with 0 variance
i2=apply(df.pse,1,function(x) sum(x==0)>0.9*nrow(dataset))
sum(i2!=0) #3516 genes with >90% zero
df.pse=df.pse[i!=0 & i2==0,]
dim(df.pse) # 8224 gene x  860pseudocell


## create data frame: sample by feature + 1column as label
data = data.frame(t(df.pse), y =sex.label)
head(data)
tail(colnames(data)) #last column is sex label

## create training and validation dataset
library(e1071)
library(caTools)

dataset=data
dim(dataset) #860 8225
colnames(dataset) #feature names
(label.column=ncol(dataset)) #5nd column is the cluster label

if(T){
  start.time=Sys.time();
  set.seed(321) 
  split = sample.split(dataset$y, SplitRatio = 0.7)
  training_set = subset(dataset, split == TRUE)
  test_set = subset(dataset, split == FALSE)
 
  
  # Feature Scaling for train and test separately
  training_set[-label.column] = scale(training_set[-label.column])
  test_set[-label.column] = scale(test_set[-label.column])
  training_set$y=factor(training_set$y)
  sum(is.na(training_set)) #make sure there is no NA in svm
  
  #https://stackoverflow.com/questions/33688421/missing-value-where-true-false-needed-with-caret
  #https://stackoverflow.com/questions/36027510/error-in-if-anyco-missing-value-where-true-false-needed
  classifier = svm(formula = y ~ .,
                   data = training_set,
                   type = 'C-classification',
                   probability = TRUE,
                   kernel = 'linear')
  #gamma=0.05,kernel = 'radial', cost=10)
  y_pred = predict(classifier, newdata = test_set[-label.column])
  #y_pred = predict(classifier, newdata = test_set[-label.column],prob=T)
  table(test_set[,label.column],y_pred)
  #classifier
  #saveRDS(classifier,'svm_classifier_2019_train34samples.rds')
  saveRDS(classifier,'svm_2021sc_pseudobulk_train80samples.rds')
  end.time=Sys.time();
  print(end.time-start.time)
}
#Time difference of 13.47783 secs
# 2hr
#https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
w<-t(classifier$coefs) %*% classifier$SV # weight vectors
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
length(w) #8934
head(w,20)

library(org.Dm.eg.db)
dim(dataset) #49 x 5
top.genes=names(w)[1:50]
x=AnnotationDbi::select(org.Dm.eg.db,keys=top.genes,
                        keytype='FLYBASE',columns=c('SYMBOL'))
x$SYMBOL  

preds = predict(classifier, test_set[,-label.column])
table(preds) 
table(preds,test_set[,label.column])

