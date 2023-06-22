
## look for genes informative about sex labels
## identity a non-redundant gene set for sex label classificaiton in single cells
library(Seurat)
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)
library(iterators)
library(doParallel)
numCores <- parallel::detectCores() # Requires library(parallel)
registerDoParallel(2) #foreach(i=1:nrow(df.c),.combine="c") %dopar% { ...}

source('src_cal_mutual_information.R')
library(infotheo)
compute_fun<-function(x,y){return(entropy(x)+entropy(y)-entropy(cbind(x,y)))}


## read in test data
dat.test=readRDS('test.rds')
dat.test; #12878 features across 11788 samples within 1 assay 
  
cell.type=unique(dat.test$annotation) #"enterocyte-like"

## normalize data
dat.test <- NormalizeData(dat.test, normalization.method = "RC", 
                          scale.factor = 1e6) #only CPM, counts per million, no log


## remove gene which express in <30% cells
df.test=dat.test@assays$RNA@data
i=Matrix::rowSums(df.test>0)
df.test=df.test[i>ncol(dat.test) * 0.3,]
dim(df.test)

# log2(CPM+1) transformation
df.test = log(df.test+1,base=2);

## discretize expression levels
bins=seq(floor(min(df.test)),ceiling(max(df.test)))
hist.data=hist(as.numeric(df.test),breaks = bins,plot=F)
hist.data$counts=log10(hist.data$counts)
hist.data$counts[is.infinite(hist.data$counts)]=0
plot(hist.data,xlab="Expression level (log2(CPM+1))",ylab="log10(Occurrences)",main=cell.type)

bins=seq(floor(min(df.test)),ceiling(max(df.test)))
x=cut(as.numeric(df.test),breaks=bins,include.lowest = TRUE)
if(sum(table(x))!=nrow(df.test)*ncol(df.test)) cat('warning',cell.type,"cut didn't include all values\n")


#use 3 as cutoff to binarize data
df.c=df.test
binary.cutoff=3
df.c[df.c<binary.cutoff]=0;
df.c[df.c>=binary.cutoff]=1;
i=Matrix::rowSums(df.c) #remove genes with all 0
j=(i==ncol(df.c) | i==0) #expr in all or none  cells
if(sum(j)!=0) df.c=df.c[!j,]

# get sex label
table(dat.test$sex)
x=as.numeric(factor(dat.test$sex))
label.test=data.frame(name=colnames(df.test),label=x)
if(sum(label.test$name!=colnames(df.c))==ncol(df.c)) cat('warning: cell number of label and df discrete matrix do not match!\n')


## compute mutual information between each gene and cluster assignment
df_info<- foreach(i=1:nrow(df.c),.combine="c") %dopar% {
  #mutinformation(df.c[i,],label.test[,2])/log(2)
  compute_fun(df.c[i,],label.test[,2])/log(2)
}
#mi=lapply(1:nrow(df.c),function(i) compute_fun(df.c[i,],label.test[,2])/log(2)) #too slow
names(df_info)=rownames(df.c) #gene.name
df_info=sort(df_info,decreasing = T)
  
## Calculate median expression of each gene within each cluster
# as most of genes are 0, median=0, use sum or mean on orignal log(CPM+1) values
x1=Matrix::rowMeans(df.test[,label.test$label==1]) 
x2=Matrix::rowMeans(df.test[,label.test$label==2])
df_expr_labels=cbind(x1,x2)

## Entropy of classification without gene expression data
H_naive=entropy(label.test$label)/log(2)
cat("Total entropy of classification", H_naive, "bits")

## Find non redundant gene set
df_labels=label.test$label
names(df_labels)=label.test$name

topn=100; #test on topn genes
genes=names(df_info[1:topn])
# as there is entropy calculation and loop, a little bit slow
df.c=as.matrix(df.c)
genes_nonredundant=find_nonredundant_gene_set(df_discrete=df.c, genes=genes,
                                              df_labels,df_info,
                                              H_naive, N_constrain=length(genes),
                                              cumulative_information_cutoff=0.999,
                                              verbose=FALSE)
print(genes_nonredundant)


## Calculate cumulative information for top N genes
cumulative_informations = calc_cumulative_informations(df_discrete=df.c, genes_nonredundant,
                                                       df_labels, N=length(genes_nonredundant))
print(cumulative_informations)
df_info_nonredundant=cumulative_informations

## Calculate information relative to total entropy
relative_cumulative_informations = df_info_nonredundant$cumulative_mutual_information / H_naive
df_info_nonredundant$relative_cumulative_informations=relative_cumulative_informations
df_info_nonredundant$mutual_information=df_info[df_info_nonredundant$symbol]


## extract non-redundant genes and plot their gene expression between males and females
data = df_expr_labels[ df_info_nonredundant$symbol,]
# Discrete or Continuous
data_discrete = data;
data_discrete[data >= binary.cutoff]=1 #binary.cutoff is 3
data_discrete[data < binary.cutoff]=0
pheatmap::pheatmap(data,cluster_rows = F,cluster_cols = F) #sum log2(CPM+1) per gene per sex

