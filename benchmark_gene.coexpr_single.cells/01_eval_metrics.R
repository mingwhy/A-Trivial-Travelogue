
library(COGENT) #for getEdgeSimilarity
library(Seurat);
library(ggplot2);library(gridExtra);
options(stringsAsFactors = F)


gene.names=readRDS('common100genes_brain_scRNA-seq.rds')
mat.list=readRDS('ncell_2023_female_1.rds')
mat.scran=mat.list$mat.scran; #normalized by scran
mat.count=mat.list$mat.count; #raw gene count matrix
dim(mat.scran)

#######
source("src_cor.funcs.R")
source('cal_bigScale2.R');
source('cal_BaCo.R')
source('cal_CDI.R')


#pearsonA <- buildPearson(df.mat.scran)
getThresholdStability <- function(th){
  dat <- splitData(mat, propShared=0.5)
  mat.out=cal_cor.mat(method=my.method,as.matrix(dat[[1]]),as.matrix(dat[[2]]))
  abs.mat.out=lapply(mat.out,abs)
  As=lapply(abs.mat.out,function(A){
    threshold <- quantile(A[upper.tri(A)], th, na.rm=TRUE)
    A <- 1*(A>=threshold); diag(A) <- 0; 
    A
  })
  #sapply(abs.mat.out,dim)
  #return(getEdgeSimilarityCorrected(As,type="expected"))
  x=getEdgeSimilarity(As,align=FALSE)
  return(x$globalSimilarity)
}

aggregateThresholdStability <- function(th, repCount=5){
  thS <- replicate(repCount, getThresholdStability(th), simplify=FALSE)
  thS <- do.call("rbind", thS); 
  thS <- apply(thS, 2, unlist)
  return(as.data.frame(cbind(thS, threshold=th)))
}

thresholds <- c(0.9,0.95,0.98,0.99,0.995,0.998,0.999)

## begin eval
out.dir='./7metric_eval/';
dir.create(out.dir)
#for(my.method in c('pearson','pseudocell','CDI','propr.rho_p','baco','sclink','bigscale')){
for(my.method in c('pearson','pseudocell','CDI','propr.rho_p','baco','sclink')){
    
  outfile=paste0(out.dir,'/out_',my.method,'.rds')
  if(file.exists(outfile)){next}
  if(any(my.method=='CDI',my.method=='baco',my.method=='bigscale')){
    mat=mat.count;
  }else{
    mat=mat.scran;
  }
  
  thresholdComparisonDF <- mclapply(thresholds,
                  aggregateThresholdStability, mc.cores=1)
  
  thresholdComparisonDF <- do.call("rbind", thresholdComparisonDF)
  
  saveRDS(thresholdComparisonDF,file=outfile)
}
parallel::stopCluster(cl)
## end eval


(files=Sys.glob('7metric_eval/out*.rds'))
result=lapply(files,function(i){
  x=readRDS(i);
  name=gsub('out_|.rds','',basename(i))
  df=x
  #df=as.data.frame(Reduce(`rbind`,x))
  colnames(df)=c('value','threshold')
  df$metric=name;
  df
})
df.all=Reduce(`rbind`,result)
names=unique(df.all$metric)
df.all$metric=factor(df.all$metric,levels=names[c(4,5,1,2,7,3,6)])
df.all$topT=1-df.all$threshold

pdf('7metric_eval.pdf',useDingbats = T,width = 12)
print( ggplot(df.all, aes(x=metric, y=value,col=metric)) +
  theme_classic()+
    facet_wrap(.~topT)+
  geom_boxplot(outlier.shape = NA)+geom_jitter(size=0.5)+
  ggtitle("Consistency of gene co-expression measures") +
  scale_x_discrete("Method") +
  scale_y_continuous("Global similarity")+
    theme(axis.text.x=element_text(angle=45,vjust=0.8,hjust=0.8))
)
dev.off()


## if you want to use scran to normalize data
if(F){
  #https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html
  library(scran)
  sce <- SingleCellExperiment(list(counts=mat),
                              colData=DataFrame(cell.type=rep(cell.type,ncol(mat)),sex=sex),
                              rowData=DataFrame(gene=rownames(mat)) )
  sce
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  sce <- logNormCounts(sce)
  log.sce=sce@assays@data$logcounts
  dim(log.sce)
  sum(mat==0) #the same number of 0
  sum(log.sce==0) #the same number of 0
  mat.scran=log.sce
}


