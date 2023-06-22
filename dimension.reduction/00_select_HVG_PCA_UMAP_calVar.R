
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);library(grid)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(scMerge)

##########################################################################
## select only male cells and cell types with >=50cells in both 3m and 24m
if(!file.exists('TMS_male.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,4)]) #3m and 24m
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  #data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & x[,3]>=50, ] #both age groups contain >=50 cells
  dim(y) #72
  #data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  pick.cell.types=as.character(y$tc);
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  unique(sce$age) #'3m','18m','21m','24m'
  sce$age=droplevels(sce$age)
  unique(sce$age)
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  
  #remove 18m data
  sce=sce[,sce$age!='18m']
  sce$age=droplevels(sce$age)
  table(sce$age)
  #3m   24m 
  #21820 26078 
  cell.meta=colData(sce)
  saveRDS(cell.meta,'TMS_male_colData.rds')
  writeH5AD(sce, 'TMS_male.h5ad')
}

##########################################################################################################
## select HVG. log
if(!file.exists('hvg2k_mat.rds')){
  sce=readH5AD('TMS_male.h5ad')
  
  # normalization
  assayNames(sce) #'X'
  assayNames(sce)<-'counts'
  assayNames(sce) # "counts"
  sce <- logNormCounts(sce,log = TRUE) 
  #logNormCounts from scuttle R package, internal call `normalizeCounts` and `librarySizeFactors` function 
  #it's log2-transformation, https://rdrr.io/github/LTLA/scuttle/man/normalizeCounts.html
  #Library sizes are converted into size factors by scaling them so that their mean across cells is unity. https://rdrr.io/bioc/scuttle/man/librarySizeFactors.html
  assayNames(sce) # "counts"    "logcounts"
  
  # select HVG
  allf <- modelGeneVar(sce) #function from scran #https://rdrr.io/github/MarioniLab/scran/man/modelGeneVar.html
  #plot(allf$mean, allf$total)
  #curve(metadata(allf)$trend(x), add=TRUE, col="dodgerblue")
  allf=allf[order(allf$bio, decreasing=TRUE),]
  pick.genes=rownames(allf)[1:2000] #select top HGV genes
  
  # extract matrix for NMF (HGV gene x ncell, log-normalized)
  input.mat=assay(sce[pick.genes,],'logcounts')
  dim(input.mat) #2000 47898
  tmp=rowMeans(input.mat)
  summary(tmp) #check if there are any gene with 0 expr across cells
  
  saveRDS(input.mat,'hvg2k_mat.rds')
}

input.mat=readRDS('hvg2k_mat.rds')

if(!file.exists('hvg2k_mat.npy')){
  # save a npy or npz version for python
  library(reticulate)
  np=import('numpy')
  x=input.mat
  #np$savez('hvg2k_mat.npz',x=input.mat)
  np$save('hvg2k_mat.npy', x)
}
#######################################################################################################################
## PCA, https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/used.cases/dimension.reduction/dimension.reduction_PCA.R
## create single cell obj to perform PCA
library(Seurat)
obj=CreateSeuratObject(input.mat)
#obj=NormalizeData(obj)
var.genes=rownames(obj)
obj <- ScaleData(obj, features =var.genes) #for PCA
ncell=ncol(obj);

outDir='select_PCA_rank/';
if(!dir.exists(outDir)){dir.create(outDir)}
dim(input.mat) # feature by cell

start.time=Sys.time()
for(rank in c(20,50,80,100,150,200,250,500,1000,1500,2000)){
  output.file=paste0(outDir,'/pca_rank',rank,'.rds');

  #obj <- RunPCA(obj, features =var.genes,npcs = 200,verbose = F)
  obj <- RunPCA(obj, features =var.genes,npcs = rank,verbose = F)
  #dim(obj@reductions$pca@cell.embeddings) #ncell x #PC
  
  # PC variance explained (https://github.com/satijalab/seurat/issues/982)
  mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "scale.data")
  mat=mat[var.genes,] 
  pca <- obj[["pca"]]
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(mat))
  total_variance; #gene.var has already been scaled, show equal to #gene
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  cat('Seurat PCA,varExplained',sum(varExplained),'\n') #0.5576383 
  saveRDS(pca,file=output.file)
}
end.time=Sys.time()
print(end.time-start.time) #rank=200, 12min. rank=1000,4hr

(files=Sys.glob('select_PCA_rank/pca_rank*.rds'))
total_variance=2000; #top2k gene, scaled matrix
out=c();
for(file in files){
  rank=gsub('pca_rank|.rds','',basename(file))
  pca=readRDS(file)
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = sum(eigValues / total_variance)
  out=rbind(out,c(rank,varExplained))
}
df=as.data.frame(out)
colnames(df)=c('rank','var.explained')
head(df)
df$rank=as.numeric(df$rank)
pca.var.explained.plot<-ggplot(df,aes(x=rank,y=var.explained,col=factor(rank)))+geom_point()+theme_classic()
# choose 1000
df[which.max(df$var.explained),]

#######################################################################################################################
## compute aging vectors per tc: Centroid(Young)-> Centroid(Old)
cell.embed=readRDS('select_PCA_rank/pca_rank1000.rds')
dim(cell.embed@feature.loadings) #2000 gene by 1000 pc 
dim(cell.embed@cell.embeddings) # ncell by 000 pc
cell.coord=t(cell.embed@cell.embeddings)

cell.meta=readRDS('TMS_male_colData.rds')
colnames(cell.coord)=rownames(cell.meta); #barcode
cell.meta$sex=droplevels(cell.meta$sex) #only male data
tc.names=sort(unique(cell.meta$tissue_cell.type))
unique(cell.meta$age) #3m and 24m
sum(colnames(cell.coord)==rownames(cell.meta)) #47898

age.vectors=lapply(tc.names,function(tc){
  young= cell.meta$tissue_cell.type==tc & cell.meta$age=='3m'
  old= cell.meta$tissue_cell.type==tc & cell.meta$age=='24m'
  young.center=Matrix::rowMeans(cell.coord[,young])
  old.center=Matrix::rowMeans(cell.coord[,old])
  age.vector=old.center- young.center
  return(age.vector)
})
df=as.data.frame(Reduce(`cbind`,age.vectors))
colnames(df)=tc.names
rownames(df)=paste('dim',1:nrow(df))
#pheatmap::pheatmap(df)

# compute cosine similarities between the aging trajectories of each cell type in each tissue
# https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
X=df
cos.sim <- function(ix){
  A = X[,ix[1]]
  B = X[,ix[2]]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}   
n <- ncol(df) 
cmb <- expand.grid(i=1:n, j=1:n) 
C <- matrix(apply(cmb,1,cos.sim),n,n)
dim(C) #cos.similarity matrix 72 tc x 72 tc
rownames(C)=colnames(C)=colnames(df)

pdf('pca_cosine.similarity.heatmap.pdf',height = 14,width = 12,useDingbats = T)
print(pca.var.explained.plot)
print(pheatmap::pheatmap(C, show_colnames = FALSE))
dev.off()


