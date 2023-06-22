
library(ggplot2);library(gridExtra)
library(dplyr)
library(ggpubr)
#install.packages('Seurat')
library(Seurat)
library(irlba) # for fast PCA on large matrix: https://github.com/bwlewis/irlba
library(RcppML) #for fast NMF
library(grDevices);library(RColorBrewer)
########################################################
## read in dataset
#BiocManager::install("zellkonverter")
library(zellkonverter)
library(SummarizedExperiment)

## read in data
inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas//TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
class(assay(inp.sce))
inp.sce

dim(assay(inp.sce))
df.expr=assay(inp.sce)
df.expr[1:10,1:10] #umi count
dim(df.expr) #22966 110824

## sample meta information
cell.meta=colData(inp.sce)
head(cell.meta)
dim(cell.meta) #110824     13
colnames(cell.meta)


table(cell.meta$age) #6 age group: 1, 3, 18, 21, 24,30month. fac: 3,18,21,24 month
table(cell.meta$tissue) #23 tissue
length(table(cell.meta$cell_ontology_class)) #120 cell types
unique(paste(cell.meta$tissue,cell.meta$cell_ontology_class)) #207 unique #consistent with `elife-62293-supp1-v2.xlsx`

tissue_cell.type=paste(cell.meta$tissue,cell.meta$cell_ontology_class,sep=':')
cell.meta$tissue_cell.type=tissue_cell.type

table(cell.meta[cell.meta$age=='21m',]$tissue_cell.type) #remove these cells
i=cell.meta$age!='21m'
cell.meta=cell.meta[i,]
df.expr=df.expr[,i]
cell.meta$age=factor(cell.meta$age,levels = c('3m','18m','24m'))
dim(cell.meta) # 110096     14
dim(df.expr) # 22966 110096

cell.meta$binary.age='old';
cell.meta[cell.meta$age=='3m',]$binary.age<-'young';
#######################################################################
##filter tc to cell.types which contain>=100 cells in both young and old
#x=as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age) %>% summarise(n=n())
#tcs=names(which(table(x[x$n>=100,]$tissue_cell.type)==2)) #76 tc
x=as.data.frame(cell.meta) %>% group_by(tissue_cell.type,age) %>% summarise(n=n())
tcs=names(which(table(x[x$n>=20,]$tissue_cell.type)==3)) #115 tc
#tcs=tcs[grep('Marrow|Heart',tcs)] # example tissue for now

df.expr=df.expr[,cell.meta$tissue_cell.type %in% tcs]
cell.meta=cell.meta[cell.meta$tissue_cell.type %in% tcs,]
dim(df.expr);dim(cell.meta) 


###############################################################################################
## perform DR (dimension reduction)
table(cell.meta$tissue_cell.type,cell.meta$binary.age)
all.tcs=sort(unique(cell.meta$tissue_cell.type))
#my.colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- brewer.pal(8,'Dark2')

pdf('PCA_high.var.genes.pdf',useDingbats = T,width = 12,height = 5)
plots=list();
out=list();mse.out=list();
for(tc in all.tcs){
  expression_matrix=df.expr[,cell.meta$tissue_cell.type==tc]
  gene.names <- rownames(expression_matrix)
  cell.names <- colnames(expression_matrix)
  cell_metadata<-cell.meta[cell.meta$tissue_cell.type==tc,]
  rownames(cell_metadata)=cell.names
  
  # gene filter
  geneFilter<-Matrix::rowSums(expression_matrix>=3)>=10
  expression_matrix <- expression_matrix[geneFilter, ]
  
  sce <- CreateSeuratObject(counts = expression_matrix)
  sce@meta.data=as.data.frame(cell_metadata)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e6) #ln(CPM+1)
  
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  #sce <- FindVariableFeatures(sce, selection.method = "mean.var.plot", mean.cutoff=c(0.5,7),dispersion.cutoff =c(0.5,Inf) )
  #VariableFeaturePlot(sce)
  var.genes <- VariableFeatures(sce)
  #length(var.genes)
  Idents(sce)=sce$age #color by age
  
  # PCA
  #sce <- ScaleData(sce, features = rownames(sce))
  sce <- ScaleData(sce, features = var.genes)
  sce <- RunPCA(sce, features =var.genes,npcs = 50)
  p1=DimPlot(sce, reduction = "pca")+ggtitle(tc)
  p2=ElbowPlot(sce,ndims = 50)+ggtitle(tc)
  # UMAP
  sce <- RunUMAP(sce, dims = 1:50) #use all 50PC
  p3=DimPlot(sce, reduction = "umap")+ggtitle(tc)
  print( grid.arrange(p1,p2,p3,ncol=3) )
  if(F){
    
    # MDS (https://satijalab.org/seurat/articles/dim_reduction_vignette.html)
    d <- dist(t(GetAssayData(sce, slot = "scale.data")))
    # Run the MDS procedure, k determines the number of dimensions
    mds <- cmdscale(d = d, k = 50)
    # cmdscale returns the cell embeddings, we first label the columns to ensure downstream consistency
    colnames(mds) <- paste0("MDS_", 1:50)
    # We will now store this as a custom dimensional reduction called 'mds'
    sce[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(sce))
    # We can now use this as you would any other dimensional reduction in all downstream functions
    DimPlot(sce, reduction = "mds", pt.size = 0.5)
    
    
    # NMF
    A = sce@assays$RNA@data[var.genes,] #ngene by ncell
    
    choose.k=500; #as you have 2000 var.genes, 500 always works
    nmf_model<-RcppML::nmf(A,k=choose.k,tol=1e-5)
    dim(nmf_model@w) #
    dim(nmf_model@h) # k x ncell
    length(nmf_model@d) #k
    nmf.h=t(nmf_model@h)
    colnames(nmf.h) <- paste0("NMF_", 1:choose.k)
    sce[["nmf"]] <- CreateDimReducObject(embeddings = nmf.h, key = "NMF_", assay = DefaultAssay(sce))
    # We can now use this as you would any other dimensional reduction in all downstream functions
    #DimPlot(sce, reduction = "nmf", pt.size = 0.5)
    
    mse=sapply(1:choose.k, function(i){
      #R_mse <- mean((A - nmf_model$w %*% Diagonal(x = nmf_model$d) %*% nmf_model$h)^2)
      R_mse <- mean((A - nmf_model$w[,1:i] %*% Diagonal(x = nmf_model$d[1:i]) %*% nmf_model$h[1:i,])^2)
      R_mse
    })
    #plot(1:choose.k,mse,main=tc)
    
    out[[tc]]=nmf_model
    mse.out[[tc]]=mse
  }
  cat('tc',tc,'is done\n')
}
#saveRDS(out,'nmf_marrow.heart.rds')
#saveRDS(mse.out,'nmf.mse_marrow.heart.rds')

dev.off()


