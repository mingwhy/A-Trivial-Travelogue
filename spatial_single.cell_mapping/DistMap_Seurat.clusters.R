
#2017-The Drosophila embryo at single-cell transcriptome resolution

#https://github.com/rajewsky-lab/distmap
# files download from: https://shiny.mdc-berlin.de/DVEX/
#library(devtools)
#install_github("rajewsky-lab/DistMap")
library(DistMap)
###################################################
## read in 4 types of data
#a matrix with genes as rows and cells as columns.
raw.data<-read.delim('dge_raw.txt.gz', 
                     header = F, sep = "\t", row.names = 1)
dim(raw.data) #8924gene by 1297 cell
raw.data[1:2,1:3]
raw.data=as.matrix(raw.data)

#normalized data
data<-read.delim('dge_normalized.txt.gz', 
                 header = T, sep = "\t", row.names = 1)
dim(data) #8924gene by 1297 cell
data[1:2,1:3]
data=as.matrix(data)

#insitu.matrix
insitu.matrix=read.delim('binarized_bdtnp.csv.gz',
                         header = TRUE, sep = ",")
dim(insitu.matrix) #3039 bin x 84 gene
insitu.matrix[1:3,1:3]
insitu.matrix=as.matrix(insitu.matrix)

#geometry, a matrix containing the cartesian coordinates of each bin in three dimensional space.
# use geometry file on github: https://github.com/rajewsky-lab/distmap 
geometry=read.delim('geometry_3039.txt',
                    header = TRUE,sep=' ')
dim(geometry) #3039 bin x 3 coords (the as insitu.matrix)
geometry[1:3,1:3]
geometry=as.matrix(geometry)

#####################################################
## compare lib size between raw.data and norm.data
cell.lib=Matrix::colSums(raw.data)
head(cell.lib)

cell.lib2=Matrix::colSums(data)
head(cell.lib2)

max(cell.lib);max(cell.lib2)

cell.lib3=Matrix::colSums(2^(data)-1)
max(cell.lib3) #I guess the normalized data is log2(x/colSum(x)*scaling.factor+1) transformed 

##########################################################################################
## Seurat generate clusters: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
library(Seurat)

head(colnames(raw.data))
head(colnames(data))
colnames(raw.data)=colnames(data)

fly.embryo <- CreateSeuratObject(counts = raw.data, 
        project = "fly.embryo", min.cells = 0, min.features = 0)
fly.embryo #8924 features across 1297 samples within 1 assay 
head(rownames(fly.embryo))
head(colnames(fly.embryo))

# fill in slot with read in normalized data
fly.embryo@assays$RNA@data=data;

#Identification of highly variable features (feature selection)
fly.embryo <- FindVariableFeatures(fly.embryo, selection.method = "vst", nfeatures = 2000)
length(fly.embryo@assays$RNA@var.features) #2000 selected features

# scaling the data as PCA requires this step 
all.genes <- rownames(fly.embryo)
fly.embryo <- ScaleData(fly.embryo, features = all.genes)
dim(fly.embryo@assays$RNA@scale.data) #8924gene by 1297cell

# Perform linear dimensional reduction
fly.embryo <- RunPCA(fly.embryo, features = VariableFeatures(object = fly.embryo))
fly.embryo@reductions$pca
fly.embryo@reductions$pca@cell.embeddings
print(fly.embryo[["pca"]], dims = 1:5, nfeatures = 5)
# to DimPlot, make sure, matrix colnames == meta.data rownames
head(Idents(fly.embryo))
DimPlot(fly.embryo, reduction = "pca")

#Run non-linear dimensional reduction: tSNE
fly.embryo <- RunTSNE(fly.embryo, dims = 1:16)
DimPlot(fly.embryo, reduction = "tsne")

# Cluster the cells
fly.embryo <- FindNeighbors(fly.embryo, dims = 1:16)
dim(fly.embryo@graphs$RNA_nn)
dim(fly.embryo@graphs$RNA_snn)

fly.embryo=FindClusters(fly.embryo,resolution = 0.5)
head(fly.embryo$RNA_snn_res.0.5)
fly.embryo=FindClusters(fly.embryo,resolution = 0.4)
head(fly.embryo$RNA_snn_res.0.4)

fly.embryo=FindClusters(fly.embryo,resolution = 0.8)
head(fly.embryo$seurat_clusters)

DimPlot(fly.embryo,reduction = 'tsne')

head(fly.embryo@meta.data)
saveRDS(fly.embryo@meta.data,file='embryo_clusters.rds')


