
library(ggplot2)
library(Seurat)
library(gridExtra)
########################################################
## read in gene chro info
df.gene.table=data.table::fread('~/Documents/single.cell_datasets/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol

##############################################################
## integrate female and male seurat object, then read in dataset
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

dat.both=readRDS('../../../single.cell_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
dat.both; #12185 genes x 17142 cells
grep("Sxl|msl-2",rownames(dat.both))

sub.dat=dat.both[grep("Sxl|msl-2",rownames(dat.both))]
dim(sub.dat) #17142 cells
sub.dat@assays$RNA@counts[1:2,1:3]
sub.dat@assays$RNA@data[1:2,1:3]
df=as.data.frame(cbind(t(as.matrix(sub.dat@assays$RNA@data)),sub.dat$sex))
head(df)
df[,1]=as.numeric(df[,1])
df[,2]=as.numeric(df[,2])
colnames(df)[3]='sex'

p1=ggplot(df,aes(x=`msl-2`,y=`Sxl`,col=sex))+
  geom_jitter(size=0.5)+theme_classic()+ggtitle('log1p normalization')

p2=ggplot(subset(df,sex=='embryoFemale'),aes(x=`msl-2`,y=`Sxl`))+
  #geom_jitter(size=0.2)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  ggtitle('Female')+theme_classic()

p3=ggplot(subset(df,sex=='embryoMale'),aes(x=`msl-2`,y=`Sxl`))+
  #geom_jitter(size=0.2)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  ggtitle('Male')+theme_classic()

pdf('sex.samples_Sxl_msl-2.pdf',useDingbats = T,width = 16,height = 5)
grid.arrange(grobs=list(p1,p2,p3),ncol=3)
dev.off()

####################################################3
DimPlot(dat.both, group.by = 'sex', reduction = "umap",
        cols = c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))

seurat.umap=dat.both@reductions$umap@cell.embeddings
dim(seurat.umap) #17142     2
seurat.umap=cbind(seurat.umap,dat.both@meta.data)

#saveRDS(seurat.umap,'integrated.sexed.samples_seurat.umap.rds')

ggplot(seurat.umap,aes(x=UMAP_1,y=UMAP_2,col=sex))+
  geom_point(pch=16,size=0.5)+scale_color_manual(values=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))+
  theme_classic()

########################################################
## use monocle3 to integrate data
library(monocle3)
expression_matrix=dat.both@assays$RNA@counts
class(expression_matrix)
gene.names <- rownames(expression_matrix)
cell.names <- colnames(expression_matrix)

cell_metadata<-dat.both@meta.data
rownames(cell_metadata)=cell.names
gene_annotation<-data.frame('gene_short_name'=gene.names)
rownames(gene_annotation)=gene.names

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds,reduction_method ='UMAP')

plot_cells(cds, label_groups_by_cluster=FALSE)

# cluster cells, each cell is assigned to a cluster and a partition
cds <- cluster_cells(cds,cluster_method='leiden')

plot_cells(cds, color_cells_by = "partition")
# Learn the trajectory graph
#cds <- learn_graph(cds)
#plot_cells(cds,color_cells_by = "partition",label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
ciliated_genes <- c("pgc","gcl",'nos','vas')
ciliated_genes <- c("pgc","gcl",'nos','vas','fs(1)h','vig2') #female, male marker
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

#https://github.com/cole-trapnell-lab/monocle3/blob/master/R/reduce_dimensions.R
S=reducedDims(cds) #extract cell coordinates on UMAP
dim(S$UMAP)  # 17142     2
head(S$UMAP)
monocle3.umap=S$UMAP;
monocle3.umap=cbind(monocle3.umap,cell_metadata)
colnames(monocle3.umap)[c(1,2)]=c('UMAP_1','UMAP_2')
saveRDS(monocle3.umap,'integrated.sexed.samples_monocle3.umap.rds')

pdf('integrated.sexed.samples_umap.pdf')
print( ggplot(seurat.umap,aes(x=UMAP_1,y=UMAP_2,col=sex))+
  geom_point(pch=16,size=0.5)+scale_color_manual(values=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))+
  theme_classic()+ggtitle('seurat') )
FeaturePlot(dat.both, features = ciliated_genes)

print( ggplot(monocle3.umap,aes(x=UMAP_1,y=UMAP_2,col=sex))+
  geom_point(pch=16,size=0.5)+scale_color_manual(values=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))+
  theme_classic()+ggtitle('monocle3')  )
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

