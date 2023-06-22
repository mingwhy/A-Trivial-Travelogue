
library(ggplot2)
library(Seurat)

########################################################
## read in gene chro info for single cell embryo
df.gene.table=data.table::fread('~/Documents/Data_fly_FCA//embryo/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
table(df.gene.table$LOCATION_ARM) #15 Y-chromosome genes

df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol


########################################################
## read in 2017 fly embryo single-cell dataset
tissue='embryo'
dat=readRDS('~/Documents/Data_fly_FCA//embryo/whole_embryo_filtered_valid.rds')

mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)
umi.mat=as.matrix(mat)

# use gene chr info to filter genes
# (some genes annotated to non A or X chr, discard these genes)
i=rownames(umi.mat) %in% df.gene.table$current_symbol
sum(i);#9705
umi.mat=umi.mat[i,]
dim(umi.mat) #9705 1297

pick.genes=df.gene.table[df.gene.table$LOCATION_ARM=='Y',]$SYMBOL
chrY.expr=Matrix::colSums(umi.mat[c(pick.genes),])
sum(chrY.expr==0) #1075
sum(chrY.expr!=0) #222

##########################################################################
## 2021 single-cell germline cell 

## read in gene chro info
df.gene.table=data.table::fread('~/Documents/Data_fly_FCA/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM) #29 chrY genes

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol
if(F){
tissue='embryoSexed'
dat1=readRDS('~/Documents/Data_fly_FCA/embryo_germline/whole_embryo.germline_female_filtered_valid.rds')
dat2=readRDS('~/Documents/Data_fly_FCA/embryo_germline/whole_embryo.germline_male_filtered_valid.rds')

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
table(dat$sex)#female 10386, male 6756   
saveRDS(dat,'integrated_germline.sexed.samples.rds')
}

dat=readRDS('integrated_germline.sexed.samples.rds')
mat <- dat@assays$RNA@counts
gene.names <- rownames(mat)
cell.names <- colnames(mat)
sex=dat$sex
cell.type='unknown';

umi.mat=as.matrix(mat)

# use gene chr info to filter genes
# (some genes annotated to non A or X chr, discard these genes)
i=rownames(umi.mat) %in% df.gene.table$current_symbol
sum(i);#12185
dim(umi.mat) #12224 17142
umi.mat=umi.mat[i,]
dim(umi.mat) #12211 17142

pick.genes=df.gene.table[df.gene.table$LOCATION_ARM=='Y',]$SYMBOL
pick.genes #29 genes
pick.genes=rownames(umi.mat) [rownames(umi.mat) %in% pick.genes]
chrY.expr=Matrix::colSums(umi.mat[pick.genes,])
sum(chrY.expr==0) #16405
sum(chrY.expr!=0) #737
dim(umi.mat)
length(dat$sex)
chrY.expr.male=Matrix::colSums(umi.mat[pick.genes,dat$sex=='embryoMale'])
chrY.expr.female=Matrix::colSums(umi.mat[pick.genes,dat$sex!='embryoMale'])
sum(chrY.expr.male==0)/length(chrY.expr.male) #0.8932
sum(chrY.expr.female==0)/length(chrY.expr.female) #0.9984

#########################################################
## some preliminary analysis with seurat
dat.both=readRDS('integrated_germline.sexed.samples.rds')
dim(dat.both) # 12224 17142
sum(colnames(dat.both)==colnames(dat)) #17142
table(dat.both$sex)
#embryoFemale   embryoMale 
#10386         6756 

dat.both<- NormalizeData(dat.both, verbose = FALSE)
dat.both <- FindVariableFeatures(dat.both, selection.method = "vst", 
                                 nfeatures = 10000, 
                                 verbose = TRUE)
#output of nfeatures=5000 or 10000 is consistent with monocle3 UMAP result (original publication use monocle3 in data analysis)
dat.both <- ScaleData(dat.both, features =rownames(dat.both))
dat.both <- RunPCA(dat.both,features=VariableFeatures(object = dat.both), npcs = 100)
dat.both <- RunUMAP(dat.both, dims = 1:100)

DimPlot(dat.both, group.by = 'sex', reduction = "umap",
        cols = c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))

seurat.umap=dat.both@reductions$umap@cell.embeddings
dim(seurat.umap) #17142     2
seurat.umap=cbind(seurat.umap,dat.both@meta.data)

#saveRDS(seurat.umap,'umap.rds')

pdf('integrated_germline.sexed.samples_umap.pdf',useDingbats = T)
ggplot(seurat.umap,aes(x=UMAP_1,y=UMAP_2,col=sex))+
  geom_point(pch=16,size=0.5)+scale_color_manual(values=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))+
  theme_classic()
dev.off()

#########################################################################
## use monocle3 to integrate data as in the original 2021 publication
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

#######################################################
## check Sxl, msl-2 expression
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

################################################################
