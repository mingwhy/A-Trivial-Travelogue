
library(ggplot2)
library(Seurat)
library(gridExtra)
########################################################
## read in gene chro info
df.gene.table=data.table::fread('~/Documents/Data_fly_FCA/embryo_germline/gene.meta_embryo.txt',header=T,sep='\t')
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

dat.both=readRDS('~/Documents/to_be_removed/sc_sex.differences/embryo_sex.cells/integrated.sexed.samples_seurat.obj.rds')
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

yup=ceiling(max(df$Sxl)*1.05)
p1=ggplot(df,aes(x=`msl-2`,y=`Sxl`,col=sex))+
  ylim(0,yup)+
  geom_jitter(size=0.5)+theme_classic()+ggtitle('log1p normalization')

p2=ggplot(subset(df,sex=='embryoFemale'),aes(x=`msl-2`,y=`Sxl`))+
  ylim(0,yup)+
  #geom_jitter(size=0.2)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  ggtitle('Female')+theme_classic()

p3=ggplot(subset(df,sex=='embryoMale'),aes(x=`msl-2`,y=`Sxl`))+
  ylim(0,yup)+
  #geom_jitter(size=0.2)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  ggtitle('Male')+theme_classic()

pdf('sex.samples_Sxl_msl-2.pdf',useDingbats = T,width = 16,height = 5)
grid.arrange(grobs=list(p1,p2,p3),ncol=3)
dev.off()

########################################################
## pseudobulk: combine 20 cells into one pseudocell
class(df)
df.list=list()
for(sex in unique(df$sex)){
  tmp=df[df$sex==sex,]
  tmp=tmp[,-3]
  ncell=floor(nrow(tmp)/20)
  df.pse=c()
  for(i in 1:(ncell+1)){
    start=20*(i-1)+1;
    end=start+20;
    if(end>nrow(tmp)){end=nrow(tmp)}
    j=Matrix::colSums(tmp[start:end,])
    df.pse=rbind(df.pse,j)
  }
  df.pse=as.data.frame(df.pse)
  df.pse$sex=sex;
  df.list[[sex]]=df.pse
}
df.pse=as.data.frame(rbind(df.list[[1]],df.list[[2]]))

max(df.pse$Sxl)
p1=ggplot(df.pse,aes(x=`msl-2`,y=`Sxl`,col=sex))+
  ylim(0,40)+
  geom_jitter(size=0.5)+theme_classic()+ggtitle('log1p normalization')


p2=ggplot(subset(df.pse,sex=='embryoFemale'),aes(x=`msl-2`,y=`Sxl`))+
  #geom_jitter(size=0.2)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  
  ggtitle('Female')+theme_classic()

p3=ggplot(subset(df.pse,sex=='embryoMale'),aes(x=`msl-2`,y=`Sxl`))+
  #geom_jitter(size=0.2)+
  ylim(0,40)+
  geom_bin_2d()+scale_fill_continuous(type = "viridis") +
  ggtitle('Male')+theme_classic()

pdf('pseudobulk_sex.samples_Sxl_msl-2.pdf',useDingbats = T,width = 16,height = 5)
grid.arrange(grobs=list(p1,p2,p3),ncol=3)
dev.off()

