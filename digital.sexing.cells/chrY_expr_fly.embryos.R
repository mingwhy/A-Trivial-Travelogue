
library(ggplot2)
library(Seurat)

########################################################
## read in gene chro info for single cell embryo
df.gene.table=data.table::fread('~/Documents/Data_fly_FCA/embryo/gene.meta_embryo.txt',header=T,sep='\t')
head(df.gene.table)
table(df.gene.table$LOCATION_ARM) #15 Y-chromosome genes

df.gene.table=df.gene.table[df.gene.table$LOCATION_ARM %in% c('2R','3R','3L','4','2L','X','Y'),]
table(df.gene.table$LOCATION_ARM)

df.gene.table=as.data.frame(df.gene.table)
rownames(df.gene.table)=df.gene.table$current_symbol


########################################################
## read in 2017 fly embryo single-cell dataset
tissue='embryo'
dat=readRDS('~/Documents/Data_fly_FCA/embryo/whole_embryo_filtered_valid.rds')

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

#/Users/mingyang/Documents/digital.sexing.cells/recipe/TSP/TSP_on_single.cell/00_embryo_integrate.data_chrY.R
dat=readRDS('~/Documents/digital.sexing.cells/recipe/TSP/TSP_on_single.cell/integrated_germline.sexed.samples.rds');
dat #12224 features across 17142 samples within 1 assay 
mat=dat@assays$RNA@counts;
gene.names <- rownames(mat)
cell.names <- colnames(mat)
sex=dat$sex
cell.type='unknown';
table(sex)
#embryoFemale   embryoMale 
#10386         6756 

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
table(sex[which(chrY.expr!=0)])

dim(umi.mat)
length(dat$sex)
chrY.expr.male=Matrix::colSums(umi.mat[pick.genes,dat$sex=='embryoMale'])
chrY.expr.female=Matrix::colSums(umi.mat[pick.genes,dat$sex!='embryoMale'])
sum(chrY.expr.male==0)/length(chrY.expr.male) #0.8932
sum(chrY.expr.female==0)/length(chrY.expr.female) #0.9984

##########################################################################
## 202X Jay data
# read in gene.meta info which contains chr info
gene.meta=data.table::fread('~/Documents/Data_Jay_fly_development/validated_15392genes.txt')
y.genes=gene.meta[gene.meta$chromosome_name=='Y',]$SYMBOL
length(y.genes) #34 y genes

# read in embryo data by pred time window
(files=Sys.glob('~/Documents/Data_Jay_fly_development/RNA_seurat_object/pred_windows/*.rds'))

out=c()
for(file in files){
  #x=readRDS('./RNA_seurat_object/pred_windows/18-20hrs_finished_processing.rds')
  x=readRDS(file)
  expr.mat=x@assays$RNA@counts
  
  dosage.genes=grep('msl|sxl',ignore.case = T,rownames(expr.mat))
  tmp=expr.mat[dosage.genes,]
  dosage.genes.ncell=Matrix::rowSums(tmp>0)
  
  tmp=expr.mat[y.genes,]
  y.cell=Matrix::colSums(tmp) #y.expr per cell
  #sum(y.cell==0)/length(y.cell)
  
  
  time=gsub('_finished_processing.rds','',basename(file))
  x=c(time, length(y.cell), sum(y.cell!=0),dosage.genes.ncell)
  print(x)
  
  out=rbind(out,x)
}

df.out=as.data.frame(out)
colnames(df.out)=c('time.window','ncell','ncell.expr.Y','ncell.expr.msl-3','ncell.expr.Sxl',
                   'ncell.expr.msl-2','ncell.expr.msl-1')
df.out$start=as.numeric(unlist(lapply(strsplit(df.out$time.window,'-'),'[',1)))
df.out=df.out[order(df.out$start),]
df.out
#data.table::fwrite(df.out,'chrY.expr_per.time.window.txt')

