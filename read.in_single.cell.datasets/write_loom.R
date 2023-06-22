
# install loomR (https://github.com/mojaveazure/loomR)
# 1) The pre-requisites for loomR include the C++ HDF5 library
#$brew install hdf5
https://github.com/mojaveazure/loomR
# 2) install loomR following: https://satijalab.org/loomr/loomr_tutorial
# Use devtools to install hdf5r and loomR from GitHub
if(F){
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}
library(loomR)

# Generate a loom object containing sample scRNA-seq data
#if(F){remotes::install_github("aertslab/SCopeLoomR")}
#library(SCopeLoomR)


##read in data and extract three pieces of info: matrix, genes, cells
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/fly.single.cell.datasets/read_txt_tsv_csv.R

## csv.gz, Load Data: data from GSE120537
raw_counts<-read.table(file="GSE120537_midgut/GSE120537_counts.csv.gz",sep=",",
                       header = TRUE, row.names = 1)
raw_counts[1:2,1:3]
dim(raw_counts) # 16986 10605
head(colnames(raw_counts))

cell.info<-read.table(file='GSE120537_midgut/GSE120537_metadata.csv.gz',
                      sep=',',header = TRUE,as.is=T)
dim(cell.info) #10605    16
head(cell.info$cell)
colnames(cell.info)

gene.names=rownames(raw_counts)

## all three pieces of information collected
## generate a sample data stored in loom format
## select two cell types and save a new loom object
pick=names(sort(table(cell.info$library)))[c(2,3)]
pick
#[1]  "gut3" "G2"  

cell.info.sub=cell.info[cell.info$library %in% pick,]
raw.sub=raw_counts[,cell.info$library %in% pick]
dim(cell.info.sub)
#[1] 297  16
dim(raw.sub)
#[1] 16986   297
length(gene.names)
#[1] 16986
sum(rownames(raw.sub)==gene.names)
#[1] 16986

cell.info.sub.list=list(
  cell.type=cell.info.sub$library,
  batch.label=cell.info.sub$batch)

#https://rdrr.io/github/mojaveazure/loomR/man/create.html
library(loomR)
create('midgut_2celltypes_10x.loom',data=raw.sub,cell.attrs=cell.info.sub.list)


## read in small loomR file
library(loomR)
ds <- connect('midgut_2celltypes_10x.loom')
mat <- t(ds$matrix[,])
gene.names <- ds$row.attrs$Gene[]
names(ds$col.attrs)
cell.names <- ds$col.attrs$CellID[]
cell.types<-ds$col.attrs$cell.type[]
batch.labels<-ds$col.attrs$batch[]

names(ds$col.attrs)
cell.info=c()
for(var in names(ds$col.attrs)){
  cell.info=cbind(cell.info,ds$col.attrs[[var]][])
}
colnames(cell.info)=names(ds$col.attrs)
head(cell.info)

ds$close_all()



dim(mat)
length(gene.names)
length(cell.names)
dim(cell.info)


