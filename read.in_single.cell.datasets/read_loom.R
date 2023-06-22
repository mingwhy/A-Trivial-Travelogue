## use SCopeLoomR to read loom file
if(F){remotes::install_github("aertslab/SCopeLoomR")}
library(SingleCellExperiment)
library(SCopeLoomR)

if(T){
# data download from https://flycellatlas.org/
loom_path <- './s_fca_biohub_leg_10x.loom';
loom <- open_loom(loom_path, mode="r+")

cell.annotation.all=get_cell_annotation(loom)
dim(cell.annotation.all) #11788   435

labels=colnames(cell.annotation.all)
tmp=cell.annotation.all[,-grep('TrackRegulonsAUC|MotifRegulonsAUC',labels)]
colnames(tmp)
table(tmp$batch_id)
sort(table(tmp$batch))
sort(table(tmp$id))  #only two samples
sort(table(tmp$batch))==sort(table(tmp$id)) #batch <=> id

write.table(tmp,file='cell.annotation.csv',quote=F,sep='\t',row.names = F) #too big

cell.annotation.df=tmp
colnames(cell.annotation.df)
head(cell.annotation.df$annotation)
table(cell.annotation.df$annotation)

genes=get_genes(loom)
length(genes) #13411 genes

raw <- get_dgem(loom)
raw[1:5,1:5]
dim(raw) #13203 gene by 37254 cell

close_loom(loom)

## filter condition
nGeneLowCutOff <- 200; nGeneHighCutOff <- Inf
nUMILowCutOff <- 400; nUMIHighCutOff <- Inf
MitoLowCutOff <- -Inf; MitoHighCutOff <- 0.30
min.cell <- 4; # keep gene expresses at >=4 cells

## filter step: 1) calculate mito% for each cell, remove cell whose mito.perc > 0.3
# 2) remove mito gene
# 3) filter cell: gene (<200), umi (<400)
# 4) filter gene:  discard genes which express <4 cells

# 1) calculate mito% for each cell, remove cell whose mito.perc > 0.3
(mito.genes <- grep(pattern = "^mt:",ignore.case = T,x = rownames(raw), value = TRUE))
# 24 mito genes
prop.mito <- Matrix::colSums(raw[mito.genes, ]) / Matrix::colSums(raw)
summary(prop.mito)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.001161 0.002257 0.003319 0.004043 0.049934 
summary(cell.annotation.df$percent_mito) #confirm

i=prop.mito>MitoHighCutOff
cat(sum(i),'cells fail MitoHighCutOff\n');
if(sum(i)!=0){
  raw=raw[,!i];
  cell.annotation.df=cell.annotation.df[,!i]
}

# 2) remove mito gene
sum(rownames(raw) %in% mito.genes) #22 mito gene
raw<-raw[!rownames(raw) %in% mito.genes,]
dim(raw) #13181gene, 37254cell
cat(sum(rownames(raw) %in% mito.genes) ,'mito genes removed\n')

# 3) filter cell: gene(<200), umi(<400)
x=Matrix::colSums(raw>0)
sum(x<nGeneLowCutOff) #0, all cells express at least 200 genes 
i=sum(x<nGeneLowCutOff) 
cat(sum(i),'cells fail nGeneLowCutOff\n');

if(sum(i) !=0 ){
  raw=raw[,!i]
  cell.annotation.df=cell.annotation.df[!i,]
}

x=Matrix::colSums(raw)
sum(x<nUMILowCutOff) #0 cells, total UMI<400
i=x<nUMILowCutOff
cat(sum(i),'cells fail nUMILowCutOff\n');

if(sum(i) !=0){
  raw=raw[,!i]
  cell.annotation.df=cell.annotation.df[!i,]
}

# 4) filter gene:  discard genes which express <4 cells
x<-Matrix::rowSums(raw>0)
sum(x==0) #0 gene
sum(x<min.cell) #483 genes, express in less than 4 cells
i=x<min.cell

cat(sum(i),'genes fail nUMILowCutOff\n');

if(sum(i)!=0){
  raw<-raw[!rownames(raw) %in% names(which(x<min.cell)),]
  dim(raw) #12732 37254
  dim(cell.annotation.df) #10686    34
}

library(Seurat)
whole<- CreateSeuratObject(counts = raw,
      min.cells=0, min.features = 0, project = "whole")
whole #12174 features across 10686 samples within 1 assay
head(rownames(whole))
head(colnames(whole))
head(rownames(cell.annotation.df))
x=sum(colnames(whole)==rownames(cell.annotation.df)) #37254
cat('check barcode match,ncell',ncol(whole),x,'matched\n')
whole@meta.data=cell.annotation.df

saveRDS(whole,'./whole_leg_filtered.rds')
}

