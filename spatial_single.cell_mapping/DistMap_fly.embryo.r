
#2017-The Drosophila embryo at single-cell transcriptome resolution

#https://github.com/rajewsky-lab/distmap
# files download from: https://shiny.mdc-berlin.de/DVEX/
#library(devtools)
#install_github("rajewsky-lab/DistMap")
library(DistMap)
###################################################
## read in 4 types of data
#a matrix with genes as rows and cells as columns.
raw.data<-read.delim('./DVEX/dge_raw.txt.gz', 
                header = F, sep = "\t", row.names = 1)
dim(raw.data) #8924gene by 1297 cell
raw.data[1:2,1:3]
raw.data=as.matrix(raw.data)

#normalized data
data<-read.delim('DVEX/dge_normalized.txt.gz', 
                header = T, sep = "\t", row.names = 1)
dim(data) #8924gene by 1297 cell
data[1:2,1:3]
data=as.matrix(data)

#insitu.matrix
insitu.matrix=read.delim('DVEX/binarized_bdtnp.csv.gz',
                         header = TRUE, sep = ",")
dim(insitu.matrix) #3039 bin x 84 gene
insitu.matrix[1:3,1:3]
insitu.matrix=as.matrix(insitu.matrix)

#geometry, a matrix containing the cartesian coordinates of each bin in three dimensional space.
# use geometry file on github: https://github.com/rajewsky-lab/distmap 
geometry=read.delim('DVEX/geometry_3039.txt',
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
max(cell.lib3)

#####################################################
# with geometry, plot one marker gene in spatial
# NovoSpaRc: https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_embryo_tutorial.ipynb
plot(geometry[,1],geometry[,2],pch=16) #x and y
plot(geometry[,2],geometry[,3],pch=16) #y and z
plot(geometry[,1],geometry[,3],pch=16) #x and z (use this)

dim(geometry)  #3039 x 3
dim(insitu.matrix) # 3039 x 84 marker gene
i.gene=81; #pick one gene
genename=colnames(insitu.matrix)[i.gene] #twi
tmp=cbind(geometry,insitu.matrix[,i.gene])
plot(tmp[,1],tmp[,2],col=tmp[,4],pch=16)
# compare with https://shiny.mdc-berlin.de/DVEX/

library(ggplot2)
tmp=as.data.frame(tmp)
unique(tmp[,4]) #binary: 0 and 1

# lateral view  
ggplot(tmp,aes(x=tmp[,1],tmp[,3],col=factor(tmp[,4])))+
  geom_point()+theme_classic()+ggtitle(paste('gene',genename))

######################################################
## check gene names, make sure insitu.matrix match data
sum(rownames(raw.data) %in% rownames(data)) #8924
dim(insitu.matrix)
sum(colnames(insitu.matrix) %in% rownames(data)) #82
colnames(insitu.matrix)[!colnames(insitu.matrix) %in% rownames(data)]
#"Blimp.1"      "E.spl.m5.HLH"
genes=rownames(data)
genes[grep('Blimp',rownames(data))]
genes[grep('E.spl.m5.HLH',rownames(data))]

genes2=colnames(insitu.matrix)
grep(c('Blimp|E.spl.m5.HLH'),genes2)
genes2[c(6,36)]
colnames(insitu.matrix)[c(6,36)]=c("Blimp-1","E(spl)m5-HLH")

# set up object
dm = new("DistMap",
         raw.data=raw.data,
         data=data,
         insitu.matrix=insitu.matrix,
         geometry=geometry)

# compute binarized single cell data, binarized.data slot is filled
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dim(dm@binarized.data) # 84 marker gene x 1297 cells

# map cells to reference atals, mcc.scores slot is filled
dm <- mapCells(dm)
dim(dm@mcc.scores) #3039 bin x 1297 cell
dm@mcc.scores[1:3,1:3]
head(apply(dm@mcc.scores,2,max))


mcc=dm@mcc.scores
dim(mcc) #3039 bins X 1297 cells
distmap.mcc=mcc;
data.table::fwrite(distmap.mcc,'distmap.mcc.txt')

# https://github.com/rajewsky-lab/distmap/blob/master/R/DistMap.R
# in 2017 paper supplements
# We employed the Mathews correlation coefficient (MCC) to weight the confusion matrices and assign a cell-bin score. 
# The MCC scores were subsequently exponentiated.

distmap.mcc=as.matrix(data.table::fread('distmap.mcc.txt'))
log.distmap.mcc=log(distmap.mcc)

pdf('distmap.mcc_1297cell_3039bin.pdf') #fig2A in 2017 paper
print(pheatmap::pheatmap(log.distmap.mcc, treeheight_row = 0, treeheight_col = ))
dev.off()

#############################################
#############################################
# Once the cells have been mapped, 
# the DistMap functions can be used to compute a vISH 
# or a gradient of a gene and visualize the expression pattern
out=computeVISH(dm, 'sna', threshold=0.75)
length(out) #3039

g=computeGeneGradient(dm, 'sna') #return a ggplot showing gradient
plot(g)


