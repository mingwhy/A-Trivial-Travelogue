
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
options(future.globals.maxSize = 3.2 * 1024^3)

#################################################################
## Load Data: data from GSE127832: 1 drop.seq file and 1 10x file

# txt.gz, 1 drop.seq files
raw_counts<-read.delim('GSE127832_wing/GSM3639558_DGE_Drop-Seq_wing.disc_rep1_Lib1.txt.gz',header = TRUE, sep = "\t", row.names = 1)
dim(raw_counts) #10941   150
raw_counts[1:3,1:3]

# tsv.gz, 10x, 1 file
raw_counts<-read.delim('GSE127832_wing/GSM3639563_DGE_10x_wing.disc_rep1.tsv.gz',
                 header = TRUE, sep = "\t", row.names = 1)
dim(raw_counts) #17559gene x 2554cell
raw_counts[1:2,1:3]
head(rownames(raw_counts))

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
