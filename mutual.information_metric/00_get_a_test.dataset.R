
library(Seurat)
options(stringsAsFactors = F)

## read in processed wholebrain data
file="../../single.cell_datasets/FCA_gut/whole_gut_filtered_valid.rds";
dat0=readRDS(file); 
dat0; #12878 features across 11788 samples

# for each cell type, the count or proportion of male, female and mix cell numbers
# remove 'mix' cell types
table(dat0@meta.data$sex)
#female   male 
#6338   5450 
dat=subset(dat0,sex!='mix')
dat # 12878 features across 11788 samples

## pick one cell type
cell.types=unique(dat$annotation)
(cell.type=cell.types[1]) #"enterocyte-like"

dat.test=dat[,dat$annotation==cell.type]
dat.test #12878 features across 985 samples 
saveRDS(dat.test,'test.rds')




