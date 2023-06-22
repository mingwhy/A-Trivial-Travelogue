library(R.matlab)

genes=readRDS('../brain_scRNA-seq/all_common_genes.rds')
#genes=readRDS('../brain_snRNA-seq/all_common_genes.rds')
length(genes$common.genes) #865

write.table(genes$common.genes,file='geneid.txt',quote=F,row.names = F,col.names = F)

files=Sys.glob('../brain_scRNA-seq/pearson_coexpr.mat_common/*.mat');
#files=Sys.glob('../brain_snRNA-seq/pearson_coexpr.mat_common/*.mat');
files
for(file in files){
  x=basename(file)
  cell=gsub('.mat','.txt',x)
  #cell=gsub('\\*','-star-',cell) #A.B*-KC_perb.mat"
  cell=gsub('\\s','_',cell) #matlab file names don't allow space
  
  dat=readMat(file)
  dim(dat$cell.type)
  #write.table(dat$cell.type,file=paste0('sc_brain_58/',cell),quote=F,row.names = F,col.names = F)
  #write.table(dat$cell.type,file=paste0('sn_brain_100/',cell),quote=F,row.names = F,col.names = F)
  #write.table(dat$cell.type,file=paste0('sn_500cells/',cell),quote=F,row.names = F,col.names = F)
  write.table(dat$cell.type,file=paste0('sc_500cells/',cell),quote=F,row.names = F,col.names = F)
  
}
#txtfiles=Sys.glob('sc_brain_58/*pearson.txt')
#txtfiles=Sys.glob('sn_500cells/*pearson.txt')
txtfiles=Sys.glob('sc_500cells/*pearson.txt')
txtfiles2=sapply(txtfiles,basename)
write.table(txtfiles2,'datasets_list.txt',quote=F,row.names = F,col.names = F)
# mv datasets_list.txt and geneid.txt to 'sn_brain_100/'

######################
