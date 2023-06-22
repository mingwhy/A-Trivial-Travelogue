

file='./common100_genes_brain_scRNA-seq/ncell_1005_male_5.rds'
mat=readRDS(file)
expr.data=mat$mat.count
dim(expr.data) #914 gene in 1005 cell
gene.names=rownames(expr.data)

source('cal_bigScale2.R')

system.time( {
  out=compute.network(expr.data,gene.names,modality='pca',model=NA,
          clustering='direct',quantile.p=0.9,speed.preset='fast')
})
#user  system elapsed 
#11.559   0.227  11.843

# 11seconds
# 11* 2cell.clusters x 2reps x 1thresholds ~ 44s
# 11* 120 x 20reps x 6thresholds ~ 2640min=44hr 
# 11*120*10*6/60/60 ~ 22hr
# 11*5*10*5/60/60 ~ 1hr
system.time( {
  out=compute.network(expr.data,gene.names,modality='pca',model=NA,
                      clustering='recursive',quantile.p=0.9,speed.preset='fast')
})
#user  system elapsed 
#75.906   3.396  84.207 
#84*10*10*5/60/60
