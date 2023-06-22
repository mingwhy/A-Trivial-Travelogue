library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra)
library(tidyverse);
library(SummarizedExperiment)
library(ggpointdensity)
library(viridis)
library(patchwork); #for plot_annotation
library(scMerge)

##########################################################################
## select only male cells and cell types with >=50cells in both 3m and 24m
if(!file.exists('TMS_male.h5ad')){
  inp.sce<-readH5AD('~/Documents/Data_mouse_aging_atlas/TMS.gene.data_final/tabula-muris-senis-facs-official-raw-obj.h5ad');       
  colData(inp.sce)$tissue_cell.type=paste(inp.sce$tissue,inp.sce$cell_ontology_class,sep=':')
  
  sub.sce=inp.sce[,inp.sce$sex=='male']
  x=table(sub.sce$tissue_cell.type,sub.sce$age)
  x=as.data.frame(x[,c(1,4)]) #3m and 24m
  colnames(x)=c('tc','age','ncell')
  x=x %>% spread(age,ncell)
  dim(x) #202
  #data.table::fwrite(x,'ncell_per.age_per.tc_raw.txt',sep='\t',quote=F)
  
  y=x[x[,2]>=50 & x[,3]>=50, ] #both age groups contain >=50 cells
  dim(y) #72
  #data.table::fwrite(y,'ncell_per.age_per.tc.txt',sep='\t',quote=F)
  
  pick.cell.types=as.character(y$tc);
  sce=sub.sce[,sub.sce$tissue_cell.type %in% pick.cell.types]
  unique(sce$age) #'3m','18m','21m','24m'
  sce$age=droplevels(sce$age)
  unique(sce$age)
  sce$age=factor(sce$age,levels=c('3m','18m','24m'))
  cell.meta=colData(sce)
  #rowData(sce)
  length(table(cell.meta$tissue_cell.type)) 
  
  #remove 18m data
  sce=sce[,sce$age!='18m']
  sce$age=droplevels(sce$age)
  table(sce$age)
  #3m   24m 
  #21820 26078 
  writeH5AD(sce, 'TMS_male.h5ad')
}

##########################################################################################################
## select HVG. log
if(!file.exists('hvg2k_mat.rds')){
  sce=readH5AD('TMS_male.h5ad')
  
  # normalization
  assayNames(sce) #'X'
  assayNames(sce)<-'counts'
  assayNames(sce) # "counts"
  sce <- logNormCounts(sce,log = TRUE) 
  #logNormCounts from scuttle R package, internal call `normalizeCounts` and `librarySizeFactors` function 
  #it's log2-transformation, https://rdrr.io/github/LTLA/scuttle/man/normalizeCounts.html
  #Library sizes are converted into size factors by scaling them so that their mean across cells is unity. https://rdrr.io/bioc/scuttle/man/librarySizeFactors.html
  assayNames(sce) # "counts"    "logcounts"
  
  # select HVG
  allf <- modelGeneVar(sce) #function from scran #https://rdrr.io/github/MarioniLab/scran/man/modelGeneVar.html
  #plot(allf$mean, allf$total)
  #curve(metadata(allf)$trend(x), add=TRUE, col="dodgerblue")
  allf=allf[order(allf$bio, decreasing=TRUE),]
  pick.genes=rownames(allf)[1:2000] #select top HGV genes
  
  # extract matrix for NMF (HGV gene x ncell, log-normalized)
  input.mat=assay(sce[pick.genes,],'logcounts')
  dim(input.mat) #2000 47898
  tmp=rowMeans(input.mat)
  summary(tmp) #check if there are any gene with 0 expr across cells
  
  saveRDS(input.mat,'hvg2k_mat.rds')
}

input.mat=readRDS('hvg2k_mat.rds')

if(!file.exists('hvg2k_mat.npy')){
  # save a npy or npz version for python
  library(reticulate)
  np=import('numpy')
  x=input.mat
  #np$savez('hvg2k_mat.npz',x=input.mat)
  np$save('hvg2k_mat.npy', x)
}

#######################################################################################################################
## NMF (call python in R, key points: `Sys.setenv` for load modules, `.->$`, as.interger() or 'xxx' to pass parameters)
library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
Sys.setenv(RETICULATE_PYTHON ='/Users/mingyang/anaconda3/bin')
py_config()
sys=import('sys')
sys$path
nimfa=import('nimfa') #open terminal, if `import nimfa` successfully, use `import sys;sys.path` to get the env set as above
np=import('numpy')

if(F){
  #https://nimfa.biolab.si/
  V = nimfa$examples$medulloblastoma$read(normalize='True') #feature by sample matrix
  dim(V) # 5893   34
  lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
  lsnmf_fit = lsnmf()
  
  cat('Iterations:' ,lsnmf_fit$n_iter,'\n')
  cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
  cat('Evar: ', lsnmf_fit$fit$evar(),'\n') #Compute the explained variance of the NMF estimate of the target matrix.
  cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
  cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')
  #https://github.com/mims-harvard/nimfa/issues/54
  
  W = lsnmf_fit$basis()
  dim(W) #5839 x 50, Basis matrix
  H = lsnmf_fit$coef()
  dim(H) #50 x 34, Mixture matrix
  
  target.mat=np$dot(W,H) #Target estimate
  dim(target.mat) #5893   34
}

#big dataset test
#LSNMF on numpy dense matrix with quality and performance measures.
V = as.matrix(input.mat) #change it into as.matrix!!! important
lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
start.time=Sys.time()
lsnmf_fit = lsnmf() 
end.time=Sys.time()
end.time-start.time; # 2000 x 47898, Time difference of 27.97791 mins for sparse matrix
# Time difference of 40.35332 secs for as.matrix() or dense matrix

cat('Iterations:' ,lsnmf_fit$n_iter,'\n')
cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
cat('Evar: ', lsnmf_fit$fit$evar(),'\n')
cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')
#https://github.com/mims-harvard/nimfa/issues/54


#Standard NMF - Divergence on scipy.sparse matrix with matrix factors estimation.
V = input.mat
lsnmf = nimfa$Nmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100),
                  update='euclidean', objective='fro')
start.time=Sys.time()
lsnmf_fit = lsnmf() 
end.time=Sys.time()
end.time-start.time; # 2000 x 47898, Time difference of 49.1204 mins
cat('Iterations:' ,lsnmf_fit$n_iter,'\n')
cat('Rss:' , lsnmf_fit$fit$rss(),'\n')
cat('Evar: ', lsnmf_fit$fit$evar(),'\n')
cat('K-L divergence:' , lsnmf_fit$distance(metric='kl'),'\n')
cat('Sparseness, W: , H: ' , unlist(lsnmf_fit$fit$sparseness()),'\n')
#https://github.com/mims-harvard/nimfa/issues/54

W = lsnmf_fit$basis()
dim(W) #5839 x 50, Basis matrix
H = lsnmf_fit$coef()
dim(H) #50 x 34, Mixture matrix

#######################################################################################################################
# https://nimfa.biolab.si/nimfa.methods.factorization.lsnmf.html
# Latent dimensionality selection
rank_cands = c(20,50,100,200)
V=as.matrix(input.mat)
dim(V) # 2000 47898
lsnmf_estrank = nimfa$Lsnmf(V, seed='random_vcol', max_iter=100) #set up function

start.time=Sys.time() # start at 1:44 pm
summary = lsnmf_estrank$estimate_rank(rank_range=as.integer(rank_cands), n_run=as.integer(2), what='all') #run function
end.time=Sys.time()

#https://nimfa.biolab.si/nimfa.methods.factorization.lsnmf.html
#parameter of Lsnmf `evar`
#Compute the explained variance of the NMF estimate of the target matrix.
#This measure can be used for comparing the ability of models for accurately 
#reproducing the original target matrix. 
#Some methods specifically aim at minimizing the RSS and maximizing the explained variance 
#while others not, which one should note when using this measure.

