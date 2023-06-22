
library(zellkonverter)
library(SingleCellExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);library(grid)
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
  cell.meta=colData(sce)
  saveRDS(cell.meta,'TMS_male_colData.rds')
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
## select rank in NMF
input.mat=readRDS('hvg2k_mat.rds')

library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/Users/mingyang/anaconda3/envs/myenv/bin/python")
Sys.setenv(RETICULATE_PYTHON ='/Users/mingyang/anaconda3/bin')
py_config()
sys=import('sys')
sys$path
nimfa=import('nimfa') #open terminal, if `import nimfa` successfully, use `import sys;sys.path` to get the env set as above
np=import('numpy')

#https://nimfa.biolab.si/
#LSNMF on numpy dense matrix with quality and performance measures.
outDir='select_NMF_rank/';
if(!dir.exists(outDir)){dir.create(outDir)}
V = as.matrix(input.mat) #important for computing speed!!
nrep=3;
start.time=Sys.time() #one round with no rep, start at 350pm, end at 510pm, 1.5hr
for(rank in c(20,50,80,100,150,200,250,500,1000,1500,2000)){
  for(irep in 1:nrep){
    output.file=paste0(outDir,'/nmf_',rank,'_rep',irep,'.rds');
    if(!file.exists(output.file)){
      #lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(50), max_iter=as.integer(100))
      lsnmf = nimfa$Lsnmf(V, seed='random_vcol', rank=as.integer(rank), max_iter=as.integer(100))
      start.time=Sys.time()
      lsnmf_fit = lsnmf() 
      end.time=Sys.time()
      cat('rank=',rank,', ')
      print(end.time-start.time); # 2000 x 47898, Time difference of 27.97791 mins
      
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
      tmp=lsnmf_fit$fit$sparseness()
      
      saveRDS(list(W=W,H=H,n_iter=lsnmf_fit$n_iter,Rss=lsnmf_fit$fit$rss(),
        Evar=lsnmf_fit$fit$evar(),kl=lsnmf_fit$distance(metric='kl'),
        W.sparseness=tmp[[1]],H.sparseness=tmp[[2]]),file=output.file)
    }
  }
}
end.time=Sys.time()

(files=Sys.glob('select_NMF_rank/nmf_*rep*.rds'))
out=c()
for(file in files){
  rank=gsub('nmf_|_rep.*.rds','',basename(file))
  rep=gsub('nmf_.*_|.rds','',basename(file))
  x=readRDS(file)
  out=rbind(out,c(rank,rep,x$Evar,x$Rss))
}
df=as.data.frame(out)
colnames(df)=c('rank','rep','evar','rss')
head(df)
df$rank=as.numeric(df$rank)
ggplot(df,aes(x=rank,y=evar,group=rep,col=factor(rank)))+geom_point()+theme_classic()
ggplot(df,aes(x=rank,y=rss,group=rep,col=factor(rank)))+geom_point()+theme_classic()

df[which.max(df$evar),]

# choose rank=250
cell.embed=readRDS('select_NMF_rank/nmf_250_rep2.rds')
dim(cell.embed$W) #2000 gene by 50 dims
dim(cell.embed$H) #50 dim by 47898 cells
cell.coord=cell.embed$H

cell.meta=readRDS('TMS_male_colData.rds')
colnames(cell.coord)=rownames(cell.meta); #barcode
cell.meta$sex=droplevels(cell.meta$sex) #only male data
tc.names=sort(unique(cell.meta$tissue_cell.type))
unique(cell.meta$age) #3m and 24m

###########################################################################
## UMAP or tSNE 2D plot of cells with 50 NMF embedding (both ages included)
dim(cell.coord) # 50 47898

if(!file.exists('nmf250_umap_tsne.rds')){
  x=scater::calculateUMAP(cell.coord) 
  df.reduce=as.data.frame(x)
  x=scater::calculateTSNE(cell.coord)
  df.tsne=as.data.frame(x)
  saveRDS(list(umap=df.reduce,tsne=df.tsne),file='nmf250_umap_tsne.rds')
}

x=readRDS('nmf250_umap_tsne.rds')
df.reduce=x$umap;
#df.reduce=x$tsne; 
colnames(df.reduce)=c('dim1','dim2')
df.reduce$tissue=cell.meta$tissue
df.reduce$tissue_cell.type=cell.meta$tissue_cell.type
df.reduce$cell.type=cell.meta$cell_ontology_class
df.reduce$age=cell.meta$age
df.reduce$cell.type=droplevels(df.reduce$cell.type)
unique(df.reduce$cell.type) #52 cell types

p0=ggplot(df.reduce,aes(x=dim1,y=dim2,col=cell.type))+
  geom_point(size=0.5)+theme_classic()+
  scale_color_viridis(name='cell.state',option='turbo',alpha=0.6,discrete = T)

df.reduce.sub=df.reduce[grep('endo|Brain',ignore.case = T,df.reduce$tissue_cell.type),]
#unique(df.reduce.sub$tissue_cell.type);unique(df.reduce.sub$cell.type)
df.reduce.sub$tissue2='non-Brain'
df.reduce.sub$tissue2[grep('Brain',ignore.case = T,df.reduce.sub$tissue)]='Brain';
df.reduce.sub$cell.identity=as.character(df.reduce.sub$cell.type)
df.reduce.sub$cell.identity[grep('endothelia',ignore.case = T,df.reduce.sub$tissue_cell.type)]='endothelia cell';

p=ggplot(df.reduce.sub,aes(x=dim1,y=dim2))+
  theme_classic(base_size = 15)+
  guides(color = guide_legend(override.aes = list(size = 8)),
         shape = guide_legend(override.aes = list(size = 8)))+
  scale_color_viridis(option='turbo',alpha=0.4,discrete = T)
  #theme(legend.position = 'bottom')
p+geom_point(aes(col=tissue_cell.type,shape=tissue2),size=4)
p+geom_point(aes(col=cell.identity,shape=tissue2),size=4)+
  scale_color_manual(values=grDevices::adjustcolor(c("#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"), alpha.f = 0.3))
p+geom_point(aes(col=age,shape=tissue2),size=4)+
  scale_color_manual(values=grDevices::adjustcolor(c("#009E73","#F0E442"),alpha=0.2))
# plot UMAP for 3m and 24m separately, maybe

pdf('umap.pdf',width = 9,useDingbats = T)
#pdf('tsne.pdf',width = 9,useDingbats = T)
print(p0+theme(legend.position = 'none'));
legend <- cowplot::get_legend(p0+ 
              guides(color = guide_legend(override.aes = list(size = 8)))+
              scale_color_viridis(name='cell state',option='turbo',discrete=T,alpha=0.6))
grid.newpage()
grid.draw(legend)
# theme(legend.position = 'bottom'))
print(p+geom_point(aes(col=tissue_cell.type,shape=tissue2),size=4))
p+geom_point(aes(col=cell.identity,shape=tissue2),size=4)+
  scale_color_manual(values=grDevices::adjustcolor(c("#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"), alpha.f = 0.3))
p+geom_point(aes(col=age,shape=tissue2),size=4)+
  scale_color_manual(values=grDevices::adjustcolor(c("#009E73","#F0E442"),alpha=0.2))
dev.off()

########################################################
## for cell.types in multiple tissues, variance in NMF distribution explained by tissue, cell type and their interaction
library('variancePartition');library(tidyverse)

cell.meta=readRDS('TMS_male_colData.rds');
cell.meta$cell.identity=cell.meta$cell_ontology_class
#unique(cell.meta[grep('B cell',cell.meta$cell.identity),]$cell.identity)
cell.meta[grep('endothelia',cell.meta$cell.identity),]$cell.identity='endothelial cell';
cell.meta[grep('T cell',cell.meta$cell.identity),]$cell.identity='T cell';
#unique(cell.meta[cell.meta$tissue=='Thymus',]$cell_ontology_class)
cell.meta[cell.meta$tissue=='Thymus',]$cell.identity='T cell';

cell.meta$cell.identity=droplevels(cell.meta$cell.identity)
tmp=cell.meta %>% as.data.frame() %>% group_by(cell.identity,tissue) %>% 
  summarise(ncell.3m=sum(age=='3m'),ncell.24m=sum(age=='24m'))
tmp

as.data.frame(sort(table(tmp$cell.identity),decreasing = T)) %>% 
  rename(cell.identity = Var1, show.up.in.n.tissue=Freq) %>% 
  filter(show.up.in.n.tissue>1) %>% select(cell.identity) ->pick.cell.types
pick.cell.types=as.character(unlist(pick.cell.types))

cell.meta.sub=cell.meta[cell.meta$cell.identity %in% pick.cell.types,]
#cell.meta.sub=cell.meta[grep('endothelial|B cell',ignore.case = T,cell.meta$tissue_cell.type),]
#cell.meta.sub=cell.meta.sub[-grep('Marrow:',cell.meta.sub$tissue_cell.type),]
unique(cell.meta.sub$tissue_cell.type); #31 tissues
unique(cell.meta.sub$tissue) #16 tissues
cell.meta.sub$cell.identity=droplevels(cell.meta.sub$cell.identity)
unique(cell.meta.sub$cell.identity) #8 

cell.coord.sub=cell.coord[,rownames(cell.meta.sub)]
dim(cell.coord.sub) #feature by sample matrix 
cell.meta.sub$tissue=droplevels(cell.meta.sub$tissue)

cell.meta.sub$age=as.numeric(gsub('m','',cell.meta.sub$age)) #for LMM to work
form <- ~ age + (1|tissue)+(1|cell.identity) 
varPart <- fitExtractVarPartModel( cell.coord.sub, form, cell.meta.sub )

dim(varPart) #50  3
head(varPart)
vp <- sortCols( varPart )
plotPercentBars( vp[1:20,] ) + xlab('Dimensions from NMF')
fig=plotVarPart(vp)

pdf('var_explained.pdf',useDingbats = T,width = 9)
grid.arrange(
  plotPercentBars( vp[1:20,] ) + xlab('Dimensions from NMF'),fig,ncol=2)
dev.off()

#######################################################################################################################
## compute aging vectors per tc: Centroid(Young)-> Centroid(Old)
age.vectors=lapply(tc.names,function(tc){
  young= cell.meta$tissue_cell.type==tc & cell.meta$age=='3m'
  old= cell.meta$tissue_cell.type==tc & cell.meta$age=='24m'
  young.center=Matrix::rowMeans(cell.coord[,young])
  old.center=Matrix::rowMeans(cell.coord[,old])
  age.vector=old.center- young.center
  return(age.vector)
})
df=as.data.frame(Reduce(`cbind`,age.vectors))
colnames(df)=tc.names
rownames(df)=paste('dim',1:nrow(df))
#pheatmap::pheatmap(df)

# compute cosine similarities between the aging trajectories of each cell type in each tissue
# https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
X=df
cos.sim <- function(ix){
  A = X[,ix[1]]
  B = X[,ix[2]]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}   
n <- ncol(df) 
cmb <- expand.grid(i=1:n, j=1:n) 
C <- matrix(apply(cmb,1,cos.sim),n,n)
dim(C) #cos.similarity matrix 72 tc x 72 tc
rownames(C)=colnames(C)=colnames(df)

pdf('nmf_cosine.similarity.heatmap.pdf',height = 14,width = 12,useDingbats = T)
print(pheatmap::pheatmap(C, show_colnames = FALSE))
dev.off()

