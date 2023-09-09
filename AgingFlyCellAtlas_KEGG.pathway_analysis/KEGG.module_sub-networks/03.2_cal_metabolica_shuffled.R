
library(ggplot2);library(gridExtra)
library(tidyverse);library(Matrix)
#https://www.preprints.org/manuscript/202101.0280/v1
source('propagationFuncs.R')

load("./files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)
load("./files/rxn2gene.RData")
dim(rxn2gene) #df,6420 obs. of  4 variables, reaction and its associated enzymes
head(rxn2gene) #3rd column: Entrez 
length(unique(rxn2gene$Reaction)) #1901

# end product of each (sub)module
mz.produced=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$isProduced
})
table(unlist(lapply(mz.produced,length))) #1559,each sub-graph only produce one met
mz.produced=unlist(mz.produced)
tail(sort(table(mz.produced))) #there are metabolites which are end product of multiple module
length(unique(mz.produced)) #1270
which(mz.produced=='cpd:C00388')
all_producers[[718]] #sub-network which produce Histamine

# read in fly-human orthologs (00_process_h5ad_AFCA.R)
orthologs=readRDS('~/Documents/Data_AgingFlyCellAtlas/human_fly_ortholog_entrezid.rds')
colnames(orthologs)
sum(orthologs$human_entrezid=='')
sum(is.na(orthologs$human_entrezid)) #1558
orthologs=orthologs[!is.na(orthologs$human_entrezid),]
dim(orthologs) # 24085    10
gene.order=unique(orthologs$human_entrezid)
length(gene.order)#11310

##########################
#overlap with metabolome measured data
load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86
head(mz.age.betas)
tmp=paste0('cpd:', mz.age.betas$KEGGid)
sum(tmp %in% mz.produced) #54 overlapped
mz.age.betas[tmp %in% mz.produced,]

keep.mz.index=which(mz.produced %in% tmp)
length(keep.mz.index) #79
sort(table(mz.produced[keep.mz.index])) #there are metabolites which are product of multiple modules

########################################################
dataset.name='fly.female'
(files=Sys.glob('shuffle_out/*.rds'))
if(length(grep('endNode',files))!=0){
  files=files[-grep('endNode',files)]
}

for(file in files){
  #dat=readRDS('log1p_female_gene.mean.expr.rds')
  dat=readRDS(file)
  out.file=gsub('.rds','_endNode.res.rds',file)
  if(file.exists(out.file)){next}

  length(dat) #27
  length(dat[[1]]) #4
  names(dat[[1]])
  # tissue;cell.type;age, 4 age groups
  head(dat[[1]][[1]])
  #replace fly gene symbols to human gene ENTREZID

  sapply(dat[[2]],length)
  head(names(dat[[2]][[1]]))
  expr.mat=matrix(0,nrow=length(orthologs$human_entrezid),ncol=27*4)
  rownames(expr.mat)=orthologs$human_entrezid
  col.names.vec=c();
  z=0;
  for(i in 1:length(dat)){
    for(j in 1:length(dat[[i]])){
      z=z+1;
      expr.vec=dat[[i]][[j]]
      expr.vec=expr.vec[names(expr.vec) %in% orthologs$original.id]
      #length(expr.vec)
      tmp=orthologs[match(names(expr.vec),orthologs$original.id),]
      #sum(tmp$original.id==names(expr.vec))
      df.tmp=data.frame(gene.id=tmp$human_entrezid,expr=expr.vec[tmp$original.id])
      #sum(table(df.tmp$gene.id)>1)
      df.tmp2=df.tmp %>% group_by(gene.id) %>% summarize(expr.value=mean(expr))
      #length(unique(df.tmp2$gene.id))==nrow(df.tmp2)
      expr.mat[df.tmp2$gene.id,z]=df.tmp2$expr.value
      #colnames(expr.mat)[z]=names(dat[[i]])[[j]]
      col.names.vec=c(col.names.vec,names(dat[[i]])[[j]])
    }
  }
  dim(expr.mat)
  colnames(expr.mat)=col.names.vec
  expr.mat[1:3,1:3]

  ##########################
  ## in `propagationFuncs.R`
  #only look at metabolites which also have metabolome data
  #i=718;
  #c=lapply(1:length(all_producers),function(i){
  endNode.res=lapply(keep.mz.index,function(i){
    cat('start',names(all_producers)[i],'\n');
    my_producer <- all_producers[[i]] # i can 1 to length(all_producers)
    nodeVals <- do.calc.rxnvals(all_producers,rxn2gene,expr.mat)
    #dim(nodeVals); #3056 108
    #length(unique(unlist(rownames(nodeVals)))) #3056, #node=(rxn+cpd)
    #nodeVals[1:10,1:2]
    #table(Matrix::rowSums(nodeVals))
    trial<-tryCatch(
      res <- metabolica(nodes.vals=nodeVals, subgraph = my_producer$produceRigraph, 
                        ininodes = my_producer$initialNodes, 
                        endnode = my_producer$isProduced, algR = "minprod", algM = "nsqrt", 
                        maxnum = my_producer$maxnum*1.5),
      error=function(c) 'error')
    #names(res) #"node.signal" "sig.dif"  
    #length(res$node.signal) # for this metabolite, node.signal across 108 samples
    #head(res$node.signal)
    #length(res$sig.dif) #113: 10.00000 29.94719 ... 
    #tail(res$sig.dif)
    # difference bewteen iterations when calculating this node signal
    if(trial=='error'){return(NA)}
    if(is.na(res$sig.dif[length(res$sig.dif)])){return(NA)}
    res$node.signal
  })
  names(endNode.res)<-names(all_producers)[keep.mz.index]
  #saveRDS(endNode.res,'endNode.res')
  saveRDS(endNode.res,out.file);
}

