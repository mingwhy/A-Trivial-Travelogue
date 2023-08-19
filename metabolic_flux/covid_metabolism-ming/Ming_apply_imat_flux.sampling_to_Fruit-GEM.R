
##############################################
## prepare
# install gurobi
# https://github.com/ruppinlab/gembox
# need to install gurobi (solver algorithm) (or Rcplex2: https://github.com/ruppinlab/Rcplex2)
# library(gurobi)
# check out install_gurobi.R
#devtools::install_github("ruppinlab/gembox")
library(gembox)

## algorithm of lee et al. BMC Syst Biol 2012 is implemented in https://github.com/ruppinlab/gembox/blob/master/R/exprs.R

#README.md from https://github.com/ruppinlab/covid_metabolism
# 1, download data
#$ cd covid_metabolism-main/data
#$pip3 install gdown #for use gdown in ./dnload.sh
#$ ./dnload.sh
# covid_metabolism_files.tar.gz (270MB), liao.RDS(3.6G), chua.RDS(1.7G) files.

# process data `collect.validation.data.R`
#BiocManager::install('GSA')    
#devtools::install_github("ImNotaGit/my.utils")
#library(my.utils)
(files=Sys.glob('../gembox-ming/my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R`, install 'httpgd', then source files.
#install.packages("httpgd")
lapply(files, source)

##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('./Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
# run GEM_convertOldStyleModel.m to get `rules` in GEM fields, used in exprs2fluxes.R
rules=readMat('rules.mat')
length(rules)

ifly=gem[[1]][,,1]
length(ifly)
lapply(ifly,dim)
names(ifly)
dim(ifly$S) # 8135mz X 11898 fluxes (how many overlapped with metabolome assay)
#genes=unlist(ifly$genes)
#length(genes) #1753 genes
ifly$rules=rules$rules
length(ifly$subSystems)
ifly$subSystems=unlist(ifly$subSystems)

##
#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns)) #11898
ifly$rxns=unlist(ifly$rxns)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)
length(table(ifly$subSystems)) #149
ifly$subSystems[grep('hist',ifly$subSystems)]
ifly$subSystems[grep('hist',ifly$subSystems,ignore.case = T)]

col.idx=grep('hist',ifly$subSystems,ignore.case = T);
tmp=ifly$S[,grep('hist',ifly$subSystems,ignore.case = T)]
row.idx=which(Matrix::rowSums(tmp)!=0)

ifly$rxns[col.idx]
ifly$rxnNames[col.idx]
ifly$metFormulas[col.idx]
ifly$grRules[col.idx]
ifly$metNames[row.idx]
tmp=ifly$S[row.idx,col.idx]
rowSums(tmp)
rowSums(ifly$S[row.idx,])
tmp=as.matrix(tmp)
rownames(tmp)=ifly$metNames[row.idx]
colnames(tmp)=ifly$rxnNames[col.idx]
tmp
#############################################
# id Cross references
#kegg id C01672: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02123e
#model folder in: https://github.com/SysBioChalmers/Fruitfly-GEM
mets.mapping=data.table::fread('./Fruitfly-GEM-main/model/metabolites.tsv')
head(mets.mapping)
mets.mapping[mets.mapping$metKEGGID=='C00388',]#Cytosol,Extracellular #ifly$compNames
tmp=ifly$S[ifly$mets %in% mets.mapping[mets.mapping$metKEGGID=='C00388',]$mets,]
col.index=which( colSums(abs(tmp))!=0 ) 
tmp[,col.index] #all rxn which involve `C00388`, not necessariliy histine metabolism
ifly$rxns[col.index] #"MAR04428" "MAR00619"
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124c
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124e
ifly$subSystems[col.index]

rxns.mapping=data.table::fread('./Fruitfly-GEM-main/model/reactions.tsv')
rxns.mapping[rxns.mapping$rxns=="MAR01442",]
##############################################
## iMAT
(files=Sys.glob('~/Documents/bioinfo_software/GEM_metabolic.models/covid_metabolism-ming/my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R` then source files
#install.packages("httpgd")
lapply(files, source)

nc <- 2L # number of CPUs/cores to use for parallelization
library(gembox)
model=ifly;
model$genes=unlist(model$genes)
length(model$genes)
length(model$rules)
model$ub=model$ub[,1]
model$lb=model$lb[,1]
model$rowub=rep(0,nrow(model$S)) #used in iMat, check with human data `data("recon1")`, all 0
model$rowlb=rep(0,nrow(model$S))

dataset.name='fly.body'
dat=data.table::fread('log1p_female_body.csv')
genes=dat[[1]]
sum(model$genes %in% genes) #1480

dat=as.data.frame(dat[,-1])
dim(dat) #14426 gene x 80 (cell.type_age)
dat=dat[,-grep('unanno',colnames(dat))]
dim(dat) #76
source('src_exprs2fluxes.R')
source('src_sample.warmup.pnts.R')
dir.create('warmup_out'); #save warmup files used in iMat
dir.create('achr_out')
dir.create('imat_out');
# run iMAT 
#for(i in 1:ncol(dat)){
i=1
  id=colnames(dat)[i]
  exprs=dat[,i]
  names(exprs)=genes
  exprs=exprs[exprs!=0]
  # in `exprs2fluxes` should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
  # gene order handled in `exprs2int` in `iMat.R`
  exprs.int<- exprs2int(model, exprs)
  #table(genes.int)
  #length(genes.int)==length(model$genes)
  
  ## step by step of imat
  expr=exprs.int; imat.pars=list(); solv.pars=list();
  #################
  # formulate iMAT model,run format.imat mannually
  #imat.model <- form.imat(model, expr, imat.pars)
  
  # formulate iMAT model (the original version)
  # expr is the output from exprs2int()
  # imat.pars: the parameters for iMAT
  pars <- get.pars("imat", imat.pars)
  
  # get intended activity level of rxns (-1/0/1) from discretized expression data
  rxns.int <- exprs2fluxes(model, expr)
  rxns.int[model$lb==0 & model$ub==0] <- 0
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)
  S <- model$S
  
  # 1. Active reactions: specify the y+ indicator variables, representing activation in the forward direction (i.e. v>flux.act)
  rxns.act <- which(rxns.int>0)
  n.act <- length(rxns.act)
  if (n.act!=0) {
    m1 <- sparseMatrix(1:n.act, rxns.act, dims=c(n.act, n.rxns))
    m2 <- Diagonal(n.act, x=(-pars$flux.act-pars$flux.bound))
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.act))), cbind(m1, m2))
  }
  
  # 2. Reversible active reactions: for those reversible ones among the active reactions, specify the extra y- indicator variables, representing activation in the backward direction (i.e. v<-flux.act)
  # thus, an reversible active reaction has both the y+ and y- indicator variables, because it can be active in either direction (but never both, i.e. 1 XOR 2)
  rxns.act.rev <- which(rxns.int>0 & model$lb<0)
  n.act.rev <- length(rxns.act.rev)
  if (n.act.rev!=0) {
    m1 <- sparseMatrix(1:n.act.rev, rxns.act.rev, dims=c(n.act.rev, ncol(S)))
    m2 <- Diagonal(n.act.rev, x=pars$flux.act+pars$flux.bound)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.act.rev))), cbind(m1, m2))
  }
  
  # 3. Inactive reactions: specify the y0 indicator variables
  # 3a. specify inactivation in the forward direction (i.e. v<flux.inact)
  rxns.inact <- which(rxns.int<0)
  n.inact <- length(rxns.inact)
  if (n.inact!=0) {
    m1 <- sparseMatrix(1:n.inact, rxns.inact, dims=c(n.inact, ncol(S)))
    m2 <- Diagonal(n.inact, x=pars$flux.bound-pars$flux.inact)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.inact))), cbind(m1, m2))
  }
  # 3b. for those reversible inactive reactions, need to further specify inactivation in the backward direction (i.e. v>-flux.inact)
  # note that a reversible inactive reaction has only one y0 indicator variable, because for these reactions we want -flux.inact<v<flux.inact (3a AND 3b)
  rxns.inact.rev <- which(rxns.int<0 & model$lb<0)
  n.inact.rev <- length(rxns.inact.rev)
  if (n.inact.rev!=0) {
    m3 <- sparseMatrix(1:n.inact.rev, rxns.inact.rev, dims=c(n.inact.rev, ncol(S)-n.inact))
    m4 <- sparseMatrix(1:n.inact.rev, match(rxns.inact.rev, rxns.inact), x=pars$flux.inact-pars$flux.bound, dims=c(n.inact.rev, n.inact))
    S <- rbind(S, cbind(m3, m4))
  }
  
  # other parameters
  n <- n.act + n.act.rev + n.inact + n.inact.rev
  rowlb <- c(model$rowlb, rep(-pars$flux.bound, n))
  rowub <- c(model$rowub, rep(pars$flux.bound, n))

  n <- ncol(S) - n.rxns
  c <- rep(c(0, 1/sum(rxns.int!=0,na.rm=TRUE)), c(n.rxns, n))
  vtype <- ifelse(c==0, "C", "I")
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  var.ind <- rep(c("v","y+","y-","y0"), c(n.rxns, n.act, n.act.rev, n.inact)) # iMAT variable type indicators (v: fluxex; y+/-/0: indicator variables)
  
  # iMAT model
  imat.model <- list(exprs.int=expr, fluxes.int=rxns.int,
                     rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                     rxns=model$rxns, mets=model$mets, csense="max", c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
  if ("irxn.ids" %in% names(model)) imat.model$irxn.ids <- model$irxn.ids # for multicellular imat models, keep the irxns.ids variable (intracellular reactions indices)
  #imat.model
  unique(imat.model$fluxes.int) #0 -1  1
  
  ###############
  # solve the iMAT model, two new things are added to imat.model
  imat.model <- run.imat(imat.model, imat.pars, solv.pars) #run.imat mode 0, only determine (de)activated reactions

  unique(imat.model$fluxes.int.imat) #9 0 1 -1
  str(imat.model$solver.out)
  imat.model$solver.out[[1]]$obj
  summary(imat.model$solver.out[[1]]$xopt)
  length(imat.model$solver.out[[1]]$xopt) #longer than original S due to reverse reaction addition
  dim(model$S) # 8135 11898
  length(imat.model$fluxes.int.imat) #11898
  table(imat.model$vtype)
  #C     I 
  #11898  2821 
  x=which(imat.model$vtype=='C') #1~11898
  summary(imat.model$solver.out[[1]]$xopt[x])
  
  ###############
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, imat.pars)
  unique(res.model$lb)
  unique(res.model$ub)
  # based on this updated model, could be further fed into matlab or python to sample fluxes.
  writeMat("test.mat", res.model=res.model) #took too long
  
  ################ 
  #begin sample.model
  #res.model <- sample.model(res.model, samp.pars)
  # sample.model() from sampling.R in gembox
  #sample.model <- function(model, pars=list()) {
  samp.pars=list(nc=6, n.sample=2)
  pars=samp.pars;
  pars <- get.pars("samp", pars)
  pars$nc
  pars$method #"achr" #achr is implemented in `achr.cpp`
  pars$steps.per.pnt #400 default value
  
  model=res.model; 
  names(names(model)) #NULL
  
  # get warmup data mannually
  pnt.file=paste0('./warmup_out/id_',id,'_warmup.pnts.rds')
  if(!file.exists(pnt.file)){
    #warmup.pnts <- sample.warmup.pnts(model, pars$n.warmup, pars$nc)
    n=pars$n.warmup; nc=pars$nc;
    n.rxns <- ncol(model$S)
    if (n<2*n.rxns) {
      n <- 2*n.rxns
      warning(sprintf("#{warmup points} should be at least 2*#{reactions}=%d.", 2*n.rxns), immediate.=TRUE)
    }
    message("Begin generating ", n, " warmup points...")
    orth.pnts.list <- get.orth.pnts(model, n, nc)
    orth.pnts=do.call(cbind,orth.pnts.list)
    dim(orth.pnts) #11898 23796
    
    rand.pnts.list <- get.rand.pnts(model, n, nc)
    rand.pnts=do.call(cbind,rand.pnts.list)
    dim(rand.pnts) #11898 23796
    
    r <- rep(runif(n), each=n.rxns)
    dim(r) <- c(n.rxns, n)
    res <- orth.pnts*r + rand.pnts*(1-r)
    dim(res) # 11898 23796
    message("Done generating warmup points.")
    warmup.pnts=res
    saveRDS(warmup.pnts,pnt.file)
  }
  warmup.pnts=readRDS(pnt.file)

  ##
  center.pnt <- rowMeans(warmup.pnts)
  init.stat <- list(center.pnt=center.pnt, prev.pnt=center.pnt, n.tot.steps=0)
  message("Will sample ", pars$n.sample, " points.")

  dim(model$S) #2766 x 3788
  dim(warmup.pnts) #3788 7576 (3788x2)

  
  achr.file=paste0('achr_out/achr_id_',id,".RDS");
  if(!file.exists(achr.file)){
    #Rcpp::sourceCpp('achr_test.cpp')
    Rcpp::sourceCpp('achr_ming.cpp')
    #achr(model, init.stat, warmup.pnts, 1, 2)
    #achr(model, state, warmupPnts, nPnts, stepsPerPnt) 
    
    #took too long
    start=Sys.time()
    res=achr(model, init.stat, warmup.pnts, 2, 400)
    end=Sys.time() # 2reps=
    end-start;
    
    #saveRDS(res, file=paste0("achr.res.",id,".ctrl.RDS"))
    saveRDS(res, file=!file.exists(achr.file))
  }
  res=readRDS(achr.file)
  
  model$sample <- list()
  model$sample$warmup.pnts <- warmup.pnts
  model$sample$pnts <- res$pnts
  model$sample$stat <- res$stat
  res.model<-model
  
  saveRDS(list(imat.model=imat.model, result.model=res.model),
    paste0('imat_out/id_,',id,'.rds'))
#}

