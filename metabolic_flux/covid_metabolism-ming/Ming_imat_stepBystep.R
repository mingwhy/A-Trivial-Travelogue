# github: https://github.com/ruppinlab/covid_metabolism
# download `covid_metabolism_files.tar.gz` from google drive link:
# https://drive.google.com/file/d/1bVPCQlDR3G8TTMx09jkN3IuGwzq9n9hi/view

(files=Sys.glob('../gembox-ming/my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R` then source files
#install.packages("httpgd")
lapply(files, source)

library(gembox)
data("recon1")
 # checked in unique(model$rowlb)
  # check code pasted:
  #data("recon1")
  #model<-recon1
  #length(model$rowub) #equal to #mz in S
  #length(model$rowlb)
  #dim(model$S)
  #unique(model$rowub)
  #unique(model$rowlb) #all 0
  
data("media")
nc <- 2L # number of CPUs/cores to use for parallelization
model <- recon1

dat <- readRDS("covid_metabolism_files_ming/GEM/data.for.gem.RDS")
names(dat)
# "exprs.int" "dflux.int"
names(dat$exprs.int)
#[1] "Vero"           "293T"           "Swab.Lieberman" "BALF"   

# run iMAT and rMTA on each dataset
# here the dataset id is passed via command line argument, such that the "imat.and.mta.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs

#Rscript imat.and.mta.R Vero 0 0.05
#Rscript imat.and.mta.R Swab.Lieberman 0 0
#x <- commandArgs(TRUE)
id = 'Swab.Lieberman'; #id <- x[1] # dataset ID
media = 0; #media <- x[2] # media
bm.lb = 0; #bm.lb <- as.numeric(x[3]) # minimal biomass requirement

exprs.int <- dat$exprs.int[[id]]
dflux.int <- dat$dflux.int[[id]]
str(exprs.int) #1905,1905
str(dflux.int)

if (media=="1") model <- set.medium(model, media$dmem.lo.glc, set.all=TRUE)
if (bm.lb!=0) {
  bm <- get.biomass.idx(model)
  model <- set.rxn.bounds(model, bm, lbs=bm.lb, relative=TRUE)
}
#bm
str(model); #a huge GEM object

# run iMAT (imat.R from gembox)
# sample.par indicate return sampled flux solutions satisfying the obj func (see iMat.R in gembox)
#ctrl <- imat(model, exprs.int$ctrl, samp.pars=list(nc=nc, n.sample=5)) # control samples
#imat <- function(model, expr, imat.pars=list(), solv.pars=list(), samp.pars=list()) {

####################################################
# step by step of imat
pick='ctrl';
#pick='trt';
if(pick=='ctrl'){expr=exprs.int$ctrl;imat.pars=list(); solv.pars=list();samp.pars=list(nc=nc, n.sample=100)
}else{expr=exprs.int$trt;imat.pars=list(); solv.pars=list();samp.pars=list(nc=nc, n.sample=100)}

# formulate iMAT model
imat.model <- form.imat(model, expr, imat.pars)
length(model$lb); #3788
length(imat.model$lb); #5120
dim(model$S) #2766 3788 
dim(imat.model$S) #4197 5120
unique(imat.model$fluxes.int) #1  0 NA -1

# solve the iMAT model, two new things are added to imat.model
# imat.model$fluxes.int.imat;
# imat.model$solver.out;
imat.model <- run.imat(imat.model, imat.pars, solv.pars) #run.imat mode 0, only determine (de)activated reactions
unique(imat.model$fluxes.int.imat) #9  0  1 -1
str(imat.model$solver.out)
imat.model$solver.out[[1]]$obj
head(imat.model$solver.out[[1]]$xopt)

# update the original metabolic model based on iMAT result
res.model <- update.model.imat(model, imat.model, imat.pars)
unique(res.model$lb)
unique(res.model$ub)

###### begin sample.model
#error message from:
#res.model <- sample.model(res.model, samp.pars)
# sample.model() from sampling.R in gembox
#sample.model <- function(model, pars=list()) {
model=res.model;pars=samp.pars;
pars <- get.pars("samp", pars)
pars$method #"achr" #achr is implemented in `achr.cpp`
names(names(model)) #NULL
pars$steps.per.pnt #400 default value

pnt.file=paste0(id,'_',pick,'_warmup.pnts.rds')
if(!file.exists(pnt.file)){
  warmup.pnts <- sample.warmup.pnts(model, pars$n.warmup, pars$nc) #10min
  #saveRDS(warmup.pnts,'Swab.Lieberman_ctrl_warmup.pnts.rds')
  saveRDS(warmup.pnts,pnt.file)
}
#warmup.pnts=readRDS('Vero_warmup.pnts.rds')
#warmup.pnts=readRDS('Swab.Lieberman_trt_warmup.pnts.rds')
warmup.pnts=readRDS(pnt.file)

center.pnt <- rowMeans(warmup.pnts)
init.stat <- list(center.pnt=center.pnt, prev.pnt=center.pnt, n.tot.steps=0)
message("Will sample ", pars$n.sample, " points.")

dim(model$rxnGeneMat) #3788 x 1905
dim(model$S) #2766 x 3788
dim(warmup.pnts) #3788 7576 (3788x2)

#Rcpp::sourceCpp('~/Downloads/gembox-ming/src/achr.cpp')
#res <- achr(model, init.stat, warmup.pnts, pars$n.sample, pars$steps.per.pnt)
#res <- achr(model, init.stat, warmup.pnts, 1, 2)

#In achr.cpp: 
#res <- achr(model, model$sample$stat, model$sample$warmup.pnts, pars$n.sample, pars$steps.per.pnt)
length(model$lb) #3788
length(model$ub) #3788

achr.file=paste0("achr.res.",id,'.',pick,".RDS");
if(!file.exists(achr.file)){
  #Rcpp::sourceCpp('achr_test.cpp')
  Rcpp::sourceCpp('achr_ming.cpp')
  #achr(model, state, warmupPnts, nPnts, stepsPerPnt) 
  #achr(model, init.stat, warmup.pnts, 1, 2)
  
  start=Sys.time()
  res=achr(model, init.stat, warmup.pnts, 1e4, 400)
  end=Sys.time() #50min ctrl, 25min trt
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
######## finish sample.model
saveRDS(list(imat.model=imat.model, result.model=res.model),paste0(pick,'.rds'))

###save imat.res.xxx.RDS, used in `dflux.R`
ctrl=readRDS('ctrl.rds')
trt=readRDS('trt.rds')
imat.res <- list(ctrl=ctrl, trt=trt)
names(imat.res$ctrl$result.model) #make sure `sample` exists
saveRDS(imat.res, file=paste0("imat.res.",id,".RDS"))

####################################################
library(data.table)
data("recon1")
nc <- 2L # number of CPUs/cores to use for parallelization
model <- recon1

# run rMTA (mta.R from gembox)
start=Sys.time()
dim(imat.res$trt$result.model$sample$pnts)
vtrt <- rowMeans(imat.res$trt$result.model$sample$pnts[, 4001:1e4])
mta.res <- rmta(model, vtrt, dflux.int, nc=nc, detail=FALSE)
end=Sys.time() 
end-start; #8min

saveRDS(mta.res, file=paste0("mta.res.",id,".RDS"))

if(FALSE){
#### rmta source code below (in gembox 'mta.R')
#rmta <- function(model, flux0, dflux, rxns="all+ctrl", ko=NULL, nc=1L, detail=TRUE, k=100, mta.pars=list(), mip.pars=list(), qp.pars=list()) {
flux0=vtrt; dflux=dflux.int; 
rxns="all+ctrl";ko=NULL; detail=FALSE; k=100;
mta.pars=list();mip.pars=list(); qp.pars=list();
  
# formulate MTA model for either direction (dflux and -dflux)
mta.model <- form.mta(model, flux0, dflux, mta.pars)
mta.model0 <- form.mta(model, flux0, -dflux, mta.pars)

if (!is.null(ko)) {
  ko <- all2idx(model, ko)
  mta.model$lb[ko] <- 0
  mta.model$ub[ko] <- 0
  mta.model0$lb[ko] <- 0
  mta.model0$ub[ko] <- 0
  model$lb[ko] <- 0
  model$ub[ko] <- 0
}
start=Sys.time()
# solve the MTA models across the rxns for both models
message("rmta(): Running MTA for dflux.")
res1 <- run.ko.screen(mta.model, rxns, run.mta, solv.pars=mip.pars, detail=detail, nc=nc)
message("rmta(): Running MTA for -dflux.")
res0 <- run.ko.screen(mta.model0, rxns, run.mta, solv.pars=mip.pars, detail=FALSE, nc=nc)
end=Sys.time();
end-start; #32s

# MOMA for dflux
start=Sys.time()
message("rmta(): Running MOMA.")
tmpf <- function(x) get.mta.score(model=mta.model, x, detail=detail)
res.moma <- moma(model, rxns, nc, flux0, obj=tmpf, solv.pars=qp.pars)
end=Sys.time();
end-start;#9min

# rMTA score
res <- data.table(id=res1$id, rxn=res1$rxn, bTS=res1$score.mta, wTS=res0$score.mta, mTS=res.moma$score.mta)
res[, rTS:=ifelse(bTS>0 & mTS>0 & wTS<0, k*mTS*(bTS-wTS), mTS)]
mta.res=list(mta.model=mta.model, result.mta=res1, result.moma=res.moma, result.rmta=res)
# := symbol in R
# https://stackoverflow.com/questions/32817780/what-is-the-r-assignment-operator-for
}

