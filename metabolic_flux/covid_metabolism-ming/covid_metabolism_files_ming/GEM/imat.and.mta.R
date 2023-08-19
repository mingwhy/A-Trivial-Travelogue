
(files=Sys.glob('../my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R` then source files
#install.packages("httpgd")
lapply(files, source)

library(gembox)
data("recon1")
data("media")
nc <- 4L # number of CPUs/cores to use for parallelization
model <- recon1

dat <- readRDS("data.for.gem.RDS")
names(dat)
# "exprs.int" "dflux.int"
names(dat$exprs.int)
#[1] "Vero"           "293T"           "Swab.Lieberman" "BALF"   

# run iMAT and rMTA on each dataset
# here the dataset id is passed via command line argument, such that the "imat.and.mta.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs

#Rscript imat.and.mta.R Vero 0 0.05
#x <- commandArgs(TRUE)
id = 'Vero'; #id <- x[1] # dataset ID
median = 0; #media <- x[2] # media
bm.lb = 0.05; #bm.lb <- as.numeric(x[3]) # minimal biomass requirement

exprs.int <- dat$exprs.int[[id]]
dflux.int <- dat$dflux.int[[id]]
str(exprs.int) #1905,1905
str(dflux.int)

if (media=="1") model <- set.medium(model, media$dmem.lo.glc, set.all=TRUE)
if (bm.lb!=0) {
  bm <- get.biomass.idx(model)
  model <- set.rxn.bounds(model, bm, lbs=bm.lb, relative=TRUE)
}
bm
str(model); #a huge GEM object

# run iMAT
# sample.par indicate return sampled flux solutions satisfying the obj func (see iMat.R in gembox)
start=Sys.time()
Rcpp::sourceCpp('achr_ming.cpp') #reload arch() which has been fixed
ctrl <- imat(model, exprs.int$ctrl, samp.pars=list(nc=nc, n.sample=2)) # control samples
#trt <- imat(model, exprs.int$trt, samp.pars=list(nc=nc, n.sample=5)) # virus-infected samples
#ctrl <- imat(model, exprs.int$ctrl, samp.pars=list(nc=nc, n.sample=5000)) # control samples
#trt <- imat(model, exprs.int$trt, samp.pars=list(nc=nc, n.sample=5000)) # virus-infected samples
imat.res <- list(ctrl=ctrl, trt=trt)
end=Sys.time()
end-start;
saveRDS(imat.res, file=paste0("imat.res.",id,".RDS"))

if(F){
# run rMTA
vtrt <- rowMeans(imat.res$trt$result.model$sample$pnts[, 2001:5000])
mta.res <- rmetal(model, vtrt, dflux.int, nc=nc, detail=FALSE)
saveRDS(mta.res, file=paste0("mta.res.",id,".RDS"))
}