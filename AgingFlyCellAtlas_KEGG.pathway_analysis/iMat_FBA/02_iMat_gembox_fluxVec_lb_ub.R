##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
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

ifly$genes=unlist(ifly$genes)
length(ifly$genes)
length(ifly$rules)
ifly$ub=ifly$ub[,1]
ifly$lb=ifly$lb[,1]
ifly$rowub=rep(0,nrow(ifly$S)) #used in iMat, check with human data `data("recon1")`, all 0
ifly$rowlb=rep(0,nrow(ifly$S))

##
#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns)) #11898
ifly$rxns=unlist(ifly$rxns)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)

#some quick access test
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
mets.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/metabolites.tsv')
head(mets.mapping)
mets.mapping[mets.mapping$metKEGGID=='C00388',]#Cytosol,Extracellular #ifly$compNames
tmp=ifly$S[ifly$mets %in% mets.mapping[mets.mapping$metKEGGID=='C00388',]$mets,]
col.index=which( colSums(abs(tmp))!=0 ) 
tmp[,col.index] #all rxn which involve `C00388`, not necessariliy histine metabolism
ifly$rxns[col.index] #"MAR04428" "MAR00619"
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124c
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124e
ifly$subSystems[col.index]

rxns.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/reactions.tsv')
rxns.mapping[rxns.mapping$rxns=="MAR01442",]

##############################################
## iMAT
#https://github.com/mingwhy/A-Trivial-Travelogue/tree/main/metabolic_flux/gembox-ming
(files=Sys.glob('./my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R` then source files
#install.packages("httpgd")
lapply(files, source)

nc <- 4L # number of CPUs/cores to use for parallelization
library(gembox);
source('src_exprs2fluxes.R')
model=ifly;

dataset.name='fly.femlae'
#dat=data.table::fread('log1p_female_gene.mean.expr.csv')
dat=readRDS('log1p_female_gene.mean.expr.rds')
length(dat)
length(dat[[1]])
sum(model$genes %in% names(dat[[1]][[1]])) #1480


dir.create('imat_out');

# run iMAT 

for(tc in names(dat)){
  #tc=names(dat)[13] #'head;adult brain perineurial glial cell'
  #tc=names(dat)[20] #"head;gamma Kenyon cell"
  names(dat[[tc]]) #4 age groups
  for(id in names(dat[[tc]]) ){
    exprs=dat[[tc]][[id]]
    id=gsub('\\/','\\-',id) #save filename can not contain '/'
    #id; tissue; cell.type; age, .rds
    out.file=paste0('imat_out/id;',id,'.rds')
    if(file.exists(out.file)){cat('file',out.file,'already exists\n');next}
    cat('start id',out.file,'\n')
  
    exprs=exprs[exprs!=0] #remove 0 (NA is encoded as 0 in `01_cal_gene.per_per.cell.type.R`)
    # in `exprs2fluxes` should be named by gene symbols as used in model$genes, or if it's unnamed and length being length(model$genes), assume it's already in the same order as model$genes
    # gene order handled in `exprs2int` in `iMat.R`
    exprs.int<- exprs2int(model, exprs)
    #table(genes.int)
    #length(genes.int)==length(model$genes)
    
    ## step by step of imat (code copied from gembox/* functions)
    expr=exprs.int; imat.pars=list(); solv.pars=list();
    #################
    # formulate iMAT model,run format.imat mannually
    #imat.model <- form.imat(model, expr, imat.pars)
    
    # formulate iMAT model (the original version)
    # expr is the output from exprs2int()
    # imat.pars: the parameters for iMAT
    pars <- get.pars("imat", imat.pars) #function imat.pars from `global.R`
    
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
    # imat.model$fluxes.int.imat: a vector in the order of the model rxns, with values 0/9/1/-1 representing 
    # a rxn being inactive(0), activity level not enforced(9), active in the forward direction(1), and active in the backward direction as determined by iMAT
    # imat.model$solver.out: if imat.pars$mode==0, then this is the output of MILP solver; 
    imat.model <- run.imat(imat.model, imat.pars, solv.pars) #run.imat mode 0, only determine (de)activated reactions
  
    unique(imat.model$fluxes.int.imat) #9 0 1 -1
    #str(imat.model$solver.out)
    #imat.model$solver.out[[1]]$obj
    #head(imat.model$solver.out[[1]]$xopt)
    summary(imat.model$solver.out[[1]]$xopt)
    
    ###############
    # update the original metabolic model based on iMAT result, model$lb and $ub are updated
    # based on imat.model$fluxes.int.imat, which was generated in run.imat
    res.model <- update.model.imat(model, imat.model, imat.pars) #update.model.imat from imat.R
    unique(res.model$lb)
    unique(res.model$ub)
    
    ###############
    saveRDS(list(imat.model=imat.model, result.model=res.model),
      out.file)
    cat('id',id,'is done\n') #file name, 2hr
  }
}

