#exprs: a vector with gene names
gembox_fluxVec<-function(exprs=exprs,model=model){
  #exprs=dat[,i]
  #names(exprs)=genes
  exprs=exprs[exprs!=0]
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

  #unique(imat.model$fluxes.int.imat) #9 0 1 -1
  #str(imat.model$solver.out)
  #imat.model$solver.out[[1]]$obj
  #head(imat.model$solver.out[[1]]$xopt)
  
########### 02.1_get_opt.flux_mat.R
  # imat.R in gembox
  # the mode 1 of imat (the fva-like approach): solve the iMAT MILP under the forced inactivaton/activation of each rxn, determine (de)activated rxns based on the resulting objective values
  # return a list(solver.out, flux.int.imat), solver.out is a data.table of the optimal iMAT objectives for all rxns, flux.int.imat is a vector in the order of the model rxns, 
  # with values 0/9/1/-1 representing a rxn being inactive, activity level not enforced, active in the forward direction, and active in the backward direction as determined by iMAT
  
  #length(x$imat.model$fluxes.int) # 11898, 1  0 -1
  #length(x$imat.model$fluxes.int.imat) # 11898, 9  1  0 -1
  #length(x$imat.model$solver.out[[1]]$xopt) #14495
  flux.vec=imat.model$solver.out[[1]]$xopt[1:length(imat.model$fluxes.int.imat)]
  #summary(flux.vec[x$imat.model$fluxes.int.imat==0]) # (-0.1,0.1)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==-1]) #(-1000,-1)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==1]) #(1,1000)
  #summary(flux.vec[x$imat.model$fluxes.int.imat==9]) #(-1000,1000)
  #df.flux.vec=cbind(df.flux.vec,flux.vec)

  ###############
  # update the original metabolic model based on iMAT result, model$lb and $ub are updated
  # based on imat.model$fluxes.int.imat, which was generated in run.imat
  #res.model <- update.model.imat(model, imat.model, imat.pars) #update.model.imat from imat.R
  #unique(res.model$lb)
  #unique(res.model$ub)

  return(flux.vec)
}
