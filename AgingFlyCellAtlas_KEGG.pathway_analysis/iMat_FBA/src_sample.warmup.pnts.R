
sample.warmup.pnts <- function(model, n, nc) {
  # sample warmup points for ACHR
  # n: the number of warmup points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  if (n<2*n.rxns) {
    n <- 2*n.rxns
    warning(sprintf("#{warmup points} should be at least 2*#{reactions}=%d.", 2*n.rxns), immediate.=TRUE)
  }
  message("Begin generating ", n, " warmup points...")
  orth.pnts <- get.orth.pnts(model, n, nc)
  rand.pnts <- get.rand.pnts(model, n, nc)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  res <- orth.pnts*r + rand.pnts*(1-r)
  message("Done generating warmup points.")
  res
}

get.orth.pnts <- function(model, n, nc) {
  # sample orthogonal points for ACHR
  # n: the number of orthogonal points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  mat <- cbind(Diagonal(n.rxns), Diagonal(n.rxns, x=-1))
  if (n<=2*n.rxns) {
    mat <- mat[, sample(2*n.rxns, n)]
  } else {
    mat <- cbind(mat[, sample(2*n.rxns)], mat[, sample(2*n.rxns, n-2*n.rxns, replace=TRUE)])
  }
  message("1. Generating orthogonal points, progress:")
  #get.opt.pnts(model, mat, nc)
  out1=pbmcapply::pbmclapply(1:(ncol(mat)/2), function(i)   
    solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"))[[1]]$xopt, 
    mc.cores=nc) 
  out2=pbmcapply::pbmclapply((1+(ncol(mat)/2)):ncol(mat), function(i) 
    solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"))[[1]]$xopt, 
    mc.cores=nc) 
  #do.call(cbind, pbmcapply::pbmclapply(1:ncol(mat), function(i) 
  #  solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"), pars=.pkg.const$lp[[.pkg.var$solver]])[[1]]$xopt, 
  #  mc.cores=nc)    )
  out=c(out1,out2)
  return(out)

}

get.rand.pnts <- function(model, n, nc) {
  # sample random points for ACHR
  # n: the number of random points
  # nc: number of cores to use

  n.rxns <- ncol(model$S)
  #cs <- runif(n.rxns*n) - 0.5
  #dim(cs) <- c(n.rxns, n)
  mat <- runif(n.rxns*n) - 0.5
  dim(mat) <- c(n.rxns, n)

  message("2. Generating random points, progress:")
  #get.opt.pnts(model, cs, nc)
  out1=pbmcapply::pbmclapply(1:(ncol(mat)/2), function(i)   
    solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"))[[1]]$xopt, 
    mc.cores=nc) 
  out2=pbmcapply::pbmclapply((1+(ncol(mat)/2)):ncol(mat), function(i) 
    solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"))[[1]]$xopt, 
    mc.cores=nc) 
  out=c(out1,out2)
  return(out)
}

get.opt.pnts <- function(model, mat, nc) {
  # a helper function to get optimal points corresponding to running LPs with objective function coefficients being the columns of mat, for ACHR

  do.call(cbind, pbmcapply::pbmclapply(1:ncol(mat), function(i) solve.model(model, csense="max", c=mat[,i]/norm(mat[,i],"2"), pars=.pkg.const$lp[[.pkg.var$solver]])[[1]]$xopt, mc.cores=nc))
}