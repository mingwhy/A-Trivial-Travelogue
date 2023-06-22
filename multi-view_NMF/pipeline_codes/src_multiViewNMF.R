#multiViewNMF( X, K, lambda, 50 );
multiViewNMF<-function( X, K, lambda, maxIter ){
# Multi-View Non-negative symmetric Matrix Factorization
#
# INPUT:
#   X: a cell which contains symmetric matrices
#   K: the number of hidden factors
#   lambda: a vector which contains the parameters for balancing the relative
#                weight among different views
#   maxiter: the maximum number of iterations
#
# OUTPUT:
#   H: a cell which contains factor matrices for all views
#   Hc: the result consensus factor matrix
#   objValue: the value of objective function
  #X=readRDS('X_for_multiViewNMF.rds')
  #names(X)#"Strength"  "Participation"
  #K = 5;
  #lambda = c(0.01, 0.05);
  #[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 );

  Xcount=length(X)
  Ncount=nrow(X[[1]])
  n=Ncount

  H = list()
  Hc = matrix(0,nrow=Ncount, ncol=K);

  # Normalize input matrix X{i}
  for(i in 1:Xcount){
    X[[i]] = X[[i]] / (sum(apply(X[[i]] ^2,1,sum)))^0.5
    #X[[i]] = X[[i]]/sqrt( sum(diag( t(X[[i]]) %*% X[[i]])) );
  }


  # Initialize H and Hc
  cat('Initialize H and Hc\n')
  iterNum = 30;
  if(Ncount > 2000) iterNum = 20

  out = SNMF(X[[1]], K, H_init=NULL, maxIter=iterNum, epsilon=1e-5);
  H[[1]]=out$H
  for(i in 2:Xcount){
    out=SNMF(X[[i]], K, H[[i-1]], maxIter=iterNum, epsilon=1e-5);
    H[[i]] = out$H
  }
  for(i in 1:Xcount){
    Hc = Hc + lambda[i]*H[[i]];
  }
  Hc = Hc/sum(lambda);

  obj_old = 0;
  for(i in 1:Xcount){
      obj_body = norm(X[[i]] - H[[i]] %*% t(H[[i]]), type='F')^2;
      obj_consensus = norm(H[[i]] - Hc, type='F')^2;
      obj_old = obj_old + obj_body + lambda[i]*obj_consensus;
  }

  ##  Update process
  cat('Begin update process\n')
  
  objValue = obj_old;
  for(iter in 1:maxIter){
    cat('begin iteration',iter,'\n')

    # Fixing Hc, minimize objective function over H
    maxIterforView = 40;
    for(i in 1:Xcount){
      H[[i]] = SNMFforView(X[[i]], Hc, H[[i]], lambda[i], maxIterforView, epsilon=1e-6);
    }
    
    # Fixing H, minimize objective function over Hc
    Hc = matrix(0,nrow=Ncount,ncol=K) 
    for(i in 1:Xcount)
      Hc = Hc + lambda[i]*H[[i]];
    
    Hc = Hc/sum(lambda);
    
    # Object function value and the relative error
    obj = 0;
    for(i in 1:Xcount){
      obj_body = norm(X[[i]] - H[[i]] %*% t(H[[i]]), type='F')^2;
      obj_consensus = norm(H[[i]] - Hc, type='F')^2;
      obj = obj + obj_body + lambda[i]*obj_consensus;
    }
    
    Delta = obj_old - obj;
    obj_old = obj;
    objValue = c(objValue, obj);
    
    if(Delta < 1e-6) break
    
  }
  # keep H, Hc, objValue
  #sapply(H,dim)
  #dim(Hc)
  #length(objValue)
  return(list(H=H,Hc=Hc,objValue=objValue))
}
