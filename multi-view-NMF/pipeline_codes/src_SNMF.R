#SNMF H{1} = SNMF(X{1}, K, H_init, iterNum, 1e-5);
#K=5;H_init=NULL;iterNum=30;epsilon=1e-05

SNMF<-function(X, K, H_init=NULL, maxIter, epsilon){
  # Symmetric Non-Negtive Matrix Factorization
  
  # INPUT:
  #       X: the adjacency matrix of a network
  #       K: the number of hidden factors
  #       maxIter: the maximal number of iterations for alternating minimization
  #       epsilon: the convergence parameter
  #
  # OUTPUT:
  #       H: the factor matirx
  #       objj: the value of objective function
  #
  
  ## Normalize the network
  #X = X/sqrt(sum(diag(t(X) %*% X)));
  X= X / (sum(apply(X^2,1,sum)))^0.5
  
  
  ## Initializaiton
  N = nrow(X)
  if (is.null(H_init)){
    H = matrix(runif(N*K),nrow=N,ncol=K)
  }else{
    H = H_init;
  }
  #H = H/sqrt(sum(diag(t(H) %*% H))) ;
  H = H / (sum(apply(H^2,1,sum)))^0.5

  
  #test = [2 0 1;-1 1 0;-3 3 0]; #in matlab
  #n = norm(test,'fro')^2 #in matlab
  #norm(rbind(c(2,0,1),c(-1,1,0),c(-3,3,0)),type='F')^2
  obj_old = norm(X - H %*% t(H), type='F')^2;
  beta = 0.5;
  objj = obj_old;
  
  ## Alternating update 
  #for matrix, * in matlab, %*% in R.
  # .* in matlab, * in R (dot product)
  #eps=pracma::eps(1); #floating-point relative accuracy
  eps=2.220446e-16;
  for(iter in 1:maxIter){  
      temp_1 = X %*% H;
      temp_2 = t(H) %*% H;
      temp_2 = H %*% temp_2;
      H = H*(1 - beta + beta*(temp_1/(temp_2+eps)));
  
      obj = norm(X - H %*% t(H), type='F')^2;
      Delta = obj_old - obj;
      obj_old = obj;
      
      objj = c(objj, obj);
      if (Delta < epsilon) break
  }
  return(list(H=H,objj=objj))
}
