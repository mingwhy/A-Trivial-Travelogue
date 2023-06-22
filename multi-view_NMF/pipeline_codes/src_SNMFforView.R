#SNMFforView  [H{i}] = SNMFforView(X[[i]], Hc, H[[i]], lambda[i], maxIterforView, 1e-6);

SNMFforView<-function(X, Hc, H, lambda, maxIter, epsilon){
  # Multi-View Non-negative symmetric Matrix Factorization for each view
  # 
  # INPUT:
  #   X: the adjacency matrix of a network
  #   Hc: initialization for consensus factor matrix 
  #   H: initialization for factor matrix of each view
  #   lambda: a vector which contains the parameters for balancing the relative
  #           weight among different views
  #   MaxIter: the maximal number of iterations for alternating minimization
  #   epsilon: the convergence parameter
  #
  # OUTPUT:
    #   H: the factor matrix
  #
  # Peizhuo Wang (wangpeizhuo_37@163.com)
  
  obj_old = norm(X - H %*% t(H), type='F')^2 + 
            lambda*norm(H - Hc, type='F')^2;
  
  # Fixing Hc, minimize objective function over H
  #eps=pracma::eps(1); #floating-point relative accuracy
  eps=2.220446e-16;
  for(iter in 1:maxIter){
      # Update rule
      temp_1 = 2* X %*% H + lambda*Hc;
      temp_2 = 2* H %*% (t(H) %*% H) + lambda*H;
      H = H * (temp_1/(temp_2+eps));
    
      # Objective function
      obj_body = norm(X - H %*% t(H), type='F')^2;
      obj_consensus = norm(H - Hc, type='F')^2;
      obj = obj_body + lambda*obj_consensus;
    
      Delta = obj_old - obj;
      obj_old = obj;
        
      if(Delta < epsilon) break
  }
  return(H)
}