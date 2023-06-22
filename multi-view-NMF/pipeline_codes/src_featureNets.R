
featureNets<-function(multiNetworks){
  network_count = length(multiNetworks);
  n=nrow(multiNetworks[[1]])
  
  Strength = matrix(0,nrow=n,ncol=n)
  temp = Strength
  A = Strength
  for(k in 1:network_count){
    
    # The edge weight is transformed using a logistic function, such that for the
    # element less than 0.3, we make it close to 0; for the element more than
    # 0.6, we make it close to 1.
    weight_max = max(multiNetworks[[k]]);
    weight_min = 0;
    matrix_weight_1 = (multiNetworks[[k]] - weight_min) / (weight_max-weight_min);
    matrix_weight_2 = 1/(1+exp(log(9999)-2*log(9999)*matrix_weight_1));
    
    matrix_weight_2[matrix_weight_1 <= 0.3] = 0;
    A = A + matrix_weight_2;
    temp = temp + matrix_weight_2^2;
    Strength = Strength + multiNetworks[[k]];
  }
  
  Participation = (network_count/(network_count-1)) * (1-(temp/(A^2)));
  Participation[is.na(Participation)]=0 #diagonal is 0
  Participation[is.nan(Participation)]=0
  diag(Participation)=0;
  
  Strength = A/network_count;
  diag(Strength)=0  #The diagonal is 0
  return(list(Strength=Strength,Participation=Participation))
}
