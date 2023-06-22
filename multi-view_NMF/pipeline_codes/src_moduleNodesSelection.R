
moduleNodesSelection<-function( Hc, xita ){
  # Assigning the module members by a soft node selection procedure
  # and then truing the modules to obtain more accurate results
  #
  # INPUT:
  #   Hc: the consensus factor matrix
  #   xita: the parameter for selecting nodes
  #
  # OUTPUT:
  #   modulesFinal: a cell which contains the final result modules
  #
  
  N = nrow(Hc);
  K = ncol(Hc);
  
  candidateModules = list();
  moduleSignal = list()
  H_mean = Matrix::colMeans(Hc);
  #A = rbind(c(4,-5,1),c(2, 3, 5),c(-9, 1, 7),c(10,2,3));
  #S = apply(A,2,sd)
  H_std = apply(Hc, 2, sd);

  for(k in 1:K){
    i=which(Hc[,k] > H_mean[k] + xita*H_std[k]);  # Z-score>=t
    if(length(i)==0){
      candidateModules[[k]]=NULL;
      moduleSignal[[k]] = NULL
    }else{
      candidateModules[[k]] =i     
      moduleSignal[[k]] = mean(Hc[i, k]);
    }
  }
  #remove NULL in a list
  candidateModules=Filter(Negate(is.null), candidateModules)
  moduleSignal=Filter(Negate(is.null), moduleSignal)

  HPI = setSimilarity( candidateModules );
  modulesFinal = candidateModules;
  
  for(roundi in 1:2){

    for(i in 1:(nrow(HPI)-1) ){
      for(j in (i+1):ncol(HPI)){
        if(HPI[i,j]>=0.3){ # merge these two modules based on mean module signal
          I = which.max(c(moduleSignal[[i]], moduleSignal[[j]]));
          if(I == 1){
            #modulesFinal[[j]] = NULL; #i is bigger
            modulesFinal[[i]]= unique(c(modulesFinal[[i]],modulesFinal[[j]]))
            modulesFinal[[j]] = ''; #i is bigger
            moduleSignal[[j]] = 0;
            HPI[j,] = rep(0,ncol(HPI));
            HPI[,j] = rep(0,nrow(HPI));
          }else{
            #modulesFinal[[i]] = NULL;
            modulesFinal[[j]]= unique(c(modulesFinal[[i]],modulesFinal[[j]]))
            modulesFinal[[i]] = '';
            moduleSignal[[i]] = 0;
            HPI[i,] = rep(0,ncol(HPI));
            HPI[,i] = rep(0,nrow(HPI));
          }
        }
      }
    }
    
    # Only modules with no less than 5 nodes are kept
    remove.list=c()
    for(i in 1:length(modulesFinal)){
      if(any(is.null(modulesFinal[[i]]),length(modulesFinal[[i]])==1)){
        remove.list=c(remove.list,i)     
      }
    }
    if(length(remove.list)!=0){
      modulesFinal=modulesFinal[-remove.list]
      moduleSignal=moduleSignal[-remove.list]
    }
    HPI=setSimilarity(modulesFinal)
  }

  for(i in 1:length(modulesFinal)){
      if(any(is.null(modulesFinal[[i]]),length(modulesFinal[[i]])<5)){
        remove.list=c(remove.list,i)     
      }
  }
  if(length(remove.list)!=0){
      modulesFinal=modulesFinal[-remove.list]
  }
  return(modulesFinal)
}


      