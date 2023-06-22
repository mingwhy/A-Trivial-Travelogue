#function [ TPR, FPR, Accuracy, MCC ] =
#evaluation(modules_final, realLabels, num_Nodes);

evaluation<-function(C_preticted, C_reference, N){
  # Compute the performance measures (TPR, FPR, Accuracy and MCC)
  #
  # INPUT: 
  #       C_preticted: a cell containing the preticted cluster results
  #       C_reference: a cell containing the ground truth
  #       N: total number of items 
  # OUTPUT: 
  #       TPR: the True Positive Rate
  #       FPR: the False Positive Rate
  #       Accuracy:
  #       MCC: the Matthews Correlation Coefficient
  #
  
  ## Confusion matrix
  K_preticted = length(C_preticted);
  K_reference = length(C_reference);
  Ncount_P = 0;
  Ncount_R = 0;
  c_num = matrix(0,nrow=K_preticted+1, ncol=K_reference+1);
  C =  vector(mode = "list", length = K_preticted)
  for(i in 1:K_preticted){
    Ncount_P = Ncount_P + length(C_preticted[[i]]);
    theRef_nonoverlap = c();
    for(j in 1:K_reference){
      C[[i]][[j]] <- intersect(as.numeric(C_preticted[[i]]),as.numeric(C_reference[[j]]));
      c_num[i, j] = length(C[[i]][[j]]);
      if (i == 1){
        Ncount_R = Ncount_R + length(C_reference[j]);
      }
      theRef_nonoverlap = unique(c(theRef_nonoverlap, C[[i]][[j]]))
    }
    c_num[i, j+1] = length(C_preticted[[i]]) - length(theRef_nonoverlap); # Background noise nodes
  }
  
  for(j in 1:K_reference){
    thePre_nonoverlap = c();
    for(i in 1:K_preticted){
      thePre_nonoverlap = unique(c(thePre_nonoverlap, C[[i]][[j]]));
    }
    c_num[i+1, j] = length(C_reference[[j]]) - length(thePre_nonoverlap); # Lost reference nodes
  }
  
  ## TP, FP, FN, TN for nodes pairs
  TP = sum(sum(c_num[1:K_preticted, 1:K_reference]*
                 (c_num[1:K_preticted, 1:K_reference]-1)/2));
  FP = 0;
  for(j in 1:K_reference){
    tempC = c_num[1:K_preticted, (j+1):(K_reference+1),drop=FALSE];
    FP = FP + sum(c_num[1:K_preticted,j]*Matrix::rowSums(tempC));
  }
  FN = 0;
  for(i in 1:K_preticted){
    tempC = c_num[(i+1):(K_preticted+1), 1:K_reference,drop=FALSE];
    FN = FN + sum(c_num[i,1:K_reference]*Matrix::colSums(tempC));
  }
  TN = N*(N-1)/2 - (TP+FN+FP);
  I = rbind(c(TP, FP),c(FN, TN));
  
  ## TPR, FPR, MCC
  TPR = TP/(TP+FN);
  FPR = FP/(FP+TN);
  Accuracy = (TP+TN)/(TP+FP+TN+FN);
  MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
 
  return(list(TPR=TPR,FPR=FPR,Accuracy=Accuracy,MCC=MCC)) 
}