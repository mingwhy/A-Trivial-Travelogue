#setSimilarity
setSimilarity<-function( modules ){
  moduleCounts = length(modules);
  HPI = matrix(0,nrow=moduleCounts,ncol=moduleCounts);
  for(i in 1:(moduleCounts-1)){
    module_i = modules[[i]];
    for(j in (i+1):moduleCounts){
      module_j = modules[[j]];
      
      HPI[i, j] = length(intersect(module_i,module_j))/min(length(module_i),length(module_j));
      
    }
  }
  return(HPI)
}