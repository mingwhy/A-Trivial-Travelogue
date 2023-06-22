#function [ pvalues_modulePerNet, FDR ] = significantModules( modules, multiNetworks, N, permu_times )

significantModules<-function( modules, multiNetworks, permu_times,numCores=2){
# Using a permutation test to assess the significance of functional modules across multiple networks.
# This allows identifying the specific conditions where each module is detected.
# 
# INPUT:
#   modules: a list of extracted modules with genes' index
#   multiNetworks: a cell containing adjacency matrices of multiple
#   networks
#   N: the total number of nodes in multiple networks
#   permu_times: the number of permutation times(default value is 1000)
#
# OUTPUT:
  #   pvalues_modulePerNet: individual p-value of each module in each network
#   FDR: Benjamin-Hochberg adjusted p-values
#
  
  ## Initialization
  moduleCounts = length(modules);
  networkCounts = length(multiNetworks);
  nodeCounts = nrow(multiNetworks[[1]]);
  x=nrow(multiNetworks[[1]])
  y=ncol(multiNetworks[[1]])
  
  ## Permutation test
  moduleDensity = matrix(0,nrow=moduleCounts,ncol= networkCounts);
  clusterQuality = matrix(0,nrow=moduleCounts,ncol=networkCounts);
  null_counts = array(0,dim=c(moduleCounts, networkCounts, permu_times));
  FDR = matrix(0,nrow=moduleCounts, ncol=networkCounts);
  
  for(m in 1:moduleCounts){
    sprintf('Module %d/%d...\n', m, moduleCounts);
    theModule = modules[[m]];
    moduleN = length(theModule)*(length(theModule)-1)/2;
    
    #obs values
    for(n in 1:networkCounts){
      theNetwork = multiNetworks[[n]];
      # Compute the individual cluster quality and density
      # x==y, full matrix format
      inDen = sum(sum(theNetwork[theModule, theModule]))/2;
      outDen = (sum(sum(theNetwork))/2 - inDen)/(nodeCounts-length(theModule));
      outDen[outDen==0] = 1e-6;
      moduleDensity[m, n] = inDen/moduleN;
      clusterQuality[m, n] = moduleDensity[m, n]/outDen;
    
        
      # Compute the random individual cluster quality and density
      # And compute individual p-value
       simu.values=parallel::mclapply(1:permu_times, FUN = function(t) {         
        randnum = sample(nodeCounts,replace = F);
        randModule = randnum[1:length(theModule)];
        # (x == y) # full matrix format
        randModule_matrix = theNetwork[randModule, randModule];
        inDen = sum(sum(randModule_matrix))/2;
        outDen = (sum(sum(theNetwork))/2 - inDen)/(nodeCounts-length(theModule));
        outDen[outDen==0]= 1e-6;
        ranModuleDensity = inDen/moduleN;
        randClusterQuality = ranModuleDensity/outDen;
        randClusterQuality
      },mc.cores=numCores)

      null_counts[m,n,]=as.numeric(unlist(simu.values)> clusterQuality[m, n]);
    }
    cat('module',m,'in network',n,'complete,',networkCounts,'nets in total\n')
  }

  #null_counts[,,2]
  #pvalues_modulePerNet = apply(null_counts,MARGIN = c(1,2),sum)/permu_times;  
  pvalues_modulePerNet = (apply(null_counts,MARGIN = c(1,2),sum)+1)/(1+permu_times);
  for(j in 1:networkCounts){
    FDR[,j] = p.adjust(pvalues_modulePerNet[,j],method='BH');
  }
  
  return(list(pvalues_modulePerNet=pvalues_modulePerNet,FDR=FDR))
}
  