#function [ pvalues_modulePerNet, FDR ] = significantModules( modules, multiNetworks, N, permu_times )

significantModules_diff.size<-function( modules, multiNetworks,permu_times,numCores=1){
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
  
  
  ## Permutation test
  moduleDensity = matrix(0,nrow=moduleCounts,ncol= networkCounts);
  clusterQuality = matrix(0,nrow=moduleCounts,ncol=networkCounts);
  null_counts = array(0,dim=c(moduleCounts, networkCounts, permu_times));
  FDR = matrix(0,nrow=moduleCounts, ncol=networkCounts);
  
  for(m in 1:moduleCounts){
    sprintf('Module %d/%d...\n', m, moduleCounts);
    theModule = modules[[m]];
    moduleN = length(theModule)*(length(theModule)-1)/2;
       
    for(n in 1:networkCounts){
      theNetwork = multiNetworks[[n]];  
      diag(theNetwork)=0;
      # obs, Compute the individual cluster quality and density
      # intersected gene set
      
      theNetwork.c=theNetwork;
      nodeCounts=nrow(theNetwork.c);      
      theModule.c=intersect(theModule,rownames(theNetwork.c))

      obsModule_matrix=theNetwork.c[theModule.c,theModule.c];
      inDen = sum(obsModule_matrix)/2;
      outDen = (sum(theNetwork.c)/2 - inDen)/(nodeCounts-length(theModule.c));
      outDen[outDen==0] = 1e-6;
      moduleDensity[m, n] = inDen/moduleN;
      clusterQuality[m, n] = moduleDensity[m, n]/outDen;
    
        
      # Compute the random individual cluster quality and density
      # And compute individual p-value
      simu.values=parallel::mclapply(1:permu_times, FUN = function(t) {       
        randnum = sample(nodeCounts,replace = F);
        randModule = randnum[1:length(theModule.c)];
        # (x == y) # full matrix format

        randModule_matrix = theNetwork.c[randModule, randModule];
        inDen = sum(randModule_matrix)/2;
        outDen = (sum(theNetwork.c)/2 - inDen)/(nodeCounts-length(theModule.c));
        outDen[outDen==0]= 1e-6;
        ranModuleDensity = inDen/moduleN;
        randClusterQuality = ranModuleDensity/outDen;   
        randClusterQuality
      },mc.cores = numCores)              
    
      null_counts[m,n,]=as.numeric(unlist(simu.values)> clusterQuality[m, n]);
    } 
    
    cat('module',m,'in network',n,'complete,',networkCounts,'nets in total,')
    #for a given module, cell.cluster by permu.i
    #null_counts[m,,]
    cat('raw.p sig in ',sum(apply(null_counts[m,,],1,function(i) sum(i)/permu_times)<0.05),'networks\n')
  }

  #module by cell type matrix
  #pvalues_modulePerNet = apply(null_counts,MARGIN = c(1,2),sum)/permu_times;
  pvalues_modulePerNet = (apply(null_counts,MARGIN = c(1,2),sum)+1)/(1+permu_times);
  for(j in 1:networkCounts){
    FDR[,j] = p.adjust(pvalues_modulePerNet[,j],method='BH');
  }
  
  return(list(pvalues_modulePerNet=pvalues_modulePerNet,FDR=FDR))
}
  