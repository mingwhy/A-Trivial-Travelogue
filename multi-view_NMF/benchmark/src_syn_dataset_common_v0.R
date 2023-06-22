syn_dataset_common <- function(alpha){
#function [ dataset, realLabels ] = syn_dataset_common(alpha, isSaveToFiles, path)
# Synthetic dataset #1
# Conserved modules have the same size and are common to a given set of
# networks.
#
# INPUT:
#   alpha: The probability of the edge connected inside a module
#   isSaveToFiles: a flag for deciding whether to store the results to
#                  files (only for sparse matrix format)
#   path: output file path (it must be denoted if 'isSaveToFiles' is true)
#
# OUTPUT:
#   dataset: the generated dataset
#   realLabels: labels indicating the conserved module to which each point is allocated.
    
    N = 500; # Number of nodes
    M = 30; # Number of networks
    dataset = list();
    C_size = 80; # Size of each cluster
    C1 = (1:C_size);
    C2 = C1 + C_size;
    C3 = C2 + C_size;
    C4 = C3 + C_size;
    C5 = C4 + C_size;
    C = list(C1, C2, C3, C4, C5);
    realLabels = C;
    M_C = length(C) * 5; # Number of networks with defined patterns
    
    p_in = alpha;
    p_out = 0.05;
    if (p_out*(N-C_size)> p_in*(C_size-1)){
        p_out = (p_in*(C_size-1)) / (N-C_size); # (N-n)*p_out < (n-1)*p_in
    }
    
    mc = 1;
    C_temp = list();
    for(m in 1:M){
        # # Background network
        W = matrix(runif(N*N,0,1),nrow=N,ncol=N) #edge weigh from runif
        WW = matrix(0,nrow=N,ncol=N);
        WW[W < p_out] = 1; # [0~p_out] prob -> edges
        
        if (m <= M_C){
            if ((m-1)%%5 == 0){
                C_temp[[mc]]=C[[mc]] #1,6,11,16,21. each 5 has 1,2,3,4,5 modules
                mc = mc + 1;
            }
            for(i in 1:length(C_temp)){
                # inside the cluster
                C_in = C_temp[[i]];
                WW1 = matrix(0,nrow=C_size,ncol=C_size);
                    
                WW1[W[C_in, C_in] < p_in] = 1; 
                WW[C_in, C_in] = WW1;
            }
        }else{
            # One random module in each of the last 5 networks
            s = sample(N,N,replace = F)
            C_in = s[1:80];
            WW1 = matrix(0,length(C_in),length(C_in))
            WW1[W[C_in, C_in] < p_in] = 1; 
            WW[C_in, C_in] = WW1;
        }
        
        # Gaussian noise, sigma=0.1 or 0.15
        WW_tril = WW;
        WW_tril[upper.tri(WW_tril,diag = F)]=0
        E = matrix(rnorm(N^2,0.25, 0.1),N,N);
        E[upper.tri(E,diag=F)] = 0 
        WW_0 = WW_tril + E; # X0+E
        WW_0[WW_tril == 1] = 0;
        WW_0[WW_0 < 0] = 0;
        WW_0[WW_0 > 1] = 1;
        WW_1 = WW_tril - E; # X1-E
        WW_1[WW_tril == 0] = 0;
        WW_1[WW_1 > 1] = 1;
        WW_1[WW_1 < 0] = 0;
        WW = WW_1 + WW_0;
    
        WW = WW - diag(diag(WW));
        WW = WW + t(WW);
        #isSymmetric(WW)
        dataset[[m]] = WW;
    }
    return(list(dataset=dataset,realLabels=realLabels))
}
