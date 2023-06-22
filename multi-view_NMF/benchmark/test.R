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
    
    N = 800; # Number of nodes
    M = 100; # Number of networks
    dataset = list();
    C_size = 30; # Size of each cluster

    C1 = (1:C_size);
    C2 = C1 + C_size;
    C3 = C2 + C_size;
    C4 = C3 + C_size;
    C5 = C4 + C_size;
    C = list(C1, C2, C3, C4, C5); #gene modules
    realLabels = C;
    
    #simulate a matrix: module by network
    r=matrix(0,nrow=length(C),ncol=M)
    Module_occur=c(1,2,5,8,10)*10; #module1 in 10nets, module2 in 20nets,... module5 in 100nets.
    for(i in 1:length(C)){
        r[i,1:Module_occur[i]]=1
    }
    # apply(r,1,sum)

    M_C = max(Module_occur); # Number of networks with defined patterns
    
    p_in = alpha;
    p_out = 0.05;
    if (p_out*(N-C_size)> p_in*(C_size-1)){
        p_out = (p_in*(C_size-1)) / (N-C_size); # (N-n)*p_out < (n-1)*p_in
    }
    
    
    for(m in 1:M){
        C_temp = list();

        # # Background network
        W = matrix(runif(N*N,0,1),nrow=N,ncol=N) #edge weigh from runif
        WW = matrix(0,nrow=N,ncol=N);
        WW[W < p_out] = 1; # [0~p_out] prob -> edges
        
        num.module=sum(r[,m])
        C_temp=C[which(r[,m]==1)];

        if (m <= M_C){            
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
            C_in = s[1:C_size];
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
