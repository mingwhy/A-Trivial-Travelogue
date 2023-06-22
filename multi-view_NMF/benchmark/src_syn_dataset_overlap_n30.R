#[multiNetworks, realLabels, lables_specific] = ming_syn_dataset_overlap(0.3, true,'simu_dat_overlap/'); 
# 15 nets with 500 noedes each, overlap type
#function [dataset, realLabels, labels_specific] = syn_dataset_overlap(alpha, isSaveToFiles, path)
syn_dataset_overlap<-function(alpha){
# Synthetic dataset #2
# Conserved modules are present only in a subset of networks and they are
# the overlapping parts of specific modules across different networks.
#
# INPUT: 
#   alpha: the probability of the edge connected inside a module
#   isSaveToFiles: a flag for deciding whether to store the results to
#                  files (only for sparse matrix format)
#   path: output file path (it must be denoted if 'isSaveToFiles' is true)
#
# OUTPUT:
#   dataset: the generated dataset
#   realLabels: labels indicating the conserved module to which each point is allocated.
#   labels_specific:lables of each specific module on each network
#

    N = 300; # Number of nodes
    M = 30; # Number of networks
    dataset = list();
    C_common_1 = 1:50; # Common part 1
    C_rest_1 = setdiff((1:N), C_common_1);

    C_common_2 = 201:230; # Common part 2
    C_rest_2 = setdiff((1:N), C_common_2);

    C_rest = intersect(C_rest_1, C_rest_2);
    realLabels = list(C_common_1, C_common_2);
    labels_specific = matrix(0,N, M);
    
    p_in = alpha;
    k = 0;
    for(m in 1:M){
        if(m <= 10){
            #s = C_rest_1[sample(1:length(C_rest_1))];
            increment=1:(m*3)
            C_specific = list(c(C_common_1,max(C_common_1)+increment));
        }else if((10 < m) && (m <= 20)){            
            s_1 = C_rest[sample(1:length(C_rest))];
            #s_1 = C_rest;
            C_specific_1 = c(C_common_1, s_1[1:(5 * (21-m))]); # 
            s_1_rest = setdiff(s_1, C_specific_1);

            s_2 = s_1_rest[sample(1:length(s_1_rest))];
            #s_2 = s_1_rest;
            C_specific_2 = c(C_common_2, s_2[1:(3 * (m-10))]) # 
            C_specific = list(C_specific_1, C_specific_2);  
        }else{
            #s = C_rest_2[sample(1:length(C_rest_2))];         
            increment=1:((m-20)*3) 
            C_specific = list(c(C_common_2,max(C_common_2)+increment));
        }
        
        p_out = 0.05;
        if( p_out*(N-60) > p_in*(60-1) ){
            p_out = (p_in*(60-1)) / (N-60); # (N-n)*p_out < (n-1)*p_in
        }
        # Background network
        W = matrix(runif(N*N,0,1),nrow=N,ncol=N) #edge weigh from runif
        WW = matrix(0,nrow=N,ncol=N);
        WW[W < p_out] = 1; # [0~p_out] prob -> edges
        k = k + 1; 
        
        # inside the module
        for(i in 1:length(C_specific)){
            labels_specific[C_specific[[i]], m] = i;
            WW1 = matrix(0,length(C_specific[[i]]),length(C_specific[[i]]));
            WW1[W[C_specific[[i]], C_specific[[i]]] < p_in] = 1;
            WW[C_specific[[i]], C_specific[[i]]] = WW1;
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
    return(list(dataset=dataset,realLabels=realLabels,labels_specific=labels_specific))
}
