
## generate simu nets
source('src_syn_dataset_common.R')
source('src_evaluation.R')
source('src_syn_dataset_overlap.R')
out=syn_dataset_common(alpha=0.4) #0.4, 0.5,0.6 ok. 15 nets with 500 noedes each, overlap type
names(out)
sapply(out$dataset,dim)

#out=syn_dataset_overlap(alpha=0.6) 

if(F){
  pdf('syn_overlap.nets.pdf',useDingbats = T)
  sapply(1:length(out$dataset),function(i)image(out$dataset[[i]]))
  dev.off()
}
realLabels=out$realLabels
labels_specific=out$labels_specific
dim(labels_specific)

simu.nets=out$dataset

## detect modules, strength and uniformity
library(R.matlab)
multiNetworks=simu.nets;
network_count=length(multiNetworks) #number of simulated matrix
gene_count=nrow(multiNetworks[[1]])

if(F){
  out=featureNets(simu.nets)
  strength=out$Strength
  participate=out$Participation
  
  diag(strength)=0;
  diag(participate)=0;
}

# data transform
if(T){ #yes transform
  multiNetworks.transform=list()
  for(k in 1:network_count){
    weight_max = max(multiNetworks[[k]]);
    weight_min = 0;
    matrix_weight_1 = (multiNetworks[[k]] - weight_min) / (weight_max-weight_min);
    matrix_weight_2 = 1/(1+exp(log(9999)-2*log(9999)*matrix_weight_1));
    
    #matrix_weight_2[matrix_weight_1 <= 0.3] = 0;
    multiNetworks.transform[[k]]=matrix_weight_2
  }
}
#multiNetworks.transform=multiNetworks;
# compute strength matrix
all.weigh=Reduce(`+`,multiNetworks.transform)
strength=all.weigh/network_count
diag(strength)=0;

# compute uniformity matrix
tmp=array(0,dim=c(network_count,gene_count,gene_count)) #number of matrix, nrow, ncol
for(i in 1:network_count){
  tmp[i,,]=multiNetworks.transform[[i]]
}
max.weigh <- apply(tmp, c(2,3), max)
dim(max.weigh)
max.weigh[1:3,1:3]

C.nets=list();
for(i in 1:network_count){
  mat=multiNetworks.transform[[i]] #all pos values
  #C.nets[[i]]=1-(mat/max.weigh)^2 #no need square
  C.nets[[i]]=1-(mat/max.weigh)
}

all.C=Reduce(`+`,C.nets)
uni=1 - all.C/network_count
diag(uni)=0;
uni[is.nan(uni)]=0;
uni[1:3,1:3]
max(uni) #0.74
min(uni) #0
uniformity=uni

library(MESS)
MESS::cmd(uniformity,strength) #0.3


## mvNMF
source('../ming_pipeline/src_SNMF.R')
source('../ming_pipeline/src_SNMFforView.R')
source('../ming_pipeline/src_multiViewNMF.R')


start=Sys.time()
X=list(strength=strength,uniformity=uniformity)
#X=list(strength=strength,information=information) #working
#X=list(strength=strength,participate=participate) #working
#X=list(participate=participate,information=information) #not working at all
lambda = c(0.01, 0.01);
maxIter=50;
#K=3
K=5;
#K=10
out=multiViewNMF( X, K, lambda, maxIter)
names(out)
end=Sys.time()
print(end-start)

## select node to become modules
Hc=out$Hc
Hc[1:10,1:3]
sum(Hc<0) #0, non-negative matrix

source('../ming_pipeline/src_moduleNodesSelection.R')
source('../ming_pipeline/src_setSimilarity.R')
xita =1; #cutoff for node seletion in a module
modules_final = moduleNodesSelection( Hc, xita ); 
length(modules_final)
sapply(modules_final,length)

## evaluate
num_Nodes=500
out=evaluation(modules_final, realLabels, num_Nodes)
out

## see what's their significance in each simulated net
source('../ming_pipeline/src_significantModules.R')
module.in.nets=significantModules(modules_final, simu.nets,
                                  num_Nodes,permu_times=100);
names(module.in.nets)
apply(module.in.nets$FDR,1,function(i) sum(i<0.05)) ##modules detected in each net


