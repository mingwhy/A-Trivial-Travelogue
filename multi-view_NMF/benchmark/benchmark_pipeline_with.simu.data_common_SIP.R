
## generate simu nets
source('src_syn_dataset_common.R')
source('src_evaluation.R')
out=syn_dataset_common(alpha=0.5) #0.4,0.5 and 0.8 all match. 
sapply(out$dataset,dim) #30nets
realLabels=out$realLabels
simu.nets=out$dataset

if(F){
pdf('syn_common.nets.pdf',useDingbats = T)
sapply(1:length(out$dataset),function(i)image(out$dataset[[i]]))
dev.off()
}

## detect modules
library(R.matlab)
source('../ming_pipeline/src_featureNets.R')
M=length(simu.nets)

C.nets=list();
for(i in 1:M){
  mat=simu.nets[[i]]
  mat[mat==1]=0.9999;
  mat[mat==0]=0.0001;
  #max(mat);
  H=-1*mat*log2(mat)  - (1-mat) * log2(1-mat)
  C=1-1/exp(abs(H))
  C.nets[[i]]=C
}
out=featureNets(simu.nets)
strength=out$Strength
participate=out$Participation

diag(strength)=0;
diag(participate)=0;

all.C=Reduce(`+`,C.nets)
#all.C[1:3,1:3]
integrate.net=matrix(0,nrow=nrow(all.C),ncol=ncol(all.C))

for(i in 1:length(simu.nets)){
  x=C.nets[[i]]
  alpha=x/all.C
  integrate.net=integrate.net+ alpha*mat
}
dim(integrate.net)
max(integrate.net)
diag(integrate.net)=0;
information=integrate.net

library(MESS)
MESS::cmd(information,strength) #0.3
MESS::cmd(participate,strength) #0.19
MESS::cmd(information,participate) #0.12 #didn't work at all


## mvNMF
source('../ming_pipeline/src_SNMF.R')
source('../ming_pipeline/src_SNMFforView.R')
source('../ming_pipeline/src_multiViewNMF.R')


start=Sys.time()

#X=list(strength=strength,information=information) #working
X=list(strength=strength,participate=participate)
#X=list(participate=participate,information=information) #not working at all
lambda = c(0.01, 0.01);
maxIter=50;
#K=3;
K=5;
#K=10;
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
xita = 1.5; #cutoff for node seletion in a module
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
apply(module.in.nets$FDR,1,function(i) sum(i<0.05))
#5 10 15 20 25
#5 10  0 20 15  6 16 k=10
#5 20 10 11 16 15 k=10
#15 25 20 k=3