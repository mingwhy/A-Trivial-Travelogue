
library(ESCO)

# in my case, there is some error when running 'makeCluster'
# I took suggestion from: https://stackoverflow.com/questions/62730783/error-in-makepsockclusternames-spec-cluster-setup-failed-3-of-3-work
# to modify my ~/.Rprofile and solved this issue.

corr.single = list()
corr.single[[1]] = diag(75)
corr.single[[1]][1:15, 1:15] <- corr.single[[1]][16:30, 16:30] <- corr.single[[1]][31:45, 31:45] <- 
  corr.single[[1]][46:60, 46:60] <- corr.single[[1]][61:75, 61:75]  <- 0.9
corr.single[[1]][1:15, 16:30] <- corr.single[[1]][16:30, 31:45] <- corr.single[[1]][31:45, 46:60] <-
  corr.single[[1]][46:60, 61:75] <- corr.single[[1]][16:30, 1:15] <- corr.single[[1]][31:45, 16:30] <- 
  corr.single[[1]][46:60, 31:45] <- corr.single[[1]][61:75, 46:60] <- 0.7
corr.single[[1]][1:15, 31:45] <- corr.single[[1]][16:30, 46:60] <- corr.single[[1]][31:45, 61:75] <-
  corr.single[[1]][31:45, 1:15] <- corr.single[[1]][46:60, 16:30] <- 
  corr.single[[1]][61:75, 31:45] <- 0.5
corr.single[[1]][1:15, 46:60] <- corr.single[[1]][16:30, 61:75] <- corr.single[[1]][46:60, 1:15] <-
  corr.single[[1]][61:75, 16:30] <- 0.3
corr.single[[1]][1:15, 61:75] <- corr.single[[1]][61:75, 1:15] <- 0.1
diag(corr.single[[1]]) = 1

#alpha = 0.1; lib.loc = 9;
sim.single <- escoSimulateSingle(nGenes = 100, nCells = 200, 
                                 #numCores=2, 
                withcorr = TRUE, corr = corr.single, verbose = TRUE)
corr.single = metadata(sim.single)$Params@corr
dim(corr.single[[1]]) #75 75

gene.order = colnames(corr.single[[1]])[hclust(dist(corr.single[[1]]))$order]
corr.single[[1]][gene.order[seq(1, 75, 15)], gene.order[seq(1, 75, 15)]]
gene.order = gene.order[c( 1:15, 16:30, 31:45, 61:75, 46:60)]
image(corr.single[[1]][gene.order, gene.order])

tmp=corr.single[[1]][gene.order, gene.order]
tmp[1:3,1:3]

gene.order = c(gene.order, setdiff(paste0('Gene', 1:100), gene.order))

# get the data
datalist = list("simulated truth" = assays(sim.single)$TrueCounts, 
                "zero-inflated" = assays(sim.single)$counts, 
                "down-sampled" = assays(sim.single)$observedcounts)
datalist$`simulated truth`[1:3,1:3] #gene by cell matrix

source('../src_estimate_K_RMT.R')
estimate_K(datalist$`simulated truth`)
estimate_K(datalist$`zero-inflated`)
estimate_K(datalist$`down-sampled`)

##########################################
## use syn_dataset_common as input
source('src_syn_dataset_common.R');

out=syn_dataset_common(alpha=0.5) #0.4,0.5 and 0.8 all match. 
sapply(out$dataset,dim) #30nets
realLabels=out$realLabels
simu.nets=out$dataset
dim(simu.nets[[1]]) #a cor.mat ngene by ngene

sim.single <- escoSimulateSingle(nGenes = 500, nCells = 200, 
                                 #numCores=2, 
              withcorr = TRUE, corr =list(simu.nets[[1]]), verbose = TRUE)

corr.single = metadata(sim.single)$Params@corr
dim(corr.single[[1]]) #75 75
corr.single[[1]][1:3,1:3]
