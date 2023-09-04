#https://www.preprints.org/manuscript/202101.0280/v1
#https://github.com/babelomics/Metabolica
source('propagationFuncs.R')
load("./files/combat.vals.LUAD.RData") 
dim(combat.vals); #372 genes X 4950 samples
load("./files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)
load("./files/rxn2gene.RData")
dim(rxn2gene) #df,6420 obs. of  4 variables, reaction and its associated enzymes
head(rxn2gene) #3rd column: Entrez 
length(unique(rxn2gene$Reaction)) #1901
combat.vals[1:3,1:3] #gene in row, encoded via Entrez
##########################
# all_producers, in human, precomputed subpathways of KEGG metabolic networks. Each subpathway produces a metabolite. 
x=lapply(1:length(all_producers),function(i){
  c(all_producers[[i]]$isProduced,all_producers[[i]]$initialNodes)
})
length(x)
x[[1]]
x[[30]]
all_producers[[30]]
sum(table(unlist(lapply(x,length)))) #1559
which(unlist(lapply(x,length))==21)
all_producers[[1060]]$initialNodes

##
mz.produced=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$isProduced
})
table(unlist(lapply(mz.produced,length))) #each sub-graph only produce one met
mz.produced=unlist(mz.produced)
which(mz.produced=='cpd:C00388')
all_producers[[718]] #sub-network which produce Histamine
length(unique(mz.produced)) #1270
length(mz.produced) #1559
# one metabolite might be generated via multiple modules

cluster.names=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$cluster
})
cluster.names=unlist(cluster.names)
cluster.names[grep('hsa00340',cluster.names)]
mz.produced[grep('hsa00340',cluster.names)]


mz.start=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$initialNodes
})
table(unlist(lapply(mz.start,length))) #each sub-graph only produce one met

##########################
## in `propagationFuncs.R`
i=718
my_producer <- all_producers[[i]] # i can 1 to length(all_producers)
nodeVals <- do.calc.rxnvals(all_producers,rxn2gene,combat.vals)
dim(nodeVals); #3056 4950
length(unique(unlist(rownames(nodeVals)))) #3056, #node=(rxn+cpd)
nodeVals[1:10,1:2]
table(Matrix::rowSums(nodeVals))

res <- metabolica(nodes.vals=nodeVals, subgraph = my_producer$produceRigraph, 
                  ininodes = my_producer$initialNodes, 
                  endnode = my_producer$isProduced, algR = "minprod", algM = "nsqrt", 
                  maxnum = my_producer$maxnum*1.5)
names(res) #"node.signal" "sig.dif"  
length(res$node.signal) # for this metabolite, node.signal across 4950 samples
tail(res$node.signal)
length(res$sig.dif) #113: 10.00000 29.94719 ... 
# difference bewteen iterations when calculating this node signal

##################################################
#overlap with metabolome measured data
load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
dim(mz.age.betas) #86
head(mz.age.betas)
tmp=paste0('cpd:', mz.age.betas$KEGGid)
sum(tmp %in% mz.produced) #54 overlapped
mz.age.betas[tmp %in% mz.produced,]

###################################################
# pathview: https://bioconductor.org/packages/release/bioc/html/pathview.html
library(pathview)
library(igraph)
load("./files/all_producers.RData")
length(all_producers) #1559, each one <=> one metabolite
names(all_producers)

mz.produced=lapply(1:length(all_producers),function(i){
  all_producers[[i]]$isProduced
})
table(unlist(lapply(mz.produced,length))) #each sub-graph only produce one met
mz.produced=unlist(mz.produced)
which(mz.produced=='cpd:C00388')
all_producers[[718]] #sub-network which produce Histamine
all_producers[[718]]$cluster #"hsa00340_SIF"
grep('hsa00340',names(all_producers)) #14 sub-networks
names(all_producers)[grep('hsa00340',names(all_producers))]
mz.start[grep('hsa00340',names(all_producers))]

grep('C00135', names(all_producers)); #C00135 L-Histidine
names(all_producers)[grep('C00135', names(all_producers))]
#"cpd:C00135-hsa00340_SIF" "cpd:C00135-hsa00410_SIF"
# one metabolite might be generated via multiple modules
i=718; cpd='C00388'; #Histamine
i=710; cpd='C00025'; #L-Glutamate;
# hsa0041, beta-Alanine metabolism,https://www.genome.jp/entry/hsa00410
#	hsa00340, Histidine metabolism,https://www.genome.jp/entry/pathway+hsa00340
pathway.id='hsa00410';
#for(i in grep(pathway.id,names(all_producers))){
for(i in grep(pathway.id,names(all_producers))){
    
  my_producer <- all_producers[[i]] # i can 1 to length(all_producers)
  cpd=gsub('cpd:','',my_producer$isProduced)
  #g<-my_producer$clusteRigraph
  g<-my_producer$produceRigraph
  node.names=V(g) #E(g)
  node.names=names(node.names)
  node.names=node.names[grep('cpd',node.names)]
  node.names
  inp.cpd=rep(1,length(node.names))
  names(inp.cpd)=gsub('cpd:','',node.names)
  my_producer$cluster #"hsa00340_SIF"
  
  # figure: kegg.native 
  pv.out <- pathview(cpd.data = inp.cpd, 
                     pathway.id =gsub('hsa','',pathway.id), species = "hsa", 
                     #out.suffix = "hsa00340_C00388.cpd");
                     out.suffix = paste0(pathway.id,"_",cpd,'.cpd'));
  
  str(pv.out)
  names(pv.out)
  head(pv.out$plot.data.gene)
  head(pv.out$plot.data.cpd)
  # figure: graphviz
  pv.out <- pathview(cpd.data = inp.cpd, 
                     pathway.id =gsub('hsa','',pathway.id), species = "hsa", 
                     kegg.native = F,
                     out.suffix = paste0(pathway.id,"_",cpd,'.cpd'));
}  



