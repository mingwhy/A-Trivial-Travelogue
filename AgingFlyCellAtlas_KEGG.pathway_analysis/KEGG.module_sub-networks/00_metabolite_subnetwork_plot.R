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
table(unlist(lapply(mz.start,length))) #each sub-graph can have multiple start mz.start

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

overlap.mz=mz.age.betas[tmp %in% mz.produced,]$KEGGid
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

for(i in overlap.mz){
    pathway.index=grep(paste0('cpd:',i),names(all_producers))
    names(all_producers)[pathway.index]
    
    if(length(pathway.index)==0){next}
    for(j in pathway.index){
      
      my_producer <- all_producers[[j]] # i can 1 to length(all_producers)
      pathway.id=gsub('hsa|_SIF','',my_producer$cluster)
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
      r=tryCatch(
        pv.out <- pathview(cpd.data = inp.cpd, 
                           pathway.id =pathway.id, species = "hsa", 
                           #out.suffix = "hsa00340_C00388.cpd");
                           out.suffix = cpd),
        error=function(c)'error')
      if(r=='error'){next }
      str(pv.out)
      names(pv.out)
      head(pv.out$plot.data.gene)
      head(pv.out$plot.data.cpd)
      # figure: graphviz
      pv.out <- pathview(cpd.data = inp.cpd, 
                         pathway.id =pathway.id, species = "hsa", 
                         kegg.native = F,
                         out.suffix = cpd);
    }
}  



