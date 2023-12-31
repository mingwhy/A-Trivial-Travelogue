library(pathview)
library(XML)
library(graph) #use graph not igraph!!!
library(KEGGgraph)
library(tidyverse)
(files=Sys.glob('pathview_Rscripts/*'))
#for(file in files){source(file)}
###########################
# get all kegg pathway names of a species
# https://support.bioconductor.org/p/109871/
if(F){
library(KEGGREST)
dme_all_path  <- keggLink("pathway", "dme") %>% 
  tibble(pathway = ., eg = sub("dme:", "", names(.)))
dme_all_pathways=unique(dme_all_path$pathway)
length(dme_all_pathways) #142

# locate metabolite's all associated pathways
kegg_id_names = c('C00388')
path_ids = map(kegg_id_names, ~keggLink('pathway', .)) %>% 
  unlist() %>%
  unique()
length(path_ids) #12
path_ids=gsub('map','dme',path_ids)
#path_ids='path:dme00340' #only look at one pathway
}

##########
# download xml files from KEGG database
tempfile='hsa04010.xml'; #file would be download and saved in this tempfile
tempfile='hsa00340.xml';
#tempfile='dme00340.xml'; #file would be download and saved in this tempfile
if(!file.exists(tempfile)){
  #KEGGgraph::retrieveKGML(pathwayid="04010", organism="hsa", destfile=tempfile,quiet=TRUE)
  #KEGGgraph::retrieveKGML(pathwayid="00340", organism="dme", destfile=tempfile,quiet=TRUE)
  KEGGgraph::retrieveKGML(pathwayid="00340", organism="hsa", destfile=tempfile,quiet=TRUE)
}

########## begin pathview
#in `pathview.R`,
#gR1=try(parseKGML2Graph2(xml.file[i], genes=F, expand=expand.node, split.group=split.group), silent=T)
#in `parseKGML2Graph2.R`
#parseKGML2Graph2 <-function (file, ...)
#{
#  pathway <- parseKGML2(file)
#  gR <- KEGGpathway2Graph2(pathway, ...)
#  return(gR)
#}

## now, `parseKGML2.R`
#file='hsa00340.xml'
file='dme00340.xml'
doc <- XML::xmlTreeParse(file, getDTD = FALSE)
r <- xmlRoot(doc)
childnames <- sapply(xmlChildren(r), xmlName)
isEntry <- childnames == "entry"
isRelation <- childnames == "relation"
isReaction <- childnames == "reaction"
kegg.pathwayinfo <- parsePathwayInfo(r)
kegg.nodes <- sapply(r[isEntry], parseEntry)
kegg.edges <- sapply(r[isRelation], parseRelation)
source('pathview_Rscripts/parseReaction2.R')
kegg.reactions <- sapply(r[isReaction], parseReaction2)
names(kegg.nodes) <- sapply(kegg.nodes, getEntryID)
pathway <- new("KEGGPathway", pathwayInfo = kegg.pathwayinfo,
               nodes = kegg.nodes, edges = kegg.edges, reactions = kegg.reactions)
#return(pathway)
length(kegg.nodes) #87
length(kegg.edges) #24
length(kegg.reactions) #16

## now, `KEGGpathway2Graph2.R`
genesOnly = FALSE;expandGenes = FALSE; split.group=FALSE; check.reaction=TRUE;

stopifnot(is(pathway, "KEGGPathway"))
if(split.group) pathway <- splitKEGGgroup(pathway)
rdata=(pathway@reactions)

if (expandGenes){
  if(check.reaction & length(rdata)>0) message("Note: ", "Gene nodes not expanded when reactions are converted to edges!")
  else pathway <- expandKEGGPathway(pathway)
}
knodes <- graph::nodes(pathway)
kedges <- graph::edges(pathway)
node.entryIDs <- getEntryID(knodes) #kegg.nodes
edge.entryIDs <- getEntryID(kedges) #each row, a pair of nodes
V <- node.entryIDs
edL <- vector("list", length = length(V))
names(edL) <- V
if (is.null(nrow(edge.entryIDs))) {
  for (i in seq(along = edL)) {
    edL[[i]] <- list()
  }
}else {
  for (i in 1:length(V)) {
    id <- node.entryIDs[i]
    hasRelation <- id == edge.entryIDs[, "Entry1ID"]
    if (!any(hasRelation)) {
      edL[[i]] <- list(edges = NULL)
    }
    else {
      entry2 <- unname(unique(edge.entryIDs[hasRelation,
                                            "Entry2ID"]))
      edL[[i]] <- list(edges = entry2)
    }
  }
}
length(V) #87
length(edL) #87, each element <=> each node
gR <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")

if(check.reaction & length(rdata)>0){
  r2e.res=reaction2edge(pathway, gR) #`reaction2edge.R`
  # `reaction2edge` is the critical function to construct gene-compound graph
  gR=r2e.res[[1]]
  kedges=r2e.res[[2]]
  knodes=r2e.res[[3]]
}
length(kedges) #36
length(knodes) #88

names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x),
                                                  collapse = "~"))
env.node <- new.env()
env.edge <- new.env()
assign("nodes", knodes, envir = env.node)
assign("edges", kedges, envir = env.edge)
nodeDataDefaults(gR, "KEGGNode") <- env.node
edgeDataDefaults(gR, "KEGGEdge") <- env.edge
if (genesOnly) {
  gR <- subGraphByNodeType(gR, "gene")
}
#return(gR)
gR

######################################################
# R package graph, https://bioconductor.org/packages/release/bioc/html/graph.html
#https://bioconductor.org/packages/release/bioc/vignettes/graph/inst/doc/GraphClass.html
#nodes(gR);edges(gR);degree(gR) #inDegree, outDegree
# locate target metabolite and extract 1-step neighbors for plot

node.data=try(node.info(gR), silent=T)
table(node.data$type)
length(node.data$kegg.names) #88 nodes

type.sel=node.data$type %in% c('compound')
node.data$kegg.names[type.sel]

i=which(node.data$kegg.names=='C00388')
names(node.data$kegg.names[i]) #'99', C00388
adj(gR,'99') #only out-neighbors of C00388

edges=edgeData(gR)
edges[grep('99',names(edges))] 

# map genes, use node.map function in pathview
gene.node.type='gene'
gene.data=NULL
plot.data.gene=node.map(gene.data, node.data,
                        node.types=gene.node.type, node.sum="sum", 
                        entrez.gnodes=TRUE)
head(plot.data.gene) #check if plot.data.gene already have labels column
#rownames(plot.data.gene); #EntryID in kegg

if(F){
  library(AnnotationDbi) #in `NAMESPACE`, to use `columns` function in `eg2id`
  plot.data.gene$labels=eg2id(as.character(plot.data.gene$kegg.names), 
                              category="SYMBOL", pkg.name="org.Hs.eg.db")[,2]
  mapped.gnodes=rownames(plot.data.gene) #EntryID in kegg
  node.data$labels[mapped.gnodes]
  plot.data.gene$labels
  #node.data$labels[mapped.gnodes]=plot.data.gene$labels
}

# map compounds, use node.map function in pathview
cpd.data='C00388'
plot.data.cpd=node.map(cpd.data, node.data, node.types="compound", node.sum="sum")

plot.data=rbind(plot.data.gene,plot.data.cpd)

##
i=which(node.data$kegg.names=='C00388')
i
target.mz.id=names(node.data$kegg.names[i]) #'99', C00388
gi=igraph::graph_from_graphnel(gR) 

out1=igraph::neighbors(gi,target.mz.id,mode='out') #gene, 1-step neighbor from target metabolite
in1=igraph::neighbors(gi,target.mz.id,mode='in') #gene, 1-step neighbor from target metabolite
vids=c(target.mz.id,out1$name,in1$name)
if(length(out1)>0){
  out2=unlist(lapply(out1$name,function(i){
    x=igraph::neighbors(gi,i,mode='out') #compound, 2-step neighbor from target metabolite
    x$name
  }))
  vids=c(vids,out2)
}
if(length(in1)>0){
  in2=unlist(lapply(in1$name,function(i){
    x=igraph::neighbors(gi,i,mode='in') #compound, 2-step neighbor from target metabolite
    x$name
  }))
  vids=c(vids,in2)
}
sgi=igraph::induced_subgraph(gi, vids=vids)
plot(sgi)
igraph::V(sgi)$label=plot.data[igraph::V(sgi)$name,]$labels
node.types=plot.data[igraph::V(sgi)$name,]$type
node.shapes=ifelse(node.types=='gene','square','circle')
       
plot(sgi, vertex.shape=node.shapes,vertex.size=10, #vertex.color = 'lightblue', edge.color="darkgreen", edge.label=links$value,
     vertex.label=igraph::V(sgi)$label, vertex.label.font=1, 
     vertex.label.cex = 1 , layout=igraph::layout.fruchterman.reingold)


