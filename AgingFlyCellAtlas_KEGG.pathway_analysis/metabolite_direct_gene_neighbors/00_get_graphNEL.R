library(pathview)
library(XML)
library(graph) #use graph not igraph!!!
library(KEGGgraph)
library(tidyverse)

##########################################################
## read in mz measurement data for trajectory comparison
load('~/Documents/aging_metabolism/batch_block_corrected_mzdata')
ls() #'dat' 'mzs'
head(dat)
dim(dat) #181  93
dat[1:3,1:8]
sort(unique(dat$AgeNum)) #4 10 24 45 69 80 days
colnames(dat)[1:8]

load('~/Documents/aging_metabolism/age-associated.metabolites.for.Ming')
ls() #"cmat",age-independent covariance between metabolites 
# "for.Ming",  mz ~ age betas 
mz.age.betas=for.Ming
mz.age.betas=mz.age.betas[mz.age.betas$KEGGid!='',]
dim(mz.age.betas) #84
###########################
# get all kegg pathway names of a species
# https://support.bioconductor.org/p/109871/
if(!file.exists('all_mapped_pathways.rds')){
  library(KEGGREST)
  if(F){
    dme_all_path  <- keggLink("pathway", "dme") %>% 
      tibble(pathway = ., eg = sub("dme:", "", names(.)))
    dme_all_pathways=unique(dme_all_path$pathway)
    length(dme_all_pathways) #142
  }
  # locate metabolite's all associated pathways
  all.pathways<-lapply(mz.age.betas$KEGGid,function(kegg_id_name){
    path_ids = map(kegg_id_name, ~keggLink('pathway', .)) %>% 
      unlist() %>%
      unique()
    length(path_ids) #12
    path_ids=gsub('map','dme',path_ids)
    #path_ids='path:dme00340' #only look at one pathway
    path_ids
  })
  names(all.pathways)<-mz.age.betas$KEGGid
  KEGGid=rep(names(all.pathways),sapply(all.pathways,length))
  pathway=unlist(all.pathways)
  df.all.pathway=data.frame(KEGGid=KEGGid,pathway=pathway)
  saveRDS(df.all.pathway,'all_mapped_pathways.rds')
}
df.all.pathway=readRDS('all_mapped_pathways.rds')
##########
# download xml files from KEGG database
if(F){
  length(unique(df.all.pathway$pathway)) #130
  for(i in unique(df.all.pathway$pathway)){
    tempfile=paste0(gsub('path:','',i),'.xml')
    #tempfile='dme00340.xml'; #file would be download and saved in this tempfile
    if(!file.exists(tempfile)){
      KEGGgraph::retrieveKGML(pathwayid="00340", organism="dme", destfile=tempfile,quiet=TRUE)
    }
  }
}
# save xml files to kgml_xmls folder
########## begin pathview
## now, `parseKGML2.R`

files=Sys.glob('kgml_xmls/*')
out.dir='granphNEL';
dir.create(out.dir)

for(file in files){
  #file='kgml_xmls/dme00340.xml'
  pathway.id=basename(file)
  out.file=paste0(out.dir,'/',pathway.id,'_gr.rds')
  if(file.exists(out.file)){next}
  
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
  if (genesOnly) {gR <- subGraphByNodeType(gR, "gene")}
  #return(gR)
  gR
  saveRDS(gR,out.file)
}

######################################################
# R package graph, https://bioconductor.org/packages/release/bioc/html/graph.html
#https://bioconductor.org/packages/release/bioc/vignettes/graph/inst/doc/GraphClass.html
#nodes(gR);edges(gR);degree(gR) #inDegree, outDegree
# locate target metabolite and extract 1-step neighbors for plot
(gr.files=Sys.glob('granphNEL/*gr.rds'))

yes.mz=c();
all.pathways=readRDS('all_mapped_pathways.rds')
pdf('test.pdf')
for(target.mz in unique(all.pathways$KEGGid)){
  cat('start ',target.mz,'\n')
  mapped.path=all.pathways[all.pathways$KEGGid==target.mz,]$pathway
  mapped.path=gsub('path:','',mapped.path)
  
  for(path in mapped.path){
    file=gr.files[grep(path,gr.files)]
    gR=readRDS(file)
    node.data=try(node.info(gR), silent=T)
       
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
                                  category="SYMBOL", pkg.name="org.Dm.eg.db")[,2]
      mapped.gnodes=rownames(plot.data.gene) #EntryID in kegg
      node.data$labels[mapped.gnodes]
      plot.data.gene$labels
      #node.data$labels[mapped.gnodes]=plot.data.gene$labels
    }
    
    # map compounds, use node.map function in pathview
    cpd.data=target.mz
    plot.data.cpd=node.map(cpd.data, node.data, node.types="compound", node.sum="sum")
    
    plot.data=rbind(plot.data.gene,plot.data.cpd)

    ##
    i=which(node.data$kegg.names==target.mz)
    if(length(i)==0){next}
    target.mz.id=names(node.data$kegg.names[i]) #'99', C00388
    gi=igraph::graph_from_graphnel(gR) 
    
    out1=igraph::neighbors(gi,target.mz.id,mode='out') #gene, 1-step neighbor from target metabolite
    in1=igraph::neighbors(gi,target.mz.id,mode='in') #gene, 1-step neighbor from target metabolite
    vids=c(target.mz.id,out1$name,in1$name)
    
    if(length(out1)>0){
      out2=unlist(lapply(out1$name,function(i){
        #node.data$kegg.names[[i]]
        if(node.data$type[[i]]=='gene'){
          x=igraph::neighbors(gi,i,mode='out') #compound, 2-step neighbor from target metabolite
          return(x$name)
        }else {return(NA)}
      }))
      vids=c(vids,out2)
    }
    if(length(in1)>0){
      in2=unlist(lapply(in1$name,function(i){
        if(node.data$type[[i]]=='gene'){
          x=igraph::neighbors(gi,i,mode='in') #compound, 2-step neighbor from target metabolite
          return(x$name)
        }else {return(NA)}
      }))
      vids=c(vids,in2)
    }
    vids=vids[!is.na(vids)]
    vids2=unlist(lapply(vids,function(i) if(any(node.data$type[[i]] %in% c('gene','compound'))) i))
    
    if(length(vids2)==1){next}
    sgi=igraph::induced_subgraph(gi, vids=vids2)
    #plot(sgi)
    igraph::V(sgi)$label=plot.data[igraph::V(sgi)$name,]$labels
    v.types=plot.data[igraph::V(sgi)$name,]$type
    v.shapes=ifelse(v.types=='gene','square','circle')
           
    plot(sgi, vertex.shape=v.shapes,vertex.size=10, #vertex.color = 'lightblue', edge.color="darkgreen", edge.label=links$value,
         vertex.label=igraph::V(sgi)$label, vertex.label.font=1, 
         vertex.label.cex = 1 , layout=igraph::layout.fruchterman.reingold,
         main=paste0(target.mz,', ',path))
    yes.mz=c(yes.mz,target.mz)
  }
}
dev.off()

unique(yes.mz) #only 4 metabolites have info
#for example:https://www.genome.jp/pathway/map01240+C01081
#C01081, no direct neighbors contain enzymes
#but if we look at GEM: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/reaction/MAR04208
