#https://robertness.github.io/2015/07/30/Simulating-a-Phosphorylation-Network-from-KEGG-Pathways-Part-1.html
#https://robertness.github.io/2015/07/31/Simulating-a-Phosphorylation-Network-from-KEGG-Pathways-Part-2.html

# https://gist.github.com/robertness/ecb79071301924f43d09
# Download MAPK pathway phosphorylations into an igraph object in R
# hsa04010: https://www.genome.jp/entry/hsa04010
library(igraph)
library(KEGGgraph, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
# install.packages("devtools")
#devtools::install_github("robertness/lucy")
#library(lucy)
## instead of installing lucy, load functions mannually
#https://github.com/robertness/lucy/blob/master/R/graph-plots.R
igraphviz <- function(g, main = NULL){
  g <- name_vertices(g)
  gnell <- g %>% igraph.to.graphNEL %>% {Rgraphviz::layoutGraph(.)}
  Rgraphviz::layoutGraph(gnell, nodeAttrs=list(label=structure(V(g)$name, names=V(g)$name)))
  if(!is.null(main)) graph::graph.par(list(graph = list(main = main))) # Add a title if one is given
  Rgraphviz::renderGraph(gnell) # Render the graph
  graph::graph.par(list(graph = list(main = ""))) # Reset graph parameters
}
#https://github.com/robertness/lucy/blob/master/R/misc-utilities.R
name_vertices <- function(g){
  if(is.null(V(g)$name)){
    V(g)$name <- paste(V(g))
  }
  g
}

if(F){
# `%T>%` in R
# https://stackoverflow.com/questions/59551512/what-does-t-function-mean-in-r#:~:text=%25T%3E%25%20is%20typicaly%20used,the%20pipeline%20we%20can't.
# The tee pipe, %T>%, is useful when a series of operations have a function that does not return any value.
g_nell <- tempfile() %T>%
    {KEGGgraph::retrieveKGML("04010", organism="hsa", destfile=., quiet=TRUE)} %>%
  {KEGGgraph::parseKGML2Graph(., expandGenes=FALSE)} 
}

tempfile='hsa04010.xml'; #file would be download and saved in this tempfile
tempfile='hsa00340.xml';
#tempfile='dme00340.xml'; #file would be download and saved in this tempfile
if(!file.exists(tempfile)){
  #KEGGgraph::retrieveKGML(pathwayid="04010", organism="hsa", destfile=tempfile,quiet=TRUE)
  #KEGGgraph::retrieveKGML(pathwayid="00340", organism="dme", destfile=tempfile,quiet=TRUE)
  KEGGgraph::retrieveKGML(pathwayid="00340", organism="hsa", destfile=tempfile,quiet=TRUE)
}
g_nell <- KEGGgraph::parseKGML2Graph(tempfile, expandGenes=FALSE)
g_nell

vertex_list <- KEGGgraph::getKEGGnodeData(g_nell) %>%
    {data.frame(
      kegg = unlist(lapply(., function(item) item@name[1])),
      label = unlist(lapply(., function(item)
        strsplit(item@graphics@name, ",")[[1]][1])), stringsAsFactors = F)}
vertex_list #genes 

g_init <- igraph.from.graphNEL(g_nell) 
V(g_init)$name <- vertex_list$kegg 
vertex_list <- dplyr::filter(vertex_list, !duplicated(kegg))

edge_list <- KEGGgraph::getKEGGedgeData(g_nell) %>%
  lapply(function(item){
    if(length(item@subtype) > 0){
      subtype_info <- item@subtype
      # KEGG uses a hierarchy of term for describing terms
      # for example, the first edge type is "activation", the second is "phosphorylation"
      # where phosphorylation is a type of activation.  The second term is more specific than
      # the first, so when it is provided, use it in lieu of the first type.
      if(length(subtype_info) > 1) {
        return(subtype_info[[2]]@name)
      } else {
        return(subtype_info$subtype@name)
      }
    } 
    NA
    }) %>%
    unlist %>%
    {cbind(get.edgelist(g_init), type = .)} %>%
    data.frame   # %>%
    #{dplyr::filter(.,type == "phosphorylation")} #no filter

table(edge_list$type)
edge_list #genes which are connected through compound

edge_list <- edge_list %>%
    as.data.frame %>%
    unique
vertex_list <- vertex_list %>%
    unique %>%
    {dplyr::filter(., !duplicated(kegg))}
g <- graph.data.frame(edge_list, directed = TRUE, vertices = vertex_list)
V(g)$kid <- V(g)$name
V(g)$name <- V(g)$label

g <- g - V(g)[igraph::degree(g) == 0]
igraphviz(g)

