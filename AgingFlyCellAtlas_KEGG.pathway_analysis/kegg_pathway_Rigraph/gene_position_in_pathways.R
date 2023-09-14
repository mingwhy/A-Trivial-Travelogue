#How upstream is a gene in any of it's pathway ?
#https://www.biostars.org/p/470973/
library(tidyverse)
library(KEGGREST)

kegg_id_names = c('C00388')

path_ids = map(kegg_id_names, ~keggLink('pathway', .)) %>% 
  unlist() %>%
  unique()
path_ids

##########
library(tidyverse)
library(KEGGREST)

kegg_id_names = c('hsa:7040')

path_ids = map(kegg_id_names, ~keggLink('pathway', .)) %>% 
  unlist() %>%
  unique()

path_data = keggGet(path_ids[1:3])

get_kegg_position = function(path_record) {
  
  gene_ids = path_record$GENE[seq(1,length(path_record$GENE),2)]
  gene_names = gsub(";.*","",path_record$GENE[seq(2,length(path_record$GENE),2)])
  
  if (is.null(gene_ids)) {
    res = tibble(gene = NA,
                 position = NA)
  }else{
    gene_posns = seq(1, 0, length.out = length(gene_ids))
    
    res = tibble(gene_id = gene_ids,
                 gene_name = gene_names,
                 position = gene_posns) %>% 
      separate_rows('gene_id', sep = ',')
  }
  
  res %>% mutate(path_id = path_record$ENTRY,
                 path_name = path_record$NAME)
  
}

kegg_positions_by_module = map_df(path_data,
                                  get_kegg_position) %>% 
  filter(!is.na(position))

# TGFB1 score in MAPK pathway (should be close to 1)
kegg_positions_by_module %>%
  filter(path_id=="hsa04010") %>% filter(gene_name == "TGFB1")
# SRF score in MAPK pathway (should be lower than TGFB1 as it is only upstream of 1 gene (FOS)) 
kegg_positions_by_module %>%
  filter(path_id=="hsa04010") %>% filter(gene_name == "SRF")
# MAX score in MAPK pathway (should be close to 0)
kegg_positions_by_module %>%
  filter(path_id=="hsa04010") %>% filter(gene_name == "MAX")
######################################################################
GetKEGGigraph <- function(KEGG_pathway_id, plot = FALSE, plot_file=NULL){
  stopifnot(is.character(KEGG_pathway_id))
  
  # retrieve pathway thanks to KEGGgraph
  tmp <- tempfile()
  res = KEGGgraph::retrieveKGML(KEGG_pathway_id, organism="hsa", destfile=tmp, method="wget", quiet=TRUE)
  mapkG <- KEGGgraph::parseKGML2Graph(res,expandGenes=TRUE, genesOnly = TRUE)
  
  outs <- sapply(KEGGgraph::edges(mapkG), length) > 0
  ins <- sapply(KEGGgraph::inEdges(mapkG), length) > 0
  ios <- outs | ins
  ## translate the KEGG IDs into Gene Symbol
  if(require(org.Hs.eg.db)) {
    ioGeneID <- KEGGgraph::translateKEGGID2GeneID(names(ios))
    nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  } else {
    nodesNames <- names(ios)
  }
  names(nodesNames) <- names(ios)
  
  mapkG_igraph = igraph::igraph.from.graphNEL(mapkG, name = TRUE, weight = TRUE,
                                              unlist.attrs = TRUE)
  mapkG_igraph = igraph::simplify(mapkG_igraph, remove.multiple = TRUE, remove.loops = TRUE,
                                  edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
  Isolated = which(igraph::degree(mapkG_igraph)==0)
  mapkG_igraph = igraph::delete.vertices(mapkG_igraph, Isolated)
  
  # Minimum spanning tree graph from pathway
  mstree = igraph::mst(mapkG_igraph)
  V(mstree)$id <- seq_len(vcount(mstree))-1
  roots <- sapply(igraph::decompose(mstree), function(x) {
    V(x)$id[ igraph::topo_sort(x)[1]+1 ] })
  
  if(plot & !is.null(plot_file)){
    # Change names to gene name for graph
    mapkG@nodes = nodesNames
    names(mapkG@edgeL) = nodesNames
    
    mapkG_igraph_gene = igraph::igraph.from.graphNEL(mapkG, name = TRUE, weight = TRUE,
                                                     unlist.attrs = TRUE)
    mapkG_igraph_gene = igraph::simplify(mapkG_igraph_gene, remove.multiple = TRUE, remove.loops = TRUE,
                                         edge.attr.comb = igraph::igraph_opt("edge.attr.comb"))
    Isolated = which(igraph::degree(mapkG_igraph_gene)==0)
    mapkG_igraph_gene = igraph::delete.vertices(mapkG_igraph_gene, Isolated)
    
    # Minimum spanning tree graph from pathway
    mstree_gene = igraph::mst(mapkG_igraph_gene)
    pdf(file.path(plot_file))
    plot(mstree_gene, layout = igraph::layout_nicely(mstree_gene),
         vertex.color= ifelse(,"red","grey"), vertex.size = 3.75,
         vertex.label.cex=0.25, edge.arrow.width=0.25, edge.arrow.size=0.25, edge.width=0.5)
    dev.off()
    
  }
  
  return(mstree)
}

upstream_score_KEGG <- function(gene_list){
  ##Get the Entrez gene IDs associated with those symbols
  EG_IDs = mget(gene_list, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
  
  ##Then get the KEGG IDs associated with those entrez genes.
  KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA)
  
  results <- foreach::foreach(KEGG_id = names(KEGG_IDs), .combine=rbind,
                              .packages=c('KEGGREST',"KEGGgraph",'igraph')) %dopar% 
    {
      # results=data.frame("upstream_score"=0,"KEGG_id"="")
      # for(KEGG_id in names(KEGG_IDs)) {
      KEGG_pathway_id = KEGG_IDs[[KEGG_id]]
      print(KEGG_id)
      if(is.na(KEGG_pathway_id[1])) {
        ret = data.frame("upstream_score"=NA, "KEGG_id"=KEGG_id)
        return(ret)
        # results = rbind(results,ret)
      } else {
        list_topo_sorted <- lapply(KEGG_pathway_id, function(id){
          mstree_sorted <- data.frame("upstream_score"=0,"KEGG_id"="")
          try({
            mstree <- GetKEGGigraph(id, plot=FALSE)
            mstree_sorted <- igraph::as_ids(igraph::topo_sort(mstree))
          }, TRUE)
          
          return(mstree_sorted)
        })
        names(list_topo_sorted) <- paste0(KEGG_id,"_",KEGG_pathway_id)
        list_topo_sorted = lapply(list_topo_sorted, function(x){
          n = 1 - (which(x==paste0("hsa:",KEGG_id))/length(x))
          if(length(n)==0) return(NA) else return(n[1])
        } )
        df_topo_sorted = as.data.frame(t(as.data.frame(list_topo_sorted,drop=F)))
        df_topo_sorted$KEGG_id = rep(KEGG_id,nrow(df_topo_sorted))
        colnames(df_topo_sorted)[1] = "upstream_score"
        if(!is.null(df_topo_sorted)) return(df_topo_sorted)
        # if(!is.null(df_topo_sorted)) results=rbind(results,df_topo_sorted)
      }
      
    }
  results$KEGG_pathway = rownames(results)
  results$KEGG_pathway[grep("_",results$KEGG_pathway,invert = T)] = NA
  results$KEGG_pathway = gsub(".*_","",results$KEGG_pathway)
  
  results$Gene = sapply(mget(results$KEGG_id, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  results$Pathway = ""
  
  l = sapply(paste0("path:hsa",results$KEGG_pathway[which(!is.na(results$KEGG_pathway))]),function(x){
    print(x)
    ret = ""
    try({
      query = keggGet(x)
      ret = query[[1]]$PATHWAY_MAP
    }, TRUE)
    return(ret)
  })
  results$Pathway[which(!is.na(results$KEGG_pathway))] = l
  return(results)
}

# Get "upstream scores" for a list of gene
library(igraph)
library(doParallel)
library(tidyverse)
library(org.Hs.eg.db)
library(KEGGREST)
library(KEGGgraph)

library(foreach)
# Choose number of cores to run parallely (takes time)
registerDoParallel(cores=6)

gene_list <- c("TGFB1", "MAX")
scores <- upstream_score_KEGG(gene_list)
print(scores)

#        upstream_score KEGG_id KEGG_pathway  Gene                                      Pathway
#X7040_04010      0.7312925    7040        04010 TGFB1                       MAPK signaling pathway
#X7040_04060      0.9438596    7040        04060 TGFB1       Cytokine-cytokine receptor interaction
#X7040_04110      0.7096774    7040        04110 TGFB1                                   Cell cycle
#X7040_04144             NA    7040        04144 TGFB1                                  Endocytosis
#X7040_04350      0.3152174    7040        04350 TGFB1                   TGF-beta signaling pathway
#                                                                        ......
#X4149_04010      0.1870748    4149        04010   MAX                       MAPK signaling pathway
#X4149_05200      0.8090129    4149        05200   MAX                           Pathways in cancer
#X4149_05222      0.9534884    4149        05222   MAX                       Small cell lung cancer


