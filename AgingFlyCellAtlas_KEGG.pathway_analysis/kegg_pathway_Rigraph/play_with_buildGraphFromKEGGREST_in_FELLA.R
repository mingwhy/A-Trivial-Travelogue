library(FELLA)
library(igraph)
library(magrittr)
library(tidyverse)
library(org.Dm.eg.db)
library(KEGGREST)

if(F){
  # all pathways, modules 
  all.pathway.graph <- buildGraphFromKEGGREST(organism = "dme");
  all.pathway.graph
  saveRDS(all.pathway.graph,'dme_all.pathway.graph.rds')
}
g=readRDS('dme_all.pathway.graph.rds')
which(V(g)$name=='dme00340')
x<-igraph::neighborhood(g, order = 1, V(g)$name %in% c('dme00340'))
x #pathway.id, enzymes, reactions.
x=igraph::neighborhood(g, order = 1, V(g)$name %in% c('R00069'))
# https://www.genome.jp/entry/R00069
plot.x=induced_subgraph(graph, vids=x[[1]]$name)
plot(plot.x)

x=igraph::neighborhood(g, order = 2, V(g)$name %in% c('C05565'))
x
# https://www.genome.jp/entry/R00069
plot.x=induced_subgraph(graph, vids=x[[1]]$name)
plot(plot.x)
# directed plot, all compounds only point to reactions,reactions point to enzymes and pathways.
# not compound->reaction, reaction -> compound  bipartition hypergraph

# get all kegg pathway names of a species
# https://support.bioconductor.org/p/109871/
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

filter_out_path=setdiff(dme_all_pathways,path_ids)
length(filter_out_path) 
########################################################
# construct graph of a pathway, `buildGraphFromKEGGREST.R` from FELLA
#graph <- buildGraphFromKEGGREST(
# organism = "dme",
#  filter.path = 
#    gsub('path:dme','',filter_out_path)) #filter.path: filter out this path
source('FELLA_Rscripts/buildGraphFromKEGGREST.R')
organism = "dme"
categories <- listCategories()
#categories #"pathway"  "module"   "enzyme"   "reaction" "compound"
categories=c("enzyme","reaction","compound")
# Data from KEGGREST
# 
# List of id-name
message("Building through KEGGREST...")

info.org <- KEGGREST::keggInfo(organism)
info.geneannot <- grep(
  "ncbi-[[:lower:]]+", 
  capture.output(cat(info.org)), 
  value = TRUE)
info.geneannot <- gsub("[[:space:]]", "", info.geneannot)
cat.geneannot <- head(
  intersect(c("ncbi-geneid", "ncbi-proteinid"), info.geneannot), 
  1
)
if (length(cat.geneannot) == 0) 
  stop(
    "Organism ", organism, " does not appear to have either ", 
    "ncbi-geneid or ncbi-proteinid. Please contact ", 
    "FELLA's maintainer.")
message(
  "Available gene annotations: ", 
  paste(info.geneannot, collapse = ", "), 
  ". Using ", cat.geneannot)

list.list <- plyr::llply(
  stats::setNames(categories, categories), 
  function(category) {
    # only pathways are organism-specific now
    # modules are filtered through genes
    if (category %in% c("pathway")) {
      ans <- KEGGREST::keggList(
        database = category, 
        organism = organism)
    } else {
      ans <- KEGGREST::keggList(database = category)
    }
    names(ans) <- sanitise(names(ans), category, organism)
    
    ans
  }, 
  .progress = "text"
)
names(list.list)
head(list.list[[3]])
length(list.list[[3]])

# Map identifiers to category 
map.category <- plyr::ldply(
  list.list, 
  function(categ) 
    data.frame(id = names(categ), stringsAsFactors = FALSE), 
  .id = "category") 
map.category <- stats::setNames(
  as.character(map.category$category), 
  map.category$id)

table(map.category)
#compound   enzyme   module  pathway reaction 
#19131     8077      479      142    11957 

# List of kegg links - essentially our edges
list.link <- plyr::alply(
  expand.grid(
    categories, 
    categories, 
    KEEP.OUT.ATTRS = FALSE, 
    stringsAsFactors = FALSE)[lower.tri(
      matrix(seq_len(length(categories)^2), nrow =length(categories))), ], 
  1, 
  function(row) {
    original <- KEGGREST::keggLink(row[1], row[2])
    df <- data.frame(
      from = sanitise(original, row[1], organism), 
      to = sanitise(names(original), row[2], organism))
    attr(df, "from") <- as.character(row[1])
    attr(df, "to") <- as.character(row[2])
    
    df
  }, 
  .progress = "text"
)
split_labes=attributes(list.link)$split_labels
attributes(list.link) <- NULL
attributes(list.link[[1]])$from;attributes(list.link[[2]])$from

# To mine mapping through enzymes (pathway_gene, module_gene relationship)
m.path_gene <- KEGGREST::keggLink(organism, "pathway") %>% 
  stats::setNames(., sanitise(names(.), "pathway", organism))
m.mod_gene <- KEGGREST::keggLink(organism, "module") %>% 
  stats::setNames(., sanitise(names(.), "module", organism))

# Gene to enzyme relationship
m.gene_enzyme <- KEGGREST::keggLink("enzyme", organism) %>% 
  sanitise(., "enzyme", organism)

# Enzyme to gene, but giving the entrez id rather than KEGG's
# Map kegg to entrez
keggGene2entrez <- KEGGREST::keggConv(cat.geneannot, organism) %>% 
  sanitise(., category = "ncbi") %>%  
  split(., names(.))

# Map kegg enzymes to entrez
m.enzyme_gene <- KEGGREST::keggLink(organism, "enzyme") %>% 
  stats::setNames(., sanitise(names(.), "enzyme", organism)) %>%
  # sanitise(., "gene", organism) %>% 
  split(., names(.), drop = TRUE) %>%
  plyr::llply(
    ., function(r) sort(as.character(unique(keggGene2entrez[r]))))

## begin to build links or edges
if(F){
  # Inferred connections (through genes)
  con.infere <- list(
    infere.con2ec(
      names(list.list$pathway),  #not add pathway to network
      "pathway",  
      m.path_gene, 
      m.gene_enzyme),
    
    infere.con2ec(
      names(list.list$module),  #not add module to network
      "module", 
      m.mod_gene, 
      m.gene_enzyme)
  )
}
# Direct connections, no inference
df.noinfere <- plyr::ldply(
  list.link, 
  function(df.piece) {
    a.from <- attr(df.piece, "from")
    a.to <- attr(df.piece, "to")
    if (a.from == "enzyme" & (a.to %in% c("module", "pathway"))) 
      return(NULL)
    
    return(df.piece)
  }, 
  .id = NULL
)
dim(df.noinfere)

if(F){
  df.infere <- plyr::ldply(
    con.infere, 
    function(df.piece) {
      a.from <- attr(df.piece, "from")
      a.to <- attr(df.piece, "to")
      
      return(df.piece)
    } 
  )
  
  matrix.adjacency <- as.matrix(rbind(
    df.noinfere, 
    df.infere
  ))
}
matrix.adjacency <- as.matrix(df.noinfere)
dim(matrix.adjacency)
head(matrix.adjacency)

message("Done.")

message("Building graph...")
g.raw <- igraph::simplify(
  graph.edgelist(matrix.adjacency, directed = TRUE)
) 

V(g.raw)$com <- match(map.category[V(g.raw)$name], categories) 

# Nodes without a kegg name are either obsolete or inexistent
g.raw <- delete.vertices(
  g.raw, 
  which(is.na(V(g.raw)$com)))

# Enzymes that cannot be inferred should be deleted 
# (not found in desired species!)
g.raw <- delete.vertices(
  g.raw, 
  which((V(g.raw)$com == 3) & !(V(g.raw)$name %in% df.infere$from)))

# Same for the modules that do no belong to the species
# Keep only those that have at least one gene associated 
# 
# Alternative: keggLink("genome", "compound") and pick only 
# those of the organism code
# The downside is that it takes around 90s, but might be 
# a safer option
# Checked in 19/10/2019 and both approaches are equivalent
org.modules <- unique(names(m.mod_gene))
g.raw <- delete.vertices(
  g.raw, 
  which((V(g.raw)$com == 2) & !(V(g.raw)$name %in% org.modules)))

# Order by category and id
g.raw <- permute.vertices(
  g.raw, 
  order(order(V(g.raw)$com, V(g.raw)$name)))

# Weighting the edges
tmp <- get.edges(g.raw, E(g.raw))
E(g.raw)$weight <- abs(V(g.raw)$com[tmp[, 1]] - V(g.raw)$com[tmp[, 2]])

# Keep only reactions in a pathway
# i.e. delete reactions that don't have any 3-weight edge
g.raw <- (setdiff(
  which(V(g.raw)$com == 4), 
  get.edges(g.raw, E(g.raw)[E(g.raw)$weight == 3])[, 1]) %>%
    delete.vertices(graph = g.raw, .))

# Keep only compounds that are reactants/products in these reactions
# i.e. delete compounds that don't have any 1-weight edge
g.raw <- (setdiff(
  which(V(g.raw)$com == 5), 
  get.edges(g.raw, E(g.raw)[E(g.raw)$weight == 1])[, 1]) %>%
    delete.vertices(graph = g.raw, .)) 

# Other filtering (remove nodes?)
if (!is.null(filter.path)) {
  names.path <- V(g.raw)[V(g.raw)$com == 1]$name
  filter.out <- lapply(
    filter.path, 
    function(p) {
      which(grepl(p, names.path))
    })
  names.out <- names.path[unique(unlist(filter.out))]
  message(paste0("Filtering ", length(names.out), " pathways."))
  g.raw <- delete.vertices(g.raw, names.out)
}

g.raw <- largestcc(g.raw)

message("Done.")

message("Pruning graph...")
# CURATE GRAPH
# We start with the graph curation
edges.split <- split(seq_len(ecount(g.raw)), E(g.raw)$weight)

message(paste0("Current weight: 1 out of 4..."))
g.curated <- subgraph.edges(
  graph = g.raw, 
  eids = edges.split[[1]], 
  delete.vertices = FALSE)

for (w in names(edges.split)[-1]) {
  current.w <- as.numeric(w)
  message(paste0("Current weight: ", w, " out of 4..."))
  
  dist.matrix <- distances(g.curated, mode = "out")
  list.edges <- edges.split[[w]]
  
  list.ends <- ends(g.raw, list.edges)
  new.edges <- dist.matrix[list.ends] > E(g.raw)$weight[list.edges]
  
  g.curated <- add.edges(
    graph = g.curated, 
    edges = t(list.ends[new.edges, ]), 
    attr = list(weight = E(g.raw)[list.edges[new.edges]]$weight))
}

# Final edge weights have to be inverted
E(g.curated)$weight <- 1/E(g.curated)$weight

tmp <- list.list
names(tmp) <- NULL
tmp <- unlist(tmp)

# Tried sorting according to number of characters 
# (take shortest names). This is weird, as some names 
# are not very known. I will leave the original order
V(g.curated)$NAME <- strsplit(tmp[V(g.curated)$name], split = "; ")
V(g.curated)$entrez <- m.enzyme_gene[V(g.curated)$name]

comment(g.curated) <- info.org
g.curated$organism <- organism

message("Done.")

keggdata.graph <- g.curated

return(keggdata.graph)

