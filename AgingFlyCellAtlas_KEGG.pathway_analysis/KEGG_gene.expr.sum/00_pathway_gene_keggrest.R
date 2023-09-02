
################################################################
#https://stackoverflow.com/questions/28724674/does-anyone-know-how-to-retrieve-list-of-cell-cycle-genes-from-kegg-in-r
## pathway id, name, gene members
library(KEGGREST)
library(org.Dm.eg.db)
library(tidyverse)
# get pathways and their entrez gene ids
dme_path_entrez  <- keggLink("pathway", "dme") %>% tibble(pathway = ., eg = sub("dme:Dmel_", "", names(.)))
head(dme_path_entrez)
dim(dme_path_entrez) #7321

# get gene symbols and ensembl ids using entrez gene ids
dme_kegg_anno <- dme_path_entrez %>%
  mutate(
    symbol = mapIds(org.Dm.eg.db, eg, "SYMBOL", "FLYBASECG"),
    ensembl = mapIds(org.Dm.eg.db, eg, "ENSEMBL", "FLYBASECG")
  )
dim(dme_kegg_anno)#7321    4

# Pathway names
dme_pathways <- keggList("pathway", "dme") %>% tibble(pathway = names(.), description = .)
head(dme_kegg_anno)
head(dme_pathways)
dme_kegg_anno$pathway=gsub('path:','',dme_kegg_anno$pathway)
KEGG_pathways <- left_join(dme_kegg_anno, dme_pathways)
dim(KEGG_pathways) #7321
head(KEGG_pathways)
x=(KEGG_pathways %>% group_by(pathway,description) %>% summarise(ngene=n()))
x[order(x$ngene),]
saveRDS(KEGG_pathways,'kegg-flygenes.rds')
