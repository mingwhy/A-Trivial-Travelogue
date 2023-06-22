#BiocManager::install("pathview")
library(pathview)
data(gse)

data(gse16873.d)
gse16873.d[, 1]
colnames(gse16873.d)[1] # "DCIS_1"
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                    species = "hsa", out.suffix = "gse16873")
