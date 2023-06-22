library(FELLA)
library(org.Dm.eg.db)
library(KEGGREST)
library(igraph)
library(magrittr)


##################################################
## read in metabolite data
df=data.table::fread('candidate.mz.txt')
mz=df$metabolite
library(readxl)
mz.info=read_excel('./2020-12-08_Tatar-60_Data.xlsx',sheet='Metabolite Information')
length(mz) #85
sum(mz %in% mz.info$`Current MS Compounds`) #85
mz.id=mz.info[mz.info$`Current MS Compounds` %in% mz,]
sum(mz.id$`KEGG ID`=='N/A') #3
mz.id=mz.id[mz.id$`KEGG ID`!='N/A',]

tmp=mz.id[grep('\\/',mz.id$`KEGG ID`),] #3 with multiple id
x=strsplit(tmp$`KEGG ID`,'\\/')
mz.id[grep('\\/',mz.id$`KEGG ID`),]$`KEGG ID`=c(x[[1]][1],x[[2]][[1]],x[[3]][[2]])

dim(mz.id) #82

##################################################
# build the KEGG graph
if(!file.exists('dem_graph.rds')){
  graph<-buildGraphFromKEGGREST(organism='dme')
  saveRDS(graph,'dem_graph.rds')
}
graph=readRDS('dem_graph.rds')

tmpdir<-'FELLA_database_diff'
if(!dir.exists(tmpdir)){
  unlink(tmpdir,recursive=T)
  buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
                     internalDir=FALSE,matrices='diffusion',
                     normality='diffusion',niter = 50)
}  
tmpdir<-'FELLA_database_hyper'
if(!dir.exists(tmpdir)){
  buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
                     internalDir=FALSE,matrices='hypergeom',
                     normality='diffusion',niter = 50)
}

########################################
## diffusion method
tmpdir<-'FELLA_database_diff'

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "diffusion")
fella.data
cat(getInfo(fella.data))

# map metabolites to the KEGG graph
analysis.mz <- defineCompounds(
   compounds = mz.id$`KEGG ID`,
   data = fella.data)

length(getInput(analysis.mz)) #74
length(getExcluded(analysis.mz)) #8
analysis.mz

## enriching using diffusion, approx='normality' is a fast option
# enrich is a wrapper function of defineCompounds, runDiffusion, etc.
# myAnalysis <- enrich(compounds = mz.id$`KEGG ID`,method='diffusion',
#  data = fella.data,approx = "normality")

analysis.mz <- runDiffusion(
   object = analysis.mz,
   data = fella.data,
   #approx = "normality")
    approx = "simulation")
analysis.mz

plot(analysis.mz,method='diffusion',data=fella.data,
     main='analysis',
     #threshold=0.0001,
     nlimit=90, #bigger number leads to more nodes. green shows your input mz
     plotLegend = T)

# top k p.scores
tab.all <- generateResultsTable(
   method = "diffusion",
   nlimit = 90,
   object = analysis.mz,
   data = fella.data)
dim(tab.all) #90 rows

myGraph <- generateResultsGraph(
  object = analysis.mz, 
  method = "diffusion", 
  #threshold = 0.1, 
  nlimit = 100,
  data = fella.data)
show(myGraph)

# details from the enzymes reported among the top k KEGG entries
tab.enzyme <- generateEnzymesTable(
  method = "diffusion",
  nlimit = 100,
  object = analysis.mz,
  data = fella.data)

exportResults(
   format = "csv",
   file = 'fella_out_diffusion.csv',
   method = "diffusion",
   object = analysis.mz,
   data = fella.data)

# empirical p-value, or p-score
length(analysis.mz@diffusion@pscores)
summary(analysis.mz@diffusion@pscores)
sum(analysis.mz@diffusion@pscores<0.0001)
#########################################################
## runHypergeom

tmpdir<-'FELLA_database_hyper/'

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "diffusion")
fella.data
cat(getInfo(fella.data))

# map metabolites to the KEGG graph
analysis.mz <- defineCompounds(
  compounds = mz.id$`KEGG ID`,
  data = fella.data)

length(getInput(analysis.mz)) #74
length(getExcluded(analysis.mz)) #8
analysis.mz

analysis.mz <- runHypergeom(
  object = analysis.mz,
  data = fella.data,
  p.adjust = 'fdr')
analysis.mz
# top pathways and affected compounds
plot(analysis.mz,method='hypergeom',data=fella.data,
     main='analysis',
     plotLegend = FALSE)

# top k p.scores
tab.all <- generateResultsTable(
  method = "hypergeom",
  object = analysis.mz,
  data = fella.data)
dim(tab.all) #100 rows

myGraph <- generateResultsGraph(
  object = analysis.mz, 
  method = "hypergeom", 
  data = fella.data)
show(myGraph)

# details from the enzymes reported among the top k KEGG entries
# not applicable to method = "hypergeom" 
#tab.enzyme <- generateEnzymesTable(method = "hypergeom",object = analysis.mz,data = fella.data)

exportResults(
  format = "csv",
  file = 'fella_out_hypergeom.csv',
  method = "hypergeom",
  object = analysis.mz,
  data = fella.data)
