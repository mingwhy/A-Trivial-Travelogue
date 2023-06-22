library(FELLA)
library(org.Dm.eg.db)
library(KEGGREST)
library(igraph)
library(magrittr)

mode='tradeoff'
#mode='insulin'
#mode='fecundity'

##################################################
## read in signaficant metabolite data
input.file=paste0(mode,'_candidate.mz.csv')
outfile1=paste0(mode,'_fella.out_diffusion.csv')
outfile2=paste0(mode,'_fella.out_hypergeom.csv')
#df=data.table::fread('tradeoff_candidate.mz.csv')
df=data.table::fread(input.file)
mz=df$metabolite

## read in assay background
library(readxl)
mz.info=read_excel('./2020-12-08_Tatar-60_Data.xlsx',sheet='Metabolite Information')
dim(mz.info) #361
sum(mz.info$`KEGG ID`=='N/A') #29 measure metabolites do not have kegg id
mz.info=mz.info[mz.info$`KEGG ID`!='N/A',]
dim(mz.info) #332
tmp=mz.info[grep('\\/',mz.info$`KEGG ID`),] #1 with multiple id
tmp #8
x=strsplit(tmp$`KEGG ID`,'\\/')
x
y=sapply(x,function(i){ i[i!='NA'][1] } )
mz.info[grep('\\/',mz.info$`KEGG ID`),]$`KEGG ID`=y
dim(mz.info)#332   6

## map sig mz to bacground
length(mz) #39
sum(mz %in% mz.info$`Current MS Compounds`) #36
mz.id=mz.info[mz.info$`Current MS Compounds` %in% mz,]
dim(mz.id) #36

##################################################
# build the KEGG graph
#https://rdrr.io/bioc/FELLA/src/R/buildDataFromGraph.R
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
                     normality='diffusion',niter = 100)
}  
tmpdir<-'FELLA_database_hyper'
if(!dir.exists(tmpdir)){
  buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
                     internalDir=FALSE,matrices='hypergeom',
                     normality='diffusion',niter = 100)
}

########################################
## diffusion method

tmpdir<-'FELLA_database_diff'

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "diffusion")
fella.data
#Categories:
#+ pathway [131]
#+ module [165]
#+ enzyme [750]
#+ reaction [5494]
#+ compound [4018]

cat(getInfo(fella.data))

# map metabolites to the KEGG graph
analysis.mz <- defineCompounds(
   compounds = mz.id$`KEGG ID`,
   compoundsBackground=mz.info$`KEGG ID`, #use targeted assay as the background
   data = fella.data)

length(getInput(analysis.mz)) #33
length(getExcluded(analysis.mz)) #3
length(analysis.mz@userinput@metabolitesbackground) #291
analysis.mz

## enriching using diffusion, approx='normality' is a fast option
# enrich is a wrapper function of defineCompounds, runDiffusion, etc.
# myAnalysis <- enrich(compounds = mz.id$`KEGG ID`,method='diffusion',
#  data = fella.data,approx = "normality")

set.seed(13579) #make sure the result is reproducible
analysis.mz <- runDiffusion(
   object = analysis.mz,
   data = fella.data,
   niter = 10000, #a small number would give non-stable results
   #approx = "normality")
    approx = "simulation")
analysis.mz

# empirical p-value, or p-score
#https://rdrr.io/bioc/FELLA/src/R/generateResultsTable.R
#https://support.bioconductor.org/p/9140352/
#https://support.bioconductor.org/p/9140352/#9140355
length(analysis.mz@diffusion@pscores) #10558
summary(analysis.mz@diffusion@pscores)
sum(analysis.mz@diffusion@pscores<0.05) #740 
sum(analysis.mz@diffusion@pscores<0.0001) #8

length(grep('dme',names(analysis.mz@diffusion@pscores))) #131 dme pathway
length(grep('^C',names(analysis.mz@diffusion@pscores))) #4018 Compounds
length(grep('^M',names(analysis.mz@diffusion@pscores))) #165 Module
length(grep('^R',names(analysis.mz@diffusion@pscores))) #5494 Reactions
length(grep('^\\d',names(analysis.mz@diffusion@pscores))) #750 Enzymes

x=analysis.mz@diffusion@pscores[grep('dme',names(analysis.mz@diffusion@pscores))]
length(x) #131 pathway
dme.p.adjust=p.adjust(x,method='BH')
min(dme.p.adjust) #0.63

x=analysis.mz@diffusion@pscores[grep('^C',names(analysis.mz@diffusion@pscores))]
length(x) #4018 compounds
c.p.adjust=p.adjust(x,method='BH')
min(c.p.adjust) #0.37

x=analysis.mz@diffusion@pscores[grep('^R',names(analysis.mz@diffusion@pscores))]
length(x) #5494 
p.adjust=p.adjust(x,method='BH')
min(p.adjust) #0.47

x=analysis.mz@diffusion@pscores[grep('^M',names(analysis.mz@diffusion@pscores))]
length(x) #165 modules
p.adjust=p.adjust(x,method='BH')
min(p.adjust) #0.94


x=analysis.mz@diffusion@pscores[grep('^\\d',names(analysis.mz@diffusion@pscores))]
length(x) #750
p.adjust=p.adjust(x,method='BH')
min(p.adjust) #0.55


exportResults(
  format = "csv",
  threshold=0.05,nlimit = 800, #top 300 & pscore<0.05
  #file = 'fella_out_diffusion.csv',
  file = outfile1,
  method = "diffusion",
  object = analysis.mz,
  data = fella.data)


if(F){
  #plot(analysis.mz,method='diffusion',data=fella.data,
  #     main='analysis',#threshold=0.0001,
  #     nlimit=90, #bigger number leads to more nodes. green shows your input mz
  #     plotLegend = T)
  
  #length(analysis.mz@diffusion@pscores) #10558
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
}



#########################################################
## runHypergeom

tmpdir<-'FELLA_database_hyper/'

fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix='hypergeom')

fella.data
cat(getInfo(fella.data))
fella.data@hypergeom #big matrix
colnames(fella.data@hypergeom@matrix)[20] #"dme00270"
sum(fella.data@hypergeom@matrix[,20]) #167 for "dme00270"

# map metabolites to the KEGG graph
dim(mz.info) #332
analysis.mz <- defineCompounds(
  compounds = mz.id$`KEGG ID`,
  compoundsBackground=mz.info$`KEGG ID`, #use targeted assay as the background
  data = fella.data)

length(getInput(analysis.mz)) #33
length(getExcluded(analysis.mz)) #3
analysis.mz
length(getBackground(analysis.mz)) #291
length(analysis.mz@userinput@metabolitesbackground) #291

analysis.mz <- runHypergeom(
  object = analysis.mz,
  data = fella.data,
  p.adjust = 'fdr')
analysis.mz

summary(analysis.mz@hypergeom@pvalues)

if(F){
  # top pathways and affected compounds
  plot(analysis.mz,method='hypergeom',data=fella.data,
       main='analysis',
       plotLegend = FALSE)
  
  # top k p.scores
  tab.all <- generateResultsTable(
    method = "hypergeom",
    object = analysis.mz,
    data = fella.data)
  dim(tab.all) #
  
  myGraph <- generateResultsGraph(
    object = analysis.mz, 
    method = "hypergeom", 
    data = fella.data)
  show(myGraph)

  # details from the enzymes reported among the top k KEGG entries
  # not applicable to method = "hypergeom" 
  #tab.enzyme <- generateEnzymesTable(method = "hypergeom",object = analysis.mz,data = fella.data)
}

exportResults(
  format = "csv",
  #file = 'fella_out_hypergeom.csv',
  file = outfile2,
  method = "hypergeom",
  object = analysis.mz,
  data = fella.data)

#dme00480  Glutathione metabolism
#dme00270  Cysteine and methionine metabolism

###################################################################
## extract mapped mz to each sig pathway
sig.path=data.table::fread('./fella_out_hypergeom.csv')
sig.path$KEGG.id

## get KEGG database
getMatrix <- function(data, method) {
  return(slot(data, method)@matrix)
}
hypergeom.matrix=getMatrix(fella.data,'hypergeom')

object=analysis.mz
metabolites.input <- getInput(object)
metabolites.input.intersect <- intersect(
  metabolites.input, 
  rownames(hypergeom.matrix))

row_comp <- which(rownames(hypergeom.matrix) %in% metabolites.input)
length(row_comp) #33 mz hits

col_comp=which(colnames(hypergeom.matrix) %in% sig.path$KEGG.id)
out=list();
for(path in col_comp){
  #hypergeom.matrix[row_comp, path]
  id=names(which(hypergeom.matrix[row_comp, path])) #matched mz hits in this pathway
  x=mz.id[mz.id$`KEGG ID` %in% id,]
  
  kegg.id=colnames(hypergeom.matrix)[path]
  kegg.name=sig.path[sig.path$KEGG.id==kegg.id,]$KEGG.name
  x$mapped.pathway.id=kegg.id
  x$mapped.pathway.name=kegg.name
  out[[kegg.id]]=x
}

df.out=Reduce(`rbind`,out)
head(df.out)
table(df.out$mapped.pathway.id)
data.table::fwrite(df.out,'sig.mz_to_pathway.csv')

###################################################################
###################################################################
## play with FELLA: https://rdrr.io/bioc/FELLA/src/R/runHypergeom.R
object=analysis.mz
data=fella.data
p.adjust='fdr'
#https://rdrr.io/bioc/FELLA/src/R/get-.R
getMatrix <- function(data, method) {
return(slot(data, method)@matrix)
}
hypergeom.matrix=getMatrix(data,'hypergeom')

which(colnames(hypergeom.matrix)=='dme00270')
which(colnames(hypergeom.matrix)=='dme00592') 
hypergeom.matrix[,20]
sum(hypergeom.matrix[,20]) #167 mz in total

metabolites.input <- getInput(object)
metabolites.input.intersect <- intersect(
  metabolites.input, 
  rownames(hypergeom.matrix))
length(metabolites.input.intersect) #33
# I didn't specify bacground, so it's 0
metabolites.background.intersect <- intersect(
  getBackground(object), 
  rownames(hypergeom.matrix))
length(metabolites.background.intersect) #0

#   hypergeom.matrix <- hypergeom.matrix[, colSums(hypergeom.matrix) > 0]
dim(hypergeom.matrix) #4018  131, mz by pathways
row_comp <- which(rownames(hypergeom.matrix) %in% metabolites.input)
length(row_comp) #33 mz hits

pvalues.path <- vector("double", length = dim(hypergeom.matrix)[2])
names(pvalues.path) <- colnames(hypergeom.matrix)
length(pvalues.path) #131 pathways

# Additional variables for the report
pathbackground <- pvalues.path
pathhits <- pvalues.path

# p-values calculation
colnames(hypergeom.matrix)[20] #"dme00270"
path=20
#for (path in seq_len(ncol(hypergeom.matrix))) {
  hypergeom.matrix[row_comp, path]
  names(which(hypergeom.matrix[row_comp, path])) #matched mz hits in this pathway
  
  sample_success <- sum(hypergeom.matrix[row_comp, path])
  pathhits[path] <- sample_success
  
  # sum over row or column (current)?
  total_success <- sum(hypergeom.matrix[, path])
  pathbackground[path] <- total_success
  
  total_failure <- dim(hypergeom.matrix)[1] - total_success
  sample_size <- length(row_comp)
  
  pvalues.path[path] <- stats::phyper(
    sample_success - 1, 
    total_success, 
    total_failure, 
    sample_size, 
    lower.tail = FALSE)
#}

######################################################################################
## play with buildDataFromGraph: https://rdrr.io/bioc/FELLA/src/R/buildDataFromGraph.R
graph=readRDS('dem_graph.rds')
#buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
#                     internalDir=FALSE,matrices='hypergeom',
#                     normality='diffusion',niter = 50)
  
keggdata.graph=graph 
E(graph)$weight <- 1/E(graph)$weight

  