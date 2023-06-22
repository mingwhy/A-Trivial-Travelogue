library(pathview)
source('/Users/ming/Documents/lnR_metabolomics/visualization/ming-kegg.species.code.R')

##################################################
## read in signaficant metabolite data
df=data.table::fread('./tradeoff_candidate.mz.csv')
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
########################################
## get all pathway names
library(igraph)
graph=readRDS('./dem_graph.rds')
nodes=V(graph)
x=unique(nodes[grep('dme',names(nodes))])
all.dme=names(x)
length(all.dme) #131

sig.pathway=data.table::fread('./tradeoff_fella.out_diffusion.csv');
sig.dme=sig.pathway[grep('dme',sig.pathway$KEGG.id),]$KEGG.id


######################################################################################
## play with buildDataFromGraph: https://rdrr.io/bioc/FELLA/src/R/buildDataFromGraph.R
graph=readRDS('dem_graph.rds')
#buildDataFromGraph(keggdata.graph=graph,databaseDir = tmpdir,
#                     internalDir=FALSE,matrices='hypergeom',
#                     normality='diffusion',niter = 50)

keggdata.graph=graph 
E(graph)$weight <- 1/E(graph)$weight

id.pathway <- which(V(graph)$com == 1) #1,2,3,4,5 5 levels network, top1,pathway
id.compound <- which(V(graph)$com == 5)#botom 5, compound
length(id.pathway) #131 pathways
length(id.compound) #4018 compounds

vcount(graph)#10558 nodes
ecount(graph) #34285 edges

x= shortest.paths(
  graph = graph, 
  v = id.compound, 
  to = id.pathway, 
  mode = "out")
dim(x) #4081 x 131
#can one compound 'walk' to a pathway in 4 steps, as compound level5, pathway level1, 5->1, 4 steps

hypergeom.matrix <- shortest.paths(
  graph = graph, 
  v = id.compound, 
  to = id.pathway, 
  mode = "out") == 4

library(Matrix)
hypergeom.matrix <- Matrix(data = hypergeom.matrix, sparse = TRUE)
# save a sparse matrix

sum(hypergeom.matrix[,20]) #167 mz in total
colnames(hypergeom.matrix)[20] #"dme00270"
dim(hypergeom.matrix) #4018compounds X 131 pathways

## intersect between mz.info and hypergeom.matrix
i=intersect(rownames(hypergeom.matrix),mz.info$`KEGG ID`)
hypergeom.matrix.sub=hypergeom.matrix[i,]
dim(mz.info) #332 mz
dim(hypergeom.matrix.sub) #291 x 131

i=Matrix::colSums(hypergeom.matrix.sub)
sum(i==0)
hypergeom.matrix.sub=hypergeom.matrix.sub[,i!=0]
dim(hypergeom.matrix.sub) #291 x 127
# as FELLA 'expand' pathways, the number of mz may be bigger than kegg database

########################################################
##in pathview, it'd automatically download all xml files
cpd.data=rep(1,length(mz.info$`KEGG ID`))
names(cpd.data)=mz.info$`KEGG ID`
cpd.data[mz.id$`KEGG ID`]=-1
#sum(cpd.data==1) #36. -1=blue. 1=yellow
sum(cpd.data==-1) #36. -1=blue. 1=yellow
if(F){
pv.out<-pathview(cpd.data=cpd.data,
                 pathway.id ='dme00270',
                 species = "dme", out.suffix = "test_cpd", 
                 #keys.align='y',
                 kegg.native=T,
                 key.pos='bottomleft',sign.pos='topright',
                 cpd.lab.offset=1)
}    
x=pv.out$plot.data.cpd
unique(x$mol.data)
barplot(1:3,col=unique(x$mol.col))
y=paste(x$mol.data,x$mol.col)
unique(y) #NA:while; assay: grey; sig.mz: yellow

demo.paths
length(sig.dme) #8
dme.i=sig.dme[2] #dme00270
if(F){
for(dme.i in sig.dme){
#for(dme.i in all.dme){
  pv.out<-pathview(cpd.data=cpd.data,
                   #pathway.id ='dme00270',
                   pathway.id =dme.i,
                   species = "dme", out.suffix = paste0(dme.i,'_cpd'), 
                   keys.align='y',
                   #kegg.native=F,
                   same.layer=F,
                   #limit=list(gene=5,cpd=2),
                   #bins=list(gene=5,cpd=2),
                   #na.col='blue',
                   cpd.lab.offset=1,
                   discrete=list(gene=T,cpd=T),
                   key.pos='bottomleft',sign.pos='topright')
}
}

########################################################
## directly read in xml, download through above plotting
library(KEGGgraph)
library(org.Dm.eg.db)
keep.kegg.id=c()
(files=Sys.glob('./png_all_dme/*xml'))
for(sfile in files){
  #sfile='./png_sig_dme/dme00230.xml'
  kegg.id=basename(sfile)
  kegg.id=strsplit(kegg.id,'.xml')[[1]][1]
  mapkKGML = sfile;
  mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
  #mapkG
  mapkpathway <- parseKGML(mapkKGML)
  
  #mapkpathway@nodes$`449`@type
  #mapkpathway@nodes$`449`@name
  cpds=sapply(names(mapkpathway@nodes),function(i){
    if(mapkpathway@nodes[[i]]@type=='compound'){
      mapkpathway@nodes[[i]]@name
    }
  })
  
  cpds=Filter(Negate(is.null), cpds)
  length(cpds) #102 mz
  if(length(cpds)!=0){
    keep.kegg.id=c(keep.kegg.id, kegg.id)
  }
}

keep.kegg.id
length(keep.kegg.id) #102

## only plot for kegg.id which has overlap with the targeted assay
# -1=blue, sig.mz. 1=yellow, present in the targeted assay
for(dme.i in keep.kegg.id){
  #for(dme.i in all.dme){
  pv.out<-pathview(cpd.data=cpd.data,
                   #pathway.id ='dme00270',
                   pathway.id =dme.i,
                   species = "dme", out.suffix = paste0(dme.i,'_cpd'), 
                   keys.align='y',
                   #kegg.native=F,
                   same.layer=F,
                   #limit=list(gene=5,cpd=2),
                   #bins=list(gene=5,cpd=2),
                   #na.col='blue',
                   cpd.lab.offset=1,
                   discrete=list(gene=T,cpd=T),
                   key.pos='bottomleft',sign.pos='topright')
}

