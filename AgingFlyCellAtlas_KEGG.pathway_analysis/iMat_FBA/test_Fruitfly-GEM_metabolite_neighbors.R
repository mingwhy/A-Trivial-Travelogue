
library(ggplot2);library(gridExtra)

##############################################
## load FruitMat GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
gem=readMat('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')

ifly=gem[[1]][,,1]

length(ifly$subSystems)==ncol(ifly$S)
ifly$subSystems=unlist(ifly$subSystems)
length(table(ifly$subSystems)) #149

#not do this: ifly$rxnNames=unlist(ifly$rxnNames) #unlist would drop NULL automatically
ifly$rxnNames=as.character(ifly$rxnNames)
length(unlist(ifly$rxns))==ncol(ifly$S) #11898
ifly$rxns=unlist(ifly$rxns)

length(unlist(ifly$mets))==nrow(ifly$S)
ifly$mets=unlist(ifly$mets) #all members non-empty, checked length
ifly$metNames=unlist(ifly$metNames)

length(unlist(ifly$grRules))==ncol(ifly$S) #FALSE
ifly$grRules=as.character(ifly$grRules)

# compartment short names and full names
as.character(ifly$comps)
as.character(ifly$compNames)
#c:Cytosol, m:Mitochondria, e:Extracellular

x=substr(ifly$mets,9,9)
table(x)
y=ifly$subSystems[x=='e']
table(y)
## test extracting one metabolite
ifly$subSystems[grep('hist',ifly$subSystems,ignore.case = T)]

col.idx=grep('hist',ifly$subSystems,ignore.case = T);
tmp=ifly$S[,grep('hist',ifly$subSystems,ignore.case = T)]
row.idx=which(Matrix::rowSums(tmp)!=0)

ifly$rxns[col.idx]
ifly$rxnNames[col.idx]
#ifly$metFormulas[col.idx]
ifly$grRules[col.idx]
ifly$metNames[row.idx]

tmp=ifly$S[row.idx,col.idx]
rowSums(tmp)
rowSums(ifly$S[row.idx,])
tmp=as.matrix(tmp)
rownames(tmp)=ifly$metNames[row.idx]
colnames(tmp)=ifly$rxnNames[col.idx]
tmp

#############################################
# id Cross references
#kegg id C01672: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02123e
#model folder in: https://github.com/SysBioChalmers/Fruitfly-GEM
mets.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/metabolites.tsv')
head(mets.mapping)
mets.mapping[mets.mapping$metKEGGID=='C00388',]#Cytosol,Extracellular #ifly$compNames
tmp=ifly$S[ifly$mets %in% mets.mapping[mets.mapping$metKEGGID=='C00388',]$mets,]
col.index=which( colSums(abs(tmp))!=0 ) 
tmp[,col.index] #all rxn which involve `C00388`, not necessariliy histine metabolism
ifly$rxns[col.index] #"MAR04428" "MAR00619"
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124c
#https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/metabolite/MAM02124e
ifly$subSystems[col.index]

rxns.mapping=data.table::fread('~/Documents/bioinfo_software/GEM_metabolic.models/Fruitfly-GEM-main/model/reactions.tsv')
rxns.mapping[rxns.mapping$rxns=="MAR01442",]

########################################################
## gene expression data
dataset.name='fly.female'
#dat=data.table::fread('log1p_female_gene.mean.expr.csv')
dat=readRDS('log1p_female_gene.mean.expr.rds')
length(dat) #27
length(dat[[1]]) #4
names(dat[[1]])
# tissue;cell.type;age, 4 age groups
head(dat[[1]][[1]])
gene.expr=dat;

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
dim(mz.age.betas) #86

share.mz=intersect(mz.age.betas$KEGGid, mets.mapping$metKEGGID)
share.mz=share.mz[share.mz!='']
length(share.mz) #61

mets.mapping[mets.mapping$metKEGGID=='C00388',] #Histamine
for(mz in share.mz){
  mz='C00388';
  cat(mz,'\n')
  # might be multiple ids due to different compartments
  ids=mets.mapping[mets.mapping$metKEGGID==mz,]$mets
  ids
  #mz.age.betas[mz.age.betas$KEGGid==mz,]
  
  bait=match(ids,ifly$mets)
  tmp=ifly$S[bait,]
  col.idx=which(Matrix::colSums(abs(tmp))!=0) #rxn target mz is involved
  #ifly$subSystems[col.idx] #there might be  "Exchange/demand reactions" and "Transport reactions"
  x=ifly$grRules[col.idx]
  col.idx=col.idx[!duplicated(x)] #rxn via the same set of genes, only keep one
  x=ifly$grRules[col.idx]
  col.idx=col.idx[x !="list(character(0))"]
  if(length(col.idx)==0){cat('Genes-associated-rxns not found\n');next}
  enzyme.genes=ifly$grRules[col.idx]
  
  tmp=ifly$S[,col.idx]
  row.idx=which(Matrix::rowSums(abs(tmp))!=0) #metabolite these traget rxn requires
  tmp=ifly$S[row.idx,col.idx]
  
  tmp=as.matrix(tmp)
  #rownames(tmp)=ifly$metNames[row.idx]
  rownames(tmp)=ifly$mets[row.idx]
  colnames(tmp)=ifly$rxnNames[col.idx]
  
  #combine ids to one metabolite and ignore stoi number, only keep 1, -1,0
  tmp1=tmp;
  
  target.mz.row.ids<-match(ids,rownames(tmp))
  
}
