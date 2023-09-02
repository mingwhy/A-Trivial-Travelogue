# GEM intro: https://www.gu.se/sites/default/files/2021-09/WangH_GOTBIN_seminar_20210928.pdf
# download Fruitfly-GEM from https://github.com/SysBioChalmers/Fruitfly-GEM
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
x=readMat('~/Documents/aging_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
x
ifly=x[[1]][,,1]
ifly$S[1:3,1:3] #sparse matrix
ifly$rxnGeneMat[1:3,1:3] #sparse matrix
dim(ifly$S) #8135 mets X 11898 rxn
dim(ifly$rxnGeneMat) #11898 rxn X 1753 genes
names(ifly)
genes=unlist(ifly$genes)
length(genes) #1753 genes
length(ifly$subSystems)#11898
length(ifly$grRules) #11898
length(ifly$rxnNames) #11898
subSystems=unlist(ifly$subSystems)
table(unlist(lapply(ifly$subSystems,length)))
length(subSystems) #11899, each reaction belong to one subSystem or pathway
length(table(subSystems)) #149 subSystems or pathways

# showcase: https://metabolicatlas.org/explore/Fruitfly-GEM/gem-browser/gene/LpR2
which(genes=='LpR2') #881
which(ifly$rxnGeneMat[,881]!=0) #4502 4503

rxn.ids<-unlist(ifly$rxns)
rxn.ids[which(ifly$rxnGeneMat[,881]!=0)]
rxn.names<-unlist(ifly$rxnNames)
rxn.names[which(ifly$rxnGeneMat[,881]!=0)]
subSystems[which(ifly$rxnGeneMat[,881]!=0)]

#
reactions=data.table::fread('~/Documents/aging_GEM_metabolic.models/Fruitfly-GEM-main/model/reactions.tsv')
metabolites=data.table::fread('~/Documents/aging_GEM_metabolic.models/Fruitfly-GEM-main/model/metabolites.tsv')
dim(reactions) #11898    15
dim(metabolites) # 8135   14
human.fly.orthologs=data.table::fread('~/Documents/aging_GEM_metabolic.models/Fruitfly-GEM-main/data/human2FruitflyOrthologs.tsv')
dim(human.fly.orthologs) # 15636     8
length(unique(human.fly.orthologs$toSymbol)) #8136

length(genes) #1753
sum(genes %in% unique(human.fly.orthologs$toSymbol)) #1665
test.genes=genes[!genes %in% unique(human.fly.orthologs$toSymbol)]
head(test.genes)
# test 'AK-3': https://flybase.org/reports/FBgn0042094.html
grep('FBgn0042094',human.fly.orthologs$toGeneId)
human.fly.orthologs[grep('FBgn0042094',human.fly.orthologs$toGeneId),]

# compare with KEGG 
x=readRDS('/Users/mingyang/Documents/git_bioinfo_homemade_tools/dataBase/KEGG.decompose/kegg-flygenes.rds')
length(x) #137 pathways (by Jul 7, 2021) https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/dataBase/KEGG.decompose
kegg.genes=unique(unlist(lapply(x,'[[',2)))
length(kegg.genes) #3263 fly genes

#########################################################
# extract genes based on mz-rxn-gene associations
library(R.matlab) #in matlab `load(Fruitfly-GEM.mat)`
x=readMat('~/Documents/aging_metabolic.models/Fruitfly-GEM-main/model/Fruitfly-GEM.mat')
x
ifly=x[[1]][,,1]
names(ifly)
ifly$S[1:3,1:3] #sparse matrix
ifly$rxnGeneMat[1:3,1:3] #sparse matrix
dim(ifly$S) #8135 mets X 11898 rxn
dim(ifly$rxnGeneMat) ## 11898 rxn X 1753 genes
genes=unlist(ifly$genes)
length(genes) #1753 genes

mz='HISTAMINE';
mz.index=grep(mz,unlist(ifly$metNames),ignore.case = T)
mz.index=mz.index[c(1,2)]
unlist(ifly$metNames[mz.index])
#ifly$S[mz.index,]
#colSums(ifly$S[mz.index,])
#which(colSums(ifly$S[mz.index,])!=0)
rxn.index=which(colSums(ifly$S[mz.index,])!=0)
rxn.index
gene.index=which(colSums(ifly$rxnGeneMat[rxn.index,])!=0)
gene.index
unlist(ifly$genes)[gene.index]

subSystems=unlist(ifly$subSystems)
length(subSystems)
subSystems[rxn.index]

rxn.index=which(subSystems=='Histidine metabolism')
gene.index=which(colSums(ifly$rxnGeneMat[rxn.index,])!=0)
gene.index
unlist(ifly$genes)[gene.index]
ifly$rxnNames[rxn.index]

load('age-associated.metabolites.for.Ming')
ls()
for.Ming$mz
metabolome=for.Ming %>% arrange(desc(beta),FDR) 
head(metabolome)
