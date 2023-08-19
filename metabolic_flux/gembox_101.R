#change to working dir: ./gembox-ming/

#############################################
# https://github.com/ruppinlab/gembox
# need to install gurobi (solver algorithm) (or Rcplex2: https://github.com/ruppinlab/Rcplex2)
# library(gurobi)
# check out install_gurobi.R
#devtools::install_github("ruppinlab/gembox")
library(gembox)
## lee et al. BMC Syst Biol 2012: https://github.com/ruppinlab/gembox/blob/master/R/exprs.R

#README.md from https://github.com/ruppinlab/covid_metabolism
# 1, download data
#$ cd covid_metabolism-main/data
#$pip3 install gdown #for use gdown in ./dnload.sh
#$ ./dnload.sh
# covid_metabolism_files.tar.gz (270MB), liao.RDS(3.6G), chua.RDS(1.7G) files.

# process data `collect.validation.data.R`
#BiocManager::install('GSA')    
#devtools::install_github("ImNotaGit/my.utils")
#library(my.utils)
(files=Sys.glob('./my.utils/R/*R')) #download `ImNotaGit/my.utils` from github
# https://github.com/ImNotaGit/my.utils
# comment out `rmote::xxx` in `utils.R`, install 'httpgd', then source files.
#install.packages("httpgd")
lapply(files, source)

# an example from Genome-scale metabolic modeling reveals SARS-CoV-2-induced metabolic changes and antiviral targets
# https://www.embopress.org/doi/full/10.15252/msb.202110260
# full code repo: https://github.com/ruppinlab/covid_metabolism
# download required data on google drive 
# `Alternatively, one can manually download the data files from this Google Drive link then decompress them into the data folder.`
# you'd have `covid_metabolism_files.tar.gz`.
# there is `de.and.gsea.res.RData` in folder `expression`.
# go to `covid_metabolism/GEM` folder, run `prepare.data.R` and `imat.and.mta.R` scripts.

there is a error in `achr.cpp` in gembox, either fix it and reload while using it via `
Rcpp::sourceCpp('./achr_ming.cpp')`.
Alternatively fixing achr.cpp (check out FIX), re-install gembox from load folder
install.packages("~/Downloads/gembox-ming/", repos = NULL, type = "source")


#############################################
## showcase

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



