options(stringsAsFactors = F)
library(readxl)
library(ggplot2)
library(gridExtra)
# There are 4 sheets in the excel file.
#sheet1: sample ID, protein total abundance, 60smaples + 2x7control=74 rows
#sheet2: controls, mz abundance, 393 mz in total, the last 32 are spike stable isotope
#sheet3: sample, mz abundance, 393 mz in total, the last 32 are spike stable isotope
#sheet4: mz information

## read in total protein (ug) for all samples
protein.abundance=read_excel('./2020-12-08_Tatar-60_Data.xlsx',sheet='Sample ID & Protein Assay');
protein.abundance=protein.abundance[,c(1,2,3)];
dim(protein.abundance) # 74 = 2*7 control + 60 samples

## read in 'control control' samples
control=read_excel('./2020-12-08_Tatar-60_Data.xlsx',sheet='Data Reproducibility')
dim(control)  # 393  23
head(control)
tag=do.call(paste,list(colnames(control),control[1,]))
tag=gsub("...\\d+ ",'',tag)
tag[4]='QC(I)#1'
tag[14]='QC(S)#1'

control=control[-1,]
colnames(control)=tag
dim(control) 

sum(protein.abundance$`Original Sample ID` %in% tag) #labels match

## read in Data
df=read_excel('./2020-12-08_Tatar-60_Data.xlsx',sheet='Data')
colnames(df)

mz=df[-c(1,2,3,4),1:3]
dim(mz) #393 metabolites 

dat=df[,-c(1,2,3)]
sample.id=as.character(dat[1,])#col info
batch.id=rep(c(1,2,3),each=20) #col info
names(batch.id)=sample.id 
genotype.id=as.character(dat[2,]) #col info
names(genotype.id)=sample.id
age=as.character(dat[3,]) #col info
names(age)=sample.id


dat=df[-c(1:4),-c(1,2,3)]
dim(dat) #393 x 60
colnames(dat)=sample.id

saveRDS(list('protein.abundance'=protein.abundance,
             'control'=control,
             'dat'=dat,
             'mz'=mz,
             'sample.id'=sample.id,
             'batch.id'=batch.id,
             'genotype.id'=genotype.id, 'age'=age),
        file='mz.dat.rds',)
