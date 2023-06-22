options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)

# 1) check QC-I and QC-S CV mean and median values
# 2) get average values of each of the 12 SILIS (stable isotope-labeled internal standards)

## read in RDS data
input=readRDS('mz.dat.rds')
names(input)
attach(input)

## make some plots for protein abundance
protein.abundance$`BCA Total Protein, ug`=as.numeric(protein.abundance$`BCA Total Protein, ug`)
tmp=protein.abundance[!is.na(protein.abundance$`BCA Total Protein, ug`),]
tmp$batch=batch.id[tmp$`Original Sample ID`]
tmp$batch=factor(tmp$batch)
tmp=tmp[order(tmp$batch),]
tmp$`Original Sample ID`=factor(tmp$`Original Sample ID`,levels=tmp$`Original Sample ID`)
ggplot(tmp,aes(x=`Original Sample ID`,y=`BCA Total Protein, ug`,col=batch))+
  geom_point(size=4)+theme_bw()+
  ggtitle('Total Protein (ug)')+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

## check CV of quality control samples
dim(control)
colnames(control)

quality.control=control[,grep('QC',colnames(control))]
tags=colnames(quality.control)
x=lapply(quality.control,as.numeric)
quality.control=as.data.frame(Reduce(`cbind`,x))
colnames(quality.control)=tags

i.cv=apply(quality.control[,1:7],1,function(i) sd(i,na.rm=T)/mean(i,na.rm=T))
summary(i.cv)
mean(i.cv,na.rm=T) #0.07140
median(i.cv,na.rm=T) #0.0647

s.cv=apply(quality.control[,8:14],1,function(i) sd(i,na.rm=T)/mean(i,na.rm=T))
summary(s.cv)
mean(s.cv,na.rm=T) #0.0806
median(s.cv,na.rm=T) #0.0740

## check CV of SILIS
silis=quality.control[is.na(control$`KEGG ID NA`),]
i.cv=apply(silis[,1:7],1,function(i) sd(i,na.rm=T)/mean(i,na.rm=T))
summary(i.cv) #0.02~0.10
mean(i.cv,na.rm=T) #0.06164
median(i.cv,na.rm=T) # 0.05934

s.cv=apply(silis[,8:14],1,function(i) sd(i,na.rm=T)/mean(i,na.rm=T))
summary(s.cv) #0.05 ~ 0.06
mean(s.cv,na.rm=T) #0.0564
median(s.cv,na.rm=T) #0.0569

## use SILIS to evaluate the data reliability for each sample
## use average SILIS per sample to evaluate sample data relibility
silis.index=which(is.na(control$`KEGG ID NA`))
length(silis.index) #32
dim(dat) #393 x 60
silis.global.mean=mean(as.numeric(unlist(dat[silis.index,])),na.rm=T)

x=apply(dat[silis.index,],2,function(i) mean(as.numeric(i),na.rm=T))
summary(abs(x-silis.global.mean)/silis.global.mean)
# 0.005 ~ 0.14, all within 15%, all samples retained

