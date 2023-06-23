# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1203-5
# https://github.com/EpigenomeClock/MouseEpigeneticClock  (original Genome Biology paper)
# https://github.com/kerepesi/MouseAgingClocks (rejuvenation Science Advances paper)
# adapted from https://github.com/kerepesi/MouseAgingClocks/blob/main/AppMouseGenomicClocks.py
## Stubbs multi-tissue clock

# read in methy data (https://github.com/kerepesi/MouseAgingClocks)
#TrainingData_Babraham_Reizel_Cannon.txt (Stubbs et al. 2017. Genome Biology
dat=as.data.frame(data.table::fread('MouseAgingClocks-main/ClockData/TrainingData_Babraham_Reizel_Cannon.txt'))
dim(dat) #17992   130 
dat[1:3,1:3]
all.ages.dat=dat;

pos=dat[,1] #chr position
dat=dat[,-1]
rownames(dat)=pos
dim(dat) #17992   129 samples

# read in clock
library(readxl)
library(stringr)
# Stubbs GB paper: Additional file 6:Table of mouse epigenetic clock sites. (XLS 94 kb)
# or download from elife-40675-supp3-v2.xlsx (Meer et al. 2018. Elife, Supplementary file 3)
YOMT=read_excel('MouseAgingClocks-main/ClockData/elife-40675-supp3-v2.xlsx', sheet='Young age multi-tissue', n_max=329)
dim(YOMT) #329 site in Stubbs clock

# reduce dat to 329 site only
pos2=paste(gsub('chr','',YOMT$Chromosome),YOMT$Position,sep=':')
sum(pos%in% pos2) #329
dat=dat[pos2,]

sum(rownames(dat)==pos2) #329, order row in dat the same order as YOMT
dim(dat) #329 site x 129 samples

############################################
# find young age samples as ground state 
dat[1:3,1:3]
sample.names=colnames(dat)
age=str_extract(sample.names,'\\d+w|\\d+wk')
sample.info=data.frame(sample.names=sample.names,age=age)
table(sample.info$age)

sample.info1=sample.info[!is.na(sample.info$age),]
sample.info2=sample.info[is.na(sample.info$age),]
# sample.info1: Reizel and Cannon samples
# sample.info2: Stubbs's own sequencing data

# handle sample.info2
library("GEOquery")
if(F){ #online
  gse <- getGEO("GSE93957") #GSE93957_series_matrix.txt.gz
  e <- gse[[1]] #62 samples
  pd <- pData(e) #62 x 41 df
  head(pd$title)
}
# read from local downloaded file
gse <- getGEO(filename='Stubbs_GSE93957_series_matrix.txt.gz')
pd <- pData(gse) #62 x 41 df
head(pd$title)

sample.info2[grep('\\.',sample.info2$sample.names),]
sample.info2$sample.new.names=gsub('\\.','',sample.info2$sample.names)

df1=as.data.frame(Reduce(`rbind`,strsplit(pd$title,'_')));
colnames(df1)=c('tissue','age','id')

sample.info2$id=unlist(lapply(strsplit(sample.info2$sample.new.names,'_'),'[',1))
dim(sample.info2) #47
sum(sample.info2$id %in% df1$id) #47
df2=merge(df1,sample.info2,by='id')

sample.info2=data.frame(sample.names=df2$sample.names,age=df2$age.x)

# update sample.info
sample.info=rbind(sample.info1,sample.info2[,c('sample.names','age')])
table(sample.info$age)
sample.info$age.num=gsub('w|wk','',sample.info$age)
table(sample.info$age.num)

######################################################
## extract week 1 samples as ground state samples
dim(sample.info) #129 samples
sample.info=sample.info[match(colnames(dat),sample.info$sample.names),]
sum(sample.info$sample.names==colnames(dat))

sample.info$age.num=as.numeric(as.character(sample.info$age.num))
dat=as.data.frame(dat)
table(sample.info$age.num)

## plot freq.distribution for each age
all.ages=sort(unique(sample.info$age.num))

pdf('Stubbs_methy_age.pdf',useDingbats = T,width = 12)
par(mfrow=c(3,4))
lapply(all.ages,function(i){
  dat.1w=dat[,as.numeric(sample.info$age.num)==i,drop=FALSE]
  mean.freq=rowMeans(dat.1w)
  hist(mean.freq,main=paste0('mean.methy.level at age week',i))
})
dev.off()

#################################
## young age
table(sample.info$age.num)
#1  3  8  9 14 20 27 28 31 41 
#25 18 10  7  8 30 16  1  2 12 
dat.1w=dat[,as.numeric(sample.info$age.num)==1]
dim(dat.1w) #329 site x 25 1w samples
dim(YOMT) #329 site

## plot site mean freq VS site coeff in the clock
mean.freq=rowMeans(dat.1w)
hist(mean.freq,main='distribution of mean.methy.levels at age week1')

coeff=YOMT$Weight
cor(mean.freq,coeff,method='spearman') #-0.10
plot(mean.freq,coeff,cex=1)
cor.test(mean.freq,coeff) #-0.10
#plot(dat.1w[,4],coeff,cex=0.5)

######################################################
## check the clock is doing prediction
dim(dat) #329 * 129 sample
pred.age=t(dat) %*% coeff
plot(sample.info$age.num,pred.age)
cor(sample.info$age.num,pred.age) #0.8752743

#https://github.com/kerepesi/MouseAgingClocks/blob/main/AppMouseGenomicClocks.py
#FT_app_reind_clock_fillna_norm=(FT_app_reind_clock_fillna - Train_reind.median()) / Train_reind.std()
#MSc=((FT_app_reind_clock_fillna_norm)*YOMT['Weight']).sum(axis=1)
#MSc=MSc[0]
#Age = np.exp (0.1207 * (MSc**2) + 1.2424 * MSc + 2.5440) - 3
#Pred=Age*(7/30.5)
site.median=apply(dat,1,median)
site.sd=apply(dat,1,sd)
x=lapply(1:nrow(dat),function(i) (dat[i,]-site.median[i])/site.sd[i] )
dat.norm=as.data.frame(Reduce(`rbind`,x))
dim(dat.norm) #329 x 129 samples

pred.age=t(dat.norm) %*% coeff
plot(sample.info$age.num,pred.age)
cor(sample.info$age.num,pred.age) #0.89
pred.age.in.weeks = exp(0.1207*(pred.age^2)+1.2424*pred.age + 2.5440) 
plot(sample.info$age.num,pred.age.in.weeks)

#pred.age2=(7/30.5)* (exp(0.1207*(pred.age^2)+1.2424*pred.age + 2.5440)-3)
#cor(sample.info$age.num,pred.age2) # 0.92
#plot(sample.info$age.num,pred.age2)
