# https://elifesciences.org/articles/40675
# https://github.com/kerepesi/MouseAgingClocks
library(readxl)
library(stringr)
library(dplyr)

# read in clock
clock=read_excel('MouseAgingClocks-main/ClockData/elife-40675-supp3-v2.xlsx', sheet='Blood', n_max=91)
dim(clock) #90 site in this clock
plot(clock$Weight)
sum(clock$Weight>0) #53

# read in Meer sample age info (https://elifesciences.org/articles/40675/figures#supp1)
dat=read_excel('Meer_elife-40675-supp1-v2.xlsx',sheet='Results on the clock sets',skip=1)
table(dat$Paper)
dat=dat[dat$Paper=='Petkovich',]
dim(dat) #259 samples 
x=table(dat$`Age, days`) 
x[order(as.numeric(names(x)))]

young=dat[as.numeric(dat$`Age, days`)<=20,]
young$ID

all.ages=dat;
all.ages$age=as.numeric(all.ages$`Age, days`)

# download raw data from GSE80672_RAW
# need to perform id conversion first, as it used ncbi RefSeq id, 
# index   mipsL3_GCCAAT.CellEF.Percentage mipsL3_GCCAAT.CellEF.Coverage
# gi|149247747|ref|NT_166280.1|:23455     92.5925925925926        27
# gi|372099109|ref|NC_000067.6|:74027683  46.6666666666667        30
# NT_166280.1, https://www.ncbi.nlm.nih.gov/nuccore/NT_166280
# NC_000067.6, https://www.ncbi.nlm.nih.gov/nuccore/372099109
# click: `Assembly: GCF_000001635.27` or `Assembly: GCF_000001635.26`
# then `Download the full sequence report` (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26)
#chromAlias=data.table::fread('GCF_000001635.27_GRCm39_assembly_report.txt',skip=28)
chromAlias=data.table::fread('GCF_000001635.26_GRCm38.p6_assembly_report.txt',skip=41)
head(chromAlias)
chromAlias$`RefSeq-Accn`

sum(clock$Chromosome %in% chromAlias$`UCSC-style-name`) #all chr exist

clock.chr=merge(clock,chromAlias,by.x='Chromosome',by.y='UCSC-style-name',all.x=T)
chr.names=paste(clock.chr$`RefSeq-Accn`,clock.chr$Position,sep=':')
length(chr.names) #90 site
clock.chr$chr.names=chr.names

################################################################################################
# read in raw epi data, check if files from all samples ids have RefSeq <=> CHR mapped.
(files=Sys.glob('Petkovich_GSE80672_RAW/*.txt.gz'))

keep.samples=list()
for(i in 1:nrow(all.ages)){
  (one.file=files[grep(paste0('_',all.ages$ID[i],'.overlap.txt.gz'),ignore.case = T,files)])
  if(length(one.file)==0){next}
  
  dat=data.table::fread(one.file) #can read txt.gz file directly
  #head(dat)
  Ref.id=str_extract(dat$index,'\\|N.+\\|')
  Ref.pos=str_extract(dat$index,'\\:\\d+$')
  Ref.id=gsub('\\|','',Ref.id)
  #Ref.pos=gsub('\\:','',Ref.pos)
  Ref.id.pos=paste(Ref.id,Ref.pos,sep='')
  cat(i,one.file,'overlapped #site',sum(chr.names %in% Ref.id.pos),'\n') #all exist
  # extract clock site mythelation level data from each sample
  dat$chr.names=Ref.id.pos
  dat=as.data.frame(dat)
  dat.keep<-dat[dat$chr.name %in% chr.names,]
  dat.keep$ID=all.ages$ID[i];
  dat.keep=as.data.frame(dat.keep);
  colnames(dat.keep)=c('index','freq','coverage','chr.names','ID')
  
  # add clock.chr, site weight
  head(clock.chr)
  dat.keep2=merge(dat.keep,clock.chr[,c('Chromosome','Position','Weight','chr.names')])
  dim(dat.keep2)
  
  # add sample age
  dat.keep3=merge(dat.keep2,all.ages[,c('ID','age')])
  keep.samples[[all.ages$ID[i]]]<-dat.keep3
}
length(keep.samples) #246
saveRDS(keep.samples,'Petkovich_246samples_90clockSites.rds')
}

###########
keep.samples=readRDS('Petkovich_246samples_90clockSites.rds')
head(keep.samples[[1]])
sapply(keep.samples,dim)
df.keep.samples=as.data.frame(Reduce(`rbind`,keep.samples))

## plot freq.distribution for each age
tmp=df.keep.samples %>% group_by(age,index) %>% summarise(mean.freq=mean(freq))
(all.ages=sort(unique(tmp$age)))
length(unique(tmp$index)) #90 site

pdf('Petkovich_methy_age.pdf',height = 12,width = 16)
par(mfrow=c(5,5))
lapply(all.ages,function(i){
  dat.1w=tmp[as.numeric(tmp$age)==i,,drop=FALSE]
  mean.freq=dat.1w$mean.freq
  #hist(mean.freq,main=paste0('mean.methy.level at day',i))
  #hist(mean.freq,main=paste0('methy.level across 90 sites\nat day',i))
  hist(mean.freq,main=paste0('methy.level across 90 sites\nat month ',round(i/30,2)))
})
dev.off()


pred.out<-sapply(keep.samples,function(i){
  MSc=sum(i$Weight*i$freq)
  #a = 0.1666
  #b = 0.4185 
  #c = - 1.712
  #Age = ((MSc - c)/a)**(1/b)
  #Pred=Age/30.5
  c(i$age[1],MSc)
})
pred.out=t(pred.out)
plot(pred.out[,1],pred.out[,2],xlab='age',ylab='pred.age')
cor(pred.out[,1],pred.out[,2])

################################################################################################
# read in raw epi data, check if files from young samples ids have RefSeq <=> CHR mapped.
(files=Sys.glob('Petkovich_GSE80672_RAW/*.txt.gz'))
keep.samples=list()
for(i in 1:nrow(young)){
  (one.file=files[grep(paste0('_',young$ID[i],'.overlap.txt.gz'),ignore.case = T,files)])
  dat=data.table::fread(one.file) #can read txt.gz file directly
  #head(dat)
  
  Ref.id=str_extract(dat$index,'\\|N.+\\|')
  Ref.pos=str_extract(dat$index,'\\:\\d+$')
  Ref.id=gsub('\\|','',Ref.id)
  #Ref.pos=gsub('\\:','',Ref.pos)
  Ref.id.pos=paste(Ref.id,Ref.pos,sep='')
  cat(i,one.file,'overlapped #site',sum(chr.names %in% Ref.id.pos),'\n') #all exist
  # extract clock site mythelation level data from each sample
  dat$chr.names=Ref.id.pos
  dat=as.data.frame(dat)
  dat.keep=dat[Ref.id.pos %in% chr.names,]
  dat.keep$ID=young$ID[i];
  
  keep.samples[[young$ID[i]]]<-dat.keep
}  
x=lapply(keep.samples,function(i) {i=as.data.frame(i);colnames(i)=c('index','freq','coverage','chr.names','ID');i})
df.samples=as.data.frame(Reduce(`rbind`,x))
head(df.samples)

# add clock.chr, site weight
head(clock.chr)
df.samples2=merge(df.samples,clock.chr[,c('Chromosome','Position','Weight','chr.names')])
head(df.samples2)

# add sample age
head(young)
sum(young$ID %in% df.samples2$ID)
young$age=as.numeric(young$`Age, days`)
df.samples3=merge(df.samples2,young[,c('ID','age')])
plot(df.samples3$freq,df.samples3$Weight)
cor(df.samples3$freq,df.samples3$Weight) #-0.5478639
saveRDS(df.samples3,'Petkovich_young_90clockSites.rds')


