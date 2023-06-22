##########
# in Loyfer, Netanel, et al. "A DNA methylation atlas of normal human cell types." Nature (2023): 1-10.
# several places mentioned: `Briefly, we maintain a 1 × N table T (N = 28,217,448 CpGs) for optimal segmentation scores across all pre- fixes. `
# where does these 28,217,448 CpGs come from? (see below)

# data can be accessed from: GSE186458
# `GSE186458_RAW.tar` is huge: 320 GB

# Goal: Use these 28,217,448 CpGs per cell type info to estimate age with aging clock developed using methylation array data.
# 1) get chr coordinate of all CpGs from 2023 atlas paper
# 2) extract aging clock methyaltion position and coefficients used in clock
# 3) use the methylation freq from each cell type to calculate age per cell type


library(ggplot2)
library(viridis)

##########
## chr coordinate of the 28,217,448 CpGs from known reference CpGs
#https://www.bioconductor.org/packages/release/bioc/vignettes/methrix/inst/doc/methrix.html
#https://compepigen.github.io/methrix_docs/index.html
#https://bioconductor.org/packages/devel/bioc/vignettes/ramwas/inst/doc/RW2_CpG_sets.html

#BiocManager::install("methrix")
#options(timeout = 1200)  # is 20 minutes enough? #https://support.bioconductor.org/p/p133518/
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
hg19_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19"))
#[1] Done. Extracted 28,217,448 CpGs from 25 contigs.
dim(hg19_cpgs$cpgs)
#[1] 28217448        5

##########
## get methylation array chip annotation
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/methylation_data_analysis/meffonym_101.R
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#[1] '0.6.0'
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
dim(ann450k) #485512     33

## get Hovarth's clock chr pos
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/wgbstools/methylclockData_101.R
library('dnaMethyAge')
data('HorvathS2013')
head(coefs)
cpgs=coefs$Probe[-1]

#pick.cpg=ann450k[cpgs,]
pick.cpg=ann450k #due to imputation, use all cpgs
x1=paste(pick.cpg$chr,pick.cpg$pos)
x2=paste(hg19_cpgs$cpgs$chr,hg19_cpgs$cpgs$start)

keep=x2 %in% x1
sum(keep) #482419, all 353 exist
hg19_cpgs_keep=hg19_cpgs$cpgs[keep,]
hg19_cpgs_keep$chr_pos=paste(hg19_cpgs_keep$chr,hg19_cpgs_keep$start)

pick.cpg$chr_pos=paste(pick.cpg$chr,pick.cpg$pos)
clock.cpg= merge(pick.cpg,hg19_cpgs_keep,by='chr_pos')
clock.cpg=clock.cpg[match(hg19_cpgs_keep$chr_pos,clock.cpg$chr_pos),]
sum(hg19_cpgs_keep$chr_pos==clock.cpg$chr_pos) #482419 or 353

##########
#https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md
if(F){
  files= Sys.glob('./beta_files/*.beta')
  files=files[-grep('hg38',files)]
  length(files) #253 files

  res<-lapply(files,function(fname){
   #fname <- PATH
   sample.id= strsplit(basename(fname),'_')[[1]][[1]]
   print(sample.id)
   N <- file.info(fname)$size
   content <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
   dim(content);
   x=content[keep,]
   methy.level=x[,1]/x[,2]
   methy.level[x[,2]<10]=NA
   return(methy.level)
  })
  df.res=Reduce(`cbind`,res)
  dim(df.res) # 482419 or 353 site x 253 samples

  sample.names=lapply(files,function(fname){
    #fname <- PATH
    sample.id= strsplit(basename(fname),'_')[[1]][[1]];
    sample.id})

  colnames(df.res)=unlist(sample.names)
  rownames(df.res)=clock.cpg$Name
  saveRDS(df.res,'253samples_hg19_cpgs.rds')
}
df.res=readRDS('253samples_hg19_cpgs.rds')

########## use R package `dnaMethyAge`
## estimate age
#https://github.com/mingwhy/bioinfo_homemade_tools/blob/main/package_101/wgbstools/methylclockData_101.R
#https://github.com/yiluyucheng/dnaMethyAge
library('dnaMethyAge')
betas=df.res[order(rownames(df.res)),]
sum(is.na(betas))

x=apply(betas,2,function(i) sum(is.na(i)))
betas2=betas[,x<=nrow(betas)*0.1] #only keep samples with less than 10% NA
dim(betas2)

horvath_age <- methyAge(betas2, clock='HorvathS2013',use_cores =4)

print(horvath_age)

##########
## get each sample or cell type annotation
library("GEOquery")
gse<-getGEO(filename='GSE186458_series_matrix.txt.gz',GSEMatrix = TRUE,getGPL = FALSE) 
df2=pData(gse) #253 samples x 64 cols
colnames(df2)

df3=df2[,c('age:ch1','cell type:ch1','tissue:ch1')]
colnames(df3)=c('age','cell.type','tissue');
df3$sample=rownames(df3)
dfc=merge(df3,horvath_age,by.x='sample',by.y='Sample')

dfc$tissue_tc=paste(dfc$tissue,dfc$cell.type,sep=':')
dfc
table(dfc$cell.type)
summary(dfc$mAge)

library(ggplot2)
  
ggplot(dfc,aes(x=reorder(dfc$cell.type, dfc$mAge, median),y=mAge))+#geom_violin()+
  geom_jitter()+theme_classic()+coord_flip()

########## use R package `methylclock`
#https://rpubs.com/jrgonzalezISGlobal/methylclock
#library(devtools)
#install_github("isglobal-brge/methylclock")
library(methylclock)
df.res=readRDS('253samples_hg19_cpgs.rds')

betas=df.res[order(rownames(df.res)),]
sum(is.na(betas))
x=apply(betas,2,function(i) sum(is.na(i)))
betas2=betas[,x<=nrow(betas)*0.1] #only keep samples with less than 10% NA
dim(betas2)

#some clocks returns NA when there are more than 80% of the required CpGs are missing as we can see when typing
TestDataset=betas2
cpgs.missing <- checkClocks(TestDataset) #Chronological DNAm age (in years)
cpgs.missing.GA <- checkClocksGA(TestDataset) #Gestational DNAm age (in weeks)
names(cpgs.missing)
cpgs.missing$Hannum

#age.example55.sel <- DNAmAge(TestDataset, clocks=c("Horvath", "Hannum","BNN"))
#age.example55.sel
if(F){
age.est <- DNAmAge(TestDataset)
saveRDS(age.est,'res_methylclock.rds')
}
age.est=readRDS('res_methylclock.rds')

dfc=merge(df3, age.est,by.x='sample',by.y='id')
dfc$tissue_tc=paste(dfc$tissue,dfc$cell.type,sep=':')
dfc
dfc$age = as.numeric(as.character(dfc$age ))
table(dfc$cell.type)
summary(dfc$Hannum)
colnames(dfc)
plot(dfc$Hannum,dfc$skinHorvath)
plot(dfc$age, dfc$skinHorvath)
cor(dfc$age, dfc$skinHorvath,use="pairwise.complete.obs")

dfc$mAge=dfc$skinHorvath
ggplot(dfc,aes(x=reorder(dfc$cell.type, dfc$mAge, median),
               y=mAge,col=dfc$age))+#geom_violin()+
  geom_jitter()+theme_classic()+coord_flip()+
  scale_color_viridis()



