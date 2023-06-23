
#BiocManager::install('RnBeads.mm10')
library(WSH)

#################################################################
## get mapped genome reference info (chr length)

input.bam='SRR5195656_sorted.bam'
# both SRR5195656_sorted.bam and SRR5195656_sorted.bam.bai need to be present!!!

which = GRanges('1:1-100')
param <- ScanBamParam(which=which,what="seq",mapqFilter=20,
                      flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
reads <- readGAlignments(input.bam,param=param) #only read in one subset of all mapped bam
reads # 1224960  alignments
GenomeInfoDb::seqlengths(reads) #chrosomoe length information, make sure it's consistent with cpgs below



#################################################################
cpgs <- unlist(rnb.get.annotation("CpG",assembly="mm10")) #43735100
pos.cpgs=start(cpgs)
chr.cpgs=seqnames(cpgs)
seqlengths(cpgs)

methy.calls=as.data.frame(data.table::fread('SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov.gz'))
head(methy.calls) 
#The coverage output looks like this
#<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
colnames(methy.calls)=c('chr','start','end','percentage','methy','nonmethy')
methy.calls$coverage=methy.calls$methy+methy.calls$nonmethy
summary(methy.calls$coverage)
methy.pick<-methy.calls[methy.calls$coverage>=10,]
dim(methy.calls) #3699270
dim(methy.pick)#1027681 

unique(methy.pick$chr) #only look at autosomal genes (mt, X, Y removed)
methy.auto.pick<-methy.pick[!is.na(as.numeric(methy.pick$chr)),]
dim(methy.auto.pick) #1017542
table(methy.auto.pick$chr)

methy.range <- GRanges(Rle(paste0('chr',methy.auto.pick$chr)),
                           IRanges(start=methy.auto.pick$start,
                                   end=methy.auto.pick$start+1))
example.GRanges<-GenomicRanges::intersect(methy.range,cpgs,ignore.strand = TRUE)

example.GRanges #718259
seqlengths(example.GRanges)

######################
source('/gscratch/csde-promislow/bioinfo_software/WSHPackage-master/R/main.R')
source('/gscratch/csde-promislow/bioinfo_software/WSHPackage-master/R/utilities.R')
source('/gscratch/csde-promislow/bioinfo_software/WSHPackage-master/R/calculate_scores.R')

set.option(mapq.filter = 20)
# from calculate_scores.R script
system.time( {
  pdr=calculate.pdr(bam.file=input.bam, anno=example.GRanges,
                  cores = 8,use.sex.chromosomes=FALSE,ignore.strand=TRUE) 
}) #20min on server
#  ssh mingy16@n2357 #check if multiple threads were implemented
summary(pdr$PDR)
tmp=pdr[!is.na(pdr$PDR),] #543971
sum(tmp$PDR>0) #381224
saveRDS(pdr,'SRR5195656_pdr.rds')

set.option(fdrp.type='FDRP')
system.time( {
fdrp.out=calculate.fdrp.score(bam.file=input.bam, anno=example.GRanges,
                  cores=8,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
})
fdrp.out
summary(fdrp.out$FDRP) #if no mapped read or mapped.read<2, return NA for one CpG site
fdrp.out[which(fdrp.out$FDRP>0),]
#     user    system   elapsed 
# 57.998    88.213 28608.796 , 8hr for 718259 site
saveRDS(fdrp.out,'SRR5195656_fdrp.rds')

set.option(fdrp.type='qFDRP')
system.time( {
pfdrp.out=calculate.qfdrp(bam.file=input.bam, anno=example.GRanges,
                          cores=8,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
})
pfdrp.out
pfdrp.out$qFDRP
#user    system   elapsed 
#58.654    78.972 28798.863 , 8hr for 718259 site
saveRDS(pfdrp.out,'SRR5195656_pfdrp.rds')
