#you may need to install dependencies ‘RnBeads’, ‘RnBeads.hg38’ first
BiocManager::install("RnBeads")
BiocManager::install("RnBeads.hg38")
BiocManager::install('RnBeads.mm10')
devtools::install_github("MPIIComputationalEpigenetics/WSHPackage")

#https://github.com/MPIIComputationalEpigenetics/WSHPackage/blob/master/vignettes/WSH.md
library(WSH)
example.bam <- system.file(file.path("extData","small_example.bam"),
                           package="WSH")
example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),
                               package="WSH")
example.rnb.set <- load.rnb.set(example.rnb.set)
example.rnb.set
example.rnb.set@sites

?set.option
set.option(coverage.threshold = 10) 
#Coverage Threshold emloyed to select the sites in a RnBSet annotation that fullfill havin a coverage
#'                   higher than this threshold
qfdrp <- rnb.calculate.qfdrp(example.rnb.set,example.bam)
dim(qfdrp) #448
dim(example.rnb.set@sites)  #448

#BiocManager::install('RnBeads.mm10')

################################################################
library(WSH)
example.bam <- system.file(file.path("extData","small_example.bam"),
package="WSH")
example.GRanges <- GRanges(Rle(rep("chr2",10)),
                          IRanges(start=c(2298361,2298554,2298732,
                                          2298743,2298787,2298792,
                                          2298827,2298884,2298915,2298921),
                                  end=c(2298361,2298554,2298732,
                                        2298743,2298787,2298792,
                                        2298827,2298884,2298915,2298921)+1))

# note, need to give per CpG site GRange object
pdr <- compute.score(bam.file=example.bam,example.GRanges,score="pdr")


source('~/Documents/aging_RRBS/WSHPackage-master/R/main.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/utilities.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/calculate_scores.R')
# from calculate_scores.R script
pdr=calculate.pdr(bam.file=example.bam, anno=example.GRanges,
                  mapq.filter = 20,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
pdr$PDR

#set.option(mapq.filter = 20)
set.option(fdrp.type='FDRP')
fdrp.out=calculate.fdrp.score(bam.file=example.bam, anno=example.GRanges,
                              cores=1,
                              use.sex.chromosomes=FALSE,ignore.strand=TRUE)
fdrp.out

set.option(fdrp.type='qFDRP')
pfdrp.out=calculate.qfdrp(bam.file=example.bam, anno=example.GRanges,
                          cores=1,
                          use.sex.chromosomes=FALSE,ignore.strand=TRUE)
plot(fdrp.out$FDRP,pfdrp.out$qFDRP)
plot(fdrp.out$FDRP,pdr$PDR)

# to use WSH, https://github.com/MPIIComputationalEpigenetics/WSHPackage/blob/master/vignettes/WSH.md
# need to index bam file first, for example 

# if you check "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/WSH/extData/small_example.bam",
# there is a 'small_example.bam.bai' file.
$ samtools view small_example.bam | less -SN

# in my case, need to sort then index
$ samtools sort SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o SRR5195656_sorted.bam
$ samtools index -b SRR5195656_sorted.bam 
$ samtools view -f 2 SRR5195656_sorted.bam | less -SN 

########################################################
## play with calculate_scores.R
bam.file=example.bam; 
anno = example.GRanges; 

cores=2
window.size=unname(get.option('window.size'))
max.reads=unname(get.option('max.reads'))
mapq.filter=unname(get.option('mapq.filter'))
coverage.threshold=unname(get.option('coverage.threshold'))
use.sex.chromosomes=FALSE
ignore.strand=TRUE
#set.option(mapq.filter = 0)
#set.option(coverage.threshold=1)

output.frame <- data.frame(chromosome=seqnames(anno),start=start(anno),end=end(anno))
bam <- BamFile(bam.file)

source('~/Documents/aging_RRBS/WSHPackage-master/R/utilities.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/calculate_scores.R')
anno_by_chr <- prepare.annotation(anno,use.sex.chromosomes=use.sex.chromosomes) #separate into different chrs
anno=anno_by_chr[[1]] #one chr at one time
#anno=sort(anno)

#calculate.pdr.by.chromosome(bam,chromosome,ignore.strand = ignore.strand)
(chromosome <- as.character(seqnames(anno))[1])
is.sex.chromosome <- grepl("X|Y|23|24",chromosome)
(start <- start(ranges(anno)[1]))
(end <- end(ranges(anno)[length(anno)]))
if(!(chromosome %in% names(scanBamHeader(bam)[[1]]))){
  if(!is.sex.chromosome){
    chromosome <- unique(na.omit(as.numeric(unlist(strsplit(chromosome,"[^0-9]+")))))
  }else{
    chromosome <- unique(gsub("chr","",chromosome))
  }
}
chromosome
which <- paste0(chromosome,":",start,"-",end)
which <- GRanges(which)
which #combine all inquiry intervals into one range on 1 chr
param <- ScanBamParam(which=which,what="seq",mapqFilter=get.option('mapq.filter'),
                      flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
reads <- readGAlignments(bam,param=param) #only read in one subset of all mapped bam
reads # 694176  alignments

#reads0<-readGAlignments(input.bam)
#reads0 #14008968 alignments 

range_reads <- GRanges(reads)
rm(reads)
newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
newStyle <- newStyle[!is.na(newStyle)]
range_reads <- renameSeqlevels(range_reads,newStyle)

# we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
range_reads #694176 query read
anno #4795 interested ranges
seqlengths(range_reads) #this information is in bam file, reference chr length(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt)
seqlengths(anno) #same as above
overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)

query <- queryHits(overlap)
query <- unique(query)
range_reads <- range_reads[query]

####################################################################################################
#### play with `calculate.pdr` (calculate_scores.R)
####### REPRESENTATION ############
#This part clasifies all reads into either discordant or concordant
range_cpgs <- ranges(anno)
starts_cpgs <- start(range_cpgs)
rm(range_cpgs)
seqs_reads <- as.character(values(range_reads)$seq)
starts_reads <- start(ranges(range_reads))

overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
overlap #queryLength: 252 / subjectLength: 10
range_reads # 252 ranges, inqury reads
match_read_cpg <- as(overlap,"list")
length(match_read_cpg) #252, each read is a element in the list, length is the #CpG this read contains
table(sapply(match_read_cpg,length))
# contain reads with>=4 CpG ??
#1  2  3  4  5 
#95 51 29 23 54 

overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
pdrs <- as.list(rep(NA,length(anno)))
#rm(anno)

# we classify each read into either discordant or concordant
classified_reads <- as.list(1:length(range_reads))
classified_reads <- lapply(classified_reads,classify.read,
                           match_read_cpg,starts_cpgs,starts_reads,seqs_reads)

#table(unlist(match_read_cpg))
head(which(unlist(match_read_cpg)>4))
classify.read(20,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)

rm(match_read_cpg)
rm(starts_reads)
rm(seqs_reads)
rm(starts_cpgs)
values(range_reads) <- DataFrame(cbind('isDiscordant'=classified_reads))
rm(classified_reads)
logger.completed()


########################################################################################################################
#### play with `calculate.fdrp.by.chromosome` and `calculate.fdrp.score` (calculate_scores.R)
####### REPRESENTATION ############
# This part converts the raw sequencing reads from the alignment into
# an representation, where only CpG positions are considered and from whom we
# can infer discordance or concordance of reads
logger.start(paste('Representation',chromosome))
range_cpgs <- ranges(anno)
starts_cpgs <- start(range_cpgs)
rm(range_cpgs)
seqs_reads <- as.character(values(range_reads)$seq)
starts_reads <- start(ranges(range_reads))
overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
match_read_cpg <- as(overlap,"list")
overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
fdrps <- as.list(rep(NA,length(anno)))
#rm(anno)
# for each read we convert the covered CpG sites into a custom representation
read_representation <- as.list(1:length(range_reads))
read_representation <- lapply(read_representation,toCpGs,
                              match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
length(match_read_cpg) #659676
table(sapply(match_read_cpg,length))
#1  2  3  4  5 
#95 51 29 23 54 
length(read_representation) #659676
table(sapply(read_representation,length))
#1  2  3  4  5 
#95 51 32 20 54 

rm(match_read_cpg)
rm(starts_reads)
rm(seqs_reads)
values(range_reads) <- DataFrame(cbind(CpG=read_representation))
rm(read_representation)
logger.completed()

## FDRP CALCULATION
# we only calculate the FDRP for the CpGs that are acutally covered by any read in the
# corresponding sample
logger.start(paste('FDRP',chromosome))
match_cpg_reads <- as(overlap,"list")
rm(overlap)

null <- lapply(match_cpg_reads,function(x){length(x)>0})
null <- unlist(null)
match_cpg_reads <- match_cpg_reads[null]
starts_cpgs <- starts_cpgs[null]
toApply <- 1:length(starts_cpgs)

fdrps_actual <- lapply(toApply,calculate.fdrp.site,
                       match_cpg_reads,range_reads,starts_cpgs)

#calculate.fdrp.site(toApply[[1]],match_cpg_reads,range_reads,starts_cpgs)

# when we do not have a read that covers this site, we set the FDRP for this site to NA
length(null)
fdrps[null] <- fdrps_actual
fdrps <- unlist(fdrps)
table(fdrps)
