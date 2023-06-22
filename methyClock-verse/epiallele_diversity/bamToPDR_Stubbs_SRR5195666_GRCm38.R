#Stubbs, Thomas M., et al. "Multi-tissue DNA methylation age predictor in mouse." Genome biology 18.1 (2017): 1-14.
#GSE93957: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93957
#GSM2465660	Lung_14wk_M0420525

######
# Obtaining sample information from GEO
library("GEOquery")
gse <- getGEO("GSE93957")
e <- gse[[1]] #62 samples
pd <- pData(e) #62 x 41 df
unique(pd$`tissue:ch1`)
# "Cortex" "Heart"  "Liver"  "Lung"  
# subset for heart and lung samples only
pd <- pd[grepl("Heart|Lung", pd$title),] #31 samples 
table(pd$`age:ch1`,pd$`tissue:ch1`)

# on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93957, click `SRA`,
# go to SRA(https://www.ncbi.nlm.nih.gov/sra?term=SRP097629), use the `Send to-> File -> Format(RunInfo)` to get SRA run ID.
srr <- data.table::fread("sra_data//SraRunInfo.csv")
srrsmall <- srr[,c("SampleName", "Run", "Experiment", "Sample","BioSample",  "avgLength", "download_path")]
colnames(srrsmall)[which(colnames(srrsmall) == "SampleName")] <- "geo_accession"

# merge the GEO phenotypic information table and the sample information from SRA 
target <- merge(pd, srrsmall, by ="geo_accession")
rownames(target) <- target$Run

#The SRA names and SRA file paths were saved to help extract the SRA files from NCBI.
write.table(target$Run, file = "sra_data/sraFiles.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(target$download_path, file = "sra_data/sraFilesPath.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)

######
# Obtaining FASTQ files from SRA
# Downloading all .sra files in the sraFilesPath.txt:
$ cd sra_data
$ cat sraFilesPath.txt | wget -i -

#Extracting .fastq files
$ brew install parallel
# install fastq-dump from: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
$ cd fastq/
$ cat ../extdata/sraFiles.txt | parallel -j 6 ~/Documents/bioinfo_software/sratoolkit.3.0.1-mac64/bin/fastq-dump -I --split-files ../sra/{}.sralite.1 
($~/Documents/bioinfo_software/sratoolkit.3.0.1-mac64/bin/fastq-dump -I --split-files ../sra_data/SRR5195666.sralite.1 )
# SRR5195666.sralite.1, 1.2GB. after extraction, two files each 9GB

#######
# to use Trim_Galore, 
# require cutadapt and fastqc (check out README.txt file)

# install cutadapt 
$ pip3 install cutadapt

# install fastqc
# java version: fastqc from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# or command line: 
$ conda install -c bioconda fastqc
$ fastqc *.fastq
# Look at QC report produced by fastqc: fastqc_report.html

$cutadapt --version
# 4.3
$fastqc -v
# FastQC v0.12.1

#######
# read trimming, based on paper method, 
# following protocol: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
# Trim_Galore_User_Guide.md

# Specialised Trimming, Mouse Epigenetic Clock trimming, use option `--clock`
# 10min
$~/Documents/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --clock --fastqc SRR5195666.sralite.1_1.fastq SRR5195666.sralite.1_2.fastq

# trim following Trim_Galore_User_Guide.md, 5min
$~/Documents/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 SRR5195666.sralite.1_1.clock_UMI.R1.fq SRR5195666.sralite.1_2.clock_UMI.R2.fq
#Writing validated paired-end Read 1 reads to SRR5195666.sralite.1_1.clock_UMI.R1_val_1.fq
#Writing validated paired-end Read 2 reads to SRR5195666.sralite.1_2.clock_UMI.R2_val_2.fq

#######
# reference setup  
# First download the mouse genome file (.fa) from ENSEMBL.
$ cd reference
$ wget ftp://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
$ gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Building the bisulfite genome indexes using Bowtie2 1.5hr 
$ ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_genome_preparation --bowtie2 ./

# genome length info: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt
#Converting Genome Coordinates From One Genome Version To Another: https://www.biostars.org/p/65558/

######
# Alignment and Mapping paired-end reads
# Using the Bismark read aligner 
#Each sample will take many hours 1.5hr
$ cd ../fastq/
$ nohup ~/Documents/bioinfo_software/Bismark-0.24.0/bismark --multicore 6 --bowtie2 --bam ../reference_GRCm38_mm10/ -1 SRR5195666.sralite.1_1.clock_UMI.R1_val_1.fq -2 SRR5195666.sralite.1_2.clock_UMI.R2_val_2.fq  2>error.log &
  
# raw reads1,2. 4G each.
# trim reads1,2. 2.6G each.
# 1.1GB, SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.bam


######
# deduplicated with UmiBam in --double_umi mode
# UmiBam: https://github.com/FelixKrueger/Umi-Grinder
# it's a perl script. download: https://github.com/FelixKrueger/Umi-Grinder

#from paper: These UMI-deduplicated BAM files were then further processed with the Bismark Methylation Extractor (default parameters) to yield Bismark coverage files.

$ nohup perl ../UmiBam.pl -p --bam --double_umi SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.bam &
# output: SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam

# view alignment in terminal(http://www.htslib.org/doc/samtools-view.html)
$samtools view -f 2 SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam | less -SN 

#######
# Post-alignment step, 20min 
# Extracting methylation calls
# $ nohup ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_methylation_extractor -p --no_overlap --comprehensive  --multicore 4 --buffer_size 5G --bedGraph --counts --gzip  SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam 2>error.log &

# Assessing the alignment
# $ for f in `cat ../extdata/sraFiles.txt`; do awk -F"\t" '$1 == "22" { print $0 }' \ 
#  $f\_1_val_1.fq_bismark_bt2_pe.bismark.cov > $f.chr22.cov; done

# $ awk -F"\t" '$1=="22" {print $0}' SRR949211.sralite.1_1_val_1_bismark_bt2_pe.bismark.cov >SRR949211.chr22.cov 

methy.calls=as.data.frame(data.table::fread('fastq_SRR5195666/SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov.gz'))
head(methy.calls) 
#The coverage output looks like this
#<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
colnames(methy.calls)=c('chr','start','end','percentage','methy','nonmethy')
methy.calls$coverage=methy.calls$methy+methy.calls$nonmethy
summary(methy.calls$coverage)
methy.pick<-methy.calls[methy.calls$coverage>=10,]
dim(methy.calls) #4129321
dim(methy.pick)#2216197

#######
# to use WSH, https://github.com/MPIIComputationalEpigenetics/WSHPackage/blob/master/vignettes/WSH.md
# need to index bam file first, for example 
# in my case, need to sort then index
$ samtools sort SRR5195666.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o SRR5195666_sorted.bam
$ samtools index -b SRR5195666_sorted.bam 
$ samtools view -f 2 SRR5195666_sorted.bam | less -SN 

if(F){
  example.bam <- system.file(file.path("extData","small_example.bam"),
                             package="WSH")
  example.GRanges <- GRanges(Rle(rep("chr2",10)),
                             IRanges(start=c(2298361,2298554,2298732,
                                             2298743,2298787,2298792,
                                             2298827,2298884,2298915,2298921),
                                     end=c(2298361,2298554,2298732,
                                           2298743,2298787,2298792,
                                           2298827,2298884,2298915,2298921)+1))
  pdr <- compute.score(bam.file=example.bam,example.GRanges,score="pdr")
  # if you check "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/WSH/extData/small_example.bam",
  # there is a 'small_example.bam.bai' file.
  $ samtools view small_example.bam | less -SN
  
  # calcualte PDR
  library(WSH)
  input.bam='./fastq/SRR5195666_sorted.bam'
  # get stubbs clock info: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1203-5/MediaObjects/13059_2017_1203_MOESM6_ESM.xls
  clock=readxl::read_excel('13059_2017_1203_MOESM6_ESM.xls',skip=3)
  clock[1:10,]
  
  #example.GRanges <- GRanges(Rle(clock$Chromosome[1:10]),IRanges(start=clock$Start[1:10],end=clock$Start[1:10]+1))
  # Error in { : task 1 failed - "subscript out of bounds"
  example.GRanges <- GRanges(Rle(rep("chr1",2)),
                             IRanges(start=c(3085153,3091103),
                                     end=c(3085153,3091103)+1))
  set.option(mapq.filter = 20)
  pdr <- compute.score(bam.file=input.bam,example.GRanges,score="pdr")
  pdr
  
  #one of the criteria implemented in the PDR, each read has to contain at least 4 CpGs to
  # be used for the classification into discordant/concordant
  library(AnnotationHub)
  ahub <- AnnotationHub()
  ahub=ahub[ahub$species=='Mus musculus',]
  str(ahub) #1649
  ahub[grep('CpG',ahub$title,ignore.case = T),] #
  
  ahub["AH6117"] #$genome: mm9, different version compared with the mm39 reference I used in mapping.
  # that's why there is a error message for `findOverlaps` in `calculate.pdr`
  cpgs <- ahub[["AH6117"]]
  cpgs #16026 ranges
}

library(WSH)
input.bam='./fastq_SRR5195666/SRR5195666_sorted.bam'

###### set up example.GRanges
## build CpG from genome with annotatr
#https://rdrr.io/bioc/annotatr/man/build_annotations.html
#https://rdrr.io/bioc/annotatr/man/annotations.html
library(annotatr)
builtin_genomes() #mm10
builtin_annotations()[grep('^mm10',builtin_annotations())]
annots_gr = build_annotations(genome = 'mm10', annotations = 'mm10_cpgs')
seqlengths(annots_gr)

levels(annots_gr@seqnames)

newStyle <- mapSeqlevels(seqlevels(annots_gr),'UCSC')
newStyle <- newStyle[!is.na(newStyle)]
newStyle
annots_gr_new=annots_gr[seqnames(annots_gr) %in% newStyle,]
annots_gr_new@seqnames=droplevels(annots_gr_new@seqnames)
levels(annots_gr_new@seqnames)

annots_gr_new@seqinfo<-seqinfo(annots_gr_new)[newStyle] #https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomeInfoDb/html/Seqinfo-class.html

#cpgs_chr=annots_gr_new
#cpgs_chr=annots_gr_new[annots_gr_new@seqnames %in% c('chr1','chr3'),]
cpgs_chr=annots_gr_new[annots_gr_new@seqnames %in% c('chr3'),] #3821 ranges
table(cpgs_chr$type)
tmp <- cpgs_chr[cpgs_chr$type=='mm10_cpg_islands',] #686 ranges 
# change ranges into one CpG site format
example.GRanges <- GRanges(seqnames(tmp),
                           IRanges(start=start(tmp),
                                   end=start(tmp)+1))
example.GRanges

###### set up example.GRanges
## or use clock site 
# get stubbs clock info: https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1203-5/MediaObjects/13059_2017_1203_MOESM6_ESM.xls
clock=readxl::read_excel('13059_2017_1203_MOESM6_ESM.xls',skip=3)
clock[1:10,]
table(clock$Chromosome)
library(dplyr)
clock<-clock %>% arrange(Chromosome,Start) #need to sort before calculating score 
example.GRanges <- GRanges(Rle(rep("chr11",31)),
                           IRanges(start=clock[clock$Chromosome==11,]$Start,
                                   end=clock[clock$Chromosome==11,]$Start+1))

example.GRanges <- GRanges(Rle(paste0('chr',clock$Chromosome)),
                           IRanges(start=clock$Start,
                                   end=clock$Start+1))
example.GRanges #329 site

###### set up example.GRanges
## use methylation mapping file, coverage>=10 site
# and intersected with genome annotated CpG

cpgs <- unlist(rnb.get.annotation("CpG",assembly="mm10")) #43735100
pos.cpgs=start(cpgs)
chr.cpgs=seqnames(cpgs)

unique(methy.pick$chr) #only look at autosomal genes (mt, X, Y removed)
methy.auto.pick<-methy.pick[!is.na(as.numeric(methy.pick$chr)),]
dim(methy.auto.pick) #1017542
table(methy.auto.pick$chr)

methy.range <- GRanges(Rle(paste0('chr',methy.auto.pick$chr)),
                       IRanges(start=methy.auto.pick$start,
                               end=methy.auto.pick$start+1))
example.GRanges<-GenomicRanges::intersect(methy.range,cpgs,ignore.strand = TRUE)

example.GRanges #1273637

################################################################################
## compute with already implemented functions
#pdr <- compute.score(bam.file=input.bam,example.GRanges,score="pdr")
#pdr


source('~/Documents/aging_RRBS/WSHPackage-master/R/main.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/utilities.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/calculate_scores.R')

set.option(mapq.filter = 20)

# from calculate_scores.R script
system.time( {
  pdr=calculate.pdr(bam.file=input.bam, anno=example.GRanges,
                    cores = 8,use.sex.chromosomes=FALSE,ignore.strand=TRUE) 
}) #20min
summary(pdr$PDR)
tmp=pdr[!is.na(pdr$PDR),] #543971
sum(tmp$PDR>0) #381224
saveRDS(pdr,'fastq_SRR5195666/SRR5195666_pdr.rds')

set.option(fdrp.type='FDRP')
system.time( {
  fdrp.out=calculate.fdrp.score(bam.file=input.bam, anno=example.GRanges,
                                cores=8,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
}) #30hr
#summary(fdrp.out$FDRP) #if no mapped read or mapped.read<2, return NA for one CpG site
#fdrp.out[which(fdrp.out$FDRP>0),]
#     user    system   elapsed 
# 57.998    88.213 28608.796 , 8hr for 718259 site
saveRDS(fdrp.out,'fastq_SRR5195666/SRR5195666_fdrp.rds')

set.option(fdrp.type='qFDRP')
system.time( {
  pfdrp.out=calculate.qfdrp(bam.file=input.bam, anno=example.GRanges,
                            cores=8,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
}) #30hr
#user    system   elapsed 
#58.654    78.972 28798.863 , 8hr for 718259 site
saveRDS(pfdrp.out,'fastq_SRR5195666/SRR5195666_pfdrp.rds')


fdrp.out=readRDS('fastq_SRR5195666/SRR5195666_fdrp.rds')
pfdrp.out=readRDS('fastq_SRR5195666/SRR5195666_pfdrp.rds')

cor(pfdrp.out$qFDRP,fdrp.out$FDRP,use = "complete") #0.96

#################################################################
## get mapped genome reference info (chr length)
which = GRanges('1:1-100')
param <- ScanBamParam(which=which,what="seq",mapqFilter=20,
                      flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
reads <- readGAlignments(input.bam,param=param) #only read in one subset of all mapped bam
reads # 1224960  alignments
GenomeInfoDb::seqlengths(reads) #chrosomoe length information

#################################################################
## play with functions
cores=2
window.size=unname(get.option('window.size'))
max.reads=unname(get.option('max.reads'))
mapq.filter=unname(get.option('mapq.filter'))
coverage.threshold=unname(get.option('coverage.threshold'))
use.sex.chromosomes=FALSE
ignore.strand=TRUE
#set.option(mapq.filter = 0)
#set.option(coverage.threshold=1)

bam.file=input.bam; anno = example.GRanges; 
output.frame <- data.frame(chromosome=seqnames(anno),start=start(anno),end=end(anno))
bam <- BamFile(bam.file)

source('~/Documents/aging_RRBS/WSHPackage-master/R/utilities.R')
source('~/Documents/aging_RRBS/WSHPackage-master/R/calculate_scores.R')
anno_by_chr <- prepare.annotation(anno,use.sex.chromosomes=use.sex.chromosomes) #separate into different chrs
length(anno_by_chr)
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
reads # 1224960  alignments

#reads0<-readGAlignments(input.bam)
#reads0 #14008968 alignments 

range_reads <- GRanges(reads)
rm(reads)
newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
newStyle <- newStyle[!is.na(newStyle)]
range_reads <- renameSeqlevels(range_reads,newStyle)

# we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
range_reads #1224960 query read
anno #686 ranges
seqlengths(range_reads) #this information is in bam file, reference chr length(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt)
seqlengths(anno) #a subset of above chr names
seqlengths(annots_gr_new) #make sure the length match, for genome ref consistency
overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
overlap
#queryLength: 1224960 / subjectLength: 686

max(queryHits(overlap)) #1224960
length(unique(queryHits(overlap))) #111591 reads mapped to the 686 ranges

query <- queryHits(overlap)
query <- unique(query)
range_reads <- range_reads[query]
length(range_reads) #111591, reads mapped to the CpG island ranges

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
range_reads #111591 reads

overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
overlap #length 1229426
#queryLength: 111591 / subjectLength: 686
table(subjectHits(overlap)) # the number of reads mapped per region (686 ranges in total)

match_read_cpg <- as(overlap,"list")
length(match_read_cpg) #111591, one read per list element, showing the hitted range
head(which(sapply(match_read_cpg,length)==2))

overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
fdrps <- as.list(rep(NA,length(anno)))

#rm(anno)
# for each read we convert the covered CpG sites into a custom representation
read_representation <- as.list(1:length(range_reads))
read_representation <- lapply(read_representation,toCpGs,
                              match_read_cpg, starts_cpgs, starts_reads,seqs_reads)

toCpGs(1,match_read_cpg, starts_cpgs, starts_reads, seqs_reads)

length(match_read_cpg) #111591
table(sapply(match_read_cpg,length)) #a read need to contain at least 2 CpG site to be called 'concordant' or not.
#1 
#111591 
length(read_representation) #111591
table(sapply(read_representation,length))
#a logical vector for each CpG site in the read
#'						where TRUE represents a methylated CpG and FALSE a non-methylated (or a sequence
#'						variation)

rm(match_read_cpg)
rm(starts_reads)
rm(seqs_reads)
values(range_reads) <- DataFrame(cbind(CpG=read_representation)) #add CpG information
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


