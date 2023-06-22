# adapted from https://github.com/MPIIComputationalEpigenetics/WSHScripts

# if you check  bam file generated via bismark
#$ samtools view -f 2 SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam | less -SN 
# $samtools view -f 2 SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam | head -3 | awk '{print NF}' -
# 16 columns
# the XM info is shown in column 14

# https://support-docs.illumina.com/SW/DRAGEN_v39/Content/SW/DRAGEN/MPipelineBisMeth_fDG.htm
#Using Bismark for Methylation Calling
#The recommended approach to methylation calling is to automate the multiple required alignment runs and 
# add the XM, XR, XG tags.

# https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
# (14) XM-tag (methylation call string)
#https://support-docs.illumina.com/SW/DRAGEN_v40/Content/SW/DRAGEN/MPipelineMethBAM_fDG.htm
# z: no CpG
# Z: yes CpG
# h: no CHH
# H: yes CHH
#####################################################################################################
#' add_bismark_meth_string.R
#' This script adds a methylation call string, as introduced with bismark
#' (http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to any bam file. Currently, this 
#' scripts operates on the mouse reference genome mm10, but can be easily extended to other 
#' reference genomes.

#' Load required libraries
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

if(F){
cmds <- commandArgs(trailingOnly=TRUE)
bam.file <- cmds[1]
store.path <- cmds[2]
store.name <- cmds[3]
cores <- as.numeric(cmds[4])
}
bam.file='./fastq_SRR5195656/SRR5195656_sorted.bam'
store.path <- './test_meth_string'
store.name <- 'SRR5195656_sorted_XM.rds'
cores <- 4

### FUNCTIONS ###
#' create.string
#' This function converts a given sequence into the bisulfite conversion representation. Only if a
#' CG is found at a CpG position, was the corresponding cytosine methylated.
#' @param index	Index of the read to be converted
#' @param lens	Lenghts of the reads
#' @param seqs	Sequences of the reads (as a string)
#' @param posis Positions of the CpGs in the corresponding read 
#' @return	The methylation call string according to the one added by bismark
create.string <- function(index,lens,seqs,posis){
  len <- lens[index]
  ret.string <- rep(".",len)
  seq.string <- seqs[index]
  pos <- posis[[index]]
  for(i in pos){
    if(substr(seq.string,i,i+1)=="CG"){
      ret.string[i] <- "Z"
    }else{
      ret.string[i] <- "z"
    }
  }
  paste0(ret.string,collapse = "")
}

### MAIN ###
bam.file <- BamFile(bam.file)
which.all <- c("1:1-195471971","2:1-182113224","3:1-160039680","4:1-156508116","5:1-151834684","6:1-149736546","7:1-145441459","8:1-129401213","9:1-124595110","10:1-130694993","11:1-122082543","12:1-120129022","13:1-120421639","14:1-124902244","15:1-104043685","16:1-98207768","17:1-94987271","18:1-90702639","19:1-61431566")
return.reads <- c()
if(!file.exists(file.path(store.path,"log"))){
	dir.create(file.path(store.path,"log"))
}
cl <- makeCluster(cores,outfile=file.path(store.path,"log","loginfo.log"))
registerDoParallel(cl)
return.reads <- foreach(which=which.all, .combine="c", .packages=c("RnBeads","GenomicAlignments","rtracklayer","Rsamtools"), .export=c("create.string","bam.file")) %dopar%{
  logger.start(substr(which,1,5))
  which <- GRanges(which)
  param <- ScanBamParam(which=which,what=scanBamWhat())
  reads <- readGAlignments(bam.file,param=param,use.names = T) #only read in `which` region from bam file
  #length(reads) #891180 reads on chr1
  seqlevels(reads) #chr information from bam (the reference genome you used in reads mapping)
  newStyle <- mapSeqlevels(seqlevels(reads),'UCSC')
  newStyle <- newStyle[!is.na(newStyle)]
  reads <- renameSeqlevels(reads,newStyle) #or keepSeqlevels
  #seqlevels(reads)
  
  # only look across genome annotated CpG site
  cpgs <- unlist(rnb.get.annotation("CpG",assembly="mm10")) #43735100
  #seqlengths(cpgs) #check for chr length consistency
  op <- findOverlaps(reads,cpgs)
  op #queryLength: 891180 / subjectLength: 43735100
  op <- as.list(op) # 891180-length list, each <=> one read, #CpG on cpgs hitted
  #table(sapply(op,length)) # 0 hit ~ 15 hit per read
  starts.cpgs <- start(cpgs) #genome location of each CpG
  starts.cpgs <- lapply(op, function(x,starts){starts[x]}, starts.cpgs)
  #head(which(sapply(op,length)>1))
  #starts.cpgs[op[[13]]] #[1] 3014928 3014974, two CpG hits for this read
  
  starts.reads <- start(reads)   #genome mapped pos for each read
  pos.cpgs <- lapply(1:length(starts.reads),function(x,starts.cpgs,starts.reads){
                      starts.cpgs[[x]]-starts.reads[[x]]+1},starts.cpgs,starts.reads)
  #x=1; starts.cpgs[[x]]-starts.reads[[x]]+1 #numeric(0)
  #x=13; starts.cpgs[[x]]-starts.reads[[x]]+1 #[1]  1 47
  
  #names(values(reads)) #sam file specification
  lengths <- values(reads)$qwidth #mapped read length
  seqs <- as.character(values(reads)$seq) #ATTCG read sequence
  xm.tag.string <- lapply(1:length(starts.reads),create.string,lengths,seqs,pos.cpgs)
  #xm.tag.string[[13]]; xm.tag.string[[826240]]
  
  values(reads) <- data.frame(values(reads),"XM"=unlist(xm.tag.string),stringsAsFactors = F)
  #names(values(reads)) #XM is added
  logger.completed()
  reads
}
length(return.reads) #13677986
names(values(return.reads)) #"XM"
values(return.reads)$XM #XM tags

reads <- GAlignmentsList(return.reads)
length(reads)
#reads[[1]]->tmp
#values(tmp)$XM

rm(return.reads)
gc()
reads <- unlist(reads)
#names(values(reads))

#export(reads,file.path(store.path,store.name))
saveRDS(reads,file.path(store.path,store.name),compress = TRUE) #no difference

