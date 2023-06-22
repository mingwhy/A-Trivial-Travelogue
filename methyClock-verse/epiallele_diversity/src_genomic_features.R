
#https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html
library(GenomicRanges) #must load GenomicRanges before plyranges on my mac
library(plyranges)
#remotes::install_github("cnobles/gintools")
library(gintools) #remove redundant chr pos
library(ggplot2)

#################################################################
# https://rockefelleruniversity.github.io/Bioconductor_Introduction/exercises/answers/GenomicFeatures_answers.html
# Genomic Features in Bioconductor exercises

library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
z=keys(TxDb.Mmusculus.UCSC.mm10.knownGene)
length(z) #24594 NCBI Gene IDs #https://support.bioconductor.org/p/9138987/

allGenes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene,single.strand.genes.only=FALSE)
length(allGenes) #24594 

geneTable <- as.data.frame(allGenes)
filtGeneTable <- geneTable[geneTable$seqnames %in% c("chr1","chr2","chr3"),]

allTranscripts <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
transcriptTable <- as.data.frame(allTranscripts)
filtTranscriptTable <- transcriptTable[transcriptTable$seqnames %in% c("chr1"),]
txLengths <- transcriptLengths(TxDb.Mmusculus.UCSC.mm10.knownGene)
toPlotTxLen <- merge(filtTranscriptTable,txLengths,all=FALSE,by="tx_name")
ggplot(toPlotTxLen,aes(x=width,y=tx_len,size=nexon))+geom_point()+scale_x_log10()+scale_y_log10()+theme_minimal()

##Plot the nucleotide content for the longest transcript (by genomic size including introns) overlapping the region (chr1:3,083,473-3,876,510.) using ggplot2
allTX <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
myRanges <- GRanges("chr1",IRanges(3083473,3876510))
filtTX <- allTX[allTX %over% myRanges]
orderedFiltTx <- filtTX[order(width(filtTX),decreasing = TRUE),]
longestTX <- orderedFiltTx[1]

allSeqs <- extractTranscriptSeqs(BSgenome.Mmusculus.UCSC.mm10,TxDb.Mmusculus.UCSC.mm10.knownGene)
filtSeq <- allSeqs[longestTX$tx_id]
alpFreq <-alphabetFrequency(filtSeq)
atcgFreq <- alpFreq[,c("A","T","C","G")]
atcgFreqDF <- data.frame(Bases=names(atcgFreq),Frequency=atcgFreq)

ggplot(atcgFreqDF,aes(y=Frequency,x=Bases,fill=Bases))+geom_bar(stat = "identity")+theme_minimal()

#############################################################################
## what's the R code for __ regions of mouse genome? 
# Load the annotation package
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(txdb)
# Load the genome sequence for the mouse genome
library(BSgenome.Mmusculus.UCSC.mm10)
mm10_seq <- BSgenome.Mmusculus.UCSC.mm10

# Get the `promoter` regions for all genes in the annotation
# The upstream and downstream parameters specify how many base pairs upstream 
# and downstream of the transcription start site (TSS) should be included in the promoter region.
promoters <- promoters(txdb, upstream = 2000, downstream = 200)
head(promoters, n = 10) # a GRanges object

# Retrieve the `exonic` ranges for all genes in the annotation
exons <- exonsBy(txdb, by="gene")
#GRList <- transcriptsBy(txdb, by = "gene") 
#GRList <- cdsBy(txdb, by = "gene") 

# Get the DNA sequence for the exonic regions
length(exons)
exon_seqs <- getSeq(mm10_seq, exons)
length(exon_seqs) #24594 genes

# Print the first 10 exonic sequences
for (i in 1:10) {
  cat(paste0("Exon ", i, " sequence:\n"))
  cat(exon_seqs[[i]], "\n")
}

# get `gene body` or gene ranges
genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene,single.strand.genes.only=FALSE)
geneTable <- as.data.frame(genes)
head(geneTable) #group_name: NCBI Gene ID

# get `intron` regions
introns <- setdiff(genes, exons) #same length of genes and exons 
intron_seqs <- getSeq(mm10_seq, introns)

# get `intergenic` regions
#https://stackoverflow.com/questions/29253412/finding-intergenic-regions
genic <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
genic <- reduce(genic, ignore.strand=T)
intergenic <- gaps(genic)
intergenic <- intergenic[strand(intergenic) == "*"] #This is important!!!

# get `LTR` and `LINE` regions
#https://bioconductor.org/packages/release/data/annotation/vignettes/UCSCRepeatMasker/inst/doc/UCSCRepeatMasker.html
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("RepeatMasker", "Mus musculus"))
# AH6075 AH6155 AH6194 AH6224
rmsk_mm10 <- ah[["AH6224"]] #only chr range locations, but no LTR or LINE, repeat type annotation
metadata(rmsk_mm10) 


#################################################################
## translate between gene id and gene symbols
(mouse_gene_ids=promoters$tx_name[1:10])
mouse_gene_ids=gsub('\\.1$','',mouse_gene_ids)
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
searchFilters(m_ensembl, "exon") #ensembl_exon_id Exon ID(s) [e.g. ENSMUSE00000097910]
searchFilters(m_ensembl, "transcript") # ensembl_transcript_id    Transcript stable ID(s) [e.g. ENSMUST00000000001]
foo <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
             filters = 'ensembl_transcript_id', values = mouse_gene_ids, mart = ensembl)
foo

#################################################################
## what's the R code for CpG islands/shelf/shore of mouse genome? 
#https://rdrr.io/bioc/annotatr/man/build_annotations.html
#https://rdrr.io/bioc/annotatr/man/annotations.html
library(annotatr)
builtin_genomes() #mm10
builtin_annotations()[grep('^mm10',builtin_annotations())]
annots_gr = build_annotations(genome = 'mm10', annotations = 'mm10_cpgs')

levels(annots_gr@seqnames) #chr names

newStyle <- mapSeqlevels(seqlevels(annots_gr),'UCSC')
newStyle <- newStyle[!is.na(newStyle)]
newStyle
annots_gr_new=annots_gr[seqnames(annots_gr) %in% newStyle,]
annots_gr_new@seqnames=droplevels(annots_gr_new@seqnames)
levels(annots_gr_new@seqnames)

annots_gr_new@seqinfo<-seqinfo(annots_gr_new)[newStyle] #https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomeInfoDb/html/Seqinfo-class.html
table(annots_gr_new$type)
summary(ranges(annots_gr_new)@width)

cpgs_ranges=annots_gr_new

######################################################
## get individual CpG genomic positions
library(RnBeads) #http://www.bioconductor.org/packages/devel/bioc/vignettes/RnBeads/inst/doc/RnBeads_Annotations.pdf
library(RnBeads.mm10)
# all CpGs in the mouse genome
all.cpgs <- unlist(rnb.get.annotation("CpG",assembly="mm10"))  
#rnb.get.annotation("promoters",assembly="mm10")
length(all.cpgs) #43735100, per site, each strand count once
#1.1 Site Annotation Tables, https://rnbeads.org/data/RnBeadsAnnotationCreator.pdf
#CpG, #CpG in a 100bp window, site centered
#GC, GC content in a 100 bp window
#CGI Relation, "Open Sea", "Shelf", "Shore" and "Island".
#Genomic coordinates are 1-based and sites are sorted based on their genomic position.
seqlevels(all.cpgs)
strand(all.cpgs)
#cpgs.pos<-paste(seqnames(all.cpgs),start(all.cpgs),sep=':')
#length(unique(cpgs.pos))
#all.cpgs$pos=cpgs.pos;

#################################################################
## what's the R code for calculating CpG density of mouse genome?
#BiocManager::install("Repitools")
library(Repitools)
library(BSgenome.Mmusculus.UCSC.mm10)
#https://rdrr.io/bioc/Repitools/man/cpgDensityCalc.html
TSSTable <- data.frame(chr = c("chr1", "chr2"), position = c(100000, 200000))
cpgDensityCalc(TSSTable, organism = Mmusculus, window = 600)
#################################################################

