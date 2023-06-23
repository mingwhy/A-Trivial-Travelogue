# http://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# create folder `methylationArrayAnalysis`, 
# download data from https://figshare.com/articles/%0AmethylAnalysisDataV3_tar_gz/4800970
# following https://f1000research.com/articles/5-1281
# set up a path to the data directory
#dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
dataDirectory='methylationArrayAnalysis/'
# list the files
list.files(dataDirectory, recursive = TRUE)

# get the 450k annotation data
packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#[1] '0.6.0'
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
dim(ann450k) #485512     33

#Infinium Methylation450K manifest column headings
#https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001567
x=ann450k$UCSC_RefGene_Group
x=x[x!='']
length(x)
strsplit(x[1],';')
x1=lapply(x,strsplit,';')
table(unlist(x1))
#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200 
#53840   29937  103241  288345  122745   89029

#Relation_to_UCSC_CpG_Island: The location of the CpG relative to the CpG island.
#Shore = 0–2 kb from island.
#Shelf = 2–4 kb from island.
#N = upstream (5’) of CpG island.
#S = downstream (3’) of CpG island.
#a cartoon figure: https://www.bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
table(ann450k$Relation_to_Island)

probe.location=getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what = "Locations")
cgi_anno = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, "Islands.UCSC")
table(cgi_anno$Relation_to_Island)

# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targets
dim(targets) #11 x 9

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet #dim: 622399 11 

# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

## Quality control
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
dim(detP) #485512     11

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")

# generate a `qcReport.pdf`
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet #dim 622399 10 , one poor quality smaple removed

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

## Normalisation
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 
# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))

## Data exploration
# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

## Filtering
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
