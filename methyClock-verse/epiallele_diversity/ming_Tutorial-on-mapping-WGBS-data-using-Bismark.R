#http://statisticalrecipes.blogspot.com/2015/05/tutorial-on-mapping-wgbs-data-using.html
#https://github.com/genomicsclass/colonCancerWGBS/tree/master/scripts
#https://htmlpreview.github.io/?https://github.com/genomicsclass/colonCancerWGBS/blob/master/scripts/createObject.html

#######
# Obtaining sample information from GEO
library("GEOquery")
gse <- getGEO("GSE46644")
e <- gse[[1]]
pd <- pData(e)
# subset for colon samples only
pd <- pd[grepl("Colon", pd$title),]

#######
#The information which connects the sample information from GEO with the SRA run id is downloaded 
# from SRA using the Send to: File button. Add the SRA ID to the end of the csv file name.
# `Send to-> File -> Format(RunInfo)`
srr <- data.table::fread("extdata/SraRunInfo.csv")
srrsmall <- srr[,c("SampleName", "Run", "Experiment", "Sample","BioSample",  "avgLength", "download_path")]
colnames(srrsmall)[which(colnames(srrsmall) == "SampleName")] <- "geo_accession"

# merge the GEO phenotypic information table and the sample information from SRA 
target <- merge(pd, srrsmall, by ="geo_accession")
rownames(target) <- target$Run

#The SRA names and SRA file paths were saved to help extract the SRA files from NCBI.
write.table(target$Run, file = "extdata/sraFiles.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(target$download_path, file = "extdata/sraFilesPath.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)

#######
# recommend only download SRR949211.sralite.1 to save time and diskspace
#Obtaining FASTQ files from SRA
#Downloading all .sra files in the sraFilesPath.txt:
# $ cd <path_to>/sra
# $ cat ../extdata/sraFilesPath.txt | wget -i -

#Extracting .fastq files
#$ brew install parallel
# install fastq-dump from: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
#$ cd fastq/
#$ cat ../extdata/sraFiles.txt | parallel -j 6 ~/Documents/bioinfo_software/sratoolkit.3.0.1-mac64/bin/fastq-dump -I --split-files ../sra/{}.sralite.1 
#($~/Documents/bioinfo_software/sratoolkit.3.0.1-mac64/bin/fastq-dump -I --split-files ../sra/SRR949210.sralite.1 )
# after extraction, two files each ~60GB

#######
# to use Trim_Galore, 
# require cutadapt and fastqc (check out README.txt file)

# install cutadapt 
# $pip3 install cutadapt

# install fastqc
# java version: fastqc from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# or command line: $conda install -c bioconda fastqc
#$ fastqc *.fastq
# Look at QC report produced by fastqc: fastqc_report.html

# $cutadapt --version
# 4.3
# $fastqc -v
# FastQC v0.12.1

#######
# Adapter and quality trimming
# install Trim_Galore from: https://www.bioinformatics.babraham.ac.uk/projects/download.html#trim_galore
# user guide: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

#$ cd fastq;
#$ fastq $~/Documents/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --fastqc SRR949211.sralite.1_1.fastq SRR949211.sralite.1_2.fastq 

#######
# reference setup 
# First download the human genome file (.fa) from ENSEMBL.
# $ cd reference
# $ wget ftp://ftp.ensembl.org/pub/release-79/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# $ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Building the bisulfite genome indexes using Bowtie2,take ~1.5hr 
# $ cd reference 
# $ ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_genome_preparation --bowtie2 ./
 
####### 
# Alignment and Mapping paired-end reads
# Using the Bismark read aligner 8hr
#Each sample will take many hours (~10-20 hrs depending on the file size)
# $ cd ../fastq/
# fastq $nohup ~/Documents/bioinfo_software/Bismark-0.24.0/bismark --multicore 6 --bowtie2 --bam ../reference/ -1 SRR949211.sralite.1_1_val_1.fq -2 SRR949211.sralite.1_2_val_2.fq  2>error.log &
# raw reads1,2. 57G each.
# trim reads1,2. 45G each.
# 20GB, SRR949211.sralite.1_1_val_1_bismark_bt2_pe.bam
# report 2k, SRR949211.sralite.1_1_val_1_bismark_bt2_PE_report.txt

#######
# Post-alignment steps (3.5hr)
# Extracting methylation calls
# $ nohup ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_methylation_extractor -p --no_overlap --comprehensive  --multicore 4 --buffer_size 5G --bedGraph --counts --gzip  SRR949211.sralite.1_1_val_1_bismark_bt2_pe.bam 2>error.log &
# Assessing the alignment
# $ for f in `cat ../extdata/sraFiles.txt`; do awk -F"\t" '$1 == "22" { print $0 }' \ 
#  $f\_1_val_1.fq_bismark_bt2_pe.bismark.cov > $f.chr22.cov; done

# $ awk -F"\t" '$1=="22" {print $0}' SRR949211.sralite.1_1_val_1_bismark_bt2_pe.bismark.cov >SRR949211.chr22.cov 



