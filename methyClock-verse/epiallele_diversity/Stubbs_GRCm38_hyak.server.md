#Stubbs, Thomas M., et al. "Multi-tissue DNA methylation age predictor in mouse." Genome biology 18.1 (2017): 1-14.
#GSE93957: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93957
#SraRunInfo.csv: SRR___ run id <=> GSMxx id on GSE 

# Obtaining sample information from GEO
#in R
#if error: https://github.com/seandavi/GEOquery/issues/110
#BiocManager::install('seandavi/GEOquery')
library("GEOquery")
if(F){ #online
  gse <- getGEO("GSE93957") #GSE93957_series_matrix.txt.gz
  e <- gse[[1]] #62 samples
  pd <- pData(e) #62 x 41 df
  head(pd$title)
}
# read from local downloaded file
gse <- getGEO(filename='Stubbs_GSE93957_series_matrix.txt.gz',getGPL = FALSE)
pd <- pData(gse) #62 x 41 df
head(pd$title)

unique(pd$`tissue:ch1`)
#"Cortex" "Heart"  "Liver"  "Lung"  
#subset for heart and lung samples only
pd <- pd[grepl("Heart|Lung", pd$title),] #31 samples 
table(pd$`age:ch1`,pd$`tissue:ch1`)

#on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93957, click `SRA`,
#go to SRA(https://www.ncbi.nlm.nih.gov/sra?term=SRP097629), use the `Send to-> File -> Format(RunInfo)` to get SRA run ID.
srr <- data.table::fread("sra_data//SraRunInfo.csv")
srrsmall <- srr[,c("SampleName", "Run", "Experiment", "Sample","BioSample",  "avgLength", "download_path")]
colnames(srrsmall)[which(colnames(srrsmall) == "SampleName")] <- "geo_accession"

#merge the GEO phenotypic information table and the sample information from SRA 
target <- merge(pd, srrsmall, by ="geo_accession")
rownames(target) <- target$Run

#The SRA names and SRA file paths were saved to help extract the SRA files from NCBI.
write.table(target$Run, file = "sra_data/sraFiles.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(target$download_path, file = "sra_data/sraFilesPath.txt", quote= FALSE,row.names = FALSE, col.names = FALSE)

# Obtaining FASTQ files from SRA
#Downloading all .sra files in the sraFilesPath.txt:
$ cd sra_data
$ cat sraFilesPath.txt | wget -i -
$ rename -- sralite.1 sra SRR*

#store data in /gscratch/scrubbed/mingy16
$ du -h ./ --max-depth=2 #check storage space
$ pwd
/gscratch/scrubbed/mingy16/GSE93957

#batch rename SRR5195656.sralite.1 to SRR5195656.sra
#rename from util-linux
`$ rename -- sralite.1 sra *fastq`

# Extracting .fastq files
## local
$ brew install parallel
#install fastq-dump from: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
$~/Documents/bioinfo_software/sratoolkit.3.0.1-mac64/bin/fastq-dump -I --split-files ../sra_data/SRR5195656.sra

## server
$ module avail | grep 'fastq'
$ module load contrib/parallel-fastq-dump/0.6.3
$ which parallel
/sw/contrib/csde/parallel/20171122/bin/parallel

$ mkdir fastq && cd $_
$ squeue -p csde #check server computing node resource
### one sample test
$ /gscratch/csde-promislow/bioinfo_software/sratoolkit.3.0.1-ubuntu64/bin/fastq-dump -I \
        --split-files ../sra_data/SRR5195656.sra 

### run batch
$ cat process_fastq.slurm 
`#!/bin/bash
#SBATCH --time=06:00:00
cd /gscratch/scrubbed/mingy16/GSE93957/fastq
cat ../sra_data/sraFiles.txt | parallel -j 20 /gscratch/csde-promislow/bioinfo_software/sratoolkit.3.0.1-ubuntu64/bin/fastq-dump \
    -I --split-files ../sra_data/{}.sra `

$ sbatch -p csde -A csde process_fastq.slurm 
#on server: 10min `parallel -j 20`
#SRR5195656.sra, 1.2GB. after extraction, two files each 9GB
$ du -h -d1 #all 31 -> 62 fastq files
525G  

# to use Trim_Galore, 
require cutadapt and fastqc (check out README.txt file)

## local
#install cutadapt 
$ pip3 install cutadapt
#install fastqc
#java version: fastqc from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#or command line: 
`$ conda install -c bioconda fastqc
$ fastqc *.fastq` #Look at QC report produced by fastqc: fastqc_report.html

$cutadapt --version
#4.3
$fastqc -v
#FastQC v0.12.1

## server
$ module avail | grep 'cutadapt'
contrib/cutadapt/1.15  # version too low, install a more recent version myself                                              
$ module avail | grep 'fastqc'
contrib/fastqc/0.11.5  

#install  on server with my own anaconda
$ srun -p build --time=3:00:00 --mem=10G --pty /bin/bash 
$ which conda
/gscratch/csde-promislow/anaconda3/condabin/conda
$ conda install -c bioconda cutadapt
$ which cutadapt
/gscratch/csde-promislow/anaconda3/bin/cutadapt
$ cutadapt --version
4.3

#Install Trim Galore 0.6.10 (https://wiki.cac.washington.edu/display/hyakusers/Hyak+installing+open+source+software)
$curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
$tar xvzf trim_galore.tar.gz
$ ./TrimGalore-0.6.10/trim_galore 

# Read trimming
#following protocol: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
#Trim_Galore_User_Guide.md

## local
#Specialised Trimming, Mouse Epigenetic Clock trimming, use option `--clock`
#trim barcode, following Trim_Galore_User_Guide.md, 10min
$~/Documents/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --clock --fastqc SRR5195656_1.fastq SRR5195656_2.fastq
#trim 3prime, following Trim_Galore_User_Guide.md, 5min
$~/Documents/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 SRR5195656_1.clock_UMI.R1.fq SRR5195656_2.clock_UMI.R2.fq

## server (https://htmlpreview.github.io/?https://github.com/genomicsclass/colonCancerWGBS/blob/master/scripts/createObject.html)
$ srun -p csde -A csde --time=6:00:00 --mem=120G --pty /bin/bash
$ cd fastq;
$ /gscratch/csde-promislow/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired \
   --clock --fastqc SRR5195656.sra_1.fastq SRR5195656.sra_2.fastq
$ rename -- .sra '' *fastq # to get SRR5195656_2.fastq file names

### trim barcode (1.5hr)
$ module load contrib/fastqc/0.11.5
$ cat process_fastq.slurm 
`#!/bin/bash
#SBATCH --time=06:00:00
cd /gscratch/scrubbed/mingy16/GSE93957/fastq
module load contrib/fastqc/0.11.5
cat ../sra_data/sraFiles.txt | parallel -j 20 /gscratch/csde-promislow/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --clock --fastqc {}\_1.fastq {}\_2.fastq`

$ sbatch -p csde -A csde process_fastq.slurm 

output SRR5195656_1.clock_UMI.R1.fq file

### trim 3prime (30min)
cat ../sra_data/sraFiles.txt | parallel -j 20 /gscratch/csde-promislow/bioinfo_software/TrimGalore-0.6.10/trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 {}\_1.clock_UMI.R1.fq {}\_2.clock_UMI.R2.fq

output SRR5195656_1.clock_UMI.R1_val_1.fq file

`if run all above sequentially, 10+11+10 runs, time cost 15min`

# reference setup  
## find out which ensembl version release <=> GRCm38 mm10
#GRCm38 mm10: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/
#go to ensembl: https://useast.ensembl.org/info/data/ftp/index.html
#click mouse, on Mouse page, `other reference assemblies GRCm38(Ensembl release 102)`, click `Go` 
#on `Mouse(GRCm38.p6)` page, you can see the version: Genome assembly: GRCm38.p6 (GCA_000001635.8) 
#Download FASTA files for genes, cDNAs, ncRNA, proteins
#a ftp finder link is opened, enter 'dna/' folder, download `Mus_musculus.GRCm38.dna.primary_assembly.fa.gz`
$ cd reference_GRCm38_mm10;
#$ wget ftp://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
$ gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

## Building the bisulfite genome indexes using Bowtie2, 1hr
#Remember that the indexing is run twice in parallel already (for the top and bottom strand separately), so e.g. '--parallel 4' will use 8 threads in total.
$ ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_genome_preparatio --help 
$ ~/Documents/bioinfo_software/Bismark-0.24.0/bismark_genome_preparation --parallel 4 --bowtie2 ./

## genome length info: https://rgd.mcw.edu/rgdweb/report/genomeInformation/genomeInformation.html?species=Mouse&mapKey=35&details=true
#Converting Genome Coordinates From One Genome Version To Another: https://www.biostars.org/p/65558/


# Alignment and Mapping paired-end reads
## local
#Using the Bismark read aligner 
#This sample will take many hours (30min)
$ cd ../fastq/
$ nohup ~/Documents/bioinfo_software/Bismark-0.24.0/bismark --multicore 6 --bowtie2 --bam ../reference_GRCm38_mm10/ -1 SRR5195666.sralite.1_1.clock_UMI.R1_val_1.fq -2 SRR5195666.sralite.1_2.clock_UMI.R2_val_2.fq  2>error.log &
   
#raw reads1,2. 4G each.
#trim reads1,2. 2.6G each.
#1.1GB, SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.bam

## server
$ module avail | grep 'bismark'
$ contrib/bismark/v0.20.0
$ module load contrib/bismark/v0.20.0

$ module avail | grep 'bowtie'
contrib/bowtie/1.2.2                                                
contrib/bowtie2/2.3.3.1                            
$ module load contrib/bowtie2/2.3.3.1  

$ module avail | grep 'samtool'
$ module load contrib/samtools/1.9 

### one sample test
$ nohup bismark --multicore 8 --bowtie2 --bam ../reference_GRCm38_mm10/ -1 ../fastq/SRR5195656_1.clock_UMI.R1_val_1.fq -2 ../fastq/SRR5195656_2.clock_UMI.R2_val_2.fq  2>error.log &

### run batch
$ cat ../sra_data/sraFiles.txt | parallel -j 20 bismark --multicore 8 --bowtie2 --bam ../reference_GRCm38_mm10/ -1 ../fastq/{}\_1.clock_UMI.R1_val_1.fq -2 ../fastq/{}\_2.clock_UMI.R2_val_2.fq 

$ cat process_mapping.slurm 
`#!/bin/bash
#SBATCH --time=4-24:00:00
#SBATCH --mem=200G
cd /gscratch/scrubbed/mingy16/GSE93957/mapping
module load contrib/bismark/v0.20.0
module load contrib/bowtie2/2.3.3.1  
module load contrib/samtools/1.9 
module load contrib/fastqc/0.11.5
cat ../sra_data/sraFiles.txt | parallel -j 2 bismark --multicore 8 --bowtie2 --bam ../reference_GRCm38_mm10/ -1 ../fastq/{}\_1.clock_UMI.R1_val_1.fq -2 ../fastq/{}\_2.clock_UMI.R2_val_2.fq   `

(as each sample uses 8 core in bowtie2 mapping, run 2 samples in parallel, if 3, hyak node failed to handle)
$ sbatch -p csde -A csde process_mapping.slurm # 17hr, 77GB
(log in via: $ ssh mingy16@n2335 to check for node usage)
( if use 3 nodes, each handle 10+10+11 files, takes 7.5 hr) 

# deduplicated with UmiBam in --double_umi mode, 2min
#UmiBam: https://github.com/FelixKrueger/Umi-Grinder
#it's a perl script. download: https://github.com/FelixKrueger/Umi-Grinder

#from paper: These UMI-deduplicated BAM files were then further processed with the Bismark Methylation Extractor (default parameters) to yield Bismark coverage files.

#if on server, load modules first!!
$ nohup perl ../UmiBam.pl -p --bam --double_umi SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.bam 2>error.log &

#server, 40min
$ cat ../sra_data/sraFiles.txt | parallel -j 11 perl ../UmiBam.pl -p --bam --double_umi {}\_1.clock_UMI.R1_val_1_bismark_bt2_pe.bam 
( if use 3 nodes, each handle 10+11+10 files, takes 20min)

#view alignment in terminal(http://www.htslib.org/doc/samtools-view.html)
$ samtools view -f 2 SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam | less -SN 

# sorting bam files in the mapping folder
$ cd mapping;
$ samtools sort SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o SRR5195656_sorted.bam
$ cat ../sra_data/sraFiles.txt | parallel -j 8 samtools sort {}\_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o {}\_sorted.bam

( if use 3 nodes, each handle 10+10+11 files, 30min) 

`if run all above sequentially, 10+11+10 runs, time cost 6hr`


# Extracting methylation calls using unsorted bam!!
https://htmlpreview.github.io/?https://github.com/genomicsclass/colonCancerWGBS/blob/master/scripts/createObject.html
#This might be the result of sorting the paired-end SAM/BAM files by chromosomal position which is not compatible with correct methylation extraction. Please use an unsorted file instead or sort the file using 'samtools sort -n' (by read name). 

$ mkdir methy_extract;
$ nohup bismark_methylation_extractor -p -no_overlap --comprehensive --multicore 4 --bedGraph --buffer_size 20G --counts --gzip ../mapping/SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam 2>methylation_extract.log &

#server, 3.5hr
$ cat ../sra_data/sraFiles.txt | parallel -j 4 bismark_methylation_extractor -p -no_overlap --comprehensive --multicore 4 --bedGraph --buffer_size 20G --counts --gzip ../mapping/{}\_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam 

( if use 3 nodes, each handle 10+10+11 files, 1.5hr) 

# use metheor (written in rust language) to calculate WSH scores
## on server (metheor only supports linux), need sorted.bam files
$ which conda
$ conda create -n 'rust_env'
$ conda activate rust_env
$ conda install -c dohlee metheor #install python packages within this environment
$ metheor --help    #should be implemented already if installed successfully
#https://github.com/dohlee/metheor 
$ metheor pdr --input SRR5195656_sorted.bam --output pdr.tsv #511MB, 1min
$ metheor lpmd --input SRR5195656_sorted.bam --output lpmd.tsv 
$ nohup metheor lpmd --input SRR5195656_sorted.bam --output lpmd.tsv --pairs ？？(The argument '--pairs <PAIRS>' requires a value but none was supplied）
$ metheor mhl --input SRR5195656_sorted.bam --output mhl.tsv #511MB, 3min
$ metheor fdrp --input SRR5195656_sorted.bam --output fdrp.tsv #511MB, 3min
$ metheor qfdrp --input SRR5195656_sorted.bam --output qfdrp.tsv #511MB, 4min

$ srun -p csde -A csde --time=4:00:00 --mem=200G --pty /bin/bash

$ conda activate rust_env
$ mkdir score && cd $_
$ metheor qfdrp --input ../mapping/SRR5195656_sorted.bam --output qfdrp_SRR5195638.tsv

$ vim run_metheor.sh
#! /bin/bash
cat ../sra_data/sraFiles.txt | parallel -j 16 metheor qfdrp --input ../mapping/{}\_sorted.bam --output qfdrp\_{}.tsv
(use nohup .. &, even if you logged off hyak, sh still run, 30min )

$ cat process_metheor.slurm 
#!/bin/bash -l
#SBATCH --time=4-24:00:00
#SBATCH --mem=200G
cd /gscratch/scrubbed/mingy16/GSE93957/score

pwd
conda info --envs
conda activate /gscratch/csde-promislow/anaconda3/envs/rust_env

module load contrib/samtools/1.9 

#cat ../sra_data/sraFiles.txt | parallel -j 16 metheor qfdrp --input ../mapping/{}\_sorted.bam --output qfdrp\_{}.tsv
#cat ../sra_data/sraFiles.txt | parallel -j 16 metheor pdr --input ../mapping/{}\_sorted.bam --output pdr\_{}.tsv
#cat ../sra_data/sraFiles.txt | parallel -j 16 metheor pm --input ../mapping/{}\_sorted.bam --output pm\_{}.tsv
#cat ../sra_data/sraFiles.txt | parallel -j 16 metheor mhl --input ../mapping/{}\_sorted.bam --output mhl\_{}.tsv
cat ../sra_data/sraFiles.txt | parallel -j 16 metheor fdrp --input ../mapping/{}\_sorted.bam --output fdrp\_{}.tsv

(with slurm, 3min for pdr, 5min for pm, 5min for mhl, 20min for fdrp)
(run sequentially, 50min)

$ wc -l *tsv
  1140890 fdrp_SRR5195656.tsv
  1140890 qfdrp_SRR5195656.tsv  
   439066 pdr_SRR5195656.tsv
   437908 mhl_SRR5195656.tsv
   194023 pm_SRR5195656.tsv
$ wc -l *SRR5195674*
  1744425 fdrp_SRR5195674.tsv
  1744425 qfdrp_SRR5195674.tsv       
   745038 pdr_SRR5195674.tsv
   742425 mhl_SRR5195674.tsv
   359419 pm_SRR5195674.tsv
   359419 me_SRR5195674.tsv
  

$ conda deactivate 

######
## mapping results summary
#Number of mapped reads from BAM file
https://www.biostars.org/p/138116/
$ samtools flagstat file.sorted.bam
$ wc -l fq file

$ cd fastq;
$ nohup wc -l *_1.clock_UMI.R1_val_1.fq >../reads_in_fastq.txt &


$ cd mapping; 
$ module load contrib/samtools/1.9 
$ ( for i in mapping/*sorted.bam ; do echo $i; echo 'hello'; done )
$ ( for i in *sorted.bam ; do echo $i; samtools flagstat $i ; done) > ../reads_in_sorted.bam.txt &
https://www.biostars.org/p/413593/


## counting number of CpGs that are methylated and unmethylated in a read of whole genome bisulfite data
https://www.biostars.org/p/384127/

#install MethylDackel: https://github.com/dpryan79/MethylDackel, https://anaconda.org/bioconda/methyldackel
#require htslib and libbigwig
$ conda install -c bioconda libbigwig
$ conda install -c bioconda htslib
$ not use this one, as the version didn't support 'perRead' conda install -c "bioconda/label/main" methyldackel 
$ conda install -c "bioconda/label/main" methyldackel
$ MethylDackel 
MethylDackel: A tool for processing bisulfite sequencing alignments.
Version: 0.5.1 (using HTSlib version 1.9)

$ srun -p csde -A csde  --time=2:00:00 --mem=30G --pty /bin/bash
$ salloc -N 1 -p csde -A csde  --time=2:00:00 --mem=10G
$ nohup MethylDackel perRead ../GSE93957/reference_GRCm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa ../GSE93957/mapping/SRR5195656_sorted.bam >test.out &

$ module load contrib/samtools/1.9 
$ samtools view -f 2 ../GSE93957/mapping/SRR5195656_sorted.bam | less -SN
$ awk '{if(NR>1)print $5}' test.out | sort | uniq -c >cpg.per.read.count.txt

$ mkdir perRead;
$ cat process_perRead.slurm 
#!/bin/bash
#SBATCH --time=4-24:00:00
#SBATCH --mem=200G
cd /gscratch/scrubbed/mingy16/Cortex_Liver/perRead
module load contrib/samtools/1.9 
cat ../sra_data/sraFiles.txt | parallel -j 8 MethylDackel perRead ../reference_GRCm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa ../mapping/{}_sorted.bam ">" ./{}_perRead.out
({} usage: https://opensource.com/article/18/5/gnu-parallel)
(time cost 5min)


## use MethylDackel to call methylation level
https://github.com/dpryan79/MethylDackel
#MethylDackel can filter reads and bases according to MAPQ and Phred score, respectively. The default minimums are MAPQ >= 10 and Phred >= 5, though these can be changed with the -q and -p options. 
$ MethylDackel extract ../reference_GRCm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa ../mapping/SRR5195656_sorted.bam -q 10 -p 5 -o SRR5195656 &

$ mkdir cpg_level && cd;
$ cat ../sra_data/sraFiles.txt | parallel -j 16 MethylDackel extract ../reference_GRCm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa ../mapping/{}_sorted.bam -q 10 -p 5 -o {} 

(time cost, only takes 10min, way faster then bismark)


## the number of CpG with methylation level and FDRP
https://github.com/FelixKrueger/Bismark/tree/master/Docs
not so sure about the mapping quality filter, Devon Ryan suggested filter out mapq<10 (https://www.biostars.org/p/155605/)

https://github.com/dohlee/metheor
-q, --min-qual: Minimum quality for a read to be considered. [default: 10]
-d, --min-depth: Minimum depth of reads covering epialleles to consider. [default: 10]

#bismark output file: The coverage output looks like this
#<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
$ awk -F"\t" '{print($5+$6)}' SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov
$ awk -F"\t" '{if(($5+$6)>=10)print}' SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov | wc -l
  1027681 CpG sites
$ wc -l fdrp_SRR5195656.tsv
  1140890 fdrp_SRR5195656.tsv


## Extract Reads From A Bam File That Fall Within A Given Region
https://www.biostars.org/p/48719/
#need sorted bam, indexed bai, and reference fasta
#[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.
$ module load contrib/samtools/1.9 
$ samtools view SRR5195656_sorted.bam "12:24498313-24498313" >test.sam
$ samtools view SRR5195656_sorted.bam "1:153890241-153890241" -T ../reference_GRCm38_mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa >test.sam

#in file `SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov` 
1725360 1       153890241       153890241       26.6666666666667        4       11

#$ wc -l test.sam 
30 test.sam #15 paired-end reads


##########################################################################################
##########################################################################################
## Assessing the alignment
$ for f in `cat ../extdata/sraFiles.txt`; do awk -F"\t" '$1 == "22" { print $0 }' $f\_1_val_1.fq_bismark_bt2_pe.bismark.cov > $f.chr22.cov; done
$ awk -F"\t" '$1=="22" {print $0}' SRR949211.sralite.1_1_val_1_bismark_bt2_pe.bismark.cov >SRR949211.chr22.cov 


# read in bismark.cov.gz in R
methy.calls=as.data.frame(data.table::fread('fastq_SRR5195656/SRR5195656.sralite.1_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bismark.cov.gz'))
head(methy.calls) 
#The coverage output looks like this
#<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
colnames(methy.calls)=c('chr','start','end','percentage','methy','nonmethy')
methy.calls$coverage=methy.calls$methy+methy.calls$nonmethy
summary(methy.calls$coverage)
methy.pick<-methy.calls[methy.calls$coverage>=10,]
dim(methy.calls) #3699270
dim(methy.pick)#1027681 

#download bismark.cov.gz files from GSE93957
$ wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE93957&format=file'


######
# WSH in R

## sort and index bam files
$ samtools sort SRR5195656_1.clock_UMI.R1_val_1_bismark_bt2_pe.UMI_deduplicated.bam -o SRR5195656_sorted.bam
$ samtools index -b SRR5195656_sorted.bam 
$ samtools view -f 2 SRR5195656_sorted.bam | less -SN 

#to use WSH, https://github.com/MPIIComputationalEpigenetics/WSHPackage/blob/master/vignettes/WSH.md
#need to sort & index bam file 

## server install WSH
$ which R
~/R-4.1.0/bin/R
$ .libPaths()
[1] "/gscratch/home/mingy16/R-4.1.0/lib64/R/library"
BiocManager::install("RnBeads")
BiocManager::install("RnBeads.hg38")
devtools::install_github("MPIIComputationalEpigenetics/WSHPackage")


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
  input.bam='./fastq/SRR5195656_sorted.bam'

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
  ahub[grep('CpG',ahub$title,ignore.case = T),] 
  ahub["AH6117"] #$genome: mm9, different version compared with the mm39 reference I used in mapping.
  # that's why there is a error message for `findOverlaps` in `calculate.pdr`
  cpgs <- ahub[["AH6117"]]
  cpgs #16026 ranges
}

library(WSH)
input.bam='./fastq_SRR5195656/SRR5195656_sorted.bam'


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

example.GRanges #718259

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
}) #15min
summary(pdr$PDR)
tmp=pdr[!is.na(pdr$PDR),] #543971
sum(tmp$PDR>0) #381224
saveRDS(pdr,'fastq_SRR5195656/SRR5195656_pdr.rds')

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
saveRDS(fdrp.out,'fastq_SRR5195656/SRR5195656_fdrp.rds')

set.option(fdrp.type='qFDRP')
system.time( {
pfdrp.out=calculate.qfdrp(bam.file=input.bam, anno=example.GRanges,
                          cores=8,use.sex.chromosomes=FALSE,ignore.strand=TRUE)
})
pfdrp.out
pfdrp.out$qFDRP
#user    system   elapsed 
#58.654    78.972 28798.863 , 8hr for 718259 site
saveRDS(pfdrp.out,'fastq_SRR5195656/SRR5195656_pfdrp.rds')

pdr=readRDS('fastq_SRR5195656/SRR5195656_pdr.rds')
fdrp.out=readRDS('fastq_SRR5195656/SRR5195656_fdrp.rds')
pfdrp.out=readRDS('fastq_SRR5195656/SRR5195656_pfdrp.rds')

cor(pfdrp.out$qFDRP,fdrp.out$FDRP,use = "complete") #0.96
cor(pfdrp.out$qFDRP,pdr$PDR,use = "complete") #0.8472071
summary(pdr$PDR)
summary(fdrp.out$FDRP)
summary(pfdrp.out$qFDRP)

#################################################################
## get mapped genome reference info (chr length)
which = GRanges('1:1-100')
param <- ScanBamParam(which=which,what="seq",mapqFilter=20,
                      flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
reads <- readGAlignments(input.bam,param=param) #only read in one subset of all mapped bam
reads # 1224960  alignments
GenomeInfoDb::seqlengths(reads) #chrosomoe length information

##############################
## play with calculate_scores.R
bam.file=input.bam; 
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

source('WSHPackage-master/R/utilities.R')
source('WSHPackage-master/R/calculate_scores.R')
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
range_reads # 694176 ranges, inqury reads
match_read_cpg <- as(overlap,"list")
length(match_read_cpg) #694176
table(sapply(match_read_cpg,length))
# contain reads with>=4 CpG ??

overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
pdrs <- as.list(rep(NA,length(anno)))
rm(anno)

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
rm(anno)
# for each read we convert the covered CpG sites into a custom representation
read_representation <- as.list(1:length(range_reads))
read_representation <- lapply(read_representation,toCpGs,
                              match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
length(match_read_cpg) #659676
table(sapply(match_read_cpg,length))
#1      2 
#652877   6799 
length(read_representation) #659676
table(unlist(read_representation))
#FALSE   TRUE 
#666366     34 

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
fdrps_actual <- lapply(toApply,calculate.fdrp.site,match_cpg_reads,range_reads,starts_cpgs)

#calculate.fdrp.site(toApply[[1]],match_cpg_reads,range_reads,starts_cpgs)

# when we do not have a read that covers this site, we set the FDRP for this site to NA
length(null)
fdrps[null] <- fdrps_actual
fdrps <- unlist(fdrps)
table(fdrps)

