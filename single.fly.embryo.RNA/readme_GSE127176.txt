

data download from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127176

Samples (54): tab files

GSE127176_series_matrix.txt.gz contain sample embryo stage information
'!Sample_data_processing "Paired-end RNA-seq reads were aligned to the reference genome (dmel release 6) using STAR (version 2.5.3a) with the parameter --outFilterMultimapNmax 1." '

########################################
use tmp.pl to gzip -d all xx.gz files

@files=glob('./*tab.gz');
for(@files){
    `gzip -d $_`;
}


all contain 17485 lines.

$cut -f1 GSM3629755_E1D_GATCAG.ReadsPerGene.out.tab >../all.gene.names.txt
17481 FBgn gene names

########################################
Download and save the transcriptome GTF file from the **Ensembl database**. This will be downloaded as a zipped file with extension .gz.

# download GTF file
$ wget http://ftp.ensembl.org/pub/release-104/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.104.chr.gtf.gz

use `get_exon.length.per.gene.R` and R packge GenomicFeatures to get non-overlapping exon total length per gene.


##################################################
check out 01_get_exon.length.per.gene.R and 02_calTPM_PCA.R

##################################################
tab file documentation

https://www.biostars.org/p/218995/

Explanation from the manual:

Counting number of reads per gene. With --quantMode GeneCounts option STAR will count number reads per gene while mapping. A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. The counts coincide with those produced by htseq-count with default parameters. This option requires annotations (GTF or GFF with –sjdbGTFfile option) used at the genome generation step, or at the mapping step. 

STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:

column 1: gene ID

column 2: counts for unstranded RNA-seq

column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)

column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

The column you need might be the 2nd, if you have non-strand specific data. Unless you have specifically designed a strand-specific protocol, usually this is not the case. But for downstram analysis, it is encouraged to run htseq-count or another read-counting program again.

$2 is for unstranded hits, but those overlapping on opposite strand of features are considered ambiguous. $3 reports hits based on the strand you have given in your gff annotation, and $4 in the reverse direction of your features in gff (for PE-data the 5'3'-direction is also considered). Refer to -s option of htseq-count

Some of the reads might have a splice-form expressed higher in the preceeding exon, this can be the reason for lesser reads aligning. Splicing is usually the process which leads to differences in mapping coverages over the read. Usually what you need would be based on the analysis you need to do - gene-specific, transcript-specific or something else altogether.

If $3>$2 or $4>$2, it means that unstranded data has lesser hits than the cases of directions, happens when reads are assigned as ambiguous. Reads overlapping on the opposite strand (of the gff feature) are considered ambiguous for unstranded but are assigned to genes for $3 or $4. This is the reason for the difference, this was refered here on STAR discussion

