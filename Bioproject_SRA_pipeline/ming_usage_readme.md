
Original pipeline from: https://github.com/sjbush/expr_atlas

My adaptations are as below:

**Script 1: make_protein_coding_transcriptome.pl**

1, install `kallisto` from https://pachterlab.github.io/kallisto/about
`brew install kallisto`, as step1 index requires it.

2, download genome reference required files, 
3, then run 1.*pl and 1*pl

file: `Rattus_norvegicus.protein_coding.fa` should be generated.

**Script 2: parse_bioproject_summary_file.pl**

install pysradb from https://github.com/saketkc/pysradb
`https://github.com/saketkc/pysradb`
test:
`pysradb metadata SRP000941`
`pysradb metadata SRP075720 --detailed`

then find `PRJNA516151` line in `identify_metadata_for_each_bioproject.sh` and 
run `pysradb metadata --detailed PRJNA516151 --saveto metadata/Rattus_norvegicus/PRJNA516151.txt`
File `PRJNA516151.txt` should be in `metadata/Rattus_norvegicus` dir.

**Script 3: parse_metadata.pl**

run 3.parse_metadata_ming.pl


**Script 4: download_fqs_and_run_kallisto.pl**

This script takes as input the sample metadata table produced by the previous script. It produces a directory of shell scripts, one per sample, which in this case are each formatted for submission to a SLURM cluster. Each script in this batch executes a series of commands per sample: downloading sequencing data (using [Aspera](https://www.ibm.com/aspera/connect/) and the [SRA toolkit](https://github.com/ncbi/sra-tools), where appropriate), pre-processing (using [fastp](https://github.com/OpenGene/fastp)), downsampling (using [seqtk](https://github.com/lh3/seqtk)) to _x_ million reads _y_ times, and finally quantifying expression using [Kallisto](https://pachterlab.github.io/kallisto/about). The output is a directory containing one subdirectory per sample ID, within which are the Kallisto output files (one per downsampled replicate).

**Script 5: parse_kallisto_output.pl**



