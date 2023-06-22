
# uniprot
from `https://www.uniprot.org/help/downloads` download
Reviewed (Swiss-Prot) FAQ    fasta
Unreviewed (TrEMBL) FAQ    fasta

$wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
$wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

Why is UniProtKB composed of 2 sections, UniProtKB/Swiss-Prot and UniProtKB/TrEMBL?
https://www.uniprot.org/help/uniprotkb_sections

What is the canonical sequence? Are all isoforms described in one entry?
https://www.uniprot.org/help/canonical_and_isoforms


$grep '>sp' /Users/mingyang/Downloads/uniprot_sprot.fasta  | wc -l
  569213

# get readme file, save as FAQ_Swiss-Prot_README.txt
$wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README

fasta header info: https://www.uniprot.org/help/fasta-headers
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion

$grep 'OS=Mus musculus' /Users/mingyang/Downloads/uniprot_sprot.fasta | wc -l
   17150
