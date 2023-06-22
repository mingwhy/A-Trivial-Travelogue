
geneID
dNdS
phylostratigraphy gene age 
connectivity (PPI)

# geneID
fly: FBgnxxx
mouse: ENSMUSG00000064372 (ensembl_gene_id)

# phylostratigraphy gene age 
database: GenOrigin
http://genorigin.chenzxlab.cn/#!/download

# dn ds
fly: dnds: data downloaded from: http://www.flydivas.info/ (https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/flygene_dNdS)
mouse: biomart, mouse_rat_biomart_to_retrieve_dnds.R

fly dnds data is already in FBgn id format, may need to update to latest FBgn xxx id.


# PPI
https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/PPI_string-biogrid
$ perl parse.ppi_link.pl 400 7227.protein.info.v11.5.txt 7227.protein.links.v11.5.txt >fly_ppi-cutoff400_link.txt
$ perl parse.ppi_link.pl 400 10090.protein.info.v11.5.txt 10090.protein.links.v11.5.txt >mouse_ppi-cutoff400_link.txt

as PPI data contain `gene symbols`, translate them into `ensembl_gene_id`,
with R script `gene.symbol_to_ensembl.id.R`.
