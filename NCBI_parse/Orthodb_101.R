####
#python package taxonkit: *实战，LCA of human, mouse, rat* 
#  $ echo Homo sapiens | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
#Homo sapiens	9606	species
#$ echo Mus Musculus | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
#Mus Musculus	10090	species
#$ echo Rattus norvegicus | taxonkit name2taxid | taxonkit lineage -i 2 -r -L
#Rattus norvegicus	10116	species

#$ echo 9606 10090 10116 | taxonkit lca
#9606 10090 10116	314146

#$echo 314146 | taxonkit lineage | taxonkit reformat | cut -f 1,3
#314146	Eukaryota;Chordata;Mammalia;;;;
####

#https://github.com/lab83bio/Cotransitions/blob/master/Utilities/procedure_Orthodb_read_tables.r
library(data.table)

oDB='odb11v0'; #download from: https://data.orthodb.org/download/
#NCBI taxonomy nodes (levels) where orthologous groups (OGs) are calculated
lev <- fread(paste0(oDB,'_levels.tab.gz'),header = F, tmpdir='.', sep = '\t',
             col.names=c('Taxid','Group','num.genes','num.orthogroups','num.orgs'))
#OrthoDB orthologous groups
OGs <- fread(paste0(oDB,'_OGs.tab.gz'),header = F, tmpdir='.', sep = '\t',
             col.names=c('Orthogroup','Taxid','Description'))
#OGs to genes correspondence
Gcs <- fread(paste0(oDB,'_OG2genes.tab.gz'), header=F, tmpdir='.', sep = '\t',
             col.names=c('Orthogroup','GeneID'))

# Select level
lev[lev$Taxid=='314146',]
#Taxid            Group num.genes num.orthogroups num.orgs
#1: 314146 Euarchontoglires   1440514           29131       70
sum(OGs$Taxid == '314146')
#[1] 29131

sOGs=OGs[Taxid == '314146',Orthogroup]
length(sOGs)
#[1] 29131  

Gs <- Gcs[Orthogroup %in% sOGs]
Gs$GenomeID=gsub(':.*','',Gs$GeneID)


#number of genes and genomes per Orthogroup
Ngenes <- Gs[,.(Total=.N),by="Orthogroup"]
Ngenomes <- Gs[,.N,by=c("Orthogroup","GenomeID")][,.N,by="Orthogroup"]
Gs_count <- merge(Ngenes,Ngenomes,by="Orthogroup")

if(F){
#number of Orthogroups per genome
Ngroups<- Gs[,.N,by=c("Orthogroup","GenomeID")][,.N,by="GenomeID"]
#write(paste0("Ortogroups: ",dim(Gs_count)[1]," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())

Ngroups$Taxid <-  gsub( '_.*','',Ngroups$GenomeID )
Ngroups <- Ngroups[order(Taxid,-N)] #order by Num orthogroups
Dup <- Ngroups[duplicated(Ngroups$Taxid)]$GenomeID # duplicated species with less orthogroups
Ngroups <- Ngroups[!(GenomeID %in% Dup)]
#write(paste0("Removed duplicated Taxid: ",Dup),stderr())
}

# one-to-one orthologs
x=grep('9606|10090|10116', Gs$GenomeID)
table(Gs[x,]$GenomeID)
#10090_0 10116_0  9606_0 
#22111   21982   20313 
Gs_3species=Gs[GenomeID %in% c('9606_0','10090_0','10116_0'),]

Ngenes <- Gs_3species[,.(Total=.N),by="Orthogroup"]
group.names=Ngenes[Total==3,Orthogroup]

oneTo1_3species=Gs_3species[ Orthogroup %in% group.names,]
table(table(oneTo1_3species$Orthogroup)) #all contain 3 gene members, 13159 groups

saveRDS(oneTo1_3species,'oneTo1_3species.rds')

# convert GeneID in Orthodb to external gene identifiers
#OGs to genes correspondence (read in R, too slow)
#genes <- fread(paste0(oDB,'_genes.tab.gz'), header=F, tmpdir='.', sep = '\t',
#             col.names=c('OrthoDB.gene.id','organism.id','protein.id','synonyms','uniprot.id','Ensembl.id','ncbi.id','description'))
#genes.xref <- fread(paste0(oDB,'_gene_xrefs.tab.gz'), header=F, tmpdir='.', sep = '\t',
#               col.names=c('OrthoDB.gene.id','external.gene.identifier','external.DB.name')

# use terminal instead:
#$grep '^9606_0|^10090_0|^10116_0' odb11v0_gene_xrefs.tab >tmp.out
#$grep '^9606_0|^10090_0|^10116_0' odb11v0_genes.tab >tmp2.out
