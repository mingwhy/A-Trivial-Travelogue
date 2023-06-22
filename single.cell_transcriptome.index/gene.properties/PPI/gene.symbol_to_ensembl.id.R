##################################################
# fly, use AnnotationDbi
if(F){
  library(AnnotationDbi)
  library(org.Dm.eg.db)
  if(F){BiocManager::install("org.Mmu.eg.db")}
  library(org.Mmu.eg.db)
  package.version('AnnotationDbi')#"1.56.2"
  package.version('org.Dm.eg.db') # "3.14.0"
  package.version('org.Mmu.eg.db')# "3.14.0"
  columns(org.Dm.eg.db)
  x=as.data.frame(data.table::fread('7227.protein.info.v11.5.txt'))
  head(x[,1:3])
  dim(x) #13932
  x$FBpp=gsub('7227\\.','',x[,1])
  
  genes=select(org.Dm.eg.db,keys=x$FBpp,keytype='FLYBASEPROT',columns=c('FLYBASE','SYMBOL','ENSEMBL'))
  dim(genes) #16020
  head(genes)       
  sum(is.na(genes$ENSEMBL)) #777 without mapping
  head(genes[is.na(genes$ENSEMBL),])
  
  genes=genes[!is.na(genes$ENSEMBL),]
  dim(genes) #15243
  # one FLYBASEPROT -> one FLYBASE -> multiple ENSEMBL
  y=genes[duplicated(genes$FLYBASEPROT),]
  y=genes[genes$FLYBASEPROT %in% y$FLYBASEPROT,]
  
  genes=genes[genes$FLYBASE==genes$ENSEMBL,]
  dim(genes)#13155
  genes[duplicated(genes$FLYBASEPROT),] #none
  genes[duplicated(genes$ENSEMBL),] #none
  
  data.table::fwrite(genes,'fly_gene.id.mapping.txt')
}

##################################################
# fly (use biomart)
## extract all mouse gene id using biomart
library(biomaRt)
fly= useEnsembl(version = 99, #this archived still contain dn ds values
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'dmelanogaster_gene_ensembl')
head(listAttributes(fly) ,20)
fly_genes <- getBM(attributes = c('ensembl_gene_id', 'description',
                                    'chromosome_name',  
                                    'start_position', 
                                    'end_position'),
                     mart = fly)

## read in string protein info
x=as.data.frame(data.table::fread('7227.protein.aliases.v11.5.txt.gz'))
head(x[,1:3])
length(unique(x$`#string_protein_id`)) #13932
x[x$`#string_protein_id`==x$`#string_protein_id`[1],]
dim(x[x$source=='Ensembl_gene',]) # 13932   
x=x[x$source=='Ensembl_gene',]
dim(x)#13932  

x0=as.data.frame(data.table::fread('7227.protein.info.v11.5.txt.gz'))
colnames(x0);
x0=x0[,c(1,2)]
head(x0)

x=as.data.frame(data.table::fread('7227.protein.aliases.v11.5.txt.gz'))
head(x[,1:3])
length(unique(x$`#string_protein_id`)) #13932
x[x$`#string_protein_id`==x$`#string_protein_id`[1],]
dim(x[x$source=='Ensembl_gene',]) # 13932   
x=x[x$source=='Ensembl_gene',]
dim(x)#13932  

x=merge(x,x0)
dim(x)#13932

sum(x$alias %in% fly_genes$ensembl_gene_id) #13433

genes=merge(x,fly_genes,by.x='alias',by.y='ensembl_gene_id')
dim(genes) #13433
head(genes)

data.table::fwrite(genes,'fly_gene.id.mapping.txt')


##################################################
# mouse, use biomart
## extract all mouse gene id using biomart
library(biomaRt)
mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')
head(listAttributes(mouse) ,20)
mouse_genes <- getBM(attributes = c('ensembl_gene_id', 'description',
                              'chromosome_name',  
                              'start_position', 
                              'end_position'),
               mart = mouse)

## read in string protein info

x0=as.data.frame(data.table::fread('10090.protein.info.v11.5.txt.gz'))
colnames(x0);
x0=x0[,c(1,2)]
head(x0)

x=as.data.frame(data.table::fread('10090.protein.aliases.v11.5.txt.gz'))
head(x[,1:3])
length(unique(x$`#string_protein_id`)) #22048
x[x$`#string_protein_id`==x$`#string_protein_id`[1],]
dim(x[x$source=='Ensembl_gene',]) # 22048   
x=x[x$source=='Ensembl_gene',]
dim(x)#22048  

x=merge(x,x0)
dim(x)#22048
sum(x$alias %in% mouse_genes$ensembl_gene_id) #21726

genes=merge(x,mouse_genes,by.x='alias',by.y='ensembl_gene_id')
dim(genes) #21726
head(genes)

data.table::fwrite(genes,'mouse_gene.id.mapping.txt')

##################################################
## id conversion for mouse
genes=data.table::fread('mouse_gene.id.mapping.txt')
x=data.table::fread('mouse_ppi-cutoff400_link.txt')
head(genes)
head(x)
colnames(x)=c('gene1','gene2','score')
sum(x$gene1 %in% genes$preferred_name)
x=x[x$gene1 %in% genes$preferred_name & x$gene2 %in% genes$preferred_name,]

dim(x) #1839074
tmp=genes[match(x$gene1,genes$preferred_name),]
sum(tmp$preferred_name==x$gene1) #1839074
x$Ensembl_gene1=tmp$alias; 

dim(x) #1839074
tmp=genes[match(x$gene2,genes$preferred_name),]
sum(tmp$preferred_name==x$gene2) #1839074
x$Ensembl_gene2=tmp$alias; 

head(x)
data.table::fwrite( x,'ensembl_mouse_ppi-cutoff400_link.txt')


##################################################
## id conversion for mouse
genes=data.table::fread('fly_gene.id.mapping.txt')
x=data.table::fread('fly_ppi-cutoff400_link.txt')
head(genes)
head(x)
colnames(x)=c('gene1','gene2','score')
sum(x$gene1 %in% genes$preferred_name)
x=x[x$gene1 %in% genes$preferred_name & x$gene2 %in% genes$preferred_name,]

dim(x) #854428
tmp=genes[match(x$gene1,genes$preferred_name),]
sum(tmp$preferred_name==x$gene1) #854428
x$Ensembl_gene1=tmp$alias; 

dim(x) #854428
tmp=genes[match(x$gene2,genes$preferred_name),]
sum(tmp$preferred_name==x$gene2) #854428
x$Ensembl_gene2=tmp$alias; 

head(x)
data.table::fwrite( x,'ensembl_fly_ppi-cutoff400_link.txt')
