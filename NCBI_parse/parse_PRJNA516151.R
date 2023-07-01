
# rat-2019-Age-Related Gene Expression Signature in Rats Demonstrate Early, Late, and Linear Transcriptional Changes from Multiple Tissues
# rat-2022-Generation and network analysis of an RNA-seq transcriptional atlas for the rat 
srr<- data.table::fread("PRJNA516151.txt")
srr$sample_accession #SRS xxx

# gene-level TPM. `Rattus_norvegicus.txt` file download from https://ora.ox.ac.uk/objects/uuid:8c347ab7-7e51-4350-91a6-d9f07500cbe8
expr=data.table::fread('Rattus_norvegicus.txt')
dim(expr) #25481  7676
i=colnames(expr) %in% srr$sample_accession
srr.expr=as.data.frame(expr)[,i] #216 samples present
srr.gene=as.data.frame(expr)[,1:8]

test=srr[match(colnames(srr.expr),srr$sample_accession),]
sum(test$sample_accession==colnames(srr.expr))
srr.info=test

colnames(srr.info)
table(srr.info$strain,srr.info$subject_id)
table(srr.info$sex)  #all males
table(srr.info$age,srr.info$tissue,srr.info$sex)

saveRDS(list(srr.expr=srr.expr,srr.info=srr.info,srr.gene=srr.gene),
        file='rat_PRJNA516151.rds')

rat=readRDS('rat_PRJNA516151.rds')
srr.info=rat$srr.info
srr.gene=rat$srr.gene

