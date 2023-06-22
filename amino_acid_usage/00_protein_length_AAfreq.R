aa.table=data.table::fread('cost_per_AA.txt')
aa.table$AA.abbreviation
#BiocManager::install('seqinr')
library(seqinr)

#https://davetang.org/muse/2013/05/09/using-the-r-seqinr-package/
seqAA=read.fasta(file='uniprot_sprot.fasta.gz',seqtype = c("AA"))
length(seqAA)   #569213
seqAA[[1]]
attr(seqAA[[1]],'Annot')
length(seqAA[[1]]) #aa length
count(seqAA[[1]],wordsize=1,alphabet=aa.table$AA.abbreviation) #aa distribution

#keep mouse protein-coding genes
annotation <- getAnnot(seqAA)
mmus.seqAA<-seqAA[grep('OS=Mus musculu',annotation)]
length(mmus.seqAA) #17150

# get gene.id
annotation <- getAnnot(mmus.seqAA)
uniprot.id<-unlist(lapply(annotation,function(x){ strsplit(x,'\\|')[[1]][[2]]}))
aa.length=unlist(lapply(mmus.seqAA,length))
aa.freq<-sapply(mmus.seqAA,function(x){
  aa.freq=count(x,wordsize=1,alphabet=aa.table$AA.abbreviation) 
  aa.freq[aa.table$AA.abbreviation]
})
aa.freq=t(aa.freq)

df=data.frame(uniprot_id=uniprot.id,aa_length=aa.length)
df.aa=cbind(df,aa.freq)
head(df.aa)
saveRDS(df.aa,'aa_id_length_AAfreq.rds')

