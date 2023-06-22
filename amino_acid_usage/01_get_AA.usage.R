
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);
library(Matrix);
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);

####################################
## reading in gene id.mapping
annotLookup=readRDS('mmus_id_ensembl2uniprot.rds')
dim(annotLookup) #44806
sum(annotLookup$mgi_symbol==annotLookup$external_gene_name) #42900
sum(annotLookup$external_gene_name==annotLookup$uniprot_gn_symbol) #42947

aa.length.freq=readRDS('aa_id_length_AAfreq.rds')
head(aa.length.freq)

df.gene<-merge(annotLookup,aa.length.freq,
                by.x='uniprot_gn_id', by.y='uniprot_id')
dim(df.gene) #13506
head(df.gene)

#######################################################
## calcualte EPCA_per_gene
aa.cost=data.table::fread('cost_per_AA.txt')
head(aa.cost)

which(colnames(df.gene)=='A')
sum(colnames(df.gene)[12:(12+19)]==aa.cost$AA.abbreviation) #20aa
numerator=apply(df.gene[,12:(12+19)] * aa.cost$Y20,1,sum)
denominator=df.gene$aa_length
EPCA_per_gene=numerator/denominator
df.gene$EPCA_per_gene=EPCA_per_gene

dim(df.gene) #13506
length(unique(df.gene$mgi_symbol)) #13454
length(unique(df.gene$external_gene_name)) #13454
df.gene[df.gene$external_gene_name %in% df.gene[duplicated(df.gene$external_gene_name),]$external_gene_name,]

df.gene=df.gene[!duplicated(df.gene$mgi_symbol),]
dim(df.gene) #13454
df.gene=as.data.frame(df.gene)
rownames(df.gene)=df.gene$mgi_symbol

#######################################################
## use expression proportion (remember to change `expr.m=expr.m[,tmp$age==___ ]`)
pick.age='3m';
#pick.age='24m';
#output_file=paste0('select.tc.h5ad_result_goterms/mouse_male_GO_',pick.age,'_libNorm2_Chaperon.rds');
#output_file

if(!file.exists(output_file)){
  if(T){
    sce=readH5AD('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/select.tc.h5ad') # 22966 31001 
    unique(sce$age) #3m 18m 24m
    assayNames(sce)<-'counts'
    sce_naive=sce[,sce$age %in% c('3m','24m')]
    table(sce_naive$age)
    assayNames(sce_naive)<-'counts'
    sce_naive<-logNormCounts(sce_naive, log=TRUE, pseudo.count=1) # use CPM, if log, the cell.lib.size range 2~3 orders
    assayNames(sce_naive)
    
    if(F){
      if(!file.exists('sizeFactors.rds')){
        library(scran)
        clusters <- quickCluster(sce)
        sce <- computePooledFactors(sce, clusters=clusters)
        summary(sizeFactors(sce))
        tmp=sizeFactors(sce)
        names(tmp)<-colnames(sce)
        saveRDS(tmp,'sizeFactors.rds')
      }
      
      preCalculated.sizeFactors=readRDS('sizeFactors.rds');
      sum(names(preCalculated.sizeFactors)==colnames(sce))
      sizeFactors(sce)<-preCalculated.sizeFactors
    }
  }
  
  if(F){
    sce.shared=readRDS('~/Documents/aging_cell.turnover/1120_TMS_male_analysis/sce_minCellCounts.rds')
    names(sce.shared) #38 tc
    unique(Matrix::colSums(assay(sce.shared[[1]],'counts'))) #check for cell.lib.size
  }
  
  colnames(colData(sce_naive))
  tc.names=unique(sce_naive$tissue_cell.type)
  length(tc.names) #72 tc
  
  mouse_tcs_TDI.list<-lapply(tc.names,function(tc){
    
    tmp=sce_naive[,sce_naive$tissue_cell.type==tc] #raw count data
    
    #expr.m=assay(tmp,'normcounts')
    expr.m=assay(tmp,'logcounts')
    overlap.genes=intersect(rownames(expr.m),df.gene$mgi_symbol)
    expr.m=expr.m[overlap.genes,]
    gene.cost=df.gene[overlap.genes,]
    
    expr.3m=expr.m[,tmp$age=='3m']
    expr.24m=expr.m[,tmp$age=='24m']
    
    coeffs.3m<-unlist(lapply(1:ncol(expr.3m),function(i)cor(expr.3m[,i],gene.cost$EPCA_per_gene,method='spearman')))
    summary(coeffs.3m)
    coeffs.24m<-unlist(lapply(1:ncol(expr.24m),function(i)cor(expr.24m[,i],gene.cost$EPCA_per_gene,method='spearman')))
    summary(coeffs.24m)
    print(as.character(tc))
    print(wilcox.test(coeffs.3m,coeffs.24m))
  })
  
    
  
  
    #EPCA_per_cell=sum(gene.cost$EPCA_per_gene * expr.m[,1] * gene.cost$aa_length)/sum(expr.m[,1] * gene.cost$aa_length)
    EPCA_per_cell<-lapply(1:ncol(expr.m),function(i){
      sum(gene.cost$EPCA_per_gene * expr.m[,i] * gene.cost$aa_length)/sum(expr.m[,i] * gene.cost$aa_length)
    })
    EPCA_per_cell=unlist(EPCA_per_cell)
    
    n.expr.gene=Matrix::colSums(expr.m>0) #control for the number of expressed genes （https://academic.oup.com/gbe/article/12/4/300/5807614?login=true
    #expr.m=expr.m[,sce_naive$age=='24m' & n.expr.gene>=100] #keep cells which expr>=100 genes
    expr.m=expr.m[, n.expr.gene>=100] #keep cells which expr>=100 genes
    
    #used.genes=intersect(gene.meta$mgi_symbol,rownames(expr.m))
    #expr.m=expr.m[used.genes,]
    cell.expr=Matrix::colSums(expr.m)
    all.genes=rownames(expr.m)
    
    tmp.go=names(input.list);
    #lapply(input.list,function(i) intersect(all.genes,i))
    go.median.expr<-lapply(tmp.go, function(go){
      overlap.genes=intersect(all.genes,input.list[[go]])
      if(length(overlap.genes)==0){return(c(0,NA))}
      #if(length(overlap.genes)<size.cutoff){return(NA)}
      
      expr.m.tmp=expr.m[overlap.genes,,drop=F]
      go.expr=Matrix::colSums(expr.m.tmp)
      #return(c(length(overlap.genes),median(go.expr/cell.expr)))  
      return(c(length(overlap.genes),mean(go.expr/cell.expr)))  
    })
    
    return(as.data.frame(Reduce(`rbind`,go.median.expr)))
  })
  tmp.go=names(input.list);
  mouse_tcs_TDI.list2<-lapply(mouse_tcs_TDI.list,function(i){
    colnames(i)=c('ngene','score');
    i$GO.id=tmp.go;
    i})
  mouse_tcs_TDI=as.data.frame(Reduce(`rbind`,mouse_tcs_TDI.list2))
  mouse_tcs_TDI$cell.type=rep(tc.names, sapply(mouse_tcs_TDI.list2,nrow) )
  
  saveRDS(mouse_tcs_TDI,output_file)
}

####################################
### see correlation
mouse_tcs_TDI=readRDS(output_file)
summary(mouse_tcs_TDI$ngene) #1~3006
length(unique(mouse_tcs_TDI$GO.id))

#mouse_tcs_TDI=mouse_tcs_TDI[mouse_tcs_TDI$ngene>=10,]

df.go.score<-mouse_tcs_TDI[,-1] %>% spread(GO.id,score)
dim(df.go.score); #cell type by go term
rownames(df.go.score)<-df.go.score$cell.type
df.go.score=df.go.score[,-1]
all.go=colnames(df.go.score)
dim(df.go.score); #39 x 5177

overlap.tc=intersect(rownames(df.go.score),cell.lifespan$`tissue: cell.type in mouse`)
#overlap.tc=overlap.tc[overlap.tc!="Brain_Non-Myeloid:neuron"]

df1=df.go.score[overlap.tc,]
df2=cell.lifespan[match(overlap.tc,cell.lifespan$`tissue: cell.type in mouse`),]

df2$duplicate='tissue-specific estimate';
df2[df2$human_tc %in% df2$human_tc[duplicated(df2$human_tc)],]$duplicate='non tissue-specific estimate'
table(df2$duplicate)

na.number=apply(df1,2,function(i) sum(is.na(i)))
na.number

sum(rownames(df1)==df2$`tissue: cell.type in mouse`) #39
df1$`tissue: cell.type in mouse`=df2$`tissue: cell.type in mouse`

df1.long=reshape2::melt(df1)
colnames(df1.long)[2]='GOid'
intersect(colnames(df1.long),colnames(df2))
df3=merge(df1.long,df2)

#######################################################################################
## use average for non-tissue specific ones
x=df3 %>% group_by(human_tc,GOid) %>% summarise(average=mean(value))
df4=x %>% spread(GOid,average)
dim(df4) #21 x 5

# at least 10 cell types terms express this GO term's gene members
i=unlist(lapply(2:ncol(df4),function(i) sum(df4[,i]!=0)>=10) )
df4=df4[,c(1,c(2:ncol(df4))[i])]

intersect(colnames(df4),colnames(df2)) #human_tc
df5=merge(df4,df2[!duplicated(df2$human_tc),],all.x=TRUE)
dim(df5) #21

Matrix=df4[,-1];
dim(Matrix) # ncell.type x 1061 GO

df5$turnover=1-exp(-1/df5$lifespan)
SampleAge=df5$lifespan

cor.coeffs.list=t(apply(Matrix,2,function(x){
  #i=cor.test(x,SampleAge,method='spearman',use='pairwise.complete.obs')
  i=cor.test(x,SampleAge,method='kendall',use='pairwise.complete.obs')
  c(i$estimate,i$p.value)
  #value=bcdcor(x,SampleAge)
  #value
}))
cor.coeffs=as.data.frame(cor.coeffs.list)
colnames(cor.coeffs)=c('Spearman.rho','Pval')
cor.coeffs$GOid=rownames(cor.coeffs)
cor.coeffs=cor.coeffs[!is.na(cor.coeffs$Pval),]
cor.coeffs$FDR=p.adjust(cor.coeffs$Pval,method='BH')

dim(cor.coeffs) #5127
cor.coeffs=merge(cor.coeffs,df.go.chaperone)
cor.coeffs

plot_go_id=cor.coeffs$GOid

plots<-lapply(1:nrow(cor.coeffs),function(i){
  goid=cor.coeffs$GOid[i]
  goterm=cor.coeffs[cor.coeffs$GOid==goid,]$GOterm
  (cor0=round(cor.coeffs[cor.coeffs$GOid==goid,]$Spearman.rho,3))
  #(pval0=round(cor.coeffs[cor.coeffs$GOid==goid,]$Pval,6))
  (pval0=round(cor.coeffs[cor.coeffs$GOid==goid,]$FDR,6))
  
  ggplot(df5,aes(x=lifespan,y=df5[,paste(goid)]))+
  #ggplot(tmp,aes(x=turnover,y=tmp[,paste(goid)]))+
    #geom_point(aes(col=human_tc,shape=duplicate),size=3)+
    #geom_point(size=3,shape=2)+
    geom_point(size=3,shape=16)+
    scale_x_log10()+ylab('GO term activity score')+
    #ggtitle(paste0('Age ',pick.age,'\ngoid,', goterm,'\nSpearman\'s rho = ',round(cor0,3),', P value=',round(pval0,6)))+
    #ggtitle(paste0(goid,', ', goterm,'\nSpearman\'s rho = ',round(cor0,3),', P.adj value=',round(pval0,6)))+
    ggtitle(paste0('Age ',pick.age,'\n',goid,', ',goterm,'\nKendall\'s tau = ',cor0,', P.adj value=',pval0))+
    
    #'\nPearson\'s cor = ',round(cor1,3),
    #', P value=',round(pval1,6)))+
    #geom_text(aes(label=tissue_cell.type),size=3,nudge_x=0.6, nudge_y=0,hjust = 1,vjust=0)+
    #geom_errorbar(aes(ymin=`0.25`, ymax=`0.75`), width=.2, position=position_dodge(0.05))+
    xlab('Cell lifespan (day)')+
    #geom_smooth(method=lm , color="black", fill=grDevices::adjustcolor( "lightgrey", alpha.f = 0.2),se=TRUE,level=0.90) + 
    scale_shape_discrete(name='Cell lifespan estimate')
    #scale_color_viridis(name='Cell type annotation from Sender and Milo (2021)',option='turbo',discrete=T)
})

plots2=lapply(plots,function(i) 
  i+theme_classic(base_size = 17)+theme(legend.position = 'none',
                                        plot.margin = unit(c(1,1,1,1), "cm"),
                                        #plot.title = element_text(size = 12, face = "bold")))
                                        plot.title = element_text(size = 14)))

#pdf(paste0('GO_chaperon_cell.lifespan_',pick.age,'_rmDups.pdf'),useDingbats = T,width = 12,height = 9)
#grid.arrange(grobs=plots2[c(3,2,4,1)],ncol=2)
#dev.off()

pdf(paste0('GO_chaperon_cell.lifespan_',pick.age,'_rmDups_onePage.pdf'),useDingbats = T,
    width = 7,height = 5,pointsize=12)
for(i in plots2){print(i)}
dev.off()






