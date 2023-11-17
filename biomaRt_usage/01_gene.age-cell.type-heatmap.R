

####################################
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater);library(scran)
library(ggplot2);library(gridExtra);
library(tidyverse);library(Matrix)
options(stringsAsFactors = F)
library(org.Mm.eg.db,verbose=F,quietly=T)
library(GO.db);
#https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN <- function(x) sprintf("%.3f", x)

####################################
## read in mouse turnover rate data
cell.lifespan=data.table::fread('../src/Dataset_S1.txt')
cell.lifespan$Ron_tc=cell.lifespan$`cell type annotation in Sender and Milo 2021`
tc.orders=cell.lifespan[order(cell.lifespan$lifespan),]$`cell type annotation in TMS`
tc.orders #39 tc

gene.meta=readRDS('gene-phylostratigraphic-age_Tautz2013.rds')
table(gene.meta$Phylostratum.rank)
bg=table(gene.meta$Phylostratum.rank)

#####################################################################################
## matrix: gene.age.bin x cell.types order by cell lifespan
if(!file.exists('gene.age_cell.type.rds')){
  
  df.mouse=readRDS('../1112_yeast/gene.expr_over_cell.types.rds')
  table(df.mouse$age)
  all.out=c()
  
  for(pick.age in c('3m','24m')){
    #pick.age='3m'
    df.one.age=df.mouse[df.mouse$age==pick.age,]
    out<-lapply(tc.orders,function(tc){
      df.one.age.one.tc=df.one.age[df.one.age$cell.type==tc,]
      #10% & >=5cell
      tmp=df.one.age.one.tc[df.one.age.one.tc$nnz_per_gene >=max(10,df.one.age.one.tc$ncell*0.2),]
      tmp1=merge(tmp,gene.meta,by.x='mgi_symbol',by.y='mgi_symbol')
      x=tmp1 %>% group_by( Phylostratum,Phylostratum.rank) %>%
              summarise(ngene=n()) %>% mutate(cell.type=tc,age=pick.age)
      x=x %>% mutate(prop=ngene/sum(x$ngene))
      x
    })
    df=as.data.frame(Reduce(`rbind`,out))
    all.out=rbind(all.out,df)
  }
  saveRDS(all.out,'gene.age_cell.type.rds')
}

###############################################
all.out=readRDS('gene.age_cell.type.rds')
dim(all.out)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
all.out.one.age=all.out[all.out$age=='3m',]
head(all.out.one.age)
#https://stackoverflow.com/questions/5890584/how-to-reshape-data-from-long-to-wide-format
# co.coeff mat
tmp=reshape(all.out.one.age[,c(2,4,6)], idvar = "cell.type", timevar = "Phylostratum.rank", direction = "wide")
mat=tmp[,-1]
rownames(mat)=tmp[,1]
colnames(mat)=gsub('prop.','',colnames(mat))

#mat=as.matrix(mat)
colnames(mat)=factor(colnames(mat),levels=as.character(sort(as.numeric(colnames(mat)))))
rownames(mat)=factor(rownames(mat),levels=tc.orders)
sum(mat==0)
mat.log10=log(mat,base=10)

pdf("test.pdf");
Heatmap(mat.log10,  name='log10(percentage)',
        col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = as.character(sort(as.numeric(colnames(mat)))),
        row_order=tc.orders,
        #cluster_rows = FALSE,cluster_columns = FALSE
)
dev.off()

p.values.over=c()
p.values.under=c()
for(i in 1:nrow(all.out.one.age)){
  x=all.out.one.age[i,]
  cell.type.bg=all.out.one.age[all.out.one.age$cell.type==x$cell.type,]
  x1=x$ngene;
  x2=sum(cell.type.bg$ngene)-x1
  y1=bg[x$Phylostratum.rank]
  y2=sum(bg)-y1
  
  res=fisher.test(matrix(c(x1,x2,y1,y2),nrow=2,byrow=T),alternative='greater')
  p.values.over=c(p.values.over,res$p.value)
  
  res=fisher.test(matrix(c(x1,x2,y1,y2),nrow=2,byrow=T),alternative='less')
  p.values.under=c(p.values.under,res$p.value)
}
all.out.one.age$pvalue.over=p.adjust(p.values.over,method='bonferroni')
all.out.one.age$pvalue.under=p.adjust(p.values.under,method='bonferroni')
summary(p.values.over)
summary(p.values.under)

# p value mat
colnames(all.out.one.age)
#tmp2=reshape(all.out.one.age[,c("Phylostratum.rank","cell.type","pvalue.over")],  idvar = "cell.type", timevar = "Phylostratum.rank", direction = "wide")
tmp2=reshape(all.out.one.age[,c("Phylostratum.rank","cell.type","pvalue.under")],  idvar = "cell.type", timevar = "Phylostratum.rank", direction = "wide")

mat2=tmp2[,-1]
rownames(mat2)=tmp2[,1]
colnames(mat2)=gsub('pvalue.','',colnames(mat2))
mat2=as.matrix(mat2)

# replace kegg compound id with metabolite name
pdf("heatmap_Pval0.01_remove0.pdf",width=10,height=6.5)
Heatmap(mat.log10,  name='log10(percentage)',
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(mat2[i,j]) & mat2[i, j] < 0.01)
            # grid.text(sprintf("%.1f", mat2[i, j]), x, y, gp = gpar(fontsize = 10))
            grid.text('*',x,y,gp = gpar(fontsize = 15))
        },
        col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        column_order = as.character(sort(as.numeric(colnames(mat)))),
        row_order=tc.orders,
    column_title = "Under-representation"
)
dev.off()

