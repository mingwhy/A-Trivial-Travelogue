
library(tidyverse)
library(grid);library(gridExtra)

plot_one_cell_hist<-function(i.cell,xs,as,cutoff, i.sex='unknown',p.value=1,cell.type=''){
  tmp=data.frame('gene'=c(names(xs),names(as)),
                 'umi'=c(xs,as),
                 'chr'=rep(c('X','A'),c(length(xs),length(as))))
  #tmp=tmp[!is.na(tmp$umi),]
  tmp=tmp[tmp$umi!=0,] #remove non-expressed genes
  nx=sum(tmp$chr=='X');
  na=sum(tmp$chr!='X')
  
  my.breaks=seq(min(tmp$umi),max(tmp$umi)+1,len=20)
  #my.breaks[length(my.breaks)]=max(tmp$umi)
  my.breaks=round(my.breaks,2)
  my.breaks=c(-1,0,1,my.breaks)
  my.breaks=my.breaks[!duplicated(my.breaks)]
  my.breaks=sort(my.breaks)
  
  x=cut(tmp$umi,breaks=my.breaks,include.lowest = F)
  tmp1=tmp;
  tmp1$bin=x
  
  #tmp2=as.data.frame(table(tmp$umi,tmp$chr))
  tmp2=as.data.frame(table(tmp1$bin,tmp$chr))
  colnames(tmp2)=c('umi','chr','ngene')
  tmp3=tmp2 %>% group_by(chr) %>% mutate(total.gene=sum(ngene))
  tmp3$gene.perc=tmp3$ngene/tmp3$total.gene
  
  g1<-ggplot(tmp3,aes(x=umi,y=gene.perc,fill=chr,group=chr))+
    geom_bar(stat='identity',position = 'dodge')+
    theme_classic()+
    xlab('log10(CPM+1)')+
    theme(axis.text.x = element_text(size=5,angle = 45,hjust=1,vjust=1))+
    ggtitle(paste(#'cell',i.cell,
                  #',sex=',i.sex,'\n',
                  'confidence >',cutoff,'\n',
                  'X.gene=',nx,'A.gene=',na,'\n'));
                  #'one tail, X>A,ks.test pvalue=',round(p.value,4)))
                  #'two-sided ks.test pvalue=',round(p.value,4)))
   return(g1)
}


#plot_one_cell(i,sex[i],xs,as,x$p.value)    
plot_one_cell<-function(i.cell,i.sex,xs,as,p.value){
  tmp=data.frame('gene'=c(names(xs),names(as)),
                 'umi'=c(xs,as),
                 'chr'=rep(c('X','A'),c(length(xs),length(as))))
  tmp=tmp[!is.na(tmp$umi),]
  tmp=tmp[tmp$umi!=0,] #remove non-expressed genes
  nx=sum(tmp$chr=='X');
  na=sum(tmp$chr!='X')
  
  my.breaks=seq(min(tmp$umi),max(tmp$umi)+1,len=10)
  #my.breaks[length(my.breaks)]=max(tmp$umi)
  #my.breaks=floor(my.breaks)
  my.breaks=round(my.breaks,2)
  #my.breaks=c(-1,0,1,my.breaks)
  my.breaks=my.breaks[!duplicated(my.breaks)]
  my.breaks=sort(my.breaks)
  
  x=cut(tmp$umi,breaks=my.breaks,include.lowest = T)
  tmp1=tmp;
  tmp1$bin=x
  
  #tmp2=as.data.frame(table(tmp$umi,tmp$chr))
  tmp2=as.data.frame(table(tmp1$bin,tmp$chr))
  colnames(tmp2)=c('umi','chr','ngene')
  tmp3=tmp2 %>% group_by(chr) %>% mutate(total.gene=sum(ngene))
  tmp3$gene.perc=tmp3$ngene/tmp3$total.gene
  
  g1<-ggplot(tmp3,aes(x=umi,y=gene.perc,fill=chr,group=chr))+
    geom_bar(stat='identity',position = 'dodge')+
    theme_classic()+
    xlab('log10(CPM+1)')+
    theme(axis.text.x = element_text(size=5,angle = 45,hjust=1,vjust=1))+
    ggtitle(paste('cell',i.cell,'from',cell.type,
                  ',sex=',i.sex,
                  '\n',
                  'X.gene=',nx,'A.gene=',na,'\n',
                  #'one tail, X>A,ks.test pvalue=',round(p.value,4)))
                  'two-sided ks.test pvalue=',round(p.value,4)))
  
  # create ECDF of data
  g2<-ggplot(tmp, aes(umi, colour = chr)) +
    stat_ecdf()+theme_classic()+
    xlab('log10(CPM+1)')+
    ggtitle(paste('cell',i.cell,'from',cell.type,
                  ',sex=',i.sex,
                  '\n',
                  'X.gene=',nx,'A.gene=',na,'\n',
                  #'one tail, X>A,ks.test pvalue=',round(p.value,4)))
                  'two-sided ks.test pvalue=',round(p.value,4)))
  
  g <- gridExtra::arrangeGrob(g1, g2, ncol=1,nrow=2)
  return(g) 
  #grid::grid.draw(g)
}


plot_one_cell2<-function(i.cell,i.sex,xs,as,p.value){
  tmp=data.frame('gene'=c(names(xs),names(as)),
                 'umi'=c(xs,as),
                 'chr'=rep(c('X','A'),c(length(xs),length(as))))
  tmp=tmp[!is.na(tmp$umi),]
  #tmp=tmp[tmp$umi!=0,]
  nx=sum(tmp$chr=='X');
  na=sum(tmp$chr!='X')
  
  tmp$log.umi=log(tmp$umi)
  
  my.breaks=seq(min(tmp$log.umi)-1,max(tmp$log.umi)+1,len=10)
  my.breaks=round(my.breaks,2)
  #my.breaks[length(my.breaks)]=max(tmp$log.umi)
  #my.breaks=c(0,1,my.breaks)
  my.breaks=my.breaks[!duplicated(my.breaks)]
  my.breaks=sort(my.breaks)
  
  x=cut(tmp$log.umi,breaks=my.breaks,include.lowest = TRUE)
  tmp1=tmp;
  tmp1$bin=x
  
  #tmp2=as.data.frame(table(tmp$umi,tmp$chr))
  tmp2=as.data.frame(table(tmp1$bin,tmp$chr))
  colnames(tmp2)=c('umi','chr','ngene')
  tmp3=tmp2 %>% group_by(chr) %>% mutate(total.gene=sum(ngene))
  tmp3$gene.perc=tmp3$ngene/tmp3$total.gene
  
  g1<-ggplot(tmp3,aes(x=umi,y=gene.perc,fill=chr,group=chr))+
    geom_bar(stat='identity',position = 'dodge')+
    theme_classic()+
    #xlab('Gene expression after normalization and transformation')+
    xlab('log(gene.expression)')+
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))+
    #ggtitle(paste('cell',i.cell,'from',cell.type,
    ggtitle(paste('cell',i.cell,'from embryo',
                  ',sex=',i.sex,
                  '\n',
                  'X.gene=',nx,'A.gene=',na,'\n',
                  #'one tail, X>A,ks.test pvalue=',round(p.value,4)))
                  'two-sided ks.test pvalue=',round(p.value,4)))  
  
  # create ECDF of data
  g2<-ggplot(tmp, aes(umi, colour = chr)) +
    stat_ecdf()+theme_classic()+
    #xlab('Gene expression after normalization and transformation')+
    xlab('gene.expression')+
    #ggtitle(paste('cell',i.cell,'from',cell.type,
    ggtitle(paste('cell',i.cell,'from embryo',
                  ',sex=',i.sex,
                  '\n',
                  'X.gene=',nx,'A.gene=',na,'\n',
                  #'one tail, X>A,ks.test pvalue=',round(p.value,4)))
                  'two-sided ks.test pvalue=',round(p.value,4)))  
  
  #g <- gridExtra::arrangeGrob(g1, g2, ncol = 2)
  g <- gridExtra::arrangeGrob(g1, g2, ncol=1,nrow=2) 
  return(g) 
  #grid::grid.draw(g)
}