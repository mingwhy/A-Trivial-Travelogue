
library(ggplot2)
module.go=readRDS('sc_useall_n10c0.1_T3/K25_modules.GO.rds')
#module.go=readRDS('sn_useall_n10c0.1_T3/K17_modules.GO.rds')
go.out=lapply(module.go,function(i){
  i@result
})

sapply(go.out,nrow)
pick.go=go.out;
sapply(pick.go,nrow)
names(pick.go)=paste('module',1:length(pick.go))

sapply(pick.go,function(i) sum(i$p.adjust<0.01) )

p.adj.cut=0.01;#0.05
(p.line=-1*log(p.adj.cut,base=10))
plots=list()
for(i in 1:length(pick.go)){
  module.name=names(pick.go)[i]
  tmp=pick.go[[i]]
  tmp$log.p.adjust= -1*log(tmp$p.adjust,base=10)
  if(nrow(tmp)>10){tmp=tmp[1:10,]}
  tmp=tmp[order(tmp$log.p.adjust),]
  tmp$Description=factor(tmp$Description,tmp$Description)
  plots[[i]]=ggplot(tmp,aes(x=Description,y=log.p.adjust,col=log.p.adjust))+
    geom_bar(aes(fill=log.p.adjust),stat='identity',width = 0.1)+
    theme_bw(base_size=15)+ylab('-log10(p.adjust)')+xlab('')+
    scale_color_distiller(name='',palette = "RdYlBu")+
    scale_fill_distiller(name='',palette = "RdYlBu")+
    #scale_color_distiller(name='',palette = "Dark2")+
    geom_point(size=5)+coord_flip()+
    geom_hline(yintercept = p.line, linetype="dashed",color = "black", size=0.3)+
    ggtitle(module.name)+
    theme(legend.position = 'none',
          axis.text=element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

#plots[4]=NULL
library(grid);library(gridExtra)

#pdf('sn_GO_enrich.pdf',useDingbats = T);
pdf('sc_GO_enrich.pdf',useDingbats = T,width = 16,height = 10)
#pdf('sn_GO_enrich.pdf',useDingbats = T,width = 16,height = 10)
#for(p in plots){print(p)}
do.call(grid.arrange,c(plots,ncol=2))
dev.off()
