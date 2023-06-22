library(ggplot2)
library(gridExtra)

dat=readRDS('./mz.filter.dat.rds')
df.mat=dat$edata;
pheno=dat$pheno

## log transform
summary(as.numeric(unlist(df.mat))) #no zero, log directly
df.mat.log=log(df.mat)

## normalize by sample
df.scaled=scale(df.mat.log,center=T,scale=T)
apply(df.scaled,2,mean) #confirm
apply(df.scaled,2,sd) #confirm

## PCA, view by batch, by genotype, by age
pca.out=prcomp(t(df.scaled)) #feature as column

pairs(pca.out$x[,1:6]) #color by genotype, age, batch, run.order
pairs(pca.out$x[,1:6],col=pheno$batch,pch=16,main='by batch')
pairs(pca.out$x[,1:6],col=factor(pheno$genotype),pch=16,main='by genotype')
pairs(pca.out$x[,1:6],col=factor(pheno$age),pch=16,main='by age')

dim(pca.out$rotation) #loading of orginal features of PC
colnames(pca.out$rotation)
dim(pca.out$x) #new coordinates for each sample
colnames(pca.out$x)
rownames(pca.out$x) #MTxx samples

pc.prop=pca.out$sdev/sum(pca.out$sdev)
sum(pc.prop)
df.pc=data.frame(PC=paste0('PC',1:length(pc.prop)),explained.var.prop=pc.prop/sum(pc.prop))
df.pc$PC=factor(df.pc$PC,levels=df.pc$PC)
pc.plot=ggplot(df.pc,aes(x=PC,y=explained.var.prop*100))+geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
pc.plot
sum(df.pc$explained.var.prop[1:2]) #20.17% 

new.coord=as.data.frame(pca.out$x[,1:3])
sum(pheno$sample==rownames(new.coord))
new.coord$age=pheno$age;new.coord$genotype=pheno$genotype;
new.coord$batch=pheno$batch
p=ggplot(new.coord,aes(x=PC1,y=PC2))+geom_point()+theme_bw(base_size = 14)+
  xlab(paste0('PC1 (',round(df.pc[1,2],4)*100,'%)'))+
  ylab(paste0('PC2 (',round(df.pc[2,2],4)*100,'%)'))
p+geom_point(aes(col=factor(batch)),size=4)
p+geom_point(aes(col=age),size=4)
p+geom_point(aes(col=genotype),size=4)
p+geom_point(aes(col=genotype,shape=age),size=4)

## estimate batch effect using ANOVA
mz=df.scaled[1,]
tmp=cbind(mz,pheno)
fit1=lm(mz~batch+age+genotype,data=tmp)
summary(fit1)
x=anova(fit1)
x$`Sum Sq`[1]/sum(x$`Sum Sq`)

batch.var=sapply(1:nrow(df.scaled),function(i){
  mz=df.scaled[i,]
  tmp=cbind(mz,pheno)
  fit1=lm(mz~batch+age+genotype,data=tmp)
  summary(fit1)
  x=anova(fit1)
  x$`Sum Sq`[1]/sum(x$`Sum Sq`)
})
summary(batch.var)
hist(batch.var)

## remove batch effect for PCA visualization only
library(sva)

batch=pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=df.scaled, batch=batch, 
                      mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

dim(combat_edata) #170mz X 60samples
dim(df.scaled)
df.scaled[1:3,1:3]
combat_edata[1:3,1:3]
saveRDS(combat_edata,'mz.filter.combat.rds')

## perform PCA on batch corrected data
pca.out2=prcomp(t(combat_edata)) #feature as column

pairs(pca.out2$x[,1:6],col=pheno$batch,pch=16,main='by batch')
pairs(pca.out2$x[,1:6],col=factor(pheno$genotype),pch=16,main='by genotype')
pairs(pca.out2$x[,1:6],col=factor(pheno$age),pch=16,main='by age')


dim(pca.out2$rotation) #loading of orginal features of PC
dim(pca.out2$x) #new coordinates for each sample
colnames(pca.out2$x) #PC
rownames(pca.out2$x) #MTxx samples

pc.prop=pca.out2$sdev/sum(pca.out2$sdev)
sum(pc.prop)
df.pc=data.frame(PC=paste0('PC',1:length(pc.prop)),explained.var.prop=pc.prop/sum(pc.prop))
df.pc$PC=factor(df.pc$PC,levels=df.pc$PC)
pc.plot=ggplot(df.pc,aes(x=PC,y=explained.var.prop*100))+geom_bar(stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
pc.plot
sum(df.pc$explained.var.prop[1:2]) #20.67% 

new.coord=as.data.frame(pca.out2$x[,1:3])
sum(pheno$sample==rownames(new.coord))
new.coord$age=pheno$age;new.coord$genotype=pheno$genotype;
new.coord$batch=pheno$batch
p=ggplot(new.coord,aes(x=PC1,y=PC2))+geom_point()+theme_bw(base_size = 14)+
  xlab(paste0('PC1 (',round(df.pc[1,2],4)*100,'%)'))+
  ylab(paste0('PC2 (',round(df.pc[2,2],4)*100,'%)'))
p+geom_point(aes(col=factor(batch)),size=4)
p+geom_point(aes(col=age),size=4)
p+geom_point(aes(col=genotype),size=4)
p+geom_point(aes(col=genotype,shape=age),size=4)

## check mz loading for PC1, PC2, PC3.
dim(pca.out2$rotation) #170x60,loading of orginal features of PC
pca.out2$rotation[,1] #loading of mz to PC1
pheatmap::pheatmap(pca.out2$rotation[,1:3],cluster_cols=F,fontsize_row=4)

#I ranked metabolites using their **absolute contribution (loading) values** to PC1~3.
x1=order(abs(pca.out2$rotation[,1]),decreasing = T)
x1=data.frame(pca.out2$rotation[x1[1:10],1])
colnames(x1)='contribution to PC1'
x1

x2=order(abs(pca.out2$rotation[,2]),decreasing = T)
x2=data.frame(pca.out2$rotation[x2[1:10],1])
colnames(x2)='contribution to PC2'
x2

x3=order(abs(pca.out2$rotation[,3]),decreasing = T)
x3=data.frame(pca.out2$rotation[x3[1:10],1])
colnames(x3)='contribution to PC3'
x3
