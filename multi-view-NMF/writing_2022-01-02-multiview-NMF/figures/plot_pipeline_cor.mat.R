########################################################
## cell cluster specific gene co-expression matrix
#https://www.schemecolor.com/yellow-purple-gradient.php
#Yellow-Purple Gradient Color Scheme
#The Yellow-Purple Gradient Color Scheme palette has 6 colors 
#which are Metallic Yellow (#FFD014), American Gold (#DFB12F), Aztec Gold (#BF914A), Blast-Off Bronze (#A07265), Razzmic Berry (#805280) and Rebecca Purple (#60339B).

Yellow2Purple=c('#FFD014','#DFB12F','#BF914A','#A07265', '#805280' ,'#60339B');
barplot(1:length(Yellow2Purple),col=Yellow2Purple)
N=7
mat=matrix(runif(N^2),N,N)
image(mat,col=Yellow2Purple)

y2p.pal<-colorRampPalette(c("#FFD014","white",'#60339B'))
y2p=y2p.pal(7)
barplot(1:length(y2p),col=y2p)
image(mat,col=y2p,axes=FALSE,ylab="", xlab="")

mat=Matrix::forceSymmetric(mat)
axis.tick=0:N
x=0:nrow(mat) 
y=0:ncol(mat)
image(x,y,mat, axes=FALSE, col=y2p,ylab="", xlab="")
for(i in 1:length(axis.tick)){
  segments(axis.tick[i],min(y),axis.tick[i],max(y))#vertical
  segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
}
box()

pdf('fig1-gene.co.mat.pdf',useDingbats = T)
N=7
for(i in 1:3){
  mat=Matrix::Matrix(runif(N^2),N,N)
  mat=as.matrix(Matrix::forceSymmetric(mat))
  axis.tick=0:N
  x=0:nrow(mat) 
  y=0:ncol(mat) #x and y should be 1 unit longer than dim(mat)
  mat1 <- apply(mat, 2, rev);
  image(y,x,mat1, axes=FALSE, col=y2p,ylab="", xlab="")
  for(i in 1:length(axis.tick)){
    segments(axis.tick[i],min(y),axis.tick[i],max(y))#vertical
    segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
  }
  box()
}
dev.off()

## change to green or blue color panel 
my.col1=RColorBrewer::brewer.pal(7,'Greens')
barplot(1:length(my.col1),col=my.col1)

my.col2=RColorBrewer::brewer.pal(7,'Blues')
barplot(1:length(my.col2),col=my.col2)

pdf('fig1-feature.mat.pdf',useDingbats = T)
N=7
for(i in 1:2){
  if(i==1){my.col=my.col1}else{my.col=my.col2}
  mat=Matrix::Matrix(runif(N^2),N,N)
  mat=as.matrix(Matrix::forceSymmetric(mat))
  axis.tick=0:N
  x=0:nrow(mat) 
  y=0:ncol(mat) #x and y should be 1 unit longer than dim(mat)
  mat1 <- apply(mat, 2, rev);
  image(y,x,mat1, axes=FALSE, col=my.col,ylab="", xlab="")
  for(i in 1:length(axis.tick)){
    segments(axis.tick[i],min(y),axis.tick[i],max(y))#vertical
    segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
  }
  box()
}
dev.off()

## change to red color panel 
pdf("fig1-node.selection.pdf",useDingbats = T)

my.col=RColorBrewer::brewer.pal(7,'Reds')
#barplot(1:length(my.col),col=my.col)
my.col[1]='white'
N=20;M=4;
set.seed(123456)
mat=as.matrix(Matrix::Matrix(rexp(N*M),N,M))
mat=mat/10
mat=t(mat)
mat1=t(apply(mat,1,function(i){
  x=abs((i-mean(i))/sd(i))
  x[x<1.3]=0
  x
}))
axis.tick=0:N
x=0:ncol(mat1) 
y=0:nrow(mat1) #x and y should be 1 unit longer than dim(mat)

image(y,x,mat, axes=FALSE, col=my.col,ylab="", xlab="")
for(i in 1:length(axis.tick)){
  segments(axis.tick[i],min(x),axis.tick[i],max(x))#vertical
  #segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
}
box()

image(y,x,mat1, axes=FALSE, col=my.col,ylab="", xlab="")

for(i in 1:length(axis.tick)){
  segments(axis.tick[i],min(x),axis.tick[i],max(x))#vertical
  #segments(min(x),axis.tick[i],max(x),axis.tick[i])#horizental
}
box()


## module by cell type
N=400;M=5;
set.seed(123456)
mat=as.matrix(Matrix::Matrix(runif(N*M),N,M))
mat1=apply(mat,2,function(i){
  x=abs((i-mean(i))/sd(i))
  x[x<1]=0
  x[x>=1]=1
  x
})
print (pheatmap::pheatmap(mat1,col=c('grey90',"#99000D"),
                   treeheight_row = 0, treeheight_col = 0)
)
dev.off()
