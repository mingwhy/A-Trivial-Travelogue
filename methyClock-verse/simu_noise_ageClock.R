

n.feature=200;
n.sample=100;
n.age.step=100; #one step one sample,age 1~100, 100 samples

# ground state (all share the same ground state, before adding noise)
set.seed(230213)
ground=runif(n.feature,min=0,max=1);
#ground=rep(0.5,n.feature);
#ground=rep(0.6,n.feature);

give.me.one.replicate<-function(){
  samples=matrix(rep(ground,n.sample),nrow=n.feature) #sample by feature matrix
  sum(samples[,1]==samples[,2])
  
  # add variation to each sample, N(0,0.01)
  x=sapply(1:n.sample,function(i) rnorm(n.feature,mean=0,sd=0.01))
  dim(x) #1000 x 100 
  samples=samples+as.matrix(x)
  
  # add noise function
  add_noise<-function(state=state,step=step){
    for(i in 1:step){
      state=state+rnorm(length(state),mean=0,sd=0.05)
      #state=state+rnorm(length(state),mean=0,sd=0.1)
      state[state>1]=1;
      state[state<0]=0;
    }
    #state[state>1]=1;
    #state[state<0]=0;
    return(state)
  }
  #dim(samples) #feature by sample matrix
  #add_noise(state=samples[,1],step=10)
  
  aged.samples=samples
  for(i in 1:ncol(samples)){
    aged.samples[,i]<-add_noise(state=samples[,i],step=i)
  }
  #dim(aged.samples)
  #par(mfrow=c(2,2))
  #lapply(c(1,10,50,100),function(i) hist(aged.samples[,i],main=i))
  return(aged.samples)
}
train.samples<-lapply(1:3,function(i) give.me.one.replicate())
test.samples<-lapply(1:3,function(i) give.me.one.replicate())

## train clock
library(glmnet)

TrainData = t(Reduce(`cbind`,train.samples)) #sample by feature matrix
TestData = t(Reduce(`cbind`,test.samples)) 

plot(ground,TrainData[1,])
plot(ground,TrainData[20,])
plot(ground,TrainData[100,])
#sum(abs(TrainData[2,]-TrainData[1,]))
#sum(abs(TrainData[10,]-TrainData[1,]))
#sum(abs(TrainData[100,]-TrainData[1,]))
TrainAge = rep(1:100,3)
TestAge = rep(1:100,3)

par(mfrow=c(2,2))
#Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainData, TrainAge, nfolds=5,alpha=0.5, family="gaussian")
fit = glmnet(TrainData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)

#Examine full model
plot(TrainAge,predict(fit,TrainData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainData,s = cv$lambda.min))

plot(TestAge,predict(fit,TestData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestData,s = cv$lambda.min))


## plot ground state feature ~ coeff
coeffs<-coef(fit,s=cv$lambda.min)
feature.coeffs=as.vector(coeffs[-1]) #remove intercept
plot(ground,feature.coeffs,cex=0.5,pch=16,xlab='methy.freq.at.age.0',ylab='regression coefficient of each site')


