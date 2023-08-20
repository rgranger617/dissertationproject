library(mvtnorm) #needed for multivariate normal
library(MCMCpack) #needed for inverse wishart
library(BayesLogit) #pollygamma sampling
library(cubature) #numerical integration
library(MASS)
library(LCMCR)
library(ellipse) #helpful for plotting ellipses
library(BLRCR)

RMSE <- function(x,N){
  sqrt(mean((x-N)^2))
}


N=2000
betacoef="moderate"
mysamples=10000
dist="mixnormal"
  
################
#Create Dataset#
################
if(dist=="normal"){
  X1 = rnorm(N,0,1)
  X2 = rnorm(N,0,1)
  X = cbind(rep(1,N),X1,X2)
}else if(dist=="chisq"){
  X1 = rchisq(N,1)
  X2 = rchisq(N,1)
  X = cbind(rep(1,N),X1,X2)
}else if(dist=="gamma"){
  X1 = rgamma(N,3,1)
  X2 = rgamma(N,1,1)
  X = cbind(rep(1,N),X1,X2)
}else if(dist=="mixnormal"){
  z=sample.int(3,N,replace=TRUE,prob=c(0.5,0.2,0.3))
  X <- matrix(rep(NA,3*N),ncol=3)
  for(i in 1:N){
    if(z[i]==1){
      X[i,]=c(1,rmvnorm(1,mean=c(0,0),
                        sigma=matrix(c(1,0,
                                       0,1),byrow=TRUE,ncol=2)))
    }else if(z[i]==2){
      X[i,]=c(1,rmvnorm(1,mean=c(-2,-2),
                        sigma=matrix(c(.5,-.35,
                                       -.35,.5),byrow=TRUE,ncol=2)))
    }else if(z[i]==3){
      X[i,]=c(1,rmvnorm(1,mean=c(2,2),
                        sigma=matrix(c(.5,.45,
                                       .45,.5),byrow=TRUE,ncol=2)))
    }
  }
}

  
sigmoid <- function(x){
  1/(1+exp(-x))
}
  
J=4 #number of lists
if(betacoef=="negative"){
  mybeta = matrix(c(-2,-1,1,
                    -2,1,-1,
                    -2,1,1,
                    -2,-1,-1),nrow=J,byrow=TRUE)
}else if(betacoef=="moderate"){
  mybeta = matrix(c(-2,-1,1,
                    -2,1,-1,
                    -2,-1,1,
                    -2,1,-1),nrow=J,byrow=TRUE)
}else if(betacoef=="positive"){
  mybeta = matrix(c(-2,-1,1,
                    -2,-1,1,
                    -2,-1,1,
                    -2,-1,1),nrow=J,byrow=TRUE)
}else if(betacoef=="independent"){
  mybeta = matrix(c(-2,0,0,
                    -2,0,0,
                    -2,0,0,
                    -2,0,0),nrow=J,byrow=TRUE)
}
  
myprobs = sigmoid(X%*%t(mybeta))
Y <- matrix(rep(NA,J*N),ncol=J)
for(j in 1:J){
  Y[,j] = rbinom(n=N,size=1,prob=myprobs[,j])
}
  
mydata = data.frame(
  "y"=Y,
  "x"=X[,-1]
)
colnames(mydata) <- gsub("\\.","",colnames(mydata))
colnames(mydata) <- gsub("X","",colnames(mydata))
  
#Remove unobserved
myobserveddata = mydata[rowSums(mydata[,1:J])>0,]
mymissingdata = mydata[rowSums(mydata[,1:J])==0,]
samplesizen=nrow(myobserveddata)
  

plot(mydata[,5:6])
plot(myobserveddata[,5:6])

###########################
### Set hyperparameters ###
###########################
  
mypriorb=rep(0,3) #prior on intercept and both beta means set to 0
mypriorB=diag(3) #identity matrix, size 6
  
###########
mypriornu0 = 3
mypriorLAMBDA0 = diag(2)
mypriorkappa0 = 1
mypriorMU0 = rep(0,2)
myaalpha = 0.25
mybalpha = 0.25
myKstar = 20
  
##########################
### Run ALgorithms #######
##########################
  
##CondBLRCR
cBLRCR = condBLRCRsolver(as.matrix(myobserveddata[,1:4]),
                          as.matrix(myobserveddata[,5:6]),mypriorb,mypriorB,
                          gradparam=0.001,maxiter=1000,prior=1)
cBLRCRMAP = cBLRCR$N
cBLRCRconverge = cBLRCR$converge
  
##SP-BLRCR
SPBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                        mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                        mypriornu0,myaalpha,mybalpha,Kstar=myKstar,samples = mysamples)
  
quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
##Norm-BLRCR
normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                        mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                        mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
quantile(normBLRCR$N,c(0.025,0.5,0.975))
  

#################################################

saveRDS(mydata,"mydataALL.rds")
saveRDS(mymissingdata,"mymissingdata.rds")
saveRDS(myobserveddata,"myobserveddata.rds")
saveRDS(SPBLRCR,"SPBLRCR.rds")
saveRDS(normBLRCR,"normBLRCR.rds")

################################################

##Create a plot##
par(mfrow=c(1,3))
plot(myobserveddata[,5:6],
     xlab="X1",ylab="X2",main="Simulated Posterior - SP-BLRCR",
     xlim=c(-4,4),ylim=c(-4,4),pch=16)
points(SPBLRCR$Xmis[9001:10000,],col="blue")
legend(-4,4,legend=c("Observed","Simulated Missing"),
       col=c("black","blue"),pch=c(16,1),bg="lightblue",text.font=4,
       text.width = 4)
plot(myobserveddata[,5:6],
     xlab="X1",ylab="X2",main="Simulated Posterior - BLRCR-Normal",
     xlim=c(-4,4),ylim=c(-4,4),pch=16)
points(normBLRCR$Xmis[9001:10000,],col="blue")
legend(-4,4,legend=c("Observed","Simulated Missing"),
       col=c("black","blue"),pch=c(16,1),bg="lightblue",text.font=4,
       text.width = 4)
plot(myobserveddata[,5:6],xlab="X1",ylab="X2",main="True Population",
     xlim=c(-4,4),ylim=c(-4,4),pch=16)
points(mymissingdata[,5:6],col="blue")
legend(-4,4,legend=c("Observed","Actual Missing"),
       col=c("black","blue"),pch=c(16,1),bg="lightblue",text.font=4,
       text.width = 4)




################################
## Some Plots ##################
################################

hist(SPBLRCR$N)
hist(normBLRCR$N)
plot(SPBLRCR$N,type="l")
plot(SPBLRCR$N[1:10000],type="l")
plot(normBLRCR$N,type="l")

plot(SPBLRCR$Xmis[99500:100000,])
plot(mydata[,5:6])
plot(myobserveddata[,5:6])

SPBLRCR$beta

SPBLRCR$MUK
SPBLRCR$SIGMAK
SPBLRCR$PIK

hist(SPBLRCR$Xmis)
plot(SPBLRCR$Xmis)
plot(mymissingdata[,5:6])
sample.int(10000,1131)
plot(SPBLRCR$Xmis[sample.int(10000,1131),])

for(i in 1:200){
  first=(i-1)*500+1
  last=first+499
  plot(SPBLRCR$Xmis[first:last,],main=paste(first))
}



library(ggplot2)
library(gridExtra)
#posteriorcovariates = data.frame(SPBLRCR$Xmis[sample.int(10000,1131),])
posteriorcovariates = data.frame(SPBLRCR$Xmis)
ggp1=ggplot(posteriorcovariates, aes(x=X1, y=X2) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(-5,5)+ylim(-5,5)+ggtitle("Posterior of Missing Covariates")+
  theme(legend.position='none') 
ggp2=ggplot(mymissingdata, aes(x=x1, y=x2) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(-5,5)+ylim(-5,5)+ggtitle("Actual Missing Coariates")+
  theme(legend.position='none') 
grid.arrange(ggp1,ggp2,ncol=2) 



library(ggplot2)
library(cowplot)
posteriorcovariates = data.frame(SPBLRCR$Xmis,"type"=rep("Posterior",mysamples))
names(posteriorcovariates)=c("x1","x2","type")
mymis=data.frame(mymissingdata[,5:6],"type"=rep("Actual Missing",nrow(mymissingdata)))
datformissing=rbind(posteriorcovariates,mymis)
ggp1=ggplot(datformissing, aes(x=x1, y=x2) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(-5,5)+ylim(-5,5)+ggtitle("Posterior of Missing Covariates")+
  facet_grid(cols=vars(type))+
  theme_minimal_grid()
ggp1
