###############################################
##### BLRCR WRAPPERS ##########################
###############################################


twonormalcovswrapper <- function(data,job,instance,N,betacoef,mysamples){
  ################
  #Create Dataset#
  ################
  covnames="twonormalcovs"
  H=2
  X1 = rnorm(N,0,1)
  X2 = rnorm(N,0,1)
  
  X = cbind(rep(1,N),X1,X2)
  
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
  samplesizen=nrow(myobserveddata)

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
  
  #Conditional Maximum Likelhiood
  cmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
             family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
             data=myobserveddata)
  cmlSE=cmlfit@extra$SE.N.hat
  cmlN=cmlfit@extra$N.hat
  cmlinterval=cmlN+c(-1.96,1.96)*cmlSE
  cmlresult = c(cmlN,cmlinterval)
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  Nboot = rep(NA,numboots)
  for(boot in 1:numboots){
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:4]),as.matrix(myobserveddata[,5:6]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
  bootstrapprobs = sigmoid(as.matrix(cbind(rep(1,samplesizen),myobserveddata[,5:6])) %*% t(matrix(coef(cmlfit),nrow=J)))
  NiHAT=1/(1-apply(1-bootstrapprobs,1,prod))
  di = NiHAT - floor(NiHAT)
  Nboot = rep(NA,1000)
  for(boot in 1:1000){
    bootdata = matrix(ncol=J+H)
    for(i in 1:samplesizen){
      for(nhat in 1:(floor(NiHAT[i])+rbinom(1,1,di[i]))){
        Yrow=t(rep(NA,J))
        for(j in 1:J){
          Yrow[1,j] = rbinom(1,1,prob=bootstrapprobs[i,j])
        }
        Xrow=myobserveddata[i,5:6]
        bootdata = rbind(bootdata,as.matrix(unname(cbind(Yrow,Xrow))))
      }     
    }
    bootdata=bootdata[-1,]
    bootdata=as.data.frame(bootdata[rowSums(bootdata[,1:4])!=0,])
    names(bootdata) <- c("y1","y2","y3","y4","x1","x2")
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  hist(Nboot)
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
  
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
  
  SPBLRCRqs = quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity 2class SP-BLRCR
  HSPBLRCR2=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=2,samples = mysamples)
  
  HSPBLRCR2qs = quantile(HSPBLRCR2$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity SP-BLRCR
  HSPBLRCR=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=20,samples = mysamples)
  
  HSPBLRCRqs = quantile(HSPBLRCR$N,c(0.025,0.5,0.975))
  
  
  ##Norm-BLRCR
  normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                        mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                        mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
  normBLRCRqs = quantile(normBLRCR$N,c(0.025,0.5,0.975))
  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:4])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:4],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinlower,myloglinest,myloglinupper)
  
  
  ##LCMCR
  mydatalcmcr = myobserveddata[,1:4]
  noentries=which(colSums(mydatalcmcr)==0) #find lists with no individuals
  mydatalcmcr = data.frame(
    as.factor(myobserveddata$y1),
    as.factor(myobserveddata$y2),
    as.factor(myobserveddata$y3),
    as.factor(myobserveddata$y4))
  if(length(noentries)>0){
    mydatalcmcr=mydatalcmcr[,-noentries] #remove lists with no individuals
  }
  sampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                   not_in_list_label = '0', K = 20, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 1)
  LCMCRposterior <- lcmCR_PostSampl(sampler, burnin = 1000,
                                    samples = 10000, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.025,0.5,0.975))

  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = 10000, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.025,0.5,0.975))

  
  ########################################
  ### Track Key Population Parameters ###
  TRUEN = N
  TRUEBETA = betacoef
  
  ########################################
  ### Summarize all values to be kept ####
  ########################################
  
  myreturnlist = list("N"=TRUEN,
                      "n"=samplesizen,
                      "covariates"=covnames,
                      "BETA"=TRUEBETA,
                      "condMLCR"= cmlresult,
                      "condMLCRboot" = cmlbootstrap,
                      "condBLRCR"=c(cBLRCRMAP,cBLRCRconverge),
                      "SPBLRCR"=SPBLRCRqs,
                      "HSPBLRCR2"= HSPBLRCR2qs,
                      "HSPBLRCR"= HSPBLRCRqs,
                      "normBLRCR"=normBLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  

  return(myreturnlist) 
}




twochisqcovswrapper <- function(data,job,instance,N,betacoef,mysamples){
  ################
  #Create Dataset#
  ################
  covnames="twochisqcovs"
  H=2
  X1 = rchisq(N,1)
  X2 = rchisq(N,1)
  
  X = cbind(rep(1,N),X1,X2)
  
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
  samplesizen=nrow(myobserveddata)
  
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
  
  #Conditional Maximum Likelhiood
  cmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
              family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
              data=myobserveddata)
  coef(cmlfit)
  cmlSE=cmlfit@extra$SE.N.hat
  cmlN=cmlfit@extra$N.hat
  cmlinterval=cmlN+c(-1.96,1.96)*cmlSE
  cmlresult = c(cmlN,cmlinterval)
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  Nboot = rep(NA,numboots)
  for(boot in 1:numboots){
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:4]),as.matrix(myobserveddata[,5:6]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
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
  
  SPBLRCRqs = quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity 2class SP-BLRCR
  HSPBLRCR2=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                           mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                           mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=2,samples = mysamples)
  
  HSPBLRCR2qs = quantile(HSPBLRCR2$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity SP-BLRCR
  HSPBLRCR=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=20,samples = mysamples)
  
  HSPBLRCRqs = quantile(HSPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Norm-BLRCR
  normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
  normBLRCRqs = quantile(normBLRCR$N,c(0.025,0.5,0.975))
  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:4])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:4],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinlower,myloglinest,myloglinupper)
  
  ##LCMCR
  mydatalcmcr = myobserveddata[,1:4]
  noentries=which(colSums(mydatalcmcr)==0) #find lists with no individuals
  mydatalcmcr = data.frame(
    as.factor(myobserveddata$y1),
    as.factor(myobserveddata$y2),
    as.factor(myobserveddata$y3),
    as.factor(myobserveddata$y4))
  if(length(noentries)>0){
    mydatalcmcr=mydatalcmcr[,-noentries] #remove lists with no individuals
  }
  sampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                   not_in_list_label = '0', K = 20, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 1)
  LCMCRposterior <- lcmCR_PostSampl(sampler, burnin = 1000,
                                    samples = 10000, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.025,0.5,0.975))
  
  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = 10000, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.025,0.5,0.975))
  
  
  ########################################
  ### Track Key Population Parameters ###
  TRUEN = N
  TRUEBETA = betacoef
  
  ########################################
  ### Summarize all values to be kept ####
  ########################################
  
  myreturnlist = list("N"=TRUEN,
                      "n"=samplesizen,
                      "covariates"=covnames,
                      "BETA"=TRUEBETA,
                      "condMLCR"= cmlresult,
                      "condMLCRboot" = cmlbootstrap,
                      "condBLRCR"=c(cBLRCRMAP,cBLRCRconverge),
                      "SPBLRCR"=SPBLRCRqs,
                      "HSPBLRCR2"= HSPBLRCR2qs,
                      "HSPBLRCR"= HSPBLRCRqs,
                      "normBLRCR"=normBLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  
  
  return(myreturnlist) 
}


twogammacovswrapper <- function(data,job,instance,N,betacoef,mysamples){
  ################
  #Create Dataset#
  ################
  covnames="twogammacovs"
  X1 = rgamma(N,3,1)
  X2 = rgamma(N,1,1)
  H=2
  X = cbind(rep(1,N),X1,X2)
  
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
  samplesizen=nrow(myobserveddata)
  
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
  
  #Conditional Maximum Likelhiood
  cmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
              family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
              data=myobserveddata)
  coef(cmlfit)
  cmlSE=cmlfit@extra$SE.N.hat
  cmlN=cmlfit@extra$N.hat
  cmlinterval=cmlN+c(-1.96,1.96)*cmlSE
  cmlresult = c(cmlN,cmlinterval)
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  Nboot = rep(NA,numboots)
  for(boot in 1:numboots){
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:4]),as.matrix(myobserveddata[,5:6]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
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
  
  SPBLRCRqs = quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity 2class SP-BLRCR
  HSPBLRCR2=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                           mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                           mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=2,samples = mysamples)
  
  HSPBLRCR2qs = quantile(HSPBLRCR2$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity SP-BLRCR
  HSPBLRCR=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=20,samples = mysamples)
  
  HSPBLRCRqs = quantile(HSPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Norm-BLRCR
  normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
  normBLRCRqs = quantile(normBLRCR$N,c(0.025,0.5,0.975))
  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:4])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:4],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinlower,myloglinest,myloglinupper)
  
  ##LCMCR
  mydatalcmcr = myobserveddata[,1:4]
  noentries=which(colSums(mydatalcmcr)==0) #find lists with no individuals
  mydatalcmcr = data.frame(
    as.factor(myobserveddata$y1),
    as.factor(myobserveddata$y2),
    as.factor(myobserveddata$y3),
    as.factor(myobserveddata$y4))
  if(length(noentries)>0){
    mydatalcmcr=mydatalcmcr[,-noentries] #remove lists with no individuals
  }
  sampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                   not_in_list_label = '0', K = 20, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 1)
  LCMCRposterior <- lcmCR_PostSampl(sampler, burnin = 1000,
                                    samples = 10000, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.025,0.5,0.975))
  
  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = 10000, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.025,0.5,0.975))
  
  
  ########################################
  ### Track Key Population Parameters ###
  TRUEN = N
  TRUEBETA = betacoef
  
  ########################################
  ### Summarize all values to be kept ####
  ########################################
  
  myreturnlist = list("N"=TRUEN,
                      "n"=samplesizen,
                      "covariates"=covnames,
                      "BETA"=TRUEBETA,
                      "condMLCR"= cmlresult,
                      "condMLCRboot" = cmlbootstrap,
                      "condBLRCR"=c(cBLRCRMAP,cBLRCRconverge),
                      "SPBLRCR"=SPBLRCRqs,
                      "HSPBLRCR2"= HSPBLRCR2qs,
                      "HSPBLRCR"= HSPBLRCRqs,
                      "normBLRCR"=normBLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  
  
  return(myreturnlist) 
}



twomixture3normcovswrapper <- function(data,job,instance,N,betacoef,mysamples){
  ################
  #Create Dataset#
  ################
  covnames="twomixture3normcovs"
  H=2
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

  #plot(X[,2:3])
  
  
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
  samplesizen=nrow(myobserveddata)
  
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
  
  #Conditional Maximum Likelhiood
  cmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
              family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
              data=myobserveddata)
  coef(cmlfit)
  cmlSE=cmlfit@extra$SE.N.hat
  cmlN=cmlfit@extra$N.hat
  cmlinterval=cmlN+c(-1.96,1.96)*cmlSE
  cmlresult = c(cmlN,cmlinterval)
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  Nboot = rep(NA,numboots)
  for(boot in 1:numboots){
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:4]),as.matrix(myobserveddata[,5:6]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
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
  
  SPBLRCRqs = quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity 2class SP-BLRCR
  HSPBLRCR2=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                           mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                           mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=2,samples = mysamples)
  
  HSPBLRCR2qs = quantile(HSPBLRCR2$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity SP-BLRCR
  HSPBLRCR=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=20,samples = mysamples)
  
  HSPBLRCRqs = quantile(HSPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Norm-BLRCR
  normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:4]),x=as.matrix(myobserveddata[,5:6]),
                          mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
  normBLRCRqs = quantile(normBLRCR$N,c(0.025,0.5,0.975))
  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:4])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:4],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinlower,myloglinest,myloglinupper)
  
  ##LCMCR
  mydatalcmcr = myobserveddata[,1:4]
  noentries=which(colSums(mydatalcmcr)==0) #find lists with no individuals
  mydatalcmcr = data.frame(
    as.factor(myobserveddata$y1),
    as.factor(myobserveddata$y2),
    as.factor(myobserveddata$y3),
    as.factor(myobserveddata$y4))
  if(length(noentries)>0){
    mydatalcmcr=mydatalcmcr[,-noentries] #remove lists with no individuals
  }
  sampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                   not_in_list_label = '0', K = 20, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 1)
  LCMCRposterior <- lcmCR_PostSampl(sampler, burnin = 1000,
                                    samples = 10000, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.025,0.5,0.975))
  
  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = 10000, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.025,0.5,0.975))
  
  
  ########################################
  ### Track Key Population Parameters ###
  TRUEN = N
  TRUEBETA = betacoef
  
  ########################################
  ### Summarize all values to be kept ####
  ########################################
  
  myreturnlist = list("N"=TRUEN,
                      "n"=samplesizen,
                      "covariates"=covnames,
                      "BETA"=TRUEBETA,
                      "condMLCR"= cmlresult,
                      "condMLCRboot" = cmlbootstrap,
                      "condBLRCR"=c(cBLRCRMAP,cBLRCRconverge),
                      "SPBLRCR"=SPBLRCRqs,
                      "HSPBLRCR2"= HSPBLRCR2qs,
                      "HSPBLRCR"= HSPBLRCRqs,
                      "normBLRCR"=normBLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  
  
  return(myreturnlist) 
}



onenormwithheterocovswrapper <- function(data,job,instance,N,betacoef,mysamples){
  ################
  #Create Dataset#
  ################
  covnames="onenormwithheterocovs"
  
  
  H=1
  X = cbind(rep(1,N),rnorm(N,0,1))
  W = t(rmultinom(N,1,c(0.65,0.35)))[,-1] #unobserved class
  
  sigmoid <- function(x){
    1/(1+exp(-x))
  }
  
  J=3 #number of lists
  if(betacoef=="heterogeneity"){
    mybeta = matrix(c(-2.5,-1.5,3.0,
                      -2.5,-1.5,3.0,
                      0.5,-1.5,-3.0),nrow=J,byrow=TRUE)
  }else{
    stop("Error, unknown betacoef specified.")
  }
  
  
  
  #Bind the unobserved to the observed
  X=cbind(X,W)
  
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
  colnames(mydata) <- gsub("V","",colnames(mydata))
  
  #Remove unobserved
  myobserveddata = mydata[rowSums(mydata[,1:J])>0,1:4]
  samplesizen=nrow(myobserveddata)
  
  
  ###########################
  ### Set hyperparameters ###
  ###########################
  
  mypriorb=rep(0,2) #prior on intercept and both beta means set to 0
  mypriorB=diag(2) #identity matrix, size 4
  
  ###########
  mypriornu0 = 3
  mypriorLAMBDA0 = as.matrix(1)
  mypriorkappa0 = 1
  mypriorMU0 = as.matrix(0)
  myaalpha = 0.25
  mybalpha = 0.25
  myKstar = 20
  
  ##########################
  ### Run ALgorithms #######
  ##########################
  
  #Conditional Maximum Likelhiood
  cmlfit=vglm(cbind(y1,y2,y3)~x1,
              family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
              data=myobserveddata)
  coef(cmlfit)
  cmlSE=cmlfit@extra$SE.N.hat
  cmlN=cmlfit@extra$N.hat
  cmlinterval=cmlN+c(-1.96,1.96)*cmlSE
  cmlresult = c(cmlN,cmlinterval)
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  Nboot = rep(NA,numboots)
  for(boot in 1:numboots){
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:J]),as.matrix(myobserveddata[,(J+1):(J+H)]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3)~x1,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    Nboot[boot]=bootcmlfit@extra$N.hat
  }
  cmlbootstrap=quantile(Nboot,c(.025,.975))
  
  ##CondBLRCR
  cBLRCR = condBLRCRsolver(as.matrix(myobserveddata[,1:3]),
                           as.matrix(myobserveddata[,4]),mypriorb,mypriorB,
                           gradparam=0.001,maxiter=1000,prior=1)
  cBLRCRMAP = cBLRCR$N
  cBLRCRconverge = cBLRCR$converge
  
  ##SP-BLRCR
  SPBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:3]),x=as.matrix(myobserveddata[,4]),
                        mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                        mypriornu0,myaalpha,mybalpha,Kstar=myKstar,samples = mysamples)
  
  SPBLRCRqs = quantile(SPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity 2class SP-BLRCR
  HSPBLRCR2=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:3]),x=as.matrix(myobserveddata[,4]),
                           mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                           mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=2,samples = mysamples)
  
  HSPBLRCR2qs = quantile(HSPBLRCR2$N,c(0.025,0.5,0.975))
  
  ##Heterogeneity SP-BLRCR
  HSPBLRCR=HSPBLRCRsolver(y=as.matrix(myobserveddata[,1:3]),x=as.matrix(myobserveddata[,4]),
                          mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=myKstar,omegasize=20,samples = mysamples)
  
  HSPBLRCRqs = quantile(HSPBLRCR$N,c(0.025,0.5,0.975))
  
  ##Norm-BLRCR
  normBLRCR=SPBLRCRsolver(y=as.matrix(myobserveddata[,1:3]),x=as.matrix(myobserveddata[,4]),
                          mypriorb,mypriorB,mypriorMU0,mypriorLAMBDA0,mypriorkappa0,
                          mypriornu0,myaalpha,mybalpha,Kstar=1,samples = mysamples)
  
  normBLRCRqs = quantile(normBLRCR$N,c(0.025,0.5,0.975))
  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:3])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:3],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinlower,myloglinest,myloglinupper)
  
  ##LCMCR
  mydatalcmcr = myobserveddata[,1:3]
  noentries=which(colSums(mydatalcmcr)==0) #find lists with no individuals
  mydatalcmcr = data.frame(
    as.factor(myobserveddata$y1),
    as.factor(myobserveddata$y2),
    as.factor(myobserveddata$y3))
  if(length(noentries)>0){
    mydatalcmcr=mydatalcmcr[,-noentries] #remove lists with no individuals
  }
  sampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                   not_in_list_label = '0', K = 20, a_alpha = 0.25, b_alpha = 0.25,
                   seed = 'auto', buffer_size = 10000, thinning = 1)
  LCMCRposterior <- lcmCR_PostSampl(sampler, burnin = 1000,
                                    samples = 10000, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.025,0.5,0.975))
  
  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = 10000, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.025,0.5,0.975))
  
  
  ########################################
  ### Track Key Population Parameters ###
  TRUEN = N
  TRUEBETA = betacoef
  
  ########################################
  ### Summarize all values to be kept ####
  ########################################
  
  myreturnlist = list("N"=TRUEN,
                      "n"=samplesizen,
                      "covariates"=covnames,
                      "BETA"=TRUEBETA,
                      "condMLCR"= cmlresult,
                      "condMLCRboot" = cmlbootstrap,
                      "condBLRCR"=c(cBLRCRMAP,cBLRCRconverge),
                      "SPBLRCR"=SPBLRCRqs,
                      "HSPBLRCR2"= HSPBLRCR2qs,
                      "HSPBLRCR"= HSPBLRCRqs,
                      "normBLRCR"=normBLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  
  
  return(myreturnlist) 
}