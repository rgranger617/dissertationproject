dependencylevels <- function(data,job,instance,N,betacoef,mysamples){
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
  
  ##CondBLRCR
  cBLRCR = condBLRCRsolver(as.matrix(myobserveddata[,1:4]),
                           as.matrix(myobserveddata[,5:6]),mypriorb,mypriorB,
                           gradparam=0.001,maxiter=1000,prior=1)
  cBLRCRMAP = cBLRCR$N
  cBLRCRconverge = cBLRCR$converge
  
  
  #parametric bootstrap for confidence intervals (see Zwane 2003)
  numboots = 1000
  NbootCML = rep(NA,numboots)
  NbootBCML = rep(NA,numboots)
  for(boot in 1:numboots){
    print(boot)
    bootdata=as.data.frame(bootstrapRCR(as.matrix(myobserveddata[,1:4]),as.matrix(myobserveddata[,5:6]),matrix(coef(cmlfit),nrow=J)))
    bootcmlfit=vglm(cbind(y1,y2,y3,y4)~x1+x2,
                    family=posbernoulli.tb(drop.b = FALSE ~ 0, parallel.t = TRUE~ 0),
                    data=bootdata)
    NbootCML[boot]=bootcmlfit@extra$N.hat

    cBLRCRboot = condBLRCRsolver(as.matrix(bootdata[,1:4]),
                             as.matrix(bootdata[,5:6]),mypriorb,mypriorB,
                             gradparam=0.001,maxiter=1000,prior=1)
    if(cBLRCRboot$converge==1){
      NbootBCML[boot]=cBLRCRboot$N
    }
  }
  cmlbootstrap=c(cmlN,quantile(NbootCML,.025),quantile(NbootCML,.975))
  bootconverge=1-mean(is.na(NbootBCML)) #percentage of bootstrap samples that converged
  Bcmlbootstrap=c(cBLRCRMAP,quantile(NbootBCML,.025,na.rm=TRUE),quantile(NbootBCML,.975,na.rm=TRUE))
  
  
  ##BLRCR-Fixed Covariate Distribution
  myformula = cbind(y1,y2,y3,y4)~x1+x2
  BLRCR=SPBLRCRsolver2(myformula,df=myobserveddata,Homega=1,
                              covmethod="fixed",
                              samples = mysamples, burnin=1000,thinning=1)
  
  BLRCRqs = quantile(BLRCR$N,c(0.5,0.025,0.975))
  
  

  
  ##Log Linear
  loglinmodel=closedpMS.t(myobserveddata[,1:4])
  besthierarchy=which.min(loglinmodel$results[,6]) #smallest BIC
  llm = closedpCI.t(myobserveddata[,1:4],mX=names(besthierarchy))
  myloglinest=llm$CI[1]
  myloglinlower=llm$CI[2]
  myloglinupper=llm$CI[3]
  myloglin = c(myloglinest,myloglinlower,myloglinupper)
  
  
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
                                    samples = mysamples, thinning = 1,
                                    output = FALSE)
  
  LCMCRqs=quantile(LCMCRposterior,c(0.5,0.025,0.975))
  
  
  ## Independent
  indsampler <- lcmCR(captures = mydatalcmcr, tabular = FALSE, in_list_label = '1',
                      not_in_list_label = '0', K = 1, a_alpha = 0.25, b_alpha = 0.25,
                      seed = 'auto', buffer_size = 10000, thinning = 1)
  indposterior <- lcmCR_PostSampl(indsampler, burnin = 1000,
                                  samples = mysamples, thinning = 1,
                                  output = FALSE)
  
  indqs=quantile(indposterior,c(0.5,0.025,0.975))
  
  
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
                      "condBLRCR"=Bcmlbootstrap,
                      "condBLRCRconvergence"=c(cBLRCRconverge,bootconverge),
                      "BLRCR"=BLRCRqs,
                      "loglin"=myloglin,
                      "LCMCR"=LCMCRqs,
                      "Independent"=indqs)
  
  
  return(myreturnlist) 
}



