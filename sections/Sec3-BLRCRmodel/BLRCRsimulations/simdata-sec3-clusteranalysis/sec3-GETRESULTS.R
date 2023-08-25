library(batchtools)

#Load Registry
reg=loadRegistry("mymodelselectsims",writeable=TRUE)

loadResult(1)

simcount=9000
mysimlist = list("N"=rep(NA,simcount),
                 "n"=rep(NA,simcount),
                 "BETA"=rep(NA,simcount),
                 "covariates"=rep(NA,simcount),
                 "condMLCR"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "condMLCRboot"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "BcondMLCRboot"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "BcondMLCRconverge"=array(rep(NA,simcount),dim=c(2,simcount)),
                 "BLRCR"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "loglin"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "LCMCR"=array(rep(NA,simcount),dim=c(3,simcount)),
                 "independent"=array(rep(NA,simcount),dim=c(3,simcount))
                 )


for(sim in 1:simcount){
  mysimlist$N[sim]=loadResult(sim)$N
  mysimlist$n[sim]=loadResult(sim)$n
  mysimlist$BETA[sim]=loadResult(sim)$BETA
  mysimlist$covariates[sim]=loadResult(sim)$covariates
  mysimlist$condMLCR[,sim]=loadResult(sim)$condMLCR
  mysimlist$condMLCRboot[,sim]=loadResult(sim)$condMLCRboot
  mysimlist$BcondMLCRboot[,sim]=loadResult(sim)$condBLRCR
  mysimlist$BcondMLCRconverge[,sim]=loadResult(sim)$condBLRCRconvergence
  mysimlist$BLRCR[,sim]=loadResult(sim)$BLRCR
  mysimlist$loglin[1,sim]=loadResult(sim)$loglin[1]
  mysimlist$loglin[2,sim]=loadResult(sim)$loglin[2]
  mysimlist$loglin[3,sim]=loadResult(sim)$loglin[3]
  mysimlist$loglin[,sim]=loadResult(sim)$loglin
  mysimlist$LCMCR[,sim]=loadResult(sim)$LCMCR
  mysimlist$independent[,sim]=loadResult(sim)$Independent
}





saveRDS(mysimlist,"sec3-simdata.RDS")

