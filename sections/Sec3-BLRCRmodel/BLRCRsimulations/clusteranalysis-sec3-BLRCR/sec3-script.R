# 1. This file estimates parameters for all performances in parallel using
# a cluster. See https://mllg.github.io/batchtools/ for setup instructions
# to run on your own machine or an appropriate managed HPC system.
library(mvtnorm) #needed for multivariate normal
library(BayesLogit) #pollygamma sampling
library(LCMCR)
library(Rcapture)
library(BLRCR)
library(batchtools)
library(data.table)
library(VGAM)


#Create Registry
reg=makeExperimentRegistry("sec3simulationregistry", 
                           packages=c('mvtnorm','BayesLogit','BLRCR','Rcapture','LCMCR',"VGAM"), 
                           source = "sec3-functions.R",
                           seed = 20231508) #YYYYMMDD
reg$cluster.functions=makeClusterFunctionsSlurm("mybatchjobs.tmpl")
#reg=loadRegistry("mymodelselectsims",writeable=TRUE)

#Add the Algorithm
addAlgorithm(name="dependencylevels", fun = dependencylevels)

#Add Data
addProblem("nodata")



#Add Experiments for Data Simulated using Log Log link
addExperiments(algo.designs=list(dependencylevels = CJ(N = 2000,
                                                    betacoef=c("negative","moderate","positive","independent"),
                                                    mysamples=10000)),repls=1000)
addExperiments(algo.designs=list(dependencylevels = CJ(N = c(200,500,1000,2000,5000,10000),
                                                       betacoef=c("moderate"),
                                                       mysamples=10000)),repls=1000)


################
summarizeExperiments(by = c("problem", "algorithm","N","mysamples","betacoef"))

submitJobs(resources = list(ppn=1, 
                            nodes=1, 
                            ncpus=1,
                            memory='16gb', 
                            walltime=86400)) #24 hours
getStatus()



submitJobs(ids=findExpired(),
           resources = list(ppn=1, 
                            nodes=1, 
                            ncpus=1,
                            memory='16gb', 
                            walltime=172800)) #48 hours

submitJobs(ids=findNotSubmitted(),
           resources = list(ppn=1, 
                            nodes=1, 
                            ncpus=1,
                            memory='16gb', 
                            walltime=36000)) #10 hours
#findExpired()
#findErrors()
#getErrorMessages(1)
#submitJobs(ids=findErrors(),resources = list(ppn=1, nodes=1, memory='16gb', walltime='02:00:00'))
#getLog(46)
#writeLines(getLog(46))
#getJobStatus(46)



#Use when jobs need chunked (puts some chunk of "jobs" into a single "job")
library(data.table)
ids = getJobTable(reg = reg)[, .(job.id, problem, algorithm)]
ids[, chunk := chunk(job.id, n.chunks = 500), by = "problem"]
ids[, chunk := .GRP, by = c("problem", "chunk")]
dcast(ids, chunk ~ problem)
submitJobs(ids=ids,
           resources = list(ppn=1, 
                            nodes=1, 
                            ncpus=1,
                            memory='16gb', 
                            walltime=86400))

#Use when jobs need chunked when some come back expired
library(data.table)
ids = getJobTable(reg = reg)[, .(job.id, problem, algorithm)]
ids=ids[findExpired(),]
ids[, chunk := chunk(job.id, n.chunks = 500), by = "problem"]
ids[, chunk := .GRP, by = c("problem", "chunk")]
dcast(ids, chunk ~ problem)
submitJobs(ids=ids,
           resources = list(ppn=1, 
                            nodes=1, 
                            ncpus=1,
                            memory='16gb', 
                            walltime=36000))

