## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9,0.8,0.7,0.6,0.5)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.95,0.85,0.75,0.65,0.55)
pprev <- 0.3
nbadgers <- 1000
print(psens)
print(pspec)
print(pprev)

#simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens)))
colnames(tests) <- c("test1", "test2", "test3", "test4", "test5")

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 1:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}

###choose number of tests to infer results from
ntests <- 1
psensmod<-psens[1:ntests]
pspecmod<-pspec[1:ntests]

## Define parameters you want to report on
params <- c("pprev", "psensmod", "pspecmod")


code  <- nimbleCode({

    for(i in 1:n) {
    for(j in 1:nT) {
      pinf[i, j] <- pprev * psensmod[j] + (1 - pprev) * (1 - pspecmod[j])
      Te[i, j] ~ dbern(pinf[i, j])
    }
  }
 
pprev ~ dunif(0,1)
  
for(j in 1:nT) {
    psensmod[j] ~ dunif(0, 1)
    pspecmod[j] ~ dunif(0, 1)
  }


})


consts <- list(n = nbadgers,
               nT = ntests)

data <- list(
  Te = array(tests[, 1:ntests],dim=c(nbadgers,ntests))
)

inits <- list(
  pprev = runif(1, 0, 1),
  psensmod = runif(ntests, 0.5, 1),
  pspecmod = runif(ntests, 0.5, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run1 <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

plot(run4$samples)
run1summary<-MCMCsummary(run1$samples)
