## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.75,0.7,0.9,0.65,0.6)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.55,0.7,0.65,0.75,0.6)
pprev <- 0.3
nbadgers <- 1000
print(psens)
print(pspec)
print(pprev)

###estimate beta parameters for priors
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

Se.beta<-estBetaParams(psens,0.1^2)
Sp.beta<-estBetaParams(pspec,0.1^2)

Se.beta
Sp.beta

#define range
p = seq(0, 1, length=100)

#create plot of Beta distribution 
plot(p, dbeta(p, 13.0625, 10.6875), type='l')

#simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)
sum(inf)

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens)))
colnames(tests) <- c("test1", "test2", "test3", "test4", "test5")

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 1:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}

############manipulate test numbers from here###############

#how many tests?
###choose number of tests to infer results from
ntests <- 3

#get table of frequencies
binop<-2^seq(0,ntests-1)
testbin<-tests[,1:ntests]%*%binop
testcounts<-tabulate(testbin+1,nbins=2^ntests)
testcounts

#create omega matrix of binary sequence
omega<-t(expand.grid(replicate(ntests, c(0,1), simplify=FALSE)))

## Define parameters you want to report on
params <- c("pprev", "psens", "pspec")

code  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
psens[1]~dunif(0.5,1)
pspec[1]~dunif(0.5,1)
for(j in 2:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(runvagueprior <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
runvaguepriorsummary<-MCMCsummary(runvagueprior$samples)
#save(run2, file = "run2.RData")

code2  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
psens[1]~dbeta(2.765625,0.921875)
pspec[1]~dbeta(2.853125,2.334375)
for(j in 2:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model2 <- nimbleModel(code2, constants = consts, data = data, inits = inits)

cModel2 <- compileNimble(model2)

config2 <- configureMCMC(model2)

rMCMC2 <- buildMCMC(config2)

cMCMC2 <- compileNimble(rMCMC2, project = model2)

system.time(runmediumprior <- runMCMC(cMCMC2, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
runmediumpriorsummary<-MCMCsummary(runmediumprior$samples)
#save(run2, file = "run2.RData")

code3  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
psens[1]~dbeta(13.3125,4.4375)
pspec[1]~dbeta(13.0625,10.6875)
for(j in 2:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model3 <- nimbleModel(code3, constants = consts, data = data, inits = inits)

cModel3 <- compileNimble(model3)

config3 <- configureMCMC(model3)

rMCMC3 <- buildMCMC(config3)

cMCMC3 <- compileNimble(rMCMC3, project = model3)

system.time(rungoodprior <- runMCMC(cMCMC3, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
rungoodpriorsummary<-MCMCsummary(rungoodprior$samples)
#save(run2, file = "run2.RData")

code4  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
psens[1]~dbeta(3.5625,0.1875)
pspec[1]~dbeta(13.3125,4.4375)
for(j in 2:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model4 <- nimbleModel(code4, constants = consts, data = data, inits = inits)

cModel4 <- compileNimble(model4)

config4 <- configureMCMC(model4)

rMCMC4 <- buildMCMC(config4)

cMCMC4 <- compileNimble(rMCMC4, project = model4)

system.time(runbadprior <- runMCMC(cMCMC4, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
runbadpriorsummary<-MCMCsummary(runbadprior$samples)
#save(run2, file = "run2.RData")

code5  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
psens[1]~dbeta(0.178125,0.009375)
pspec[1]~dbeta(2.765625,0.921875)
for(j in 2:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model5 <- nimbleModel(code5, constants = consts, data = data, inits = inits)

cModel5 <- compileNimble(model5)

config5 <- configureMCMC(model5)

rMCMC5 <- buildMCMC(config5)

cMCMC5 <- compileNimble(rMCMC5, project = model5)

system.time(runmediumbadprior <- runMCMC(cMCMC5, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
runmediumbadpriorsummary<-MCMCsummary(runmediumbadprior$samples)
#save(run2, file = "run2.RData")


#pdf("prior manipulations 20221116.pdf")
postmeans.prev<-c(runvaguepriorsummary$mean[1],runmediumpriorsummary$mean[1],rungoodpriorsummary$mean[1],runmediumbadpriorsummary$mean[1],runbadpriorsummary$mean[1])
#postsds.prev<-c(run1summary$sd[1],run2summary$sd[1],run3summary$sd[1],run4summary$sd[1],run5summary$sd[1])
postmeans.Se1<-c(runvaguepriorsummary$mean[2],runmediumpriorsummary$mean[2],rungoodpriorsummary$mean[2],runmediumbadpriorsummary$mean[2],runbadpriorsummary$mean[2])
#postsds.Se1<-c(run1summary$sd[2],run2summary$sd[2],run3summary$sd[2],run4summary$sd[2],run5summary$sd[2])
postmeans.Sp1<-c(runvaguepriorsummary$mean[5],runmediumpriorsummary$mean[5],rungoodpriorsummary$mean[5],runmediumbadpriorsummary$mean[5],runbadpriorsummary$mean[5])
#postsds.Se1<-c(run1summary$sd[3],run2summary$sd[4],run3summary$sd[5],run4summary$sd[6],run5summary$sd[7])
plot(postmeans.prev~seq(1:5),ylim=c(0,1),xlim=c(0.5,5.5),xlab="Prior information",ylab="Inferred parameter",pch=16,cex=2,type="n",xaxt="n")
axis(1, at=1:5, labels=c("Vague","Weak/Good","Strong/Good","Weak/Bad","Strong/Bad"))

width.adj<-5

#prevalence densities
dd<-density(c(runvagueprior$samples$chain1[,1],runvagueprior$samples$chain2[,1],runvagueprior$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1),rev(dd$y/(width.adj*max(dd$y))+1)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(runmediumprior$samples$chain1[,1],runmediumprior$samples$chain2[,1],runmediumprior$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2),rev(dd$y/(width.adj*max(dd$y))+2)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(rungoodprior$samples$chain1[,1],rungoodprior$samples$chain2[,1],rungoodprior$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3),rev(dd$y/(width.adj*max(dd$y))+3)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(runmediumbadprior$samples$chain1[,1],runmediumbadprior$samples$chain2[,1],runmediumbadprior$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4),rev(dd$y/(width.adj*max(dd$y))+4)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(runbadprior$samples$chain1[,1],runbadprior$samples$chain2[,1],runbadprior$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5),rev(dd$y/(width.adj*max(dd$y))+5)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
points(postmeans.prev~c((seq(1:5))),pch=16,cex=2)

#Se1 densities
xadjustment<--0.2
dd<-density(c(runvagueprior$samples$chain1[,2],runvagueprior$samples$chain2[,2],runvagueprior$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(runmediumprior$samples$chain1[,2],runmediumprior$samples$chain2[,2],runmediumprior$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(rungoodprior$samples$chain1[,2],rungoodprior$samples$chain2[,2],rungoodprior$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(runmediumbadprior$samples$chain1[,2],runmediumbadprior$samples$chain2[,2],runmediumbadprior$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(runbadprior$samples$chain1[,2],runbadprior$samples$chain2[,2],runbadprior$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
points(postmeans.Se1~c((seq(1:5)+xadjustment)),pch=16,cex=2)

#Sp1 densities
xadjustment<-0.2
dd<-density(c(runvagueprior$samples$chain1[,5],runvagueprior$samples$chain2[,5],runvagueprior$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(runmediumprior$samples$chain1[,5],runmediumprior$samples$chain2[,5],runmediumprior$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(rungoodprior$samples$chain1[,5],rungoodprior$samples$chain2[,5],rungoodprior$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(runmediumbadprior$samples$chain1[,5],runmediumbadprior$samples$chain2[,5],runmediumbadprior$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(runbadprior$samples$chain1[,5],runbadprior$samples$chain2[,5],runbadprior$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
points(postmeans.Sp1~c((seq(1:5)+xadjustment)),pch=16,cex=2)

abline(0.3,0,col=rgb(1,0,0,0.5),lwd=2,lty=2)
abline(0.75,0,col=rgb(0,1,0,0.5),lwd=2,lty=2)
abline(0.55,0,col=rgb(0,0,1,0.5),lwd=2,lty=2)

#dev.off()


