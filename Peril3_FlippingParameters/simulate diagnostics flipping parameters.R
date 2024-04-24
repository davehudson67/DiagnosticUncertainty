## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9,0.8,0.7,0.7,0.7)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75,0.65,0.55,0.55,0.55)
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
omega<-t(array(expand.grid(replicate(ntests, c(0,1), simplify=FALSE))))

## Define parameters you want to report on
params <- c("pprev", "psens", "pspec")

code  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev)*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
  
 
pprev ~ dunif(0,1)
  
for(j in 1:nT) {
    psens[j] ~ dunif(0, 1)
    pspec[j] ~ dunif(0, 1)
  }


})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega=omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(ntests, 0, 1),
  pspec = runif(ntests, 0, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run3.5 <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
run3.5summary<-MCMCsummary(run3.5$samples)
#save(run2, file = "run2.RData")

###only do this for the flipping parameters problem
run1summary<-run3.1summary
run2summary<-run3.2summary
run3summary<-run3.3summary
run4summary<-run3.4summary
run5summary<-run3.5summary

run1<-run3.1
run2<-run3.2
run3<-run3.3
run4<-run3.4
run5<-run3.5
#####

pdf("flipping parameters 20220907.pdf")
postmeans.prev<-c(run1summary$mean[1],run2summary$mean[1],run3summary$mean[1],run4summary$mean[1],run5summary$mean[1])
postsds.prev<-c(run1summary$sd[1],run2summary$sd[1],run3summary$sd[1],run4summary$sd[1],run5summary$sd[1])
postmeans.Se1<-c(run1summary$mean[2],run2summary$mean[2],run3summary$mean[2],run4summary$mean[2],run5summary$mean[2])
postsds.Se1<-c(run1summary$sd[2],run2summary$sd[2],run3summary$sd[2],run4summary$sd[2],run5summary$sd[2])
postmeans.Sp1<-c(run1summary$mean[5],run2summary$mean[5],run3summary$mean[5],run4summary$mean[5],run5summary$mean[5])
postsds.Se1<-c(run1summary$sd[5],run2summary$sd[5],run3summary$sd[5],run4summary$sd[5],run5summary$sd[5])
plot(postmeans.prev~seq(1:5),ylim=c(0,1),xlim=c(0.5,5.5),xlab="Realisation",ylab="Inferred parameter",pch=16,cex=2)
width.adj<-5

#prevalence densities
dd<-density(c(run1$samples$chain1[,1],run1$samples$chain2[,1],run1$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1),rev(dd$y/(width.adj*max(dd$y))+1)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run2$samples$chain1[,1],run2$samples$chain2[,1],run2$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2),rev(dd$y/(width.adj*max(dd$y))+2)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run3$samples$chain1[,1],run3$samples$chain2[,1],run3$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3),rev(dd$y/(width.adj*max(dd$y))+3)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run4$samples$chain1[,1],run4$samples$chain2[,1],run4$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4),rev(dd$y/(width.adj*max(dd$y))+4)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run5$samples$chain1[,1],run5$samples$chain2[,1],run5$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5),rev(dd$y/(width.adj*max(dd$y))+5)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
points(postmeans.prev~c((seq(1:5))),pch=16,cex=2)

#Se1 densities
xadjustment<--0.2
dd<-density(c(run1$samples$chain1[,2],run1$samples$chain2[,2],run1$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run2$samples$chain1[,2],run2$samples$chain2[,2],run2$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run3$samples$chain1[,2],run3$samples$chain2[,2],run3$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run4$samples$chain1[,2],run4$samples$chain2[,2],run4$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run5$samples$chain1[,2],run5$samples$chain2[,2],run5$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
points(postmeans.Se1~c((seq(1:5)+xadjustment)),pch=16,cex=2)

#Sp1 densities
xadjustment<-0.2
dd<-density(c(run1$samples$chain1[,5],run1$samples$chain2[,5],run1$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(run2$samples$chain1[,5],run2$samples$chain2[,5],run2$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(run3$samples$chain1[,5],run3$samples$chain2[,5],run3$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(run4$samples$chain1[,5],run4$samples$chain2[,5],run4$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(run5$samples$chain1[,5],run5$samples$chain2[,5],run5$samples$chain3[,5]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
points(postmeans.Sp1~c((seq(1:5)+xadjustment)),pch=16,cex=2)

abline(0.3,0,col=rgb(1,0,0,0.5),lwd=2,lty=2)
abline(0.9,0,col=rgb(0,1,0,0.5),lwd=2,lty=2)
abline(0.75,0,col=rgb(0,0,1,0.5),lwd=2,lty=2)

dev.off()


