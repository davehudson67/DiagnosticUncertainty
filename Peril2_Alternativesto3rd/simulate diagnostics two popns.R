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
npopns<-2
pprev <- c(0.3,0.45)
nbadgers <- 1000
print(psens)
print(pspec)
print(pprev)

#simulate infection status of the badgers
inf <- array(0, dim = c(nbadgers, npopns))
for(p in 1:npopns){
inf[, p] <- rbinom(nbadgers, 1, pprev[p])
}

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens), npopns))
colnames(tests) <- c("test1", "test2", "test3", "test4", "test5")

#for each badger and each test, simulate test outcome
for(p in 1:npopns){
for(i in 1:nbadgers){
  for(j in 1:length(psens)){
    tests[i, j, p] <- rbinom(1, 1, ifelse(inf[i, p] == 1, psens[j], 1-pspec[j]))}}
}
############manipulate test numbers from here###############

#how many tests?
###choose number of tests to infer results from
ntests <- 2

#get table of frequencies
binop<-2^seq(0,ntests-1)
testbin<-array(0,dim=c(nbadgers,npopns))
testcounts<-array(0,dim=c(2^ntests,npopns))
for(p in 1:npopns){
testbin[,p]<-tests[,1:ntests,p]%*%binop
testcounts[,p]<-tabulate(testbin[,p]+1,nbins=2^ntests)
}
testcounts

#create omega matrix of binary sequence
#create omega matrix of binary sequence
omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
omega <- as.matrix(omega)
omega <- t(omega)

## Define parameters you want to report on
params <- c("pprev", "psens", "pspec")

code  <- nimbleCode({

    for(p in 1:nP){
    for(i in 1:n) {
	pinf[i,p] <- pprev[p]*prod(omega[1:nT,i]*psens[1:nT]+(1-omega[1:nT,i])*(1-psens[1:nT]))+
		    (1-pprev[p])*prod(omega[1:nT,i]*(1-pspec[1:nT])+(1-omega[1:nT,i])*(pspec[1:nT]))
	}
    Te[1:n,p] ~ dmulti(pinf[1:n,p],N)
    }
  
for(p in 1:npopns){ 
pprev[p] ~ dunif(0,1)
  }
  
for(j in 1:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }


})

consts <- list(n = 2^ntests,
               nT = ntests, N = nbadgers, omega=omega, nP = npopns)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(npopns, 0, 1),
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run2_0.3_0.45 <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

#plot(run3$samples)
run2_0.3_0.45summary<-MCMCsummary(run2_0.3_0.45$samples)
run2_0.3_0.45summary


###only do this for the flipping parameters problem
run1summary<-run2_0.3_0.5summary
run2summary<-run2_0.3_0.45summary
run3summary<-run2_0.3_0.4summary
run4summary<-run2_0.3_0.35summary
run5summary<-run2_0.3_0.3summary

run1<-run2_0.3_0.5
run2<-run2_0.3_0.45
run3<-run2_0.3_0.4
run4<-run2_0.3_0.35
run5<-run2_0.3_0.3
#####

pdf("two populations 20220908.pdf")
postmeans.prev1<-c(run1summary$mean[1],run2summary$mean[1],run3summary$mean[1],run4summary$mean[1],run5summary$mean[1])
postsds.prev1<-c(run1summary$sd[1],run2summary$sd[1],run3summary$sd[1],run4summary$sd[1],run5summary$sd[1])
postmeans.prev2<-c(run1summary$mean[2],run2summary$mean[2],run3summary$mean[2],run4summary$mean[2],run5summary$mean[2])
postsds.prev2<-c(run1summary$sd[2],run2summary$sd[2],run3summary$sd[2],run4summary$sd[2],run5summary$sd[2])
postmeans.Se1<-c(run1summary$mean[3],run2summary$mean[3],run3summary$mean[3],run4summary$mean[3],run5summary$mean[3])
postsds.Se1<-c(run1summary$sd[3],run2summary$sd[3],run3summary$sd[3],run4summary$sd[3],run5summary$sd[3])
postmeans.Sp1<-c(run1summary$mean[5],run2summary$mean[5],run3summary$mean[5],run4summary$mean[5],run5summary$mean[5])
postsds.Se1<-c(run1summary$sd[5],run2summary$sd[5],run3summary$sd[5],run4summary$sd[5],run5summary$sd[5])
plot(postmeans.prev1~seq(1:5),ylim=c(0,1),xlim=c(0.5,5.5),xlab="Prevalence in second population",ylab="Inferred parameter",pch=16,cex=2,type="n",xaxt="n")
axis(1, at=1:5, labels=c("0.5","0.45","0.4","0.35","0.3"))
width.adj<-5

#prevalence densities pop1
xadjustment<--0.1
dd<-density(c(run1$samples$chain1[,1],run1$samples$chain2[,1],run1$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run2$samples$chain1[,1],run2$samples$chain2[,1],run2$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run3$samples$chain1[,1],run3$samples$chain2[,1],run3$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run4$samples$chain1[,1],run4$samples$chain2[,1],run4$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run5$samples$chain1[,1],run5$samples$chain2[,1],run5$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
points(postmeans.prev1~c((seq(1:5)+xadjustment)),pch=16,cex=2)

#prevalence densities pop2
xadjustment<-0.1
dd<-density(c(run1$samples$chain1[,2],run1$samples$chain2[,2],run1$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run2$samples$chain1[,2],run2$samples$chain2[,2],run2$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run3$samples$chain1[,2],run3$samples$chain2[,2],run3$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run4$samples$chain1[,2],run4$samples$chain2[,2],run4$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run5$samples$chain1[,2],run5$samples$chain2[,2],run5$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+5+xadjustment),rev(dd$y/(width.adj*max(dd$y))+5+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
points(postmeans.prev2~c((seq(1:5)+xadjustment)),pch=16,cex=2)



#Se1 densities
xadjustment<--0.2
dd<-density(c(run1$samples$chain1[,3],run1$samples$chain2[,3],run1$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run2$samples$chain1[,3],run2$samples$chain2[,3],run2$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run3$samples$chain1[,3],run3$samples$chain2[,3],run3$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+3+xadjustment),rev(dd$y/(width.adj*max(dd$y))+3+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run4$samples$chain1[,3],run4$samples$chain2[,3],run4$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+4+xadjustment),rev(dd$y/(width.adj*max(dd$y))+4+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run5$samples$chain1[,3],run5$samples$chain2[,3],run5$samples$chain3[,3]))
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
abline(0.55,-0.05,col=rgb(1,0,0,0.5),lwd=2,lty=2)

dev.off()


