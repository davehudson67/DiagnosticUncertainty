## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.8,0.6,0.75,0.85,0.55)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75,0.65,0.55,0.85,0.7)
pprev <- 0.3
nbadgers <- 1000
print(psens)
print(pspec)
print(pprev)

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

#bounds on covariance plus
lower.Se<-(psens[1]-1)*(1-psens[2])+0.001
upper.Se<-min(psens[1],psens[2])-psens[1]*psens[2]-0.001
#bounds on covariance minus
lower.Sp<-(pspec[1]-1)*(1-pspec[2])
upper.Sp<-min(pspec[1],pspec[2])-pspec[1]*pspec[2]

###test the covariance bounds
p00inf<-(1-psens[1])*(1-psens[2])+lower.Se
p10inf<-psens[1]*(1-psens[2])-upper.Se
p01inf<-(1-psens[1])*psens[2]-upper.Se
p11inf<-psens[1]*psens[2]+lower.Se

p00uninf<-pspec[1]*pspec[2]+lower.Sp
p10uninf<-(1-pspec[1])*pspec[2]-upper.Sp
p01uninf<-pspec[1]*(1-pspec[2])-upper.Sp
p11uninf<-(1-pspec[1])*(1-pspec[2])+lower.Sp


gammaSe<-0.05
gammaSp<-0.1
p00inf<-(1-psens[1])*(1-psens[2])+gammaSe
p10inf<-psens[1]*(1-psens[2])-gammaSe
p01inf<-(1-psens[1])*psens[2]-gammaSe
p11inf<-psens[1]*psens[2]+gammaSe

p00uninf<-pspec[1]*pspec[2]+gammaSp
p10uninf<-(1-pspec[1])*pspec[2]-gammaSp
p01uninf<-pspec[1]*(1-pspec[2])-gammaSp
p11uninf<-(1-pspec[1])*(1-pspec[2])+gammaSp


####work out how to manipulate the test outcomes with covariance terms!!!
twotests<-matrix(c(0,0,1,0,0,1,1,1),nrow=4,byrow=T)
outcomes<-array(0,dim=c(length(inf),4))
for(i in 1:length(inf)){
    if(inf[i] == 1){outcomes[i,]<-rmultinom(1, size=1, prob=c(p00inf,p10inf,p01inf,p11inf))} else {outcomes[i,]<-rmultinom(1, size=1, prob=c(p00uninf,p10uninf,p01uninf,p11uninf))}
	tests[i,1:2]<-twotests[which(outcomes[i,]==1),]
}

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 3:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}


############manipulate test numbers from here###############

#how many tests?
###choose number of tests to infer results from
ntests <- 4

#get table of frequencies
binop<-2^seq(0,ntests-1)
testbin<-tests[,1:ntests]%*%binop
testcounts<-tabulate(testbin+1,nbins=2^ntests)
testcounts

#create omega matrix of binary sequence
omega<-t(array(expand.grid(replicate(ntests, c(0,1), simplify=FALSE))))

## Define parameters you want to report on
params <- c("pprev", "psens", "pspec","gammaSe", "gammaSp")  # 

code  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*((prod(omega[1:2,i]*psens[1:2]+(1-omega[1:2,i])*(1-psens[1:2]))+gammaSe*((-1)^(omega[1,i]+omega[2,i])))*
			(prod(omega[3:nT,i]*psens[3:nT]+(1-omega[3:nT,i])*(1-psens[3:nT]))))+
		    (1-pprev)*((prod(omega[1:2,i]*(1-pspec[1:2])+(1-omega[1:2,i])*(pspec[1:2]))+gammaSp*((-1)^(omega[1,i]+omega[2,i])))*
			(prod(omega[3:nT,i]*(1-pspec[3:nT])+(1-omega[3:nT,i])*(pspec[3:nT]))))
	}
    Te[1:n] ~ dmulti(pinf[1:n],N)
 

pprev ~ dunif(0,1)
gammaSe~dunif(((psens[1]-1)*(1-psens[2])+0.001),(min(psens[1],psens[2])-psens[1]*psens[2]-0.001))
gammaSp~dunif(((pspec[1]-1)*(1-pspec[2])+0.001),(min(pspec[1],pspec[2])-pspec[1]*pspec[2]-0.001))

for(j in 1:nT) {
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
  pspec = runif(ntests, 0.5, 1),
  gammaSe = runif(1,-0.01,0.01),
  gammaSp = runif(1,-0.01,0.01)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model,monitors=params)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run4cov <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run4covsummary<-MCMCsummary(run4cov$samples)
run4covsummary

plot(run3cov$samples)
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

pdf("covariance models 20221101.pdf")

par(mfrow=c(1,2))
postmeans.prev<-c(run4covsummary$mean[3],run3covsummary$mean[3])
postsds.prev<-c(run1summary$sd[1],run2summary$sd[1],run3summary$sd[1],run4summary$sd[1],run5summary$sd[1])
postmeans.Se1<-c(run4covsummary$mean[4],run3covsummary$mean[4])
postsds.Se1<-c(run1summary$sd[2],run2summary$sd[2],run3summary$sd[2],run4summary$sd[2],run5summary$sd[2])
postmeans.Sp1<-c(run4covsummary$mean[8],run3covsummary$mean[7])
postsds.Se1<-c(run1summary$sd[3],run2summary$sd[4],run3summary$sd[5],run4summary$sd[6],run5summary$sd[7])
postmeans.gammaSe<-c(run4covsummary$mean[1],run3covsummary$mean[1])
postmeans.gammaSp<-c(run4covsummary$mean[2],run3covsummary$mean[2])
plot(postmeans.prev~seq(1:2),ylim=c(0,1),xlim=c(0.5,2.5),xlab="Number of tests",ylab="Inferred parameter",pch=16,cex=2,xaxt="n")
axis(1, at=1:2, labels=c("4","3"))
abline(0.3,0,lwd=2,lty=2,col=rgb(1,0,0,0.5))
abline(0.8,0,lwd=2,lty=2,col=rgb(0,1,0,0.5))
abline(0.75,0,lwd=2,lty=2,col=rgb(0,0,1,0.5))
width.adj<-5

#prevalence densities
dd<-density(c(run4cov$samples$chain1[,3],run4cov$samples$chain2[,3],run4cov$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1),rev(dd$y/(width.adj*max(dd$y))+1)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
dd<-density(c(run3cov$samples$chain1[,3],run3cov$samples$chain2[,3],run3cov$samples$chain3[,3]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2),rev(dd$y/(width.adj*max(dd$y))+2)),c(dd$x,rev(dd$x)),col=rgb(1,0,0,0.5))
points(postmeans.prev~c((seq(1:2))),pch=16,cex=2)

#Se1 densities
xadjustment<--0.2
dd<-density(c(run4cov$samples$chain1[,4],run4cov$samples$chain2[,4],run4cov$samples$chain3[,4]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
dd<-density(c(run3cov$samples$chain1[,4],run3cov$samples$chain2[,4],run3cov$samples$chain3[,4]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,0,0.5))
points(postmeans.Se1~c((seq(1:2)+xadjustment)),pch=16,cex=2)

#Sp1 densities
xadjustment<-0.2
dd<-density(c(run4cov$samples$chain1[,8],run4cov$samples$chain2[,8],run4cov$samples$chain3[,8]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
dd<-density(c(run3cov$samples$chain1[,7],run3cov$samples$chain2[,7],run3cov$samples$chain3[,7]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,0,1,0.5))
points(postmeans.Sp1~c((seq(1:2)+xadjustment)),pch=16,cex=2)

plot(postmeans.gammaSe~seq(1:2),ylim=c(-0.2,0.25),xlim=c(0.5,2.5),xlab="Number of tests",ylab="Inferred parameter",pch=16,cex=2,xaxt="n",type="n")
axis(1, at=1:2, labels=c("4","3"))
width.adj<-5
abline(0.05,0,lwd=2,lty=2,col=rgb(1,0.5,0,0.5))
abline(0.1,0,lwd=2,lty=2,col=rgb(0,1,1))

#gammaSe densities
xadjustment<--0.2
dd<-density(c(run4cov$samples$chain1[,1],run4cov$samples$chain2[,1],run4cov$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0.5,0,0.5))
dd<-density(c(run3cov$samples$chain1[,1],run3cov$samples$chain2[,1],run3cov$samples$chain3[,1]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(1,0.5,0,0.5))
points(postmeans.gammaSe~c((seq(1:2)+xadjustment)),pch=16,cex=2)

#gammaSp densities
xadjustment<-0.2
dd<-density(c(run4cov$samples$chain1[,2],run4cov$samples$chain2[,2],run4cov$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+1+xadjustment),rev(dd$y/(width.adj*max(dd$y))+1+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,1,0.5))
dd<-density(c(run3cov$samples$chain1[,2],run3cov$samples$chain2[,2],run3cov$samples$chain3[,2]))
polygon(c((-(dd$y/(width.adj*max(dd$y)))+2+xadjustment),rev(dd$y/(width.adj*max(dd$y))+2+xadjustment)),c(dd$x,rev(dd$x)),col=rgb(0,1,1,0.5))
points(postmeans.gammaSp~c((seq(1:2)+xadjustment)),pch=16,cex=2)

dev.off()


