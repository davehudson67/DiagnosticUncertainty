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
omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
omega <- as.matrix(omega)
omega <- t(omega)

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
  psens = runif(ntests, 0.5, 1),
  pspec = runif(ntests, 0.5, 1)
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

###plotting internal correlations
plot(c(run3.5$samples$chain1[,2],run3.5$samples$chain2[,2],run3.5$samples$chain3[,2])~c(run3.5$samples$chain1[,3],run3.5$samples$chain2[,3],run3.5$samples$chain3[,3]))
check.array<-run3.5$samples$chain1[1:500,]
colnames(check.array)<-c("Prev","Sens1","Sens2","Sens3","Spec1","Spec2","Spec3")
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- (cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 4#/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

pdf("internal correlations.pdf")
pairs(check.array, upper.panel = panel.cor,
      gap=0, row1attop=TRUE,pch=16,col=rgb(0,0,0,0.2))
dev.off()