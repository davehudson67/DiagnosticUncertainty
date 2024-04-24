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
params <- c("pprev", "psens", "pspec","gammaSe", "gammaSp")  # 

code  <- nimbleCode({

    for(i in 1:n) {
	pinf[i] <- pprev*((prod(omega[1:2,i]*psens[1:2]+(1-omega[1:2,i])*(1-psens[1:2]))+gammaSe*((-1)^(omega[1,i]+omega[2,i])))*
			(omega[3,i]*psens[3]+(1-omega[3,i])*(1-psens[3])))+
		    (1-pprev)*((prod(omega[1:2,i]*(1-pspec[1:2])+(1-omega[1:2,i])*(pspec[1:2]))+gammaSp*((-1)^(omega[1,i]+omega[2,i])))*
			(omega[3,i]*(1-pspec[3])+(1-omega[3,i])*(pspec[3])))
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

system.time(run3cov <- runMCMC(cMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run3cov, "Peril4_NonIndependenceDiagTests/3Tests.rds")

################################################################################

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
omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
omega <- as.matrix(omega)
omega <- t(omega)

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

saveRDS(run4cov, "Peril4_NonIndependenceDiagTests/4Tests.rds")

run3 <- readRDS("Peril4_NonIndependenceDiagTests/3Tests.rds")
run4 <- readRDS("Peril4_NonIndependenceDiagTests/4Tests.rds")

r3_sum <- run3$summary
r3_sum <- as.matrix(r3_sum$all.chains[, 1])
r3_samp <- as.matrix(run3$samples)

r3_samp <- r3_samp %>%
  as.data.frame() %>%
  rename(ggammaSE = 1, ggammaSP = 2, pprev = 3, psens = 4, psens2 = 5, psens3 = 6, pspec = 7, pspec2 = 8, pspec3 = 9) %>%
  pivot_longer(cols = c(ggammaSE, ggammaSP, pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 3) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.05, 0.1, 0.3, 0.9, 0.75), times = 714000/5)) %>%
  mutate(Mean = rep(c(r3_sum[1], r3_sum[2], r3_sum[3], r3_sum[4], r3_sum[7]), times = 714000/5))

r4_sum <- run4$summary
r4_sum <- as.matrix(r4_sum$all.chains[, 1])
r4_samp <- as.matrix(run4$samples)

r4_samp <- r4_samp %>%
  as.data.frame() %>%
  rename(ggammaSE = 1, ggammaSP = 2, pprev = 3, psens = 4, psens2 = 5, psens3 = 6, psens4 = 7, pspec = 8, pspec2 = 9, pspec3 = 10, pspec4 = 11) %>%
  pivot_longer(cols = c(ggammaSE, ggammaSP, pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 4) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.05, 0.1, 0.3, 0.9, 0.75), times = 714000/5)) %>%
  mutate(Mean = rep(c(r4_sum[1], r4_sum[2], r4_sum[3], r4_sum[4], r4_sum[8]), times = 714000/5))


data <- rbind(r3_samp, r4_samp)
data$NumberOfTests <- as.factor(data$NumberOfTests)
#data <- data %>%
#  mutate(Inferred_parameter = rep(c("Prevalence", "Sensitivity", "Specificity"), times = 2142000/3))

ggplot(data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin(position = "dodge") +
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 2) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(~ NumberOfTests) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  labs(x = "Number of Tests", fill = "Inferred parameter", colour = "Inferred parameter") +
  scale_fill_discrete(labels = c("Co-dependency Sensitivity", "Co-dependency Specificity", "Prevalence", "Sensitivity", "Specificity")) +
  scale_colour_discrete(labels = c("Co-dependency Sensitivity", "Co-dependency Specificity", "Prevalence", "Sensitivity", "Specificity"))



