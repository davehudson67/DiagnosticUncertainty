## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(tidyverse)
library(MCMCvis)
library(coda)
library(mcmcplots)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.75,0.7,0.9,0.65,0.6)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.55,0.7,0.65,0.75,0.6)
pprev <- 0.3
nbadgers <- 500
print(psens)
print(pspec)
print(pprev)

niter <- 100000
nburnin <- 24000

###estimate beta parameters for priors
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

Se.beta <- estBetaParams(psens, 0.1^2)
Sp.beta <- estBetaParams(pspec, 0.1^2)

Se.beta
Sp.beta

#define range
p = seq(0, 1, length = 100)

#create plots of Beta distribution 
plot(p, dbeta(p, 2.765625, 0.921875), type = 'l') #medium sens (0.75)
plot(p, dbeta(p, 2.853125, 2.334375), type='l') #medium spec (0.55)

plot(p, dbeta(p, 13.3125, 4.4375), type = 'l') #good sens
plot(p, dbeta(p, 13.0625, 10.6875), type = 'l') #good spec

plot(p, dbeta(p, 3.5625, 0.1875), type = 'l') #bad sens
plot(p, dbeta(p, 13.3125, 4.4375), type = 'l') #bad spec

plot(p, dbeta(p, 0.178125, 0.009375), type = 'l') #medbad sens
plot(p, dbeta(p, 2.765625, 0.921875), type = 'l') #medbad spec


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
ntests <- 5

#get table of frequencies
binop <- 2^seq(0, ntests - 1)
testbin <- tests[, 1:ntests] %*% binop
testcounts <- tabulate(testbin + 1, nbins = 2^ntests)
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
  
psens[1] ~ dunif(0,1)
pspec[1] ~ dunif(0,1)

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
                           niter = niter, 
                           nburnin = nburnin, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

runvaguepriorsummary <- MCMCsummary(runvagueprior$samples)

saveRDS(runvagueprior, "Peril7_PriorInfo/VaguePrior.rds")

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
                           niter = niter, 
                           nburnin = nburnin, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

runmediumpriorsummary<-MCMCsummary(runmediumprior$samples)
saveRDS(runmediumprior, "Peril7_PriorInfo/MediumPrior.rds")

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
                           niter = niter, 
                           nburnin = nburnin, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

rungoodpriorsummary <- MCMCsummary(rungoodprior$samples)

saveRDS(rungoodprior, "Peril7_PriorInfo/GoodPrior.rds")

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
                           niter = niter, 
                           nburnin = nburnin, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

runbadpriorsummary<-MCMCsummary(runbadprior$samples)
saveRDS(runbadprior, "Peril7_PriorInfo/BadPrior.rds")

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
                           niter = niter, 
                           nburnin = nburnin, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

runmediumbadpriorsummary<-MCMCsummary(runmediumbadprior$samples)

saveRDS(runmediumbadprior, "Peril7_PriorInfo/MediumBadPrior.rds")

parms <- c("pprev", "psens[1]", "pspec[1]")
library(mcmcplots)

## Load samples
runbadprior <- readRDS("Peril7_PriorInfo/BadPrior.rds")
mcmcplot(runbadprior$samples, parms = parms)
rungoodprior <- readRDS("Peril7_PriorInfo/GoodPrior.rds")
runmediumbadprior <- readRDS("Peril7_PriorInfo/MediumBadPrior.rds")
runmediumprior <- readRDS("Peril7_PriorInfo/MediumPrior.rds")
runvagueprior <- readRDS("Peril7_PriorInfo/VaguePrior.rds")

## combine data
runbadprior_samples <- as.matrix(runbadprior$samples) %>%
  as.data.frame() %>%
  select(1, 2, 7) %>%
  mutate(Prior = "Strong/Bad") %>%
  rename(Prevalence = pprev, Sensitivity = "psens[1]", Specificity = "pspec[1]", Prior = Prior) %>%
  pivot_longer(cols = c(Prevalence, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.75, 0.55), times = 8784000/3))

rungoodprior_samples <- as.matrix(rungoodprior$samples) %>%
  as.data.frame() %>%
  select(1, 2, 7) %>%
  mutate(Prior = "Strong/Good") %>%
  rename(Prevalence = pprev, Sensitivity = "psens[1]", Specificity = "pspec[1]", Prior = Prior) %>%
  pivot_longer(cols = c(Prevalence, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.75, 0.55), times = 8784000/3))

runvagueprior_samples <- as.matrix(runvagueprior$samples) %>%
  as.data.frame() %>%
  select(1, 2, 7) %>%
  mutate(Prior = "Vague") %>%
  rename(Prevalence = pprev, Sensitivity = "psens[1]", Specificity = "pspec[1]", Prior = Prior) %>%
  pivot_longer(cols = c(Prevalence, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.75, 0.55), times = 8784000/3))

runmediumbadprior_samples <- as.matrix(runmediumbadprior$samples) %>%
  as.data.frame() %>%
  select(1, 2, 7) %>%
  mutate(Prior = "Weak/Bad") %>%
  rename(Prevalence = pprev, Sensitivity = "psens[1]", Specificity = "pspec[1]", Prior = Prior) %>%
  pivot_longer(cols = c(Prevalence, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.75, 0.55), times = 8784000/3))

runmediumprior_samples <- as.matrix(runmediumprior$samples) %>%
  as.data.frame() %>%
  select(1, 2, 7) %>%
  mutate(Prior = "Weak/Good") %>%
  rename(Prevalence = pprev, Sensitivity = "psens[1]", Specificity = "pspec[1]", Prior = Prior) %>%
  pivot_longer(cols = c(Prevalence, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.75, 0.55), times = 8784000/3))

plot_data <- rbind(runbadprior_samples, rungoodprior_samples, runmediumbadprior_samples, runmediumprior_samples, runvagueprior_samples)

# Plotting
ggplot(plot_data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin() +
  stat_summary(fun = "mean", geom = "point") +
  labs(x = "Prior Information", y = "Inferred parameter") +
  facet_grid( ~ Prior) + 
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15))
    

