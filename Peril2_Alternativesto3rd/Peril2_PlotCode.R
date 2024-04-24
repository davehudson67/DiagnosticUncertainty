## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)

## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9, 0.9) #c(0.9,0.8,0.7,0.7,0.7)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75, 0.75) #c(0.75,0.65,0.55,0.55,0.55)
npopns<-2
pprev <- c(0.3, 0.5)
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
colnames(tests) <- c("test1", "test2")

#for each badger and each test, simulate test outcome
for(p in 1:npopns){
  for(i in 1:nbadgers){
    for(j in 1:length(psens)){
      tests[i, j, p] <- rbinom(1, 1, ifelse(inf[i, p] == 1, psens[j], 1 - pspec[j]))}}
}

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
      pinf[i, p] <- pprev[p] * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT]))+
        (1 - pprev[p]) * prod(omega[1:nT, i] * (1 - pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
    }
    Te[1:n, p] ~ dmulti(pinf[1:n, p], N)
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

system.time(run <- runMCMC(cMCMC, 
                                     niter = 50000, 
                                     nburnin = 2400, 
                                     nchains = 3, 
                                     progressBar = TRUE, 
                                     summary = TRUE, 
                                     samplesAsCodaMCMC = TRUE, 
                                     thin = 1))

#saveRDS(run, "Peril2_Alternativesto3rd/run_030.rds")
#saveRDS(run, "Peril2_Alternativesto3rd/run_035.rds")
#saveRDS(run, "Peril2_Alternativesto3rd/run_040.rds")
#saveRDS(run, "Peril2_Alternativesto3rd/run_045.rds")
saveRDS(run, "Peril2_Alternativesto3rd/run_050.rds")

r50 <- readRDS("Peril2_Alternativesto3rd/run_050.rds")
r45 <- readRDS("Peril2_Alternativesto3rd/run_045.rds")
r40 <- readRDS("Peril2_Alternativesto3rd/run_040.rds")
r35 <- readRDS("Peril2_Alternativesto3rd/run_035.rds")
r30 <- readRDS("Peril2_Alternativesto3rd/run_030.rds")

r50_sum <- as.matrix(r50$summary$all.chains[, 1])
r50_samp <- as.matrix(r50$samples)
r50_samp <- r50_samp %>%
  as.data.frame() %>%
  rename(PrevalenceA = 1, PrevalenceB = 2, Sensitivity = 3, Sensitivity2 = 4, Specificity = 5, Specificity2 = 6) %>%
  select(PrevalenceA, PrevalenceB, Sensitivity , Specificity) %>%
  pivot_longer(cols = c(PrevalenceA, PrevalenceB, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.5, 0.9, 0.75), times = 571200/4)) %>%
  mutate(Mean = rep(c(r50_sum[1], r50_sum[2], r50_sum[3], r50_sum[4]), times = 571200/4)) %>%
  mutate(PrevalenceB_truth = 0.5)

r45_sum <- as.matrix(r45$summary$all.chains[, 1])
r45_samp <- as.matrix(r45$samples)
r45_samp <- r45_samp %>%
  as.data.frame() %>%
  rename(PrevalenceA = 1, PrevalenceB = 2, Sensitivity = 3, Sensitivity2 = 4, Specificity = 5, Specificity2 = 6) %>%
  select(PrevalenceA, PrevalenceB, Sensitivity , Specificity) %>%
  pivot_longer(cols = c(PrevalenceA, PrevalenceB, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.45, 0.9, 0.75), times = 571200/4)) %>%
  mutate(Mean = rep(c(r45_sum[1], r45_sum[2], r45_sum[3], r45_sum[4]), times = 571200/4)) %>%
  mutate(PrevalenceB_truth = 0.45)

r40_sum <- as.matrix(r40$summary$all.chains[, 1])
r40_samp <- as.matrix(r40$samples)
r40_samp <- r40_samp %>%
  as.data.frame() %>%
  rename(PrevalenceA = 1, PrevalenceB = 2, Sensitivity = 3, Sensitivity2 = 4, Specificity = 5, Specificity2 = 6) %>%
  select(PrevalenceA, PrevalenceB, Sensitivity , Specificity) %>%
  pivot_longer(cols = c(PrevalenceA, PrevalenceB, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.4, 0.9, 0.75), times = 571200/4)) %>%
  mutate(Mean = rep(c(r40_sum[1], r40_sum[2], r40_sum[3], r40_sum[4]), times = 571200/4)) %>%
  mutate(PrevalenceB_truth = 0.4)

r35_sum <- as.matrix(r35$summary$all.chains[, 1])
r35_samp <- as.matrix(r35$samples)
r35_samp <- r35_samp %>%
  as.data.frame() %>%
  rename(PrevalenceA = 1, PrevalenceB = 2, Sensitivity = 3, Sensitivity2 = 4, Specificity = 5, Specificity2 = 6) %>%
  select(PrevalenceA, PrevalenceB, Sensitivity , Specificity) %>%
  pivot_longer(cols = c(PrevalenceA, PrevalenceB, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.35, 0.9, 0.75), times = 571200/4)) %>%
  mutate(Mean = rep(c(r35_sum[1], r35_sum[2], r35_sum[3], r35_sum[4]), times = 571200/4)) %>%
  mutate(PrevalenceB_truth = 0.35)

r30_sum <- as.matrix(r30$summary$all.chains[, 1])
r30_samp <- as.matrix(r30$samples)
r30_samp <- r30_samp %>%
  as.data.frame() %>%
  rename(PrevalenceA = 1, PrevalenceB = 2, Sensitivity = 3, Sensitivity2 = 4, Specificity = 5, Specificity2 = 6) %>%
  select(PrevalenceA, PrevalenceB, Sensitivity , Specificity) %>%
  pivot_longer(cols = c(PrevalenceA, PrevalenceB, Sensitivity, Specificity),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Truth = rep(c(0.3, 0.3, 0.9, 0.75), times = 571200/4)) %>%
  mutate(Mean = rep(c(r30_sum[1], r30_sum[2], r30_sum[3], r30_sum[4]), times = 571200/4)) %>%
  mutate(PrevalenceB_truth = 0.3)

data <- rbind(r50_samp, r45_samp, r40_samp, r35_samp, r30_samp)
saveRDS(data, "Peril2_Alternativesto3rd/plotdata.rds")
data <- readRDS("Peril2_Alternativesto3rd/plotdata.rds")

ggplot(data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin(position = "dodge") +
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 2) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(~ PrevalenceB_truth) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  labs(x = "Prevalence Population B", fill = "Inferred parameter", colour = "Inferred parameter")







