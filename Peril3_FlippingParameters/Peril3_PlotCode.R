## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)

#set.seed(11)
## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9, 0.9, 0.9)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75, 0.75, 0.75)
pprev <- 0.3
nbadgers <- 1000
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)  

#simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens)))
colnames(tests) <- c("test1", "test2", "test3")

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 1:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}

############manipulate test numbers from here###############

#how many tests?
###choose number of tests to infer results from
ntests <- 3

#get table of frequencies
binop <- 2^seq(0, ntests-1)
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
    pinf[i] <- pprev * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT])) +
      (1 - pprev) * prod(omega[1:nT, i] * (1 - pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
  }
  
  Te[1:n] ~ dmulti(pinf[1:n], N)
  
  
  pprev ~ dunif(0,1)
  
  for(j in 1:nT) {
    psens[j] ~ dunif(0, 1)
    pspec[j] ~ dunif(0, 1)
  }
  
})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega = omega)

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

system.time(run <- runMCMC(cMCMC, 
                           niter = 500000, 
                           nburnin = 24000, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, file = "Peril3_FlippingParameters/run1.rds")
saveRDS(run, file = "Peril3_FlippingParameters/run2.rds")
saveRDS(run, file = "Peril3_FlippingParameters/run3.rds")
saveRDS(run, file = "Peril3_FlippingParameters/run4.rds")
saveRDS(run, file = "Peril3_FlippingParameters/run5.rds")

run1 <- readRDS("Peril3_FlippingParameters/run1.rds")
run2 <- readRDS("Peril3_FlippingParameters/run2.rds")
run3 <- readRDS("Peril3_FlippingParameters/run3.rds")
run4 <- readRDS("Peril3_FlippingParameters/run4.rds")
run5 <- readRDS("Peril3_FlippingParameters/run5.rds")

r1_sum <- run1$summary
r1_sum <- as.matrix(r1_sum$all.chains[, 1])
r1_samp <- as.matrix(run1$samples)

r1_samp <- r1_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Realisation = 1) %>%
  select(Inferred_parameter, Estimate, Realisation) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r1_sum[1], r1_sum[2], r1_sum[5]), times = 4284000/3))

r2_sum <- run2$summary
r2_sum <- as.matrix(r2_sum$all.chains[, 1])
r2_samp <- as.matrix(run2$samples)

r2_samp <- r2_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Realisation = 2) %>%
  select(Inferred_parameter, Estimate, Realisation) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r2_sum[1], r2_sum[2], r2_sum[5]), times = 4284000/3))

r3_sum <- run3$summary
r3_sum <- as.matrix(r3_sum$all.chains[, 1])
r3_samp <- as.matrix(run3$samples)

r3_samp <- r3_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Realisation = 3) %>%
  select(Inferred_parameter, Estimate, Realisation) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r3_sum[1], r3_sum[2], r3_sum[5]), times = 4284000/3))

r4_sum <- run4$summary
r4_sum <- as.matrix(r4_sum$all.chains[, 1])
r4_samp <- as.matrix(run4$samples)

r4_samp <- r4_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Realisation = 4) %>%
  select(Inferred_parameter, Estimate, Realisation) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r4_sum[1], r4_sum[2], r4_sum[5]), times = 4284000/3))

r5_sum <- run5$summary
r5_sum <- as.matrix(r5_sum$all.chains[, 1])
r5_samp <- as.matrix(run5$samples)

r5_samp <- r5_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(Realisation = 5) %>%
  select(Inferred_parameter, Estimate, Realisation) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r5_sum[1], r5_sum[2], r5_sum[5]), times = 4284000/3))


data <- rbind(r1_samp, r2_samp, r3_samp, r4_samp, r5_samp)
data$Realisation <- as.factor(data$Realisation)
data <- data %>%
  mutate(Inferred_parameter = rep(c("Prevalence", "Sensitivity", "Specificity"), times = 21420000/3))

ggplot(data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin(position = "dodge") +
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 2) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(~ Realisation) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  labs(x = "Realisation", fill = "Inferred parameter", colour = "Inferred parameter")

##########################################################################################################

#set.seed(11)
## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9, 0.9, 0.1)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75, 0.75, 0.05)
pprev <- 0.3
nbadgers <- 1000
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)  

#simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens)))
colnames(tests) <- c("test1", "test2", "test3")

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 1:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}

############manipulate test numbers from here###############

#how many tests?
###choose number of tests to infer results from
ntests <- 3

#get table of frequencies
binop <- 2^seq(0, ntests-1)
testbin <- tests[, 1:ntests] %*% binop
testcounts <- tabulate(testbin + 1, nbins = 2^ntests)
testcounts

#create omega matrix of binary sequence
omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
omega <- as.matrix(omega)
omega <- t(omega)

## Define parameters you want to report on
params <- c("pprev", "psens", "pspec")

code <- nimbleCode({
  
  for(i in 1:n) {
    pinf[i] <- pprev * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT])) +
      (1 - pprev) * prod(omega[1:nT, i] * (1 - pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
  }
  
  Te[1:n] ~ dmulti(pinf[1:n], N)
  
  pprev ~ dunif(0,1)
  
  # Specify different priors for each psens
  psens[1] ~ dunif(0.5, 1)
  psens[2] ~ dunif(0, 1)
  psens[3] ~ dunif(0, 0.5)
  
  for(j in 1:nT) {
    pspec[j] ~ dunif(0, 1)
  }
  
})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega = omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = c(runif(1, 0.5, 1), runif(1, 0, 1), runif(1, 0, 0.5)),
  pspec = runif(ntests, 0, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run <- runMCMC(cMCMC, 
                           niter = 500000, 
                           nburnin = 24000, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, file = "Peril3_FlippingParameters/run1_crap3test.rds")
plot(run$samples)
## 2 test scenario

#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.9, 0.9)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.75, 0.75)
pprev <- 0.3
nbadgers <- 1000
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)  

#simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)

#set up empty array of test outcomes
tests <- array(0, dim = c(nbadgers, length(psens)))
colnames(tests) <- c("test1", "test2")

#for each badger and each test, simulate test outcome
for(i in 1:length(inf)){
  for(j in 1:length(psens)){
    tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}
#how many tests?
###choose number of tests to infer results from
ntests <- 2

#get table of frequencies
binop <- 2^seq(0, ntests-1)
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
    pinf[i] <- pprev * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT])) +
      (1 - pprev) * prod(omega[1:nT, i] * (1 - pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
  }
  
  Te[1:n] ~ dmulti(pinf[1:n], N)
  
  
  pprev ~ dunif(0,1)
  
  for(j in 1:nT) {
    psens[j] ~ dunif(0, 1)
    pspec[j] ~ dunif(0, 1)
  }
  
})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega = omega)

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

system.time(run <- runMCMC(cMCMC, 
                           niter = 500000, 
                           nburnin = 24000, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, file = "Peril3_FlippingParameters/run2_2test.rds")

run1 <- readRDS("Peril3_FlippingParameters/run1_crap3test.rds")
run2 <- readRDS("Peril3_FlippingParameters/run2_2test.rds")

r1_sum <- run1$summary
r1_sum <- as.matrix(r1_sum$all.chains[, 1])
r1_samp <- as.matrix(run1$samples)

r1_samp <- r1_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, psens3 = 4, pspec = 5, pspec2 = 6, pspec3 = 7) %>%
  pivot_longer(cols = c(pprev, psens, pspec, psens3, pspec3),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberofTests = 3) %>%
  select(Inferred_parameter, Estimate, NumberofTests) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75, 0.1, 0.05), times = 7140000/5)) %>%
  mutate(Mean = rep(c(r1_sum[1], r1_sum[2], r1_sum[5], r1_sum[3], r1_sum[6]), times = 7140000/5))

r2_sum <- run2$summary
r2_sum <- as.matrix(r2_sum$all.chains[, 1])
r2_samp <- as.matrix(run2$samples)

r2_samp <- r2_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, pspec = 4, pspec2 = 5) %>%
  select(pprev, psens, pspec) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberofTests = 2) %>%
  select(Inferred_parameter, Estimate, NumberofTests) %>%
  mutate(Truth = rep(c(0.3, 0.9, 0.75), times = 4284000/3)) %>%
  mutate(Mean = rep(c(r2_sum[1], r2_sum[2], r2_sum[4]), times = 4284000/3))


data <- rbind(r1_samp, r2_samp)
data$NumberofTests <- as.factor(data$NumberofTests)
#data <- data %>%
#  mutate(Inferred_parameter = rep(c("Prevalence", "Sensitivity", "Specificity"), times = 4284000/3))

ggplot(data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin(position = "dodge") +
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 2) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(~ NumberofTests) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  labs(x = "Number of Tests", fill = "Inferred parameter", colour = "Inferred parameter") +
  scale_fill_discrete(labels = c("Prevalence", "Sensitivity test A", "Sensitivity test C", "Specificity test A", "Specificity test C")) +
  scale_colour_discrete(labels = c("Prevalence", "Sensitivity test A", "Sensitivity test C", "Specificity test A", "Specificity test C"))

