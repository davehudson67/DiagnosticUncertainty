## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)

set.seed(11)
## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
#psens <- round(runif(6, 0.2, 0.9), 2)
psens <- c(0.85, 0.8, 0.7, 0.6, 0.5)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.70, 0.85, 0.75, 0.65, 0.55)
pprev <- 0.3
nbadgers <- 1600
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)

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
ntests <- 1

#get table of frequencies
binop <- 2^seq(0, ntests-1)
testbin <- tests[, 1:ntests] #%*% binop
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
    pinf[i] <- pprev * (omega[1, i] * psens + (1 - omega[1, i]) * (1 - psens)) +
      (1 - pprev) * (omega[1, i] * (1 - pspec) + (1 - omega[1, i]) * (pspec))
  }
  
  Te[1:n] ~ dmulti(pinf[1:n], N)
  
  
  pprev ~ dunif(0, 1)
  psens ~ dunif(0.5, 1)
  pspec ~ dunif(0.5, 1)
  
})

consts <- list(n = length(testcounts),
               N = sum(testcounts), omega = omega)

data <- list(
  Te = testcounts)

inits <- list(
  pprev = runif(1, 0, 1),
  psens = runif(1, 0.5, 1),
  pspec = runif(1, 0.5, 1)
)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)

config <- configureMCMC(model)

rMCMC <- buildMCMC(config)

cMCMC <- compileNimble(rMCMC, project = model)

system.time(run <- runMCMC(cMCMC, 
                           niter = 80000, 
                           nburnin = 12000, 
                           nchains = 3, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, file = "Peril1_InsufficientTests/run1.rds")

### Now for multiple tests ###

#how many tests?
###choose number of tests to infer results from
ntests <- 5

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
  
  
  pprev ~ dunif(0, 0.7)
  
  for(j in 1:nT) {
    psens[j] ~ dunif(0.5, 1)
    pspec[j] ~ dunif(0.5, 1)
  }
  
})

consts <- list(n = length(testcounts),
               nT = ntests, N = sum(testcounts), omega = omega)

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

system.time(run <- runMCMC(cMCMC, 
                              niter = 80000, 
                              nburnin = 12000, 
                              nchains = 3, 
                              progressBar = TRUE, 
                              summary = TRUE, 
                              samplesAsCodaMCMC = TRUE, 
                              thin = 1))

saveRDS(run, file = "Peril1_InsufficientTests/run5.rds")

run1 <- readRDS("Peril1_InsufficientTests/run1.rds")
run2 <- readRDS("Peril1_InsufficientTests/run2.rds")
run3 <- readRDS("Peril1_InsufficientTests/run3.rds")
run4 <- readRDS("Peril1_InsufficientTests/run4.rds")
run5 <- readRDS("Peril1_InsufficientTests/run5.rds")

r1_sum <- run1$summary
r1_sum <- as.matrix(r1_sum$all.chains[, 1])
r1_samp <- as.matrix(run1$samples)

r1_samp <- r1_samp %>%
  as.data.frame() %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 1) %>%
  mutate(Truth = rep(c(0.3, 0.85, 0.7), times = 428400/3)) %>%
  mutate(Mean = rep(c(r1_sum[1], r1_sum[2], r1_sum[3]), times = 428400/3))

r2_sum <- run2$summary
r2_sum <- as.matrix(r2_sum$all.chains[, 1])
r2_samp <- as.matrix(run2$samples)

r2_samp <- r2_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, pspec = 4, pspec2 = 5) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 2) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.3, 0.85, 0.7), times = 612000/3)) %>%
  mutate(Mean = rep(c(r2_sum[1], r2_sum[2], r2_sum[3]), times = 612000/3))

r3_sum <- run3$summary
r3_sum <- as.matrix(r3_sum$all.chains[, 1])
r3_samp <- as.matrix(run3$samples)

r3_samp <- r3_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, pspec = 4, pspec2 = 5) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 3) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.3, 0.85, 0.7), times = 612000/3)) %>%
  mutate(Mean = rep(c(r3_sum[1], r3_sum[2], r3_sum[3]), times = 612000/3))

r4_sum <- run4$summary
r4_sum <- as.matrix(r4_sum$all.chains[, 1])
r4_samp <- as.matrix(run4$samples)

r4_samp <- r4_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, pspec = 4, pspec2 = 5) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 4) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.3, 0.85, 0.7), times = 612000/3)) %>%
  mutate(Mean = rep(c(r4_sum[1], r4_sum[2], r4_sum[3]), times = 612000/3))

r5_sum <- run5$summary
r5_sum <- as.matrix(r5_sum$all.chains[, 1])
r5_samp <- as.matrix(run5$samples)

r5_samp <- r5_samp %>%
  as.data.frame() %>%
  rename(pprev = 1, psens = 2, psens2 = 3, pspec = 4, pspec2 = 5) %>%
  pivot_longer(cols = c(pprev, psens, pspec),  values_to = "Estimate", names_to = "Inferred_parameter" ) %>%
  mutate(NumberOfTests = 5) %>%
  select(Inferred_parameter, Estimate, NumberOfTests) %>%
  mutate(Truth = rep(c(0.3, 0.85, 0.7), times = 612000/3)) %>%
  mutate(Mean = rep(c(r5_sum[1], r5_sum[2], r5_sum[3]), times = 612000/3))

data <- rbind(r1_samp, r2_samp, r3_samp, r4_samp, r5_samp)
data$NumberOfTests <- as.factor(data$NumberOfTests)
data <- data %>%
  mutate(Inferred_parameter = rep(c("Prevalence", "Sensitivity", "Specificity"), times = 2876400/3))

saveRDS(data, "Peril1_InsufficientTests/PlotData.rds")
data <- readRDS("Peril1_InsufficientTests/PlotData.rds")

ggplot(data, aes(x = Inferred_parameter, y = Estimate, fill = Inferred_parameter)) +
  geom_violin(position = "dodge") +
  geom_hline(aes(yintercept = Truth, colour = Inferred_parameter), linetype = 2) +
  stat_summary(fun = "mean", geom = "point") +
  facet_grid(~ NumberOfTests) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  labs(x = "Number of Tests", fill = "Inferred parameter", colour = "Inferred parameter")

