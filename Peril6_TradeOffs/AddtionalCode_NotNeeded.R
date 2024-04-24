## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)

set.seed(11)
## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
psens <- c(0.5, 0.55, 0.6, 0.65, 0.7)
pspec <- c(0.55, 0.65, 0.7)
pprev <- 0.3
nbadgers <- 1000
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)

## simulate infection status of the badgers
inf <- rbinom(nbadgers, 1, pprev)

## prepare for multiple replications
n_iterations <- 20
# array to store the summaries
summaries <- vector("list", n_iterations)

for(iteration in 1:n_iterations) {
  # Simulate test outcomes
  tests <- array(0, dim = c(nbadgers, length(psens)))
  colnames(tests) <- c("test1", "test2", "test3")
  
  for(i in 1:nbadgers) {
    for(j in 1:length(psens)) {
      tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))
    }
  }
  
  # Frequency table of test outcomes
  ntests <- 3
  binop <- 2^seq(0, ntests-1)
  testbin <- tests[, 1:ntests] %*% binop
  testcounts <- tabulate(testbin + 1, nbins = 2^ntests)
  
  # Omega matrix
  omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
  omega <- as.matrix(omega)
  omega <- t(omega)
  
  # NIMBLE model
  code <- nimbleCode({
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
  
  consts <- list(n = length(testcounts), nT = ntests, N = sum(testcounts), omega = omega)
  data <- list(Te = testcounts)
  inits <- list(pprev = runif(1, 0, 0.7), psens = runif(ntests, 0.5, 1), pspec = runif(ntests, 0.5, 1))
  
  model <- nimbleModel(code, constants = consts, data = data, inits = inits)
  cModel <- compileNimble(model)
  config <- configureMCMC(model)
  rMCMC <- buildMCMC(config)
  cMCMC <- compileNimble(rMCMC, project = model)
  
  run <- runMCMC(cMCMC, 
                 niter = 80000, 
                 nburnin = 12000, 
                 nchains = 3, 
                 progressBar = TRUE,
                 summary = TRUE)
  summaries[[iteration]] <- run$summary$all.chains
}


## extract the means
means_df <- do.call(rbind, lapply(summaries, function(summary){
  data.frame(pprev = summary[1, "Mean"],
             psens1 = summary[2, "Mean"],
             psens2 = summary[3, "Mean"],
             psens3 = summary[4, "Mean"],
             pspec1 = summary[5, "Mean"],
             pspec2 = summary[6, "Mean"],
             pspec3 = summary[7, "Mean"])
}))

ggpairs(means_df)
