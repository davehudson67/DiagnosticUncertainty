## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)
library(GGally)
library(patchwork)

set.seed(11)
## set test sensitivities and specificities
## for three different tests and prevalences
## for different populations
psens <- c(0.85, 0.8, 0.7)
pspec <- c(0.70, 0.85, 0.75)
pprev <- 0.3
nbadgers <- 1000
true_values <- c(pprev, psens[1], pspec[1])
true_parameter <- c("Prevalence", "Sensitivity", "Specificity")
truth <- data.frame(true_parameter, true_values)

n_iterations <- 30
ntests <- 3

# Initialize summaryinference and summaryinfsd inside the loop over iterations
summaryinfMean <- matrix(0, nrow = n_iterations, ncol = 2 * ntests + 1)
summaryinfSD <- matrix(0, nrow = n_iterations, ncol = 2 * ntests + 1)

## Begin loop over iterations
for(iter in 1:n_iterations) {
  
  #simulate infection status of the badgers
  inf <- rbinom(nbadgers, 1, pprev)

  #set up empty array of test outcomes
  tests <- array(0, dim = c(nbadgers, length(psens)))
  colnames(tests) <- c("test1", "test2", "test3")
  
  #for each badger and each test, simulate test outcome
  for(i in 1:length(inf)){
    for(j in 1:length(psens)){
      tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1-pspec[j]))}}
  
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
    
    
    pprev ~ dunif(0, 0.5)
    
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
    pprev = runif(1, 0, 0.7),
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
  
  runsummary <- MCMCsummary(run$samples)
  
  summaryinfMean[iter, ] <- runsummary$mean
  summaryinfSD[iter, ] <- runsummary$sd
  
  } # end loop over sensitivities

saveRDS(summaryinfMean, "Peril6_TradeOffs/SummaryInfMean2.rds")
saveRDS(summaryinfSD, "Peril6_TradeOffs/SummaryInfSD2.rds")



summaryinfMean <- as.data.frame(readRDS("Peril6_TradeOffs/SummaryInfMean2.rds"))
colnames(summaryinfMean) <- c("Prev", "Se1", "Se2", "Se3", "Sp1", "Sp2", "Sp3")

# fig 8
ggpairs(summaryinfMean, aes(alpha = 0.2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
  #theme(axis.ticks = element_blank())

dev.off()
