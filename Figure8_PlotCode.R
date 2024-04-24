## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)
library(GGally)

## set up initial parameters, number of tests and iterations
n_simulations <- 50

psens <- c(0.75, 0.55, 0.7, 0.5, 0.5)
#pspec <- round(runif(6, 0.6, 0.9), 2)
pspec <- c(0.55, 0.85, 0.75, 0.65, 0.55)
pprev <- 0.3
nbadgers <- 500
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

################################################################################

summaries <- list()

for (k in 1:n_simulations){
  
  code  <- nimbleCode({
      
      for(i in 1:n) {
        
        pinf[i] <- pprev * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT])) +
          (1 - pprev) * prod(omega[1:nT, i] * (1-pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
      }
      
      for(k in 1:nT){
        psens[k] ~ dunif(0.5, 1)
        pspec[k] ~ dunif(0.5, 1)
      }
      
      pprev ~ dunif(0, 1)
      
      Te[1:n] ~ dmulti(pinf[1:n], N)
      
    })
    
    consts <- list(n = length(testcounts),
                   nT = ntests, N = sum(testcounts), omega = omega)
    
    data <- list(
      Te = testcounts)
    
    inits <- list(
      pprev = runif(1, 0.5, 1),
      psens = runif(ntests, 0.5, 1),
      pspec = runif(ntests, 0, 1)
    )
    
    model <- nimbleModel(code, constants = consts, data = data, inits = inits)
    
    cModel <- compileNimble(model)
    
    config <- configureMCMC(model)
    
    rMCMC <- buildMCMC(config)
    
    cMCMC <- compileNimble(rMCMC, project = model)
    
    system.time(run <- runMCMC(cMCMC, 
                               niter = 10000, 
                               nburnin = 2400, 
                               nchains = 3, 
                               progressBar = TRUE, 
                               summary = TRUE, 
                               samplesAsCodaMCMC = TRUE, 
                               thin = 1))

    # Extract and store the summary for this repetition
    runsummary <- MCMCsummary(run$samples)
    summaries[[k]] <- runsummary
}

# Empty list to store parameter estimates from all simulations
all_estimates <- list()

# Combine estimates from all simulations into a single data frame
all_estimates_df <- do.call(rbind, summaries) %>%
  mutate(Sim = rep(1:n_simulations, each = 7)) %>%
  mutate(Parameter = rep(c("Prev", "Sens[1]", "Sens[2]", "Sens[3]", "Spec[1]", "Spec[2]", "Spec[3]"), times = n_simulations)) %>%
  select(mean, Sim, Parameter) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

ggpairs(all_estimates_df[-1]) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


