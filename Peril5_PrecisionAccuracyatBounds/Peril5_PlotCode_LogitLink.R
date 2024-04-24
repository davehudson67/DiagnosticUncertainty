## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)

## set up initial parameters, number of tests and iterations
ntests <- 5
psensvec <- c(0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)
n_iterations <- 20

## initialize a 3D array to store the summary statistics
summaryinference_all <- array(0, dim = c(n_iterations, length(psensvec), 4 * ntests + 2))
summaryinfsd_all <- array(0, dim = c(n_iterations, length(psensvec), 4 * ntests + 2))


## Begin loop over iterations
for(iter in 1:n_iterations) {
  
  # Initialize summaryinference and summaryinfsd inside the loop over iterations
  summaryinference <- matrix(0, nrow = length(psensvec), ncol = 4 * ntests + 2)
  summaryinfsd <- matrix(0, nrow = length(psensvec), ncol = 4 * ntests + 2)
  
  ## Begin loop over sensitivities
  for(q in 1:length(psensvec)){
  
    psens <- c(psensvec[q], 0.8, 0.7, 0.65, 0.6)
    pspec <- c(0.75, 0.7, 0.65, 0.6, 0.55)
    pprev <- 0.3
    nbadgers <- 1000
    print(psens)
    print(pspec)
    print(pprev)
    
    #simulate infection status of the badgers
    inf <- rbinom(nbadgers, 1, pprev)
    
    #set up empty array of test outcomes
    tests <- array(0, dim = c(nbadgers, length(psens)))
    colnames(tests) <- c("test1", "test2", "test3", "test4", "test5")
    
    #for each badger and each test, simulate test outcome
    for(i in 1:length(inf)){
      for(j in 1:length(psens)){
        tests[i, j] <- rbinom(1, 1, ifelse(inf[i] == 1, psens[j], 1 - pspec[j]))}}
    
    ############manipulate test numbers from here###############
    
    #get table of frequencies
    binop <- 2^seq(0, ntests - 1)
    testbin <- tests[, 1:ntests] %*% binop
    testcounts <- tabulate(testbin + 1, nbins = 2^ntests)
    testcounts
    
    #create omega matrix of binary sequence
    omega <- expand.grid(replicate(ntests, c(0, 1), simplify = FALSE))
    omega <- as.matrix(omega)
    omega <- t(omega)
    
    code  <- nimbleCode({
      
      for(i in 1:n) {
        
        pinf[i] <- pprev * prod(omega[1:nT, i] * psens[1:nT] + (1 - omega[1:nT, i]) * (1 - psens[1:nT])) +
          (1 - pprev) * prod(omega[1:nT, i] * (1-pspec[1:nT]) + (1 - omega[1:nT, i]) * (pspec[1:nT]))
      }
      
      for(k in 1:nT){
        psens[k] <- 1 / (1 + 1 / exp(psenslogit[k]))
        pspec[k] <- 1 / (1 + 1 / exp(pspeclogit[k]))
      }
      
      pprev <- 1 / (1 + 1/exp(pprevlogit))
      
      Te[1:n] ~ dmulti(pinf[1:n], N)
      
      
      pprevlogit ~ dnorm(0, sd = 100)
      
      for(j in 1:nT) {
        psenslogit[j] ~ dexp(1)
        pspeclogit[j] ~ dexp(1)
      }
      
      
    })
    
    consts <- list(n = length(testcounts),
                   nT = ntests, N = sum(testcounts), omega=omega)
    
    data <- list(
      Te = testcounts)
    
    inits <- list(
      pprevlogit = runif(1, 0, 1),
      psenslogit = rexp(ntests, 1),
      pspeclogit = rexp(ntests, 1)
    )
    
    model <- nimbleModel(code, constants = consts, data = data, inits = inits)
    
    cModel <- compileNimble(model)
    
    config <- configureMCMC(model)
    
    config$addMonitors(c("psens", "pspec", "pprev"))
    
    rMCMC <- buildMCMC(config)
    
    cMCMC <- compileNimble(rMCMC, project = model)
    
    system.time(runlogit <- runMCMC(cMCMC, 
                                    niter = 50000, 
                                    nburnin = 2400, 
                                    nchains = 3, 
                                    progressBar = TRUE, 
                                    summary = TRUE, 
                                    samplesAsCodaMCMC = TRUE, 
                                    thin = 1))
    
    runlogitsummary <- MCMCsummary(runlogit$samples)

    summaryinference[q, ] <- runlogitsummary$mean
    summaryinfsd[q, ] <- runlogitsummary$sd
  
  } # end loop over sensitivities
  
  summaryinference_all[iter, , ] <- summaryinference
  summaryinfsd_all[iter, , ] <- summaryinfsd  

} # end loop over iterations

## Save data
saveRDS(summaryinference_all, "Peril5_PrecisionAccuracyatBounds/SummaryInference_ALL.rds")
saveRDS(summaryinfsd_all, "Peril5_PrecisionAccuracyatBounds/SummaryInferenceSD_ALL.rds")


# Create an empty data frame to hold the results
plot_data <- data.frame(iteration = integer(),
                        true_sens = double(),
                        inferred_sens = double(),
                        posterior_sd = double())

# Loop through the 3D arrays to populate the data frame
for (iter in 1:dim(summaryinference_all)[1]) {
  for (q in 1:dim(summaryinference_all)[2]) {
    true_sens = psensvec[q]
    
    # Replace x, y with the indices where inferred_sens and posterior_sd are stored
    inferred_sens = summaryinference_all[iter, q, 3]
    posterior_sd = summaryinfsd_all[iter, q, 3]
    
    plot_data <- rbind(plot_data, data.frame(iteration = iter, true_sens = true_sens, 
                                             inferred_sens = inferred_sens, posterior_sd = posterior_sd))
  }
}


# Plotting inferred vs true sensitivity
ggplot(plot_data, aes(x = true_sens, y = inferred_sens)) +
  geom_point(alpha = 0.5) +
  xlab("True Sensitivity") +
  ylab("Inferred Sensitivity") +
  ggtitle("Inferred vs True Sensitivity Across Iterations") +
  geom_abline(slope = 1,
              intercept = 0,
              color="blue")


# Plotting posterior sd vs true sensitivity
ggplot(plot_data, aes(x = true_sens, y = posterior_sd)) +
  geom_point(alpha = 0.5) +
  xlab("True Sensitivity") +
  ylab("Posterior Standard Deviation") +
  ggtitle("Posterior SD vs True Sensitivity Across Iterations")
  
