## detection of disease in the absence of a gold standard

## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(patchwork)

## set up initial parameters, number of tests and iterations
ntests <- 5
psensvec <- c(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99)
n_iterations <- 20

## initialize a 3D array to store the summary statistics
summaryinference_all <- array(0, dim = c(n_iterations, length(psensvec), 2 * ntests + 1))
summaryinfsd_all <- array(0, dim = c(n_iterations, length(psensvec), 2 * ntests + 1))
summaryinfmode_all <- array(0, dim = c(n_iterations, length(psensvec), 2 * ntests + 1))

## function to calculate mode
modef <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
  
## Begin loop over iterations
for(iter in 1:n_iterations) {
  
  # Initialize summaryinference and summaryinfsd inside the loop over iterations
  summaryinference <- matrix(0, nrow = length(psensvec), ncol = 2 * ntests + 1)
  summaryinfsd <- matrix(0, nrow = length(psensvec), ncol = 2 * ntests + 1)
  summaryinfmode <- matrix(0, nrow = length(psensvec), ncol = 2 * ntests + 1)
  
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
        psens[k] ~ dunif(0, 1)
        pspec[k] ~ dunif(0.5, 1)
      }
      
      pprev ~ dunif(0, 0.5)
      
      Te[1:n] ~ dmulti(pinf[1:n], N)
      
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
                                    niter = 50000, 
                                    nburnin = 2400, 
                                    nchains = 3, 
                                    progressBar = TRUE, 
                                    summary = TRUE, 
                                    samplesAsCodaMCMC = TRUE, 
                                    thin = 1))
    
    runsummary <- MCMCsummary(run$samples)
    runM <- as.matrix(run$samples)
    
    summaryinference[q, ] <- runsummary$mean
    summaryinfsd[q, ] <- runsummary$sd
    summaryinfmode[q, ] <- modef(runM[, 2])
    
  } # end loop over sensitivities
  
  summaryinference_all[iter, , ] <- summaryinference
  summaryinfsd_all[iter, , ] <- summaryinfsd
  summaryinfmode_all[iter, , ] <- summaryinfmode
  
} # end loop over iterations

## Save data
saveRDS(summaryinference_all, "Peril5_PrecisionAccuracyatBounds/SummaryInference_ALL_IDLink.rds")
saveRDS(summaryinfsd_all, "Peril5_PrecisionAccuracyatBounds/SummaryInferenceSD_ALL_IDLink.rds")
saveRDS(summaryinfmode_all, "Peril5_PrecisionAccuracyatBounds/SummaryInferenceMode_ALL_IDLink.rds")

summaryinference_all <- readRDS("Peril5_PrecisionAccuracyatBounds/SummaryInference_ALL_IDLink.rds")
summaryinfsd_all <- readRDS("Peril5_PrecisionAccuracyatBounds/SummaryInferenceSD_ALL_IDLink.rds")
summaryinfmode_all <- readRDS("Peril5_PrecisionAccuracyatBounds/SummaryInferenceMode_ALL_IDLink.rds")

# Create an empty data frame to hold the results
plot_data <- data.frame(iteration = integer(),
                        true_sens = double(),
                        inferred_sens = double(),
                        posterior_sd = double(),
                        posterior_mode = double())

# Loop through the 3D arrays to populate the data frame
for (iter in 1:dim(summaryinference_all)[1]) {
  for (q in 1:dim(summaryinference_all)[2]) {
    true_sens = psensvec[q]
    
    inferred_sens = summaryinference_all[iter, q, 2]
    posterior_sd = summaryinfsd_all[iter, q, 2]
    posterior_mode = summaryinfmode_all[iter, q, 2]
    
    plot_data <- rbind(plot_data, data.frame(iteration = iter, true_sens = true_sens, 
                                             inferred_sens = inferred_sens, posterior_sd = posterior_sd,
                                             posterior_mode = posterior_mode))
  }
}

#A Plotting inferred vs true sensitivity
means <- ggplot(plot_data, aes(x = true_sens, y = inferred_sens)) +
  geom_point(alpha = 0.2, size = 5) +
  xlab("") +
  ylab("Inferred Sensitivity") +
  #ggtitle("Inferred vs True Sensitivity Across Iterations") +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  geom_abline(slope = 1,
              intercept = 0,
              color="blue") +
  labs(title = "a")


#B Plotting posterior sd vs true sensitivity
sd <- ggplot(plot_data, aes(x = true_sens, y = posterior_sd)) +
  geom_point(alpha = 0.2, size = 5) +
  xlab("") +
  ylab("Within sample precision") +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  labs(title = "b")
  #ggtitle("Posterior SD vs True Sensitivity Across Iterations")

#C
sdMode <- plot_data %>%
  group_by(true_sens) %>%
  mutate(sdMode = sd(inferred_sens)) %>%
  ungroup() %>%
  ggplot(aes(x = true_sens, y = sdMode)) +
  geom_point() +
  labs(x = "",
       y = "Among sample precision")  +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  labs(title = "c")

# D
mean_diffs <- plot_data %>%
  mutate(diff = inferred_sens - true_sens) %>%  # Use absolute value
  group_by(true_sens) %>%
  summarise(mean_diff = mean(diff))

# Create the plot of mean differences (among-sample accuracy)
Amongmeans <- ggplot(mean_diffs, aes(x = true_sens, y = mean_diff)) +
  geom_point() +
  labs(x = "",
       y = "Among sample accuracy")  +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  ylim(-0.1, 0.1) +
  geom_abline(slope = 0, intercept = 0, colour = "black", linetype = "dashed") +
  labs (title = "d")


Fig6 <- (means|sd) /(sdMode|Amongmeans)

wrap_elements(panel = Fig6) +
  labs(tag = "True Sensitivity") +
  theme(
    plot.tag = element_text(size = 15),
    plot.tag.position = "bottom"
  )
