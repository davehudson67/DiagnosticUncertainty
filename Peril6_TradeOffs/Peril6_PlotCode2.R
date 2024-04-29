## load libraries
library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)
library(GGally)
library(patchwork)

## set up initial parameters, number of tests and iterations
ntests <- 5
psensvec <- c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
pspecvec <- c(0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55)
n_iterations <- 20

## initialize a 3D array to store the summary statistics
summaryinference_all <- array(0, dim = c(n_iterations, length(psensvec), 2 * ntests + 1))
summaryinfsd_all <- array(0, dim = c(n_iterations, length(psensvec), 2 * ntests + 1))


## Begin loop over iterations
for(iter in 1:n_iterations) {
  
  # Initialize summaryinference and summaryinfsd inside the loop over iterations
  summaryinference <- matrix(0, nrow = length(psensvec), ncol = 2 * ntests + 1)
  summaryinfsd <- matrix(0, nrow = length(psensvec), ncol = 2 * ntests + 1)
  
  ## Begin loop over sensitivities
  for(q in 1:length(psensvec)){
    
    psens <- c(psensvec[q], 0.8, 0.7, 0.65, 0.6)
    pspec <- c(pspecvec[q], 0.75, 0.7, 0.65, 0.6)
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
        pspec[k] ~ dunif(0, 1)
      }
      
      pprev ~ dunif(0, 0.5)
      
      Te[1:n] ~ dmulti(pinf[1:n], N)
      
    })
    
    consts <- list(n = length(testcounts),
                   nT = ntests, N = sum(testcounts), omega = omega)
    
    data <- list(
      Te = testcounts)
    
    inits <- list(
      pprev = runif(1, 0, 0.5),
      psens = runif(ntests, 0.5, 1),
      pspec = runif(ntests, 0.5, 1)
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
    
    runsummary <- MCMCsummary(run$samples)
    
    summaryinference[q, ] <- runsummary$mean
    summaryinfsd[q, ] <- runsummary$sd
    
  } # end loop over sensitivities
  
  summaryinference_all[iter, , ] <- summaryinference
  summaryinfsd_all[iter, , ] <- summaryinfsd  
  
} # end loop over iterations

## Save data
saveRDS(summaryinference_all, "Peril6_TradeOffs/SummaryInference_ALL2.rds")
saveRDS(summaryinfsd_all, "Peril6_TradeOffs/SummaryInfSD_ALL2.rds")

## Save last run
saveRDS(run, "Peril6_TradeOffs/FinalRunSamples2.rds")

summaryinference_all <- readRDS("Peril6_TradeOffs/SummaryInference_ALL2.rds")

# Create an empty data frame to hold the results
plot_data <- data.frame(iteration = integer(),
                        true_sens = double(),
                        true_spec = double(),
                        true_prev = double(),
                        inferred_sens = double(),
                        inferred_spec = double(),
                        inferred_prev = double())

# Loop through the 3D arrays to populate the data frame
for (iter in 1:dim(summaryinference_all)[1]) {
  for (q in 1:dim(summaryinference_all)[2]) {
    true_sens = psensvec[q]
    true_spec = pspecvec[q]
    true_prev = 0.3
    
    # Replace x, y with the indices where inferred_sens and posterior_sd are stored
    inferred_sens = summaryinference_all[iter, q, 2]
    inferred_spec = summaryinference_all[iter, q, 7]
    inferred_prev = summaryinference_all[iter, q, 1]

    plot_data <- rbind(plot_data, data.frame(iteration = iter, true_sens = true_sens,  true_spec = true_spec,
                                             true_prev = true_prev,
                                             inferred_sens = inferred_sens, inferred_prev = inferred_prev,
                                             inferred_spec = inferred_spec))
  }
}

# Create an empty data frame to hold the results
plot_dataSD <- data.frame(iteration = integer(),
                        true_sens = double(),
                        true_spec = double(),
                        inferred_sensSD = double(),
                        inferred_specSD = double(),
                        inferred_prevSD = double())

# Loop through the 3D arrays to populate the data frame
for (iter in 1:dim(summaryinfsd_all)[1]) {
  for (q in 1:dim(summaryinfsd_all)[2]) {
    true_sens = psensvec[q]
    true_spec = pspecvec[q]
    
    # inferred_sens and posterior_sd are stored
    inferred_sensSD = summaryinfsd_all[iter, q, 2]
    inferred_specSD = summaryinfsd_all[iter, q, 7]
    inferred_prevSD = summaryinfsd_all[iter, q, 1]
    
    plot_dataSD <- rbind(plot_dataSD, data.frame(iteration = iter, true_sens = true_sens,  true_spec = true_spec,
                                             inferred_sensSD = inferred_sensSD, inferred_prevSD = inferred_prevSD,
                                             inferred_specSD = inferred_specSD))
  }
}

# Plotting inferred SD
aSD <- plot_dataSD %>%
  group_by(true_sens) %>%
  mutate(Asd = inferred_sensSD) %>%
  ungroup() %>%
  ggplot(aes(x = true_sens, y = Asd)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample precision") +
  #ggtitle("Inferred vs True Sensitivity Across Iterations") +
  #geom_abline(slope = 1,
  #            intercept = 0,
  #            color="blue") +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10)) +
  ylim(0, 0.1)
aSD
# Plotting inferred specificity vs true sensitivity
bSD <- plot_dataSD %>%
  group_by(true_sens) %>%
  mutate(Asd = inferred_specSD) %>%
  ungroup() %>%
  ggplot(aes(x = true_sens, y = Asd)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample precision") +
  #ggtitle("Inferred Specificity vs True Sensitivity Across Iterations") +
  #xlim(0, 5) +
  ylim(0, 0.1) +
#  geom_abline(slope = -1,
#              intercept = 1.5,
#              color="blue")  +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10))
  
# Plotting inferred specificity vs true sensitivity
cSD <- plot_dataSD %>%
  group_by(true_sens) %>%
  mutate(Asd = inferred_prevSD) %>%
  ungroup() %>%
  ggplot(aes(x = true_sens, y = Asd)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample precision") +
  #ggtitle("Inferred vs True Sensitivity Across Iterations") +
#  geom_hline(yintercept = 0.3,
#             color="blue") +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10)) +
  ylim(0, 0.1)

#aSD | bSD | cSD

# Plotting inferred vs true sensitivity
a <- ggplot(plot_data, aes(x = true_sens, y = inferred_sens)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Sensitivity") +
  #ggtitle("Inferred vs True Sensitivity Across Iterations") +
  geom_abline(slope = 1,
              intercept = 0,
              color="blue") +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  ylim(0.5, 1)

# Plotting inferred specificity vs true sensitivity
b <- ggplot(plot_data, aes(x = true_sens, y = inferred_spec)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Specificity") +
  #ggtitle("Inferred Specificity vs True Sensitivity Across Iterations") +
  #xlim(0, 5) +
  ylim(0.5, 1) +
  geom_abline(slope = -1,
              intercept = 1.5,
              color="blue")  +
  theme_bw() +
  theme(text = element_text(size = 13))

# Plotting inferred specificity vs true sensitivity
c <- ggplot(plot_data, aes(x = true_sens, y = inferred_prev)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Prevalence") +
  #ggtitle("Inferred vs True Sensitivity Across Iterations") +
  geom_hline(yintercept = 0.3,
              color="blue") +
  theme_bw() +
  theme(text = element_text(size = 13)) +
  ylim(0, 0.5)

#a/b/c

# Plotting mean inferred - truth
plot_data1 <- plot_data %>%
  group_by(true_sens) %>%
  mutate(DiffSens = inferred_sens - true_sens) %>%
  mutate(DiffSpec = inferred_spec - true_spec) %>%
  mutate(DiffPrev = inferred_prev - 0.3)

a1 <- ggplot(plot_data1, aes(x = true_sens, y = DiffSens)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample accuracy") +
  ylim(-0.5, 0.2) +
  geom_abline(slope = 0,
              intercept = 0,
              color="black") +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10))

b1 <- ggplot(plot_data1, aes(x = true_sens, y = DiffSpec)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample accuracy") +
  geom_abline(slope = 0,
              intercept = 0,
              color="black") +
  ylim(-0.5, 0.2) +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10))
  

c1 <- ggplot(plot_data1, aes(x = true_sens, y = DiffPrev)) +
  geom_point(alpha = 0.2) +
  xlab("") +
  ylab("Within sample accuracy") +
  geom_abline(slope = 0,
              intercept = 0,
              color="black") +
  ylim(-0.5, 0.2) +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.title.y = element_text(size = 10))
  
Fig7 <- (a/b/c)|(a1/b1/c1)|(aSD/bSD/cSD)

wrap_elements(panel = Fig7) +
  labs(tag = "True Sensitivity") +
  theme(
    plot.tag = element_text(size = 13),
    plot.tag.position = "bottom"
  )


## Pairs plot
run <- readRDS("Peril6_TradeOffs/FinalRunSamples.rds")

rundf <- as.matrix(run$samples, chain = TRUE) %>%
  as.data.frame()# %>%
  #select(c(1:5, 8:10)) %>%
  #filter(CHAIN == 1)

rundf$CHAIN <- as.factor(rundf$CHAIN)

# thin samples
rundf_filtered <- rundf[seq(1, nrow(rundf), 200), ]

# fig 8
ggpairs(rundf_filtered, columns = 2:8, aes(alpha = 0.2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15), axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank())

# fig 9
ggpairs(rundf_filtered, columns = 2:8, aes(colour = CHAIN, alpha = 0.1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 15), axis.text.x = element_blank(), axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank())
