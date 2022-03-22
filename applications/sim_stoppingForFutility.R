# ------------------------------------------------------------------------------
# Stopping for futility: Example
# ------------------------------------------------------------------------------

# Function to simulate stopping for futility

sim.stopfutil <- function(n.min, n.max, stepsize, ES.pop, futilitythreshold = 0.01, boundary=c(1/10, 10), alternative="two.sided", prior.mu=0, prior.var=1, iter = 1000, mcmc.setting = list(chains=2, burnin=0, iter=5000)){
  
  # Simulate data for n.max
  g1 <- rnorm(n.max)
  g2 <- rnorm(n.max)-ES.pop
  
  # Conduct intermediate design analyses and save results
  
  stopfutil <- list()
  n.temp <- n.min
  BF.temp <- 1
  futil.temp <- 1
  i <- 1
  
  BF <- switch(alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  while((BF.temp > boundary[1]) && (BF.temp < boundary[2]) && (n.temp <= n.max) && (futil.temp >= futilitythreshold)) {
    
    # Compute current BF
    tval1 <- unname(t.test(g1[1:n.temp], g2[1:n.temp], var.equal = TRUE)$statistic)
    BF.temp <- BF(tval = tval1, n1=n.temp, n2 = n.temp, prior.mu=prior.mu, 
              prior.var=prior.var)
    
    # Do forecast if BF is not outside thresholds already
    forecast <- NULL
    if(((BF.temp > boundary[1]) || (BF.temp < boundary[2])) && (n.temp < n.max)) {
      forecast <- BF.forecast(group1_s1 = g1[1:n.temp], group2_s1 = g2[1:n.temp],
                              group.n_s2 = (n.max-n.temp), 
                              forecastmodel = "combined", 
                              alternative = alternative, 
                              prior.mu = prior.mu, prior.var = prior.var, 
                              iter = iter, mcmc.setting = mcmc.setting)
    }
    
    # Save forecast results to list
    
    stopfutil[[i]] <- list(BF.stage = BF.temp,
                           forecast.stage = forecast)
    
    # Update markers
    n.temp <- n.temp + stepsize
    futil.temp <- min(sum(forecast$BFdist_stage_2 < boundary[1] | forecast$BFdist_stage_2 > boundary[2])/iter, 1)
    # cat(paste0(i, "\n"))
    i <- i+1
  }
  
  return(stopfutil)
}

# Conduct analysis many times for a design analysis of the futility stopping design

SIM_stopfutil_05 <- vector(mode = "list", length = 1000)
for(j in 1:1000){
  cat(paste0(j, "\n"))
  SIM_stopfutil_05[[j]] <- sim.stopfutil(n.min = 10, n.max = 50, stepsize = 5, ES.pop = 0.5)
}

save(SIM_stopfutil_05, file = "SIM_stopfutil_05.RData")

# Evaluation of simulation results

finalBFs <- rep(NA, 1000)
finalNs <- rep(NA, 1000)
reasonEnd <- rep(NA, 1000)
for(i in 1:1000){
  finalBFs[i] <- SIM_stopfutil_05[[i]][[length(SIM_stopfutil_05[[i]])]]$BF.stage
  finalNs[i] <- seq(10, 50, by = 5)[length(SIM_stopfutil_05[[i]])]
  if(finalBFs[i] < 1/10){
    reasonEnd[i] <- "H0"
  } else if (finalBFs[i] > 10) {
    reasonEnd[i] <- "H1"
  } else if (finalNs[i] < 50){
    reasonEnd[i] <- "futility"
  } else if (finalNs[i] == 50){
    reasonEnd[i] <- "maxN"
  }
}

table(cut(finalBFs, breaks = c(0, 1/10, 10, Inf))) # How many end at which threshold
barplot(table(as.factor(finalNs))/1000) # Distribution of sample sizes
pie(table(reasonEnd)) # What are the reasons for stopping
mean(finalNs)

simDat <- data.frame(finalBFs = finalBFs,
                     finalNs = as.factor(finalNs),
                     reasonEnd = as.factor(reasonEnd),
                     finalBFCat = cut(finalBFs, breaks = c(0, 1/10, 10, Inf)))

ggplot(simDat, aes(x = finalNs, fill = reasonEnd)) +
  geom_bar()


BFDAmaxN <- BFDA::BFDA.sim(expected.ES = 0.5, type = "t.between", prior = list("normal", list(prior.mean = 0, prior.variance = 1)), n.min = 10, n.max = 50, stepsize = 5, design = "sequential", boundary = c(1/10, 10), alternative = "two.sided")

BFDA.analyze(BFDAmaxN, boundary=c(1/10, 10))

layerplot_BFforecast(SIM_stopfutil_05[[2]][[4]]$forecast.stage, thresholds = c(1/10,10))
