# ==============================================================================
# Stopping for futility: Function to create example and design analysis
# ==============================================================================

source("ttest_2sample_normalprior/forecast_2sample_t.R")

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
    if(((BF.temp > boundary[1]) && (BF.temp < boundary[2])) && (n.temp < n.max)) {
      forecast <- BF.forecast(group1_s1 = g1[1:n.temp], group2_s1 = g2[1:n.temp],
                              group.n_s2 = (n.max-n.temp), 
                              forecastmodel = "combined", 
                              alternative = alternative, 
                              prior.mu = prior.mu, prior.var = prior.var, 
                              iter = iter, mcmc.setting = mcmc.setting)
    }
    
    # Save forecast results to list
    
    stopfutil[[i]] <- list(BF.temp = BF.temp,
                           BF.forecast = forecast$BFdist_stage_2,
                           n.temp = n.temp,
                           n.forecast = forecast$settings$group.n_s2,
                           futil.temp = min(sum(forecast$BFdist_stage_2 < boundary[1] | forecast$BFdist_stage_2 > boundary[2])/iter, 1))
    
    # Update markers
    n.temp <- n.temp + stepsize
    futil.temp <- min(sum(forecast$BFdist_stage_2 < boundary[1] | forecast$BFdist_stage_2 > boundary[2])/iter, 1)
    # cat(paste0(i, "\n"))
    i <- i+1
  }
  
  return(stopfutil)
}

designAnalysis.stopfutil <- function(n.min, n.max, stepsize, ES.pop, futilitythreshold = 0.01, boundary=c(1/10, 10), alternative="two.sided", prior.mu=0, prior.var=1, iter.forecast = 1000, iter.DA = 1000, mcmc.setting = list(chains=2, burnin=0, iter=5000)){
  
  res <- vector(mode = "list", length = iter.DA)
  
  for(j in 1:iter.DA){
    cat(paste0(j, "\n"))
    res[[j]] <- sim.stopfutil(n.min = n.min, n.max = n.max, 
                              stepsize = stepsize, ES.pop = ES.pop,
                              futilitythreshold = futilitythreshold,
                              boundary = boundary, alternative = alternative,
                              prior.mu=prior.mu, prior.var=prior.var, 
                              iter = iter.forecast, mcmc.setting = mcmc.setting)
  }
  
  return(list(res.DA = res,
              settings = list(n.min = n.min,
                              n.max = n.max,
                              stepsize = stepsize,
                              ES.pop = ES.pop,
                              futilitythreshold = futilitythreshold,
                              boundary = boundary,
                              alternative = alternative,
                              prior.mu = prior.mu,
                              prior.var = prior.var)))
  
  
}