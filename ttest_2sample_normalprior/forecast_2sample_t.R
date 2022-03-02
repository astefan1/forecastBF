# Continuous BFDA for 2-sample t-test with normal prior: Combined BF distribution

# source("ttest_2sample_normalprior/posterior_2sample_t.R")
source("ttest_2sample_normalprior/bf_2sample_t.R")
source("ttest_2sample_normalprior/mcmc_2sample_t.R")

BF.forecast<- function(group1_s1, group2_s1, group.n_s2, forecastmodel = "combined", alternative = "two.sided", prior.mu = 0, prior.var = 1, iter = 1000, mcmc.setting = list()){
  
  # Which BF to compute?
  BF <- switch(alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  # MCMC settings for new data
  if(is.null(mcmc.setting$chains)) mcmc.setting$chains <- 10
  if(is.null(mcmc.setting$iter)) mcmc.setting$iter <- 5000
  if(is.null(mcmc.setting$burnin)) mcmc.setting$burnin <- 1000
  
  # t value stage 1
  tval1 <- unname(t.test(group1_s1, group2_s1, var.equal = TRUE)$statistic)
  
  # Bayes factor stage 1
  n1_s1 <- length(group1_s1)
  n2_s1 <- length(group2_s1)
  BF1 <- BF(tval = tval1, n1=n1_s1, n2 = n2_s1, prior.mu=prior.mu, prior.var=prior.var)
  
  # Posterior model probability H1 (used to calculate number of samples)
  prob.H1 <- switch(forecastmodel,
                    "combined" = BF1 / (BF1 + 1),
                    "H0" = 0,
                    "H1" = 1)
  
  # Draw parameters from posterior distribution
  if(prob.H1 > 0){
    pardraws1 <- mcmc_2sample_t(mydata = c(n1_s1, n2_s1, mean(group1_s1), mean(group2_s1), var(group1_s1), var(group2_s1)),
                                priorpar = c(prior.mu, prior.var),
                                iter = mcmc.setting$iter, 
                                burnin = mcmc.setting$burnin, 
                                chains = mcmc.setting$chains,
                                alternative = alternative)
    pardraws1_comb <- apply(pardraws1, 1, rbind)
    ind1 <- sample(1:nrow(pardraws1_comb), round(prob.H1 * iter))
    param1 <- pardraws1_comb[ind1,]
  } else {
    param1 <- NULL
  }
  
  if(1-prob.H1 > 0){
    param0 <- drawH0_2sample_t(iter-round(prob.H1 * iter), group1_s1, group2_s1)
  } else {
    param0 <- matrix(NA, ncol = 2, nrow=0)
  }
  
  allparam <- rbind(param1, cbind(param0, rep(0, iter-round(prob.H1 * iter))))
  colnames(allparam) <- c("mu", "sigma", "delta")
  
  # Generate new data from posterior
  newdat <- array(NA, dim=c(group.n_s2, 2, iter))
  newdat[, 1, ] <- matrix(rnorm(group.n_s2*iter, 
                                mean = allparam[, "mu"]+(allparam[, "sigma"]*allparam[, "delta"])/2, 
                                sd=allparam[, "sigma"]), 
                          ncol=iter, byrow=T)
  newdat[, 2, ] <- matrix(rnorm(group.n_s2*iter, 
                                mean = allparam[, "mu"]-(allparam[, "sigma"]*allparam[, "delta"])/2, 
                                sd=allparam[, "sigma"]), 
                          ncol=iter, byrow=T)
  
  # Combine new data with stage 1 data
  combidat <- array(NA, dim = c(group.n_s2+max(n1_s1, n2_s1), 2, iter))
  combidat[1:n1_s1, 1, ] <- group1_s1
  combidat[1:n2_s1, 2, ] <- group2_s1
  combidat[(max(n1_s1, n2_s1)+1):(group.n_s2+max(n1_s1, n2_s1)), , ] <- newdat
  
  # Obtain combined t-values
  tvalall <- apply(combidat, 3, function(x) t.test(x[,1], x[,2], var.equal = TRUE)$statistic)
  
  # Obtain BF distribution 
  BFall <- BF(tvalall, n1 = n1_s1 + group.n_s2, n2 = n1_s1 + group.n_s2, prior.mu = prior.mu, prior.var = prior.var)
  
  # Return all results
  return(list(BF_stage_1 = BF1,
              BFdist_stage_2 = BFall,
              params = allparam,
              data = list(data1_s1 = group1_s1,
                          data2_s1 = group2_s1,
                          combidat = combidat),
              settings = list(forecastmodel = forecastmodel,
                              alternative = alternative,
                              group.n_s2 = group.n_s2,
                              prior.mu = prior.mu,
                              prior.var = prior.var)))
}

# g1 <- rnorm(100)
# g2 <- rnorm(100)-0.3
# myforecast <- BF.forecast(g1, g2, 200, forecastmodel = "combined", alternative="two.sided", prior.mu = 0, prior.var = 1)

# myforecast$BF_stage_1
# hist(myforecast$ES)
# #table(myforecast$ES)
# hist(log(myforecast$BFdist_stage_2))
