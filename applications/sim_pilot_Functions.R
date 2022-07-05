# ==============================================================================
# Simulation: Internal pilot study
# ==============================================================================

source("ttest_2sample_normalprior/forecast_2sample_t.R")
source("ttest_2sample_normalprior/mcmc_2sample_t.R")

#' @param n.pilot Sample size in the pilot study
#' @param n.main Sample size in the main study (can be a vector of length >1)
#' @param ES.pop Population effect size

sim.pilot <- function(n.pilot, n.main, ES.pop, alternative = "two.sided", prior.mu=0, prior.var=1, iter = 1000, forecastmodel = "combined", mcmc.setting = list(chains=2, burnin=0, iter=5000)){
  
  # Simulate data for n.pilot
  g1 <- rnorm(n.pilot+max(n.main))
  g2 <- rnorm(n.pilot+max(n.main))-ES.pop
  
  # Run BF forecast based on pilot data
  forecast <- BF.forecast(group1_s1=g1[1:n.pilot],
                          group2_s1=g2[1:n.pilot], 
                          group.n_s2=max(n.main), 
                          forecastmodel = forecastmodel, 
                          alternative = alternative, 
                          prior.mu = prior.mu, 
                          prior.var = prior.var, 
                          iter = iter, 
                          mcmc.setting = mcmc.setting)
  
  BFs <- NULL
  
  # Get BFs from all stages of the forecast
  if(length(n.main) > 1){
    
    # Choose Bayes factor that should be calculated depending on test direction
    BF <- switch(alternative,
                 "two.sided" = BF10_norm,
                 "greater" = BFplus0_norm,
                 "smaller" = BFmin0_norm)
    
    # Get data from forecast object
    datlist <- forecast$data$combidat
    
    # At what sample sizes to compute the BF
    stepindx <- (n.pilot+n.main)[-length(n.main)]
    
    # Calculate updating BFs
    tvals <- apply(datlist, 3, function(M) unname(sapply(stepindx, 
                                                        function(x) t.test(M[1:x, 1], M[1:x, 2], var.equal = TRUE)$statistic)))
    
    tvals <- as.matrix(tvals)
    
    dat_update <- array(NA, dim = c(nrow(tvals), ncol(tvals), 3))
    dat_update[, , 1] <- tvals
    dat_update[, , 2] <- matrix(cumsum(!is.na(datlist[,1,1]))[stepindx], nrow=nrow(tvals), ncol=ncol(tvals), byrow=FALSE)
    dat_update[, , 3] <- matrix(cumsum(!is.na(datlist[,2,1]))[stepindx], nrow=nrow(tvals), ncol=ncol(tvals), byrow=FALSE)
    
    BFs <- apply(dat_update, 2, function(M) apply(M, 1,
                                                  function(x) tryCatch(
                                                    BF(tval = x[1],
                                                       n1=x[2], 
                                                       n2 = x[3], 
                                                       prior.mu=prior.mu, 
                                                       prior.var=prior.var),
                                                                 
                                                    error = function(e) {
                                                      print(paste0("tval = ", x[1], ", n1 = ", x[2], ", n2 = ", x[3], ", prior.mu = ", prior.mu, ", prior.var = ", prior.var)
                                                            )})))
                                                                       
  }
  
  BFs <- cbind(t(BFs), forecast$BFdist_stage_2)
  
  return(list(group1_s1 = g1,
              group2_s1 = g2,
              n.pilot = n.pilot,
              n.main = n.main,
              BF_s1 = forecast$BF_stage_1,
              BFs = BFs,
              ES.pop = ES.pop,
              alternative = alternative, 
              prior.mu=prior.mu, 
              prior.var=prior.var))
  
}

designAnalysis.pilot <- function(n.pilot, n.main, ES.pop, alternative = "two.sided", prior.mu=0, prior.var=1, iter = 1000, forecastmodel = "combined", iter.DA = 1000, mcmc.setting = list(chains=2, burnin=0, iter=5000)){
  
  res <- vector(mode = "list", length = iter.DA)
  
  for(j in 1:iter.DA){
    cat(paste0(j, "\n"))
    res[[j]] <- sim.pilot(n.pilot, 
                          n.main, 
                          ES.pop, 
                          alternative = alternative, 
                          prior.mu=prior.mu, 
                          prior.var=prior.var, 
                          iter = iter, 
                          forecastmodel = forecastmodel,
                          mcmc.setting = mcmc.setting)
  }
  
  return(list(res.DA = res,
              settings = list(n.pilot=n.pilot,
                              n.main=n.main,
                              ES.pop=ES.pop, 
                              alternative = alternative, 
                              prior.mu=prior.mu, 
                              prior.var=prior.var, 
                              forecastmodel = forecastmodel,
                              iter = iter, 
                              mcmc.setting = mcmc.setting)))
}

