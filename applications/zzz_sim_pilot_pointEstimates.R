# ==============================================================================
# Internal pilot study: Bayes factor forecast based on point estimates
# ==============================================================================

source("ttest_2sample_normalprior/bf_2sample_t.R")

#'@param group1_s1 Observed data group 1
#'@param group2_s1 Observed data group 2
#'@param group.n_s2 Number of planned observations (how many observations do you look into the future)
#'@param alternative Direction of alternative hypothesis ("two.sided", "greater", "less")
#'@param prior.mu Mean of the prior distribution on delta
#'@param prior.var Variance of the prior distribution on delta
#'@param iter Number of Monte Carlo iterations for the forecast 

BF.forecast.point<- function(group1_s1, group2_s1, group.n_s2, forecastmodel = "combined", alternative = "two.sided", prior.mu = 0, prior.var = 1, iter = 1000){
  
  # Which BF to compute?
  BF <- switch(alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  # t value stage 1
  ttestres <- t.test(group1_s1, group2_s1, var.equal = TRUE)
  tval1 <- unname(ttestres$statistic)
  
  # effect size stage 1
  ES <- tval1 * sqrt(2/length(group1_s1))
  
  # Bayes factor stage 1
  n1_s1 <- length(group1_s1)
  n2_s1 <- length(group2_s1)
  BF1 <- BF(tval = tval1, n1=n1_s1, n2 = n2_s1, prior.mu=prior.mu, prior.var=prior.var)
  
  # Posterior model probability H1 (used to calculate number of samples)
  prob.H1 <- switch(forecastmodel,
                    "combined" = BF1 / (BF1 + 1),
                    "H0" = 0,
                    "H1" = 1)
  
  # Create new data based on point estimates
  if(prob.H1 > 0){
    param1 <- rep(ES, round(prob.H1 * iter))
  } else {
    param1 <- NULL
  }
  
  if(1-prob.H1 > 0){
    param0 <- rep(0, iter-round(prob.H1 * iter))
  } else {
    param0 <- NULL
  }
  
  allparam <- c(param1, param0)
  
  newdat <- array(NA, dim=c(group.n_s2, 2, iter))
  newdat[, 1, ] <- matrix(rnorm(group.n_s2*iter, 
                                mean = mean(ttestres$estimate)+(ttestres$stderr/sqrt(2/n1_s1))*allparam/2, 
                                sd=ttestres$stderr/sqrt(2/n1_s1)), # pooled sd estimator as sd
                          ncol=iter, byrow=T)
  newdat[, 2, ] <- matrix(rnorm(group.n_s2*iter, 
                                mean = mean(ttestres$estimate)-(ttestres$stderr/sqrt(2/n1_s1))*allparam/2, 
                                sd = ttestres$stderr/sqrt(2/n1_s1)), 
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
              params = ES,
              data = list(data1_s1 = group1_s1,
                          data2_s1 = group2_s1,
                          combidat = combidat),
              settings = list(forecastmodel = forecastmodel,
                              alternative = alternative,
                              group.n_s2 = group.n_s2,
                              prior.mu = prior.mu,
                              prior.var = prior.var)))
  
  
}

#' @param n.pilot Sample size in the pilot study
#' @param n.main Sample size in the main study (can be a vector of length >1)
#' @param ES.pop Population effect size
#' @param alternative Direction of alternative hypothesis ("two-sided", "greater", "less")
#' @param prior.mu Mean of the prior distribution on delta
#' @param prior.var Variance of the prior distribution on delta
#' @param iter Number of Monte Carlo iterations for the forecast 
#' @param forecastmodel Under which model should the forecast be computed? ("combined" = model-averaged, "H1" = alternative, "H0" = null)

sim.pilot.point <- function(n.pilot, n.main, ES.pop, alternative = "two.sided", prior.mu=0, prior.var=1, iter = 1000, forecastmodel = "combined", mcmc.setting = list(chains=2, burnin=0, iter=5000)){
  
  # Simulate data for n.pilot
  g1 <- rnorm(n.pilot+max(n.main))
  g2 <- rnorm(n.pilot+max(n.main))-ES.pop
  
  # Run BF forecast based on pilot data
  forecast <- BF.forecast.point(group1_s1=g1[1:n.pilot],
                                group2_s1=g2[1:n.pilot], 
                                group.n_s2=max(n.main), 
                                forecastmodel = forecastmodel, 
                                alternative = alternative, 
                                prior.mu = prior.mu, 
                                prior.var = prior.var, 
                                iter = iter)
                          
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

#' @param n.pilot Sample size in the pilot study
#' @param n.main Sample size in the main study (can be a vector of length >1)
#' @param ES.pop Population effect size
#' @param alternative Direction of alternative hypothesis ("two-sided", "greater", "less")
#' @param prior.mu Mean of the prior distribution on delta
#' @param prior.var Variance of the prior distribution on delta
#' @param iter Number of Monte Carlo iterations for the forecast 
#' @param forecastmodel Under which model should the forecast be computed? ("combined" = model-averaged, "H1" = alternative, "H0" = null)
#' @param iter.DA Number of Monte Carlo iterations for the design analysis
#' @param mcmc.setting List determining the settings of the MCMC sampling of posterior under H1 ($chains = number of chains, $iter = number of iterations in the MCMC algorithm, $burnin = number of burn-in samples)

designAnalysis.pilot.point <- function(n.pilot, n.main, ES.pop, alternative = "two.sided", prior.mu=0, prior.var=1, iter = 1000, forecastmodel = "combined", iter.DA = 1000){
  
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
                          forecastmodel = forecastmodel)
  }
  
  return(list(res.DA = res,
              settings = list(n.pilot=n.pilot,
                              n.main=n.main,
                              ES.pop=ES.pop, 
                              alternative = alternative, 
                              prior.mu=prior.mu, 
                              prior.var=prior.var, 
                              forecastmodel = forecastmodel,
                              iter = iter)))
}


