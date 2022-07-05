# ==============================================================================
# Functions to sample from the joint posterior distribution 
# p(mu, sigma, delta | y) under H1 and p(mu, sigma | y) under H0 for an 
# independent-samples t-test with a normal prior on delta. For H1, we
# use MCMC sampling. For H0, we sample sigma from the marginal distribution, 
# then mu from the conditional distribution.
# ==============================================================================

########################### H1 sampling ########################################

# Draw samples from the joint posterior distribution for delta and sigma (and mu)
# in the independent-samples t-test
#'@param mydata Vector describing the data consisting of 6 elements: Sample size group 1 and 2, means group 1 and 2, variances group 1 and 2
#'@param priorpar Vector describing the prior on delta consisting of 2 elements: mean of the prior, variance of the prior
#'@param iter Number of iterations in the MCMC algorithm (samples drawn per chain)
#'@param chains Number of chains in the MCMC algorithm
#'@param burnin Number of burn-in samples in the MCMC algorithm
#'@param alternative Direction of the test ("two.sided", "greater", "less") 

# Example vectors:
# mydata <- c(n1, n2, m1, m2, v1, v2)
# priorpar <- c(prior.mu, prior.var)

mcmc_2sample_t <- function(mydata, priorpar, iter, burnin, chains, alternative){
  
  nparam <- 3
  
  # Define restrictions on delta prior (for directional hypotheses)
  lim.low <- ifelse(alternative == "greater", 0, -Inf)
  lim.hi <- ifelse(alternative == "less", 0, Inf)
  
  # Define likelihood and prior
  likelihood <- function(x, param){
    -mean(x[1:2])*(log(2*pi)+2*log(param[2]))-(1/(2*param[2]^2)) * ((x[2]-1)*x[6] + x[2]*(x[4]-(param[1]-0.5*param[2]*param[3]))^2 + (x[1]-1)*x[5] + x[1]*(x[3]-(param[1]+0.5*param[2]*param[3]))^2)
  }
  prior <- function(param, priorpar){
    -2*log(param[2]) + LaplacesDemon::dtrunc(param[3], spec = "norm", a = lim.low, b = lim.hi, mean = priorpar[1], sd = sqrt(priorpar[2]), log = T)
  }
  
  # Initialize results object
  store <- array(NA, dim = c(nparam, chains, iter))
  store[1, , 1] <- rnorm(chains, mean(mydata[3:4]), mean(mydata[5:6]))
  store[2, , 1] <- rgamma(chains, 1, 1)
  store[3, , 1] <- LaplacesDemon::rtrunc(chains, spec = "norm", 0, 1, a=lim.low, b = lim.hi)
  
  # Run MCMC sampling procedure
  for(i in 2:iter){
    for(j in 1:chains){
      tmp1 <- sum(likelihood(mydata, store[, j, i-1])) + prior(store[, j, i-1], priorpar)
      proposal <- store[, j, i-1] + rnorm(nparam, 0, 0.2)
      proposal[3] <- LaplacesDemon::interval(proposal[3], a = lim.low, b = lim.hi)
      proposal[2] <- LaplacesDemon::interval(proposal[2], a = 0, b = Inf)
      tmp2 <- sum(likelihood(mydata, proposal))+prior(proposal, priorpar)
      if(!is.finite(tmp2)) tmp2 <- -Inf
      if(exp(tmp2-tmp1) > runif(1, 0, 1)){
        store[, j, i] <- proposal
      } else {
        store[, j, i] <- store[, j, i-1]
      }
    }
  }
  
  dimnames(store) <- list(c("mu", "sigma", "delta"), paste("c", 1:chains), paste("i", 1:iter))
  store <- store[, , (burnin+1):iter] # delete burnin samples
  
  return(store)
  
}

########################### H0 sampling ########################################

# Draw samples for sigma (and mu) in the independent-samples t-test 
#'@param nsamp: Number of samples drawn from posterior
#'@param x1: Data vector containing observations in group 1
#'@param x2: Data vector containing observations in group 2

drawH0_2sample_t <- function(nsamp, x1, x2){
  
  data <- c(x1, x2)
  
  # Draw sigma from marginal posterior
  n <- length(data)
  alpha <- (n-1)/2
  beta <- (n-1)/2 * sd(data)^2 
  sigma2 <- pscl::rigamma(nsamp, alpha, beta)
  
  # Draw mu from conditional
  xbar <- mean(data)
  mu <- rnorm(nsamp, xbar, sqrt(sigma2/n))
  
  out <- matrix(c(mu, sqrt(sigma2)), ncol = 2, byrow=F)
  colnames(out) <- c("mu", "sigma")
  
  return(out)
  
}


# Test that marginal posterior for delta from mcmc sampling is the same as 
# marginal posterior based on t-value
#
# source("ttest_2sample_normalprior/posterior_2sample_t.R")
# g1 <- rnorm(100, 11, 1)
# g2 <- rnorm(100, 12, 1)
# mydata <- c(100, 100, mean(g1), mean(g2), var(g1), var(g2))
# priorpar <- c(0, 1)
# mcmcsamp <- mcmc_2sample_t(mydata, priorpar, iter=5000, burnin=1000, chains=10, alternative="less")
# tval <- t.test(g1, g2, var.equal = TRUE)$statistic
# hist(mcmcsamp[3, , ], freq=FALSE, breaks = 40)
# x <- seq(-1.4, -0.4, by = 0.01)
# yval <- posterior(x, tval = tval, n1 = 100, n2=100, prior.mu = 0, prior.var = 1)
# points(x, yval, type="l")

