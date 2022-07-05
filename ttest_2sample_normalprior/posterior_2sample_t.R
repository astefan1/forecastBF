# ==============================================================================
# Posterior distribution for delta in an independent-samples t-test with a 
# normal prior
# ==============================================================================

# Posterior distribution
#' @param delta Parameter value of delta
#' @param tval t-value in the data
#' @param n1 Sample size group 1
#' @param prior.mu Mean of the prior distribution on delta
#' @param prior.var Variance of the prior distribution on delta

posterior <- function(delta, tval, n1, n2, prior.mu, prior.var){
  df <- n1 + n2 - 2
  n_eff <- 1/(1/n1 + 1/n2)
  scale <- sqrt(1+n_eff*prior.var)
  prior.sd <- sqrt(prior.var)
  ml1 <- 1/scale * dt(x = tval/scale, 
                      df = df, 
                      ncp = sqrt(n_eff)/scale * prior.mu)
  dt(x = tval, df = df, ncp = sqrt(n_eff) * delta) * dnorm(x = delta, mean = prior.mu, sd = prior.sd) / ml1 
}

# Example
# x <- seq(-2, 3, by=0.01)
# # yval <- posterior(delta=x, tval=tval, n1=n1, n2=n2, prior.mu=prior.mu,  prior.var=prior.var)

# Importance sampling from the posterior on delta (prior as proposal distribution)
#' @param n.sample Number of samples drawn from he posterior
#' @param tval t-value in the data
#' @param n1 Sample size group 1
#' @param prior.mu Mean of the prior distribution on delta
#' @param prior.var Variance of the prior distribution on delta
#' @param alternative Direction of alternative hypothesis ("two.sided", "greater", "less")

sample_posterior <- function(n.sample, tval, n1, n2, prior.mu, prior.var, alternative = "two.sided"){
  
  df <- n1 + n2 - 2
  prior.sd <- sqrt(prior.var)
  n_eff <- 1/(1/n1 + 1/n2)
  
  # Define sampling bounds
  lim.low <- ifelse(alternative == "greater", 0, -Inf)
  lim.hi <- ifelse(alternative == "less", 0, Inf)
  
  # samples from proposal distribution
  deltas <- LaplacesDemon::rtrunc(n = max(n.sample*2, 1000), spec = "norm", a = lim.low, b = lim.hi, mean = prior.mu, sd = prior.sd)
  # compute weights
  ws <- dt(x = tval, df = df, ncp = sqrt(n_eff) * deltas)/sum(dt(x = tval, df = df, ncp = sqrt(n_eff) * deltas))
  # draw from posterior
  draws <- sample(deltas, size=n.sample, replace = TRUE, prob = ws)
  
  return(draws)
}

# Example

# tval = 1
# n1=20
# n2=20
# prior.mu=1
# prior.var=1.5
# 
# samp  <- sample_posterior(n.sample = 10000, alternative="two.sided", tval=tval, n1=n1, n2=n2, prior.mu=prior.mu, prior.var=prior.var)

