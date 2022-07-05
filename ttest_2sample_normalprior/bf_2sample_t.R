# ==============================================================================
# Functions to compute Bayes factors for independent-samples t-tests with 
# normal priors
# ==============================================================================

# Function to compute the Bayes factor for a two-sided hypothesis
#' @param tval t value calculated based on sample
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @param prior.mu Mean of the prior distribution on delta under H1
#' @param prior.var Variance of the prior distribution on delta under H1

BF10_norm <- function(tval, n1, n2, prior.mu=0, prior.var=1){
  
  n_eff <- 1/(1/n1 + 1/n2)
  invscale <- 1/sqrt(1+n_eff*prior.var)
  df <- n1 + n2 -2
  
  ml1 <- invscale * dt(x = tval * invscale, 
                      df = df, 
                      ncp = sqrt(n_eff)* invscale * prior.mu)
  ml0 <- dt(x = tval, df = df, ncp = 0)
  
  return(ml1/ml0)
}

# Function to compute the Bayes factor for a negative directional hypothesis
#' @param tval t value calculated based on sample
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @param prior.mu Mean of the prior distribution on delta under H1
#' @param prior.var Variance of the prior distribution on delta under H1

BFmin0_norm <- function(tval, n1, n2, prior.mu=0, prior.var=1){
  
  n_eff <- 1/(1/n1 + 1/n2)
  scale <- sqrt(1+n_eff*prior.var)
  df <- n1 + n2 - 2
  prior.sd <- sqrt(prior.var)
  
  ml1 <- 1/scale * dt(x = tval/scale, 
                      df = df, 
                      ncp = sqrt(n_eff)/scale * prior.mu)
  ml0 <- dt(x = tval, df = df, ncp = 0)
  
  priorAreaSmaller0 <- pnorm(0, mean=prior.mu, sd=prior.sd)
  
  posterior <- function(delta){
    dt(x = tval, df = df, ncp = sqrt(n_eff) * delta) * dnorm(x = delta, mean = prior.mu, sd = prior.sd) / ml1 
  }
  
  postAreaSmaller0 <- integrate(posterior, lower=-Inf, upper=0)$value
  
  return(postAreaSmaller0 / priorAreaSmaller0 * ml1 / ml0)
}

BFmin0_norm <- Vectorize(BFmin0_norm, vectorize.args = "tval")

# Function to compute the Bayes factor for a positive directional hypothesis
#' @param tval t value calculated based on sample
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @param prior.mu Mean of the prior distribution on delta under H1
#' @param prior.var Variance of the prior distribution on delta under H1

BFplus0_norm <- function(tval, n1, n2, prior.mu=0, prior.var=1){
  
  n_eff <- 1/(1/n1 + 1/n2)
  scale <- sqrt(1+n_eff*prior.var)
  df <- n1 + n2 -2
  prior.sd <- sqrt(prior.var)
  
  ml1 <- 1/scale * dt(x = tval/scale, 
                      df = df, 
                      ncp = sqrt(n_eff)/scale * prior.mu)
  ml0 <- dt(x = tval, df = df, ncp = 0)
  
  priorAreaLarger0 <- pnorm(0, mean=prior.mu, sd=prior.sd, lower.tail = FALSE)
  
  posterior <- function(delta){
    dt(x = tval, df = df, ncp = sqrt(n_eff) * delta) * dnorm(x = delta, mean = prior.mu, sd = prior.sd) / ml1  
  }
  
  postAreaLarger0 <- tryCatch(1-integrate(posterior, lower=-Inf, upper=0)$value,
                              error = function(e) {
                                if(grepl("non-finite function value", e)){
                                  return(NA)
                                }
                                })
  
  BF <- postAreaLarger0 / priorAreaLarger0 * ml1 / ml0
  BF[is.na(BF)] <- ifelse(tval > 0, Inf, 0)
  
  return(BF)
}

BFplus0_norm <- Vectorize(BFplus0_norm, vectorize.args = "tval")

# Example

# BFmin0_norm(seq(-1, 1, length.out = 1000), 30, 30, prior.mu=0.5, prior.var=2)
