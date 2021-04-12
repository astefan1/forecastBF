# BFDA for independent-samples t-test

source("bf_2sample_t.R")

BFDA_2sample_t <- function(pop.ES, alternative, group.n, prior.mu, prior.var, iter = 1000){
  
  # Which BF to compute?
  BF <- switch(alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  # Generate t-values from design prior
  if(length(pop.ES) == 1){
    tvals <- rt(n = iter, df = 2*group.n-2, ncp = pop.ES / sqrt(2/group.n))
    pop.ES <- rep(pop.ES, iter)
  } else {
    stopifnot("Number of iterations does not equal number of samples from design prior." = length(pop.ES) == iter)
    tvals <- sapply(pop.ES, function(x) rt(n = 1, df = 2*group.n-2, ncp = x / sqrt(2/group.n)))
  }
  
  # Calculate BFs from t-values
  bfdist <- BF(tvals, n1 = group.n, n2 = group.n, prior.mu = prior.mu, prior.var = prior.var)

  return(data.frame(ES = pop.ES,
                    t.value = tvals,
                    BF = bfdist))
  
}


# Takes ~ 0.2 seconds
# a <- Sys.time()
# bfda1 <- BFDA_2sample_t(rnorm(1000), alternative="two.sided", group.n = 300, prior.mu = 0.5, prior.var = 1, iter = 1000)
# Sys.time() - a

# Takes ~ 22 seconds
# a <- Sys.time()  
# bfda2 <- BFDA::BFDA.sim(expected.ES = 0.5,
#                         type="t.between",
#                         prior=list("normal", list(prior.mean = 0.5, prior.variance = 1)),
#                         n.max = 300,
#                         design = "fixed.n",
#                         B = 1000,
#                         alternative = "two.sided",
#                         verbose = FALSE,
#                         seed = NULL)
# Sys.time() - a
