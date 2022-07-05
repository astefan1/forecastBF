source("applications/sim_pilot_Functions.R")

checkThreshold <- function(x, boundary, which){
  switch(which,
         upper = x > boundary,
         lower = x < boundary,
         both = x > boundary[2] | x < boundary[1])
}

# SIM_names <- dir(".", pattern = "SIM_pilot_")
# lapply(SIM_names, load, .GlobalEnv)

#' Design analysis for pilot study design
#' @param simResult Result from sim.pilot()
#' @param boundary Decision boundary: What is considered strong evidence? Provide vector of length 2
#' @param power Envisioned probability of strong evidence (for H0 or H1)
#' @param abandon What should happen if intermediate design analysis shows probability to reach strong evidence at n.max < power? abandon = TRUE = stop with n.pilot; abandon = FALSE = use largest SS
#' @param nmin What should happen if BF at n.pilot is larger than boundary? If nmin < n.pilot, stop at n.pilot; else sample nmin

analyzePilot <- function(simResult, boundary = c(1/10, 10), power = 0.8, abandon = TRUE, nmin = 0, whichThreshold = "both"){
  
  sim <- simResult[[1]]
  strongEv <- matrix(NA, nrow = length(sim), ncol = length(simResult$settings$n.main))
  strongEvSS <- rep(NA, length(sim))
  
  for(i in 1:length(sim)){
    
    # How many trajectories show strong evidence at each possible sample size?
    strongEv[i,] <- apply(sim[[i]]$BFs, 2, function(x) {sum(checkThreshold(x, boundary, whichThreshold))})
    
    # Which sample size is the first to have sufficient power?
    if(checkThreshold(sim[[i]]$BF_s1, boundary, whichThreshold)){
      strongEvSS[i] <- sim[[i]]$n.pilot
    } else if (!any(strongEv[i,] > power*length(sim))){
      strongEvSS[i] <- NA
    } else {
      strongEvSS[i] <- sim[[i]]$n.pilot+sim[[i]]$n.main[which(strongEv[i,] > power*length(sim))[1]]
    }
  }
  
  # Choose sample size
  SSFinal <- strongEvSS
  SSFinal[SSFinal < nmin] <- nmin
  if(abandon == TRUE){
    SSFinal[is.na(SSFinal)] <- sim[[1]]$n.pilot
  } else {
    SSFinal[is.na(SSFinal)] <- sim[[1]]$n.pilot+max(sim[[1]]$n.main)
  }
  
  # Calculate final BF
  BF <- switch(sim[[1]]$alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  tFinal <- unname(sapply(c(1:length(sim)), function(x) t.test(sim[[x]]$group1_s1[1:SSFinal[x]], 
                                                               sim[[x]]$group2_s1[1:SSFinal[x]], 
                                                               var.equal = TRUE)$statistic))
  BFFinal <- sapply(c(1:length(sim)), function(x) BF(tFinal[x],
                                                     n1 = SSFinal[x],
                                                     n2 = SSFinal[x],
                                                     prior.mu = sim[[x]]$prior.mu,
                                                     prior.var = sim[[x]]$prior.var))
  
  # Calculate final Cohen's d
  
  dFinal <- tFinal*sqrt(2/SSFinal)
  SEDFinal <- sqrt((2*SSFinal)/(SSFinal^2) + (dFinal^2) / (4*SSFinal))
  
  # Calculate maximum power
  maxPower <- apply(strongEv, 1, max)/length(sim)
  
  return(data.frame(maxPower, SSFinal, BFFinal, dFinal, SEDFinal))
  
} 

metaPilot <- function(resPilot){
  
  metafor::rma(yi = resPilot$dFinal, sei = resPilot$SEDFinal, method = "FE")
  
}
