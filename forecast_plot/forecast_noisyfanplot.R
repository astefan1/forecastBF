# Noisy fanplot with step size definition. The smaller the step size the longer
# it takes to plot. To show some noise in the trajectories, the default setting 
# works decently (roughly 5 steps if group n stage 2 >= 5). If stepsize >= group
# n at stage 2, the plot is identical to the non-noisy plot (just takes slightly
# longer).

source("ttest_2sample_normalprior/forecast_2sample_t.R")
source("ttest_2sample_normalprior/mcmc_2sample_t.R")
# 
# g1 <- rnorm(20)
# g2 <- rnorm(20) - 0.1
# resultBFforecast <- BF.forecast(g1, g2, 50, forecastmodel = "combined", alternative="greater", prior.mu = 0, prior.var = 1)

noisyfanplot_BFforecast <- function(resultBFforecast, thresholds, fancolor=NULL, stepsize = NULL, showdensity = TRUE){
  
  # Choose Bayes factor that should be calculated depending on test direction
  BF <- switch(resultBFforecast$settings$alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  # Get data from forecast object
  datlist <- resultBFforecast$data$combidat
  
  # Get stage 1 and total sample size
  n_s1 <- dim(datlist)[1]-resultBFforecast$settings$group.n_s2
  n_total <- dim(datlist)[1]
  
  # Split data in stage 1 data (existing dataset) and stage 2 data (simulated data)
  dat_s1 <- datlist[1:n_s1, , ]
  dat_s2 <- datlist[(n_s1+1):n_total, , ]
  
  # Get indices for updating steps
  if(is.null(stepsize)){
    stepindx <- unique(round(seq(n_s1, n_total, length.out = 5)))
  } else {
    stepindx <- unique(c(round(seq(n_s1, n_total, by = stepsize)), n_total))
  }
  
  # Calculate updating t-values
  tvals_s1 <- unname(sapply(2:n_s1, function(x) t.test(dat_s1[1:x, 1, 1], dat_s1[1:x, 2, 1], var.equal = TRUE)$statistic))
  tvals_s2 <- apply(datlist, 3, function(M) unname(sapply(stepindx, 
                                                          function(x) t.test(M[1:x, 1], M[1:x, 2], var.equal = TRUE)$statistic)))
  
  # Calculate updating Bayes factors
  dat_update_s1 <- data.frame(tvals_s1 = tvals_s1,
                              n1 = cumsum(!is.na(dat_s1[,1, 1]))[-1],
                              n2 = cumsum(!is.na(dat_s1[,2, 1]))[-1])
  BFs_s1 <- apply(dat_update_s1, 1, function(x) BF(tval = x[1],
                                                   n1=x[2], 
                                                   n2 = x[3], 
                                                   prior.mu=resultBFforecast$settings$prior.mu, 
                                                   prior.var=resultBFforecast$settings$prior.var))
  
  dat_update_s2 <- array(NA, dim = c(nrow(tvals_s2), ncol(tvals_s2), 3))
  dat_update_s2[, , 1] <- tvals_s2
  dat_update_s2[, , 2] <- matrix(cumsum(!is.na(datlist[,1,1]))[stepindx], nrow=nrow(tvals_s2), ncol=ncol(tvals_s2), byrow=FALSE)
  dat_update_s2[, , 3] <- matrix(cumsum(!is.na(datlist[,2,1]))[stepindx], nrow=nrow(tvals_s2), ncol=ncol(tvals_s2), byrow=FALSE)
  
  BFs_s2 <- apply(dat_update_s2, 2, function(M) apply(M, 1,
                                                      function(x) BF(tval = x[1],
                                                                     n1=x[2], 
                                                                     n2 = x[3], 
                                                                     prior.mu=resultBFforecast$settings$prior.mu, 
                                                                     prior.var=resultBFforecast$settings$prior.var)))
  
  # Calculate necessary numbers for plot
  
  # Set fan color
  if(is.null(fancolor)){
    fancolor <- c("#88CCEE", "#CC6677")
  } else if(length(fancolor) == 1) {
    fancolor <- rep(fancolor, 2)
  }
  
  # y axis limits
  ylimvals <- c(min(log(BFs_s1),
                    log(BFs_s2),
                    -0.7*max(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)-1),
                max(log(BFs_s1),
                    log(BFs_s2),
                    log(resultBFforecast$BFdist_stage_2),
                    -0.7*min(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)+1))
  
  # BF distribution: density, width
  scaleDens <- (pretty(c(0, n_total))[2] - pretty(c(0, n_total))[1]) * 0.5
  
  if(showdensity){
    BFforecastDens <- ks::kde.boundary(log(resultBFforecast$BFdist_stage_2), 
                                       xmin = min(log(resultBFforecast$BFdist_stage_2), -0.3 * max(log(resultBFforecast$BFdist_stage_2))), 
                                       xmax = max(log(resultBFforecast$BFdist_stage_2), -0.3 * min(log(resultBFforecast$BFdist_stage_2))), 
                                       h = ks::hpi(x = log(resultBFforecast$BFdist_stage_2)) * 2)
    
    # distance between two points on y axis
    scaleY <- pretty(ylimvals)[2]-pretty(ylimvals)[1]
    
    # what percentage of BFs is above/below/between thresholds
    abovethresh <- round(sum(resultBFforecast$BFdist_stage_2 >= thresholds[2]) / length(resultBFforecast$BFdist_stage_2) * 100, 1)
    belowthresh <- round(sum(resultBFforecast$BFdist_stage_2 <= thresholds[1]) / length(resultBFforecast$BFdist_stage_2) * 100, 1)
    betweenthresh <- round(100-abovethresh-belowthresh, 1)
    
    # x axis location of density labels
    xtext <-  n_total + 0.5*scaleDens + max(BFforecastDens$estimate)/max(BFforecastDens$estimate) * scaleDens
    
  }
  
  # number of draws from H1
  prob.H1 <- switch(resultBFforecast$settings$forecastmodel,
                    "combined" = resultBFforecast$BF_stage_1 / (resultBFforecast$BF_stage_1 + 1),
                    "H0" = 0,
                    "H1" = 1)
  nH1 <- round(length(resultBFforecast$BFdist_stage_2)*prob.H1)
  
  # Plot the plot
  
  par(mar = c(5,7,2,ifelse(showdensity,9,2)))
  
  # BF trajectory existing sample
  plot(x = 1:n_s1, y = c(log(c(1,BFs_s1))),
       type = "l", ylim = ylimvals, xlim = c(0, n_total), 
       xlab = "Sample size per group", bty="l", yaxt="n", ylab="", 
       cex.lab = 1.5)
  
  # axis and label
  axis(2, at=pretty(ylimvals), labels = format(exp(pretty(ylimvals)), 
                                               digits=2, width=5), las = 1)
  mtext("(Projected) Bayes factor", 2, line=5, cex=1.5)
  
  # trajectories for BF projection
  for(i in 1:length(resultBFforecast$BFdist_stage_2)){
    trajectorycolor <- ifelse(i <= nH1, fancolor[1], fancolor[2])
    points(c(n_s1, stepindx), log(c(tail(BFs_s1, 1), BFs_s2[, i])), type="l", 
           col = scales::alpha(trajectorycolor, 0.3))
  }
  
  # density
  if(showdensity){
    polygon(x = c(rep(n_total + 0.5*scaleDens, length(BFforecastDens$eval.points)), 
                  n_total + 0.5*scaleDens + rev(BFforecastDens$estimate/max(BFforecastDens$estimate) * scaleDens)), 
            y = c(BFforecastDens$eval.points, rev(BFforecastDens$eval.points)), 
            xpd = NA, col = "grey")
    text(x = xtext, y = log(thresholds[2]) + 0.5 * scaleY, label = paste(abovethresh, "%"), xpd=NA, adj = c(0,0))
    text(x = xtext, y = 0, label = paste(betweenthresh, "%"), xpd=NA, adj = c(0,0.5))
    text(x = xtext, y = log(thresholds[1]) - 0.5 * scaleY, label = paste(belowthresh, "%"), xpd=NA, adj = c(0,1))
    segments(x0 = n_total + 0.5*scaleDens, y0 = ylimvals[1], y1 = ylimvals[2], xpd=NA)
    
  }
  
  # thresholds for compelling evidence
  segments(x0 = 0, y0 = log(thresholds), x1 = n_total + ifelse(showdensity, 2*scaleDens, 0), lty = "dashed", xpd = NA)
  
  
}

# noisyfanplot_BFforecast(resultBFforecast = resultBFforecast, thresholds = c(1/10, 10), stepsize = 10, showdensity = T)
