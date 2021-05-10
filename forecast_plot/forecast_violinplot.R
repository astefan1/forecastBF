# 2-part density plot showing BF distribution under H0 and H1. The width of the 
# two distributions represents the posterior probability of H0 and H1 at the point of
# the current sample.

# source("ttest_2sample_normalprior/forecast_2sample_t.R")
# g1 <- rnorm(10)
# g2 <- rnorm(10) - 0
# resultBFforecast <- BF.forecast(g1, g2, 50, forecastmodel = "combined", alternative="greater", prior.mu = 0, prior.var = 1)


forecast_violinplot <- function(resultBFforecast, thresholds, distcolors=NULL){
  
  # Posterior model probabilities
  prob.H1 <- switch(resultBFforecast$settings$forecastmodel,
                    "combined" = resultBFforecast$BF_stage_1 / (resultBFforecast$BF_stage_1 + 1),
                    "H0" = 0,
                    "H1" = 1)
  
  # Number of samples from H0 and H1
  iter <- length(resultBFforecast$BFdist_stage_2)
  nH1 <- round(iter*prob.H1)
  nH0 <- iter - nH1
  
  # Y axis limits
  ylimvals <- c(min(log(resultBFforecast$BFdist_stage_2),
                    -0.7*max(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)-1),
                max(log(resultBFforecast$BFdist_stage_2),
                    -0.7*min(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)+1))
  
  if(prob.H1 > 0){
    BFDensH1 <- density(log(resultBFforecast$BFdist_stage_2[1:nH1]))
    scaleDensH1 <- prob.H1
  }
  
  if(prob.H1 < 1){
    BFDensH0 <- density(log(resultBFforecast$BFdist_stage_2[(nH1+1):iter]))
    scaleDensH0 <- 1-prob.H1
  }
  
  
  # Set colors
  if(is.null(distcolors)){
    distcolors <- c("#88CCEE", "#CC6677")
  } else if(length(distcolors) == 1) {
    distcolors <- rep(fancolor, 2)
  }

  # Initialize Plot
  
  par(mar = c(2, 10, 3, 10))
  
  plot(c(0,0), 
       ylimvals, 
       type="l", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", xlim = c(-1, 1), 
       ylim = ylimvals)
  
  # Density H1
  polygon(x = c(rep(0, length(BFDensH1$x)),
                rev(BFDensH1$y/max(BFDensH1$y) * scaleDensH1)),
          y = c(BFDensH1$x, rev(BFDensH1$x)),
          xpd = NA, col = distcolors[1])

  # Density H0
  polygon(x = -c(rep(0, length(BFDensH0$x)),
                rev(BFDensH0$y/max(BFDensH0$y) * scaleDensH0)),
          y = c(BFDensH0$x, rev(BFDensH0$x)),
          xpd = NA, col = distcolors[2])
  
  # Thresholds
  segments(x0 = -1, y0 = c(log(thresholds), 0), x1 = 1, lty = "dashed", xpd = NA)
  text(x = 1, y = log(thresholds[1]), labels = bquote("BF"[10]*" = "*.(round(thresholds[1], 2))), xpd=NA, pos = 4)
  text(x = 1, y = 0, labels = bquote("BF"[10]*" = 1"), xpd = NA, pos = 4)
  text(x = 1, y = log(thresholds[2]), labels = bquote("BF"[10]*" = "*.(round(thresholds[2], 2))), xpd=NA, pos = 4)
  
  text(x = -1, y = log(thresholds[1]), labels = bquote("BF"[10]*" = "*.(round(thresholds[1], 2))), xpd=NA, pos = 2)
  text(x = -1, y = 0, labels = bquote("BF"[10]*" = 1"), xpd = NA, pos = 2)
  text(x = -1, y = log(thresholds[2]), labels = bquote("BF"[10]*" = "*.(round(thresholds[2], 2))), xpd=NA, pos = 2)
  
  # Posterior model probabilities
  mtext(paste("p(H0 | data) = ", round(1-prob.H1, 2)), at = -0.6, side = 3, cex = 1.3, font = 2)
  mtext(paste("p(H1 | data) = ", round(prob.H1, 2)), at = 0.6, side = 3, cex = 1.3, font = 2)
  
  # Add brackets for evidence categories
  tick1 <- -round((ylimvals[2]-log(thresholds[2]))/(ylimvals[2]-ylimvals[1]), 2)
  tick2 <- -round((ylimvals[2]-log(thresholds[1]))/(ylimvals[2]-ylimvals[1]), 2)
  pBrackets::brackets(x1 = 1.5, y1 = ylimvals[2], x2 = 1.5, y2 =  ylimvals[1], type = 4, xpd = NA, ticks = c(tick1, tick2), h = 0.1)
  pBrackets::brackets(x1 = -1.5, y1 = ylimvals[1], x2 = -1.5, y2 =  ylimvals[2], type = 4, xpd = NA, ticks = c(-(tick1+1), -(tick2+1)), h = 0.1)
  
  # Add percentages
  pH1.H0 <- round(sum(resultBFforecast$BFdist_stage_2[(nH1+1):iter] > thresholds[2])/nH0*100, 1)
  pH0.H0 <- round(sum(resultBFforecast$BFdist_stage_2[(nH1+1):iter] < thresholds[1])/nH0*100, 1)
  pinc.H0 <- round(100-pH1.H0-pH0.H0, 1)
  
  pH1.H1 <- round(sum(resultBFforecast$BFdist_stage_2[1:nH1] > thresholds[2])/nH1*100, 1)
  pH0.H1 <- round(sum(resultBFforecast$BFdist_stage_2[1:nH1] < thresholds[1])/nH1*100, 1)
  pinc.H1 <- round(100-pH1.H1-pH0.H1, 1)
  
  text(x = 1.7, y = ylimvals[2]-0.5*(ylimvals[2]-log(thresholds[2])), paste(pH1.H1, "%"), xpd = NA, pos = 4)
  text(x = 1.7, y = 0, paste(pinc.H1, "%"), xpd = NA, pos = 4)
  text(x = 1.7, y = ylimvals[1]-0.5*(ylimvals[1]-log(thresholds[1])), paste(pH0.H1, "%"), xpd = NA, pos = 4)
  
  text(x = -1.7, y = ylimvals[2]-0.5*(ylimvals[2]-log(thresholds[2])), paste(pH1.H0, "%"), xpd = NA, pos = 2)
  text(x = -1.7, y = 0, paste(pinc.H0, "%"), xpd = NA, pos = 2)
  text(x = -1.7, y = ylimvals[1]-0.5*(ylimvals[1]-log(thresholds[1])), paste(pH0.H0, "%"), xpd = NA, pos = 2)
  
}

# forecast_violinplot(resultBFforecast, thresholds = c(1/10, 10))
