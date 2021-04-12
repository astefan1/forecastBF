# Fanplot to plot the development of the Bayes factor up to the current point
# and the predictions if a sample of size N is added

source("ttest_2sample_normalprior/forecast_2sample_t.R")

g1 <- rnorm(100)
g2 <- rnorm(100) - 0.3
resultBFforecast <- BF.forecast(g1, g2, 70, forecastmodel = "H0", alternative="two.sided", prior.mu = 0, prior.var = 1)

fanplot_BFforecast <- function(resultBFforecast, thresholds, fancolor="grey"){
  
  # Choose Bayes factor that should be calculated
  BF <- switch(resultBFforecast$settings$alternative,
               "two.sided" = BF10_norm,
               "greater" = BFplus0_norm,
               "smaller" = BFmin0_norm)
  
  # Sequential updating with stage 1 data
  
  # Get data from forecast object into proper format
  datlist <- resultBFforecast$data
  dat_s1 <- matrix(NA, ncol=2, nrow=max(length(datlist$data1_s1), 
                                        length(datlist$data2_s1)))
  dat_s1[1:length(datlist$data1_s1),1] <- datlist$data1_s1
  dat_s1[1:length(datlist$data2_s1),2] <- datlist$data2_s1
  
  # Calculate updating t-values
  tvals_s1 <- unname(sapply(2:nrow(dat_s1), 
                            function(x) t.test(dat_s1[1:x, 1], dat_s1[1:x, 2], var.equal = TRUE)$statistic))
  
  # Calculate updating Bayes factors
  dat_update <- data.frame(tvals_s1 = tvals_s1,
                           n1 = cumsum(!is.na(dat_s1[,1]))[-1],
                           n2 = cumsum(!is.na(dat_s1[,2]))[-1])
  BFs_s1 <- apply(dat_update, 1, function(x) BF(tval = x[1], 
                                                n1=x[2], 
                                                n2 = x[3], 
                                                prior.mu=resultBFforecast$settings$prior.mu, 
                                                prior.var=resultBFforecast$settings$prior.var))
  
  # Calculate necessary numbers for plot
  n_total <- nrow(dat_update) + resultBFforecast$settings$group.n_s2 + 1
  ylimvals <- c(min(log(BFs_s1),
                    log(resultBFforecast$BFdist_stage_2),
                    -0.7*max(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)-1),
                max(log(BFs_s1),
                    log(resultBFforecast$BFdist_stage_2),
                    -0.7*min(log(resultBFforecast$BFdist_stage_2)),
                    log(thresholds)+1))
  BFforecastDens <- ks::kde.boundary(log(resultBFforecast$BFdist_stage_2), 
                                     xmin = min(log(resultBFforecast$BFdist_stage_2), -0.3 * max(log(resultBFforecast$BFdist_stage_2))), 
                                     xmax = max(log(resultBFforecast$BFdist_stage_2), -0.3 * min(log(resultBFforecast$BFdist_stage_2))), 
                                     h = hpi(x = log(resultBFforecast$BFdist_stage_2)) * 2)
  scaleDens <- (pretty(c(0, n_total))[2] - pretty(c(0, n_total))[1]) * 0.5
  scaleY <- pretty(ylimvals)[2]-pretty(ylimvals)[1]
  abovethresh <- round(sum(resultBFforecast$BFdist_stage_2 >= thresholds[2]) / length(resultBFforecast$BFdist_stage_2) * 100, 1)
  belowthresh <- round(sum(resultBFforecast$BFdist_stage_2 <= thresholds[1]) / length(resultBFforecast$BFdist_stage_2) * 100, 1)
  betweenthresh <- round(100-abovethresh-belowthresh, 1)
  xtext <-  n_total + 0.5*scaleDens + max(BFforecastDens$estimate)/max(BFforecastDens$estimate) * scaleDens
  
  # Plot the plot
  par(mar = c(5,7,2,9))
  plot(x = 1:n_total, 
       y = c(log(c(1,BFs_s1)), rep(NA, resultBFforecast$settings$group.n_s2)),
       type = "l",
       ylim = ylimvals,
       xlab = "Sample size per group",
       bty="l",
       yaxt="n",
       ylab="",
       cex.lab = 1.5)
  axis(2, at=pretty(ylimvals), labels = format(exp(pretty(ylimvals)), digits=2, width=5), las = 1)
  mtext("(Expected) Bayes factor", 2, line=5, cex=1.5)
  segments(x0=nrow(dat_update)+1, x1 = n_total, y0 = log(resultBFforecast$BF_stage_1), y1 = log(resultBFforecast$BFdist_stage_2), col = scales::alpha(fancolor, 0.3))
  polygon(x = c(rep(n_total + 0.5*scaleDens, length(BFforecastDens$eval.points)), 
                n_total + 0.5*scaleDens + rev(BFforecastDens$estimate/max(BFforecastDens$estimate) * scaleDens)), 
          y = c(BFforecastDens$eval.points, rev(BFforecastDens$eval.points)), 
          xpd = NA, col = fancolor)
  segments(x0 = 0, y0 = log(thresholds), x1 = n_total + 2*scaleDens, lty = "dashed", xpd = NA)
  text(x = xtext, y = log(thresholds[2]) + 0.5 * scaleY, label = paste(abovethresh, "%"), xpd=NA, adj = c(0,0))
  text(x = xtext, y = 0, label = paste(betweenthresh, "%"), xpd=NA, adj = c(0,0.5))
  text(x = xtext, y = log(thresholds[1]) - 0.5 * scaleY, label = paste(belowthresh, "%"), xpd=NA, adj = c(0,1))
  segments(x0 = n_total + 0.5*scaleDens, y0 = ylimvals[1], y1 = ylimvals[2], xpd=NA)
}

fanplot_BFforecast(resultBFforecast = resultBFforecast, thresholds = c(1/10, 10), fancolor = "darkgrey")

