# ==============================================================================
# Stopping for futility: Stopping scenarios illustration
# ==============================================================================

# Function to draw a simple plot with predictions at any decision stage during 
# the sequential sampling process

#'@param SIM_stopfutil A result of the designAnalysis.stopfutil() function
#'@param ID Which Monte Carlo iteration of the design analysis?
#'@param stage What stage in the sequential sampling process

simplefan <- function(SIM_stopfutil, ID = 1, stage = 1){
  
  SIM <- SIM_stopfutil[[1]][[ID]]
  
  sim <- SIM[[stage]]
  
  if(!is.null(sim$n.forecast)){
    qnt_BFs_s2 <- quantile(sim$BF.forecast, probs = c(0.05, seq(0.1, 0.9, 0.1), 0.95))
  } else{
    qnt_BFs_s2 <- NA
  }
  
  # y axis limits
  ylimvals <- c(min(log(sim$BF.temp),
                    log(qnt_BFs_s2),
                    log(1/10)-1, na.rm = TRUE),
                max(log(sim$BF.temp),
                    log(qnt_BFs_s2),
                    log(10)+1, na.rm = TRUE))
  
  scaleDens <- (pretty(c(0, 100))[2] - pretty(c(0, 100))[1]) * 0.5
  
  if(!is.null(sim$n.forecast)){
    BFforecastDens <- ks::kde.boundary(log(sim$BF.forecast), 
                                       xmin = min(log(sim$BF.forecast), -0.3 * max(log(sim$BF.forecast))), 
                                       xmax = max(log(sim$BF.forecast), -0.3 * min(log(sim$BF.forecast))), 
                                       h = ks::hpi(x = log(sim$BF.forecast)) * 2)
    abovethresh <- round(sum(sim$BF.forecast >= 10) / length(sim$BF.forecast) * 100, 1)
    belowthresh <- round(sum(sim$BF.forecast <= 1/10) / length(sim$BF.forecast) * 100, 1)
    betweenthresh <- round(100-abovethresh-belowthresh, 1)
    xtext <-  101 + 0.5*scaleDens + max(BFforecastDens$estimate)/max(BFforecastDens$estimate) * scaleDens
    
  }
  # BF distribution
    
  scaleY <- pretty(ylimvals)[2]-pretty(ylimvals)[1]
    
  BFs_s1 <- sapply(c(1:stage), function(x) SIM[[x]]$BF.temp)
  ns_s1 <- sapply(c(1:stage), function(x) SIM[[x]]$n.temp)
  
  par(mar = c(5,7,2,9))
  plot(x = c(0, ns_s1), 
       y = c(0, log(BFs_s1)),
       type = "l",
       ylim = ylimvals,
       xlim = c(0, 100), 
       xlab = "Sample size per Group",
       bty="l",
       yaxt="n",
       ylab="",
       cex.lab = 1.5, cex.axis=1.5)
  axis(2, at=log(c(1/1000, 1/50, 1/10, 1, 10, 50, 1000)), labels = c("1/1000", "1/50", "1/10", "1", "10", "50", "1000"), las = 1, cex.axis=1.5)
  mtext("(Projected) Bayes factor", 2, line=5, cex=1.5)
  
  if(!is.null(sim$n.forecast)){
  # Plot fan
  polygon(x = c(sim$n.temp, 100, 100, sim$n.temp), 
          y = unname(log(c(sim$BF.temp, qnt_BFs_s2[1], qnt_BFs_s2[11], sim$BF.temp))),
          col = scales::alpha("grey12", 0.1), border = NA)
  polygon(x = c(sim$n.temp, 100, 100, sim$n.temp), 
          y = unname(log(c(sim$BF.temp, qnt_BFs_s2[2], qnt_BFs_s2[10], sim$BF.temp))),
          col = scales::alpha("grey12", 0.1), border = NA)
  polygon(x = c(sim$n.temp, 100, 100, sim$n.temp), 
          y = unname(log(c(sim$BF.temp, qnt_BFs_s2[3], qnt_BFs_s2[9], sim$BF.temp))),
          col = scales::alpha("grey12", 0.1), border = NA)
  polygon(x = c(sim$n.temp, 100, 100, sim$n.temp), 
          y = unname(log(c(sim$BF.temp, qnt_BFs_s2[4], qnt_BFs_s2[8], sim$BF.temp))),
          col = scales::alpha("grey12", 0.1), border = NA)
  polygon(x = c(sim$n.temp, 100, 100, sim$n.temp), 
          y = unname(log(c(sim$BF.temp, qnt_BFs_s2[5], qnt_BFs_s2[7], sim$BF.temp))),
          col = scales::alpha("grey12", 0.1), border = NA)
  segments(x0 = sim$n.temp, y0 = log(sim$BF.temp), x1 = 100, y1 = log(qnt_BFs_s2[6]),
           col = scales::alpha("grey12", 0.5))
  
  # Plot density
  polygon(x = c(rep(101 + 0.5*scaleDens, length(BFforecastDens$eval.points)), 
                101 + 0.5*scaleDens + rev(BFforecastDens$estimate/max(BFforecastDens$estimate) * scaleDens)), 
          y = c(BFforecastDens$eval.points, rev(BFforecastDens$eval.points)), 
          xpd = NA, col = scales::alpha("grey12", 0.2))
  text(x = xtext, y = log(10) + 0.5 * scaleY, label = paste(abovethresh, "%"), xpd=NA, adj = c(0,0), cex = 1.5)
  text(x = xtext, y = 0, label = paste(betweenthresh, "%"), xpd=NA, adj = c(0,0.5), cex = 1.5)
  text(x = xtext, y = log(1/10) - 0.5 * scaleY, label = paste(belowthresh, "%"), xpd=NA, adj = c(0,1), cex = 1.5)
  segments(x0 = 101 + 0.5*scaleDens, y0 = ylimvals[1], y1 = ylimvals[2], xpd=NA)
  }
  segments(x0 = 0, y0 = log(c(1/10, 10)), x1 = 101 + 2*scaleDens, lty = "dashed", xpd = NA)
  
}

# Plot the figures

load("SIM_stopfutil_01.RData")
simplefan(SIM_stopfutil_01, ID = 1, stage = 8) # Stop for futility
simplefan(SIM_stopfutil_01, ID = 10, stage = 10) # Stop at Max N
simplefan(SIM_stopfutil_01, ID = 16, stage = 8) # Stop at  upper bound
simplefan(SIM_stopfutil_01, ID = 29, stage = 9) # Stop at lower bound

