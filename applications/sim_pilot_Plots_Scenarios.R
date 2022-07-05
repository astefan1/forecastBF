# ==============================================================================
# Pilot study scenarios illustration
# ==============================================================================

load("SIM_pilot_05.RData")
load("SIM_pilot_00.RData")
source("applications/SIM_pilot_Analyze.R")
simplefan_pilot <- function(SIM_pilot, ID = 1, n.main = 30){
  
  sim <- SIM_pilot[[1]][[ID]]
  which.n.main <- which(sim$n.main == n.main)
  
  qnt_BFs_s2 <- apply(sim$BFs, 2, quantile, probs = c(0.05, seq(0.1, 0.9, 0.1), 0.95))
  
  # y axis limits
  ylimvals <- c(min(log(qnt_BFs_s2),
                    log(1/10)-1, na.rm = TRUE),
                max(log(qnt_BFs_s2),
                    log(10)+1, na.rm = TRUE))
  
  if(n.main == 0){
    ylimvals <- c(log(1/10)-1, log(10)+1)
  }
  
  scaleDens <- (pretty(c(0, 200))[2] - pretty(c(0, 200))[1]) * 0.5
  
  if(n.main > 0){
    BFforecastDens <- ks::kde.boundary(log(sim$BFs[, which.n.main]),
                                       xmin = min(log(sim$BFs[, which.n.main]), -0.3 * max(log(sim$BFs[, which.n.main]))), 
                                       xmax = max(log(sim$BFs[, which.n.main]), -0.3 * min(log(sim$BFs[, which.n.main]))), 
                                       h = ks::hpi(x = log(sim$BFs[, which.n.main])) * 2)  
    
    abovethresh <- round(sum(sim$BFs[, which.n.main] >= 10) / length(sim$BFs[, which.n.main]) * 100, 1)
    belowthresh <- round(sum(sim$BFs[, which.n.main] <= 1/10) / length(sim$BFs[, which.n.main]) * 100, 1)
    betweenthresh <- round(100-abovethresh-belowthresh, 1)
    xtext <-  201 + 0.5*scaleDens + max(BFforecastDens$estimate)/max(BFforecastDens$estimate) * scaleDens
  }
    # BF distribution
    
    scaleY <- pretty(ylimvals)[2]-pretty(ylimvals)[1]
    
    t_s1 <- sapply(c(2:sim$n.pilot), function(x) t.test(sim$group1_s1[1:x], sim$group2_s1[1:x], var.equal = TRUE)$statistic)
    BFs_s1 <- sapply(seq_along(t_s1), function(x) BFplus0_norm(tval=t_s1[x], n1=x+1, n2=x+1, prior.mu = sim$prior.mu, prior.var = sim$prior.var))
    ns_s1 <- c(2:sim$n.pilot)
    
    par(mar = c(8,9,2,12))
    plot(x = c(0, ns_s1), 
         y = c(0, log(BFs_s1)),
         type = "l",
         ylim = ylimvals,
         xlim = c(0, 200), 
         xlab = "",
         bty="l",
         yaxt="n",
         ylab="",
         cex.lab = 3, cex.axis=3, lwd=2, mgp=c(0,3,1))
    axis(2, at=log(c(1/1000, 1/50, 1/10, 1, 10, 50, 1000)), labels = c("1/1000", "1/50", "1/10", "1", "10", "50", "1000"), las = 1, cex.axis=3)
    mtext("(Projected) Bayes factor", 2, line=6.5, cex=3)
    mtext("Sample Size per Group", 1, line=5.7, cex=3)
    
    
  if(n.main > 0){
    # Plot fan
    polygon(x = c(sim$n.pilot, sim$n.pilot+n.main, sim$n.pilot+n.main, sim$n.pilot), 
            y = unname(log(c(sim$BF_s1, qnt_BFs_s2[1, which.n.main], qnt_BFs_s2[11, which.n.main], sim$BF_s1))),
            col = scales::alpha("grey12", 0.1), border = NA)
    polygon(x = c(sim$n.pilot, sim$n.pilot+n.main, sim$n.pilot+n.main, sim$n.pilot), 
            y = unname(log(c(sim$BF_s1, qnt_BFs_s2[2, which.n.main], qnt_BFs_s2[10, which.n.main], sim$BF_s1))),
            col = scales::alpha("grey12", 0.1), border = NA)
    polygon(x = c(sim$n.pilot, sim$n.pilot+n.main, sim$n.pilot+n.main, sim$n.pilot), 
            y = unname(log(c(sim$BF_s1, qnt_BFs_s2[3, which.n.main], qnt_BFs_s2[9, which.n.main], sim$BF_s1))),
            col = scales::alpha("grey12", 0.1), border = NA)
    polygon(x = c(sim$n.pilot, sim$n.pilot+n.main, sim$n.pilot+n.main, sim$n.pilot), 
            y = unname(log(c(sim$BF_s1, qnt_BFs_s2[4, which.n.main], qnt_BFs_s2[8, which.n.main], sim$BF_s1))),
            col = scales::alpha("grey12", 0.1), border = NA)
    polygon(x = c(sim$n.pilot, sim$n.pilot+n.main, sim$n.pilot+n.main, sim$n.pilot), 
            y = unname(log(c(sim$BF_s1, qnt_BFs_s2[5, which.n.main], qnt_BFs_s2[7, which.n.main], sim$BF_s1))),
            col = scales::alpha("grey12", 0.1), border = NA)
    segments(x0 = sim$n.pilot, y0 = log(sim$BF_s1), x1 = sim$n.pilot+n.main, y1 = log(qnt_BFs_s2[6, which.n.main]),
             col = scales::alpha("grey12", 0.5))
    
    # Plot density
    polygon(x = c(rep(201 + 0.5*scaleDens, length(BFforecastDens$eval.points)), 
                  201 + 0.5*scaleDens + rev(BFforecastDens$estimate/max(BFforecastDens$estimate) * scaleDens)), 
            y = c(BFforecastDens$eval.points, rev(BFforecastDens$eval.points)), 
            xpd = NA, col = scales::alpha("grey12", 0.2))
    text(x = xtext, y = log(10) + 0.5 * scaleY, label = paste(abovethresh, "%"), xpd=NA, adj = c(0,0), cex = 3)
    text(x = xtext, y = 0, label = paste(betweenthresh, "%"), xpd=NA, adj = c(0,0.5), cex = 3)
    text(x = xtext, y = log(1/10) - 0.5 * scaleY, label = paste(belowthresh, "%"), xpd=NA, adj = c(0,1), cex = 3)
    segments(x0 = 201 + 0.5*scaleDens, y0 = ylimvals[1], y1 = ylimvals[2], xpd=NA)
    
  }
    
  segments(x0 = 0, y0 = log(c(1/10, 10)), x1 = 201 + 2*scaleDens, lty = "dashed", xpd = NA)
  
}

simplefan_pilot(SIM_pilot_05, ID=43, n.main = 0)
simplefan_pilot(SIM_pilot_00, ID=122, n.main = 0)
simplefan_pilot(SIM_pilot_05, ID=43, n.main = 0)
simplefan_pilot(SIM_pilot_05, ID=931, n.main = 0)
simplefan_pilot(SIM_pilot_05, ID=931, n.main = 30)
simplefan_pilot(SIM_pilot_05, ID=931, n.main = 80)
simplefan_pilot(SIM_pilot_05, ID=931, n.main = 130)
simplefan_pilot(SIM_pilot_05, ID=931, n.main = 180)
simplefan_pilot(SIM_pilot_05, ID=3, n.main = 0)
simplefan_pilot(SIM_pilot_05, ID=3, n.main = 30)
simplefan_pilot(SIM_pilot_05, ID=3, n.main = 80)
simplefan_pilot(SIM_pilot_05, ID=3, n.main = 130)
simplefan_pilot(SIM_pilot_05, ID=3, n.main = 180)

