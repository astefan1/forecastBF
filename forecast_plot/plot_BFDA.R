library(ggplot2)
source("ttest_2sample_normalprior/bfda_fixed_2sample_t.R")

set.seed(42)
bfda1 <- BFDA_2sample_t(0.5, alternative="two.sided", group.n = 50, prior.mu = 0, prior.var = 0.5, iter = 100000)

dat <- data.frame(bfda1)
dat$logBF <- log(dat$BF)

credInt <- unname(quantile(dat$logBF, probs = c(0.05, 0.95)))
densBF <- density(dat$logBF)

plot(densBF$x, densBF$y, type="l", bty="n",  yaxt="n", xaxt="n", xlab="", ylab="", xlim = log(c(1/100, 10000)))
polygon(x = c(densBF$x[densBF$x > credInt[1] & densBF$x < credInt[2]], rev(credInt)),
        y = c(densBF$y[densBF$x > credInt[1] & densBF$x < credInt[2]], 0, 0),
        col = "grey90", border=NA)
polygon(x = c(densBF$x[densBF$x > log(10)], max(densBF$x), log(10)),
        y = c(densBF$y[densBF$x > log(10)], 0, 0),
        col = "#6BD7AF", border="#6BD7AF", density = 10, angle = 45)
segments(x0=densBF$x[which.max(densBF$y)], y0=0, y1=100, lwd=2.5, col="darkgrey", lty="dotdash")
points(densBF$x, densBF$y, type="l", lwd=1.5)
axis(1, at = log(c(1/100, 1/10, 1/3, 1, 3, 10, 100, 1000, 10000)),
     labels = c("1/100", "1/10", "1/3", "1", "3", "10", "100", "1000", "10000"),
     line = -0.8, cex.axis = 1.5, lwd=0, lwd.ticks = 1)
abline(h=0)

plot.new()
legend(legend="", x=0.5, y=0.5, lty = "dotdash", col = "darkgrey", bty="n", cex=3, lwd=3)

plot.new()
legend(legend="", x=0.5, y=0.5, fill = "grey90", bty="n", cex=3, border=NA)

plot.new()
legend(legend="", x=0.5, y=0.5, bty="n", cex=3, border="#6BD7AF", density = 10, angle = 45, fill = "#6BD7AF")
