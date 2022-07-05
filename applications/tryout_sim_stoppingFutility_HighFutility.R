SIM_stopfutil_HIGHFUTIL <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.3,
                                             alternative = "greater",
                                             futilitythreshold = 0.2)

simresults <- matrix(NA, nrow = 1000, ncol = 3)
colnames(simresults) <- c("finalBF", "finalN", "reasonEnd")
simresults <- as.data.frame(simresults)
sim <- SIM_stopfutil_HIGHFUTIL[[1]]

for(j in 1:1000){
  simresults[j, "finalBF"] <- sim[[j]][[length(sim[[j]])]]$BF.temp
  simresults[j, "finalN"] <- sim[[j]][[length(sim[[j]])]]$n.temp
  if(simresults[j, "finalBF"] < 1/10){
    simresults[j, "reasonEnd"] <- "H0"
  } else if (simresults[j, "finalBF"] > 10) {
    simresults[j, "reasonEnd"] <- "H1"
  } else if (simresults[j, "finalN"] < 100){
    simresults[j, "reasonEnd"] <- "futility"
  } else if (simresults[j, "finalN"] == 100){
    simresults[j, "reasonEnd"] <- "maxN"
  }
}

simresults$reasonEnd <- factor(simresults$reasonEnd, levels = c("H0", "maxN", "futility", "H1"))
load("SIM_summary_stopfutil.RData")
SIM03 <- SIM_summary[SIM_summary$ES.pop == 0.3 & SIM_summary$design == "futil",]

par(mfrow=c(1,2))
pie(table(SIM03$reasonEnd))
title("futility threshold = 1%")

pie(table(simresults$reasonEnd))
title("futility threshold = 20%")

library(ggplot2)
ggplot(SIM03, aes(x = as.factor(finalN), fill = as.factor(reasonEnd))) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#CC6677", "darkgrey", "lightgrey", "#7AA6DD")) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

ggplot(simresults, aes(x = as.factor(finalN), fill = as.factor(reasonEnd))) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#CC6677", "darkgrey", "lightgrey", "#7AA6DD")) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

