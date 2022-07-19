# ==============================================================================
# Evaluation of simulation results: Stopping for futility
# ==============================================================================

############################# Data wrangling ###################################

rm(list = ls())

# Extract simulation results from directory
SIM_names <- dir(".", pattern = "SIM_stopfutil_")
SIM_numbers <- gsub("SIM_stopfutil_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

# Load all simulation results
lapply(SIM_names, load, .GlobalEnv) # This may take a few seconds

# Initiate results matrix
simresults <- matrix(NA, nrow = length(SIM_stopfutil_00[[1]]), 
                     ncol = length(SIM_names)*3)
colnames(simresults) <- paste0(c("finalBF_", "finalN_", "reasonEnd_"), 
                               rep(SIM_numbers, each = 3))
simresults <- as.data.frame(simresults)

for(i in 1:(length(SIM_names))){
  
  sim <- get(paste0("SIM_stopfutil_", SIM_numbers[i]))[[1]]
  SIM_bounds <- get(paste0("SIM_stopfutil_", SIM_numbers[i]))[[2]]$boundary
  SIM_n.max <- get(paste0("SIM_stopfutil_", SIM_numbers[i]))[[2]]$n.max
  for(j in 1:length(sim)){
    simresults[j, paste0("finalBF_", SIM_numbers[i])] <- sim[[j]][[length(sim[[j]])]]$BF.temp
    simresults[j, paste0("finalN_", SIM_numbers[i])] <- sim[[j]][[length(sim[[j]])]]$n.temp
    if(simresults[j, paste0("finalBF_", SIM_numbers[i])] < SIM_bounds[1]){
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "H0"
    } else if (simresults[j, paste0("finalBF_", SIM_numbers[i])] > SIM_bounds[2]) {
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "H1"
    } else if (simresults[j, paste0("finalN_", SIM_numbers[i])] < SIM_n.max){
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "futility"
    } else if (simresults[j, paste0("finalN_", SIM_numbers[i])] == SIM_n.max){
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "maxN"
    }
  }
}

####### Stacked barplot: Design Analysis (all population effect sizes) #########

simresults_long <- data.frame(finalBF = unlist(simresults[, grep("finalBF", colnames(simresults))]),
                              finalN = unlist(simresults[, grep("finalN", colnames(simresults))]),
                              reasonEnd = unlist(simresults[, grep("reasonEnd", colnames(simresults))]),
                              ES.pop = rep(seq(0, 1, by = 0.1), each = 1000))

simresults_long$reasonEnd <- factor(simresults_long$reasonEnd, levels = c("H0", "maxN", "futility", "H1"))
simresults_long$ES.pop <- as.factor(simresults_long$ES.pop)
simresults_long$design <- "futil"

SIM_stopfutil_summary <- simresults_long
save(SIM_stopfutil_summary, file = "SIM_summary_stopfutil.RData")


library(ggplot2)

ggplot(simresults_long, aes(x = ES.pop, fill = reasonEnd)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#CC6677", "darkgrey", "lightgrey", "#7AA6DD"), labels = c("H0", "max-N", "futility", "H1")) +
  theme(text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "top") +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(breaks = seq(0, 1000, by = 250), labels = seq(0, 100, by = 25)) 

############# Stacked barplots: Distribution of sample sizes ###################

simresults[, grep("reasonEnd_", colnames(simresults))] <- lapply(simresults[, grep("reasonEnd_", colnames(simresults))],
                                                                 factor, levels = c("H0", "maxN", "futility", "H1"))

simresults[, "reasonEnd_00"] <- factor(simresults[, "reasonEnd_00"], levels = c("H1", "futility", "maxN", "H0"))

ggplot(simresults, aes(x = as.factor(finalN_00), fill = as.factor(reasonEnd_00))) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#7AA6DD", "lightgrey", "darkgrey", "#CC6677"), labels = c("H1", "futility", "max-N", "H0")) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

ggplot(simresults, aes(x = as.factor(finalN_03), fill = as.factor(reasonEnd_03))) +
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

ggplot(simresults, aes(x = as.factor(finalN_05), fill = as.factor(reasonEnd_05))) +
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

##################### Bias in Futility Stopping ################################

source("ttest_2sample_normalprior/forecast_2sample_t.R")

tval <- rep(NA, nrow(simresults_long))
for(i in 1:nrow(simresults_long)){ # figure out t-values that belong to BFs
  tval[i] <- optimize(function(x) (BFplus0_norm(x, 
                                                n1=simresults_long$finalN[i], 
                                                n2=simresults_long$finalN[i], 
                                                prior.mu = 0, prior.var = 1)-simresults_long$finalBF[i]
                                   )^2,
                      interval = c(-10, 10))$minimum
}

simresults_long$tval <- tval

simresults_long$dFinal <- simresults_long$tval*sqrt(2/simresults_long$finalN)
simresults_long$SEDFinal <- sqrt((2*simresults_long$finalN)/(simresults_long$finalN^2) + (simresults_long$dFinal^2) / (4*simresults_long$finalN))

ES.pop <- seq(0,1,by=0.1)
ES.noBias <- rep(NA, length(ES.pop))
ES.Bias <- rep(NA, length(ES.pop))
for(i in 1:length(ES.pop)){
  
  sim <- simresults_long[simresults_long$ES.pop == ES.pop[i],]
  ES.noBias[i] <- coef(metafor::rma(yi = sim$dFinal, 
                               sei = sim$SEDFinal, 
                               method = "FE"))
  ES.Bias[i] <- coef(metafor::rma(yi = sim[sim$reasonEnd %in% c("H0", "H1"),]$dFinal, 
                             sei = sim[sim$reasonEnd %in% c("H0", "H1"),]$SEDFinal, 
                             method = "FE"))
  
}

# Plot estimates for reporting everything
par(mar=c(5,5.5,2,2))
plot(ES.pop, ES.noBias, pch = 19, bty="l", xlab = "Population Effect Size", 
     ylab = "", cex.lab = 2, cex.axis = 1.5, cex=1.5, 
     ylim = c(-0.2,1.1), las=1)
mtext("Meta-Analytic Estimate", 2, cex=2, line=3.5)

abline(a=0, b=1, col = "grey")

# Add estimates for selective reporting
points(ES.pop, ES.Bias, pch = 19, bty="l", xlab = "Population Effect Size", ylab = "Meta-Analytic Estimate", cex.lab = 2, cex.axis = 1.5, cex=1.5, ylim = c(-0.2, 1), col = "red")

