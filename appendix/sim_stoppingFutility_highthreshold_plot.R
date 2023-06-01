# ==============================================================================
# Evaluation of simulation results: Stopping for futility with high threshold
# ==============================================================================

############################# Data wrangling ###################################

rm(list = ls())

# Extract simulation results from directory
SIM_names <- dir(".", pattern = "SIM_stopfutil_")
SIM_numbers <- gsub("SIM_stopfutil_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-20)

# Load all simulation results
lapply(SIM_names, load, .GlobalEnv) # This may take a few seconds

# Initiate results matrix
simresults <- matrix(NA, nrow = length(SIM_stopfutil_00_highthreshold[[1]]), 
                     ncol = length(SIM_names)*3)
colnames(simresults) <- paste0(c("finalBF_", "finalN_", "reasonEnd_"), 
                               rep(SIM_numbers, each = 3))
simresults <- as.data.frame(simresults)

for(i in 1:(length(SIM_names))){
  
  sim <- get(paste0("SIM_stopfutil_", SIM_numbers[i], "_highthreshold"))[[1]]
  SIM_bounds <- get(paste0("SIM_stopfutil_", SIM_numbers[i], "_highthreshold"))[[2]]$boundary
  SIM_n.max <- get(paste0("SIM_stopfutil_", SIM_numbers[i], "_highthreshold"))[[2]]$n.max
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
simresults[, "reasonEnd_03"] <- factor(simresults[, "reasonEnd_03"], levels = rev(c("H1", "futility", "maxN", "H0")), ordered = TRUE)
simresults[, "reasonEnd_05"] <- factor(simresults[, "reasonEnd_05"], levels = rev(c("H1", "futility", "maxN", "H0")), ordered = TRUE)


ggplot(simresults, aes(x = as.factor(finalN_00), fill = factor(reasonEnd_00))) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#7AA6DD", "lightgrey", "darkgrey", "#CC6677"), 
                    labels = c("H1", "futility", "max-N", "H0")) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

ggplot(simresults, aes(x = as.factor(finalN_03), fill = reasonEnd_03)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "#CC6677", "#7AA6DD"), 
                    #labels = c("H1", "futility", "max-N", "H0"),
                    breaks = c("futility", "maxN", "H0", "H1")) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

ggplot(simresults, aes(x = as.factor(finalN_05), fill = reasonEnd_05)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("lightgrey", "darkgrey", "#CC6677", "#7AA6DD"), 
                    #labels = c("H1", "futility", "max-N", "H0"),
                    breaks = c("futility", "maxN", "H0", "H1"),
                    drop = FALSE) +
  theme(text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = "top") +
  labs(x = "Sample Size per Group",
       y = "% Ending At") +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50), labels = seq(0, 25, by = 5)) 

