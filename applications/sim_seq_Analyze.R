# ==============================================================================
# Evaluation of simulation results: sequential design with maxN (no futility 
# stopping)
# ==============================================================================

############################# Data wrangling ###################################

# Extract simulation results from directory
SIM_names <- dir(".", pattern = "SIM_seqMaxN_")
SIM_numbers <- gsub("SIM_seqMaxN_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

# Load all simulation results
lapply(SIM_names, load, .GlobalEnv)

# Initiate results matrix
simresults <- matrix(NA, nrow = length(unique(SIM_seqMaxN_00$sim$id)), 
                     ncol = length(SIM_names)*3)
colnames(simresults) <- paste0(c("finalBF_", "finalN_", "reasonEnd_"), 
                               rep(SIM_numbers, each = 3))
simresults <- as.data.frame(simresults)

for(i in 1:(length(SIM_names))){
  
  sim <- get(paste0("SIM_seqMaxN_", SIM_numbers[i]))$sim
  SIM_bounds <- get(paste0("SIM_seqMaxN_", SIM_numbers[i]))$settings$boundary
  SIM_n.max <- get(paste0("SIM_seqMaxN_", SIM_numbers[i]))$settings$n.max
  
  simresults[ , paste0("finalBF_", SIM_numbers[i])] <- unname(tapply(exp(sim$logBF), as.factor(sim$id), tail, n = 1))
  simresults[ , paste0("finalN_", SIM_numbers[i])] <- unname(tapply(sim$n, as.factor(sim$id), tail, n = 1))
  for(j in 1:nrow(simresults)){
    if(simresults[j, paste0("finalBF_", SIM_numbers[i])] < SIM_bounds[1]){
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "H0"
    } else if (simresults[j, paste0("finalBF_", SIM_numbers[i])] > SIM_bounds[2]) {
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "H1"
    } else if (simresults[j, paste0("finalN_", SIM_numbers[i])] == SIM_n.max){
      simresults[j, paste0("reasonEnd_", SIM_numbers[i])] <- "maxN"
    }
  }
}

simresults_long <- data.frame(finalBF = unlist(simresults[, grep("finalBF", colnames(simresults))]),
                              finalN = unlist(simresults[, grep("finalN", colnames(simresults))]),
                              reasonEnd = unlist(simresults[, grep("reasonEnd", colnames(simresults))]),
                              ES.pop = rep(seq(0, 1, by = 0.1), each = 1000))

simresults_long$reasonEnd <- factor(simresults_long$reasonEnd, levels = c("H0", "maxN", "futility", "H1"))
simresults_long$ES.pop <- as.factor(simresults_long$ES.pop)
simresults_long$design <- "seq"

SIM_seqMaxN_summary <- simresults_long
save(SIM_seqMaxN_summary, file = "SIM_summary_seqMaxN.RData")

####### Stacked barplot: Design Analysis (all population effect sizes) #########

library(ggplot2)

ggplot(simresults_long, aes(x = ES.pop, fill = reasonEnd)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c("#CC6677", "darkgrey", "#7AA6DD")) +
  theme(text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "top") +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(breaks = seq(0, 1000, by = 250), labels = seq(0, 100, by = 25))

