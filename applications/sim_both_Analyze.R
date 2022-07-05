# ==============================================================================
# Evaluation of simulation results: Comparing futility stopping to simple maxN
# ==============================================================================

rm(list = ls())

############################# Data wrangling ###################################

load("SIM_summary_seqMaxN.RData")
load("SIM_summary_stopfutil.RData")

# Collate data sets
SIM_summary <- rbind(SIM_seqMaxN_summary, SIM_stopfutil_summary)

# Get matrix with mean sample sizes
meanN <- as.data.frame(tapply(SIM_summary$finalN, SIM_summary[, c("ES.pop", "design")], mean))
meanN$ES.pop <- unique(SIM_summary$ES.pop)
meanN_long <- reshape2::melt(data = meanN, id.vars = "ES.pop", measure.vars = c("futil", "seq"))

############### Plot: Comparison mean sample sizes #############################

library(ggplot2)

ggplot(data = meanN_long, aes(x = ES.pop, group = variable, y = value)) +
  geom_line(aes(lty = variable)) +
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.8, 0.9)) +
  labs(x = "Population Effect Size",
       y = "Expected Sample Size") +
  scale_linetype(labels = c("Stopping for futility", "sequential MaxN"))
