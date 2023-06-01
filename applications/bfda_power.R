# ==============================================================================
# Regular BFDA: Probability of obtaining strong evidence with specific (fixed)
# sample sizes  -- Overview for coparison & explanation pilot study design
# ==============================================================================

rm(list = ls())

############################ Conduct BFDAs #####################################

source("ttest_2sample_normalprior/bfda_fixed_2sample_t.R")

groupN <- seq(20, 200, by=10)
ES <- seq(0, 1, by=0.1)

res <- array(NA, dim = c(1000, length(groupN), length(ES)))

for(i in seq_along(groupN)){
  for(j in seq_along(ES)){
    tmp <- BFDA_2sample_t(pop.ES = ES[j], alternative = "greater", group.n = groupN[i], prior.mu = 0, prior.var = 1)
    res[, i, j] <- tmp$BF
  }
  print(i)
}

power <- array(dim=c(length(groupN), length(ES)))

for(i in seq_along(groupN)){
  for(j in seq_along(ES)){
    power[i, j] <- sum(res[, i, j] > 6 | res[, i, j] < 1/6)/1000
  }
}

rownames(power) <- groupN
colnames(power) <- ES

################## Plots: Is power > 80% for a specific sample size ############

library(reshape2)
powerDat <- reshape2::melt(as.data.frame(power))
powerDat$groupN <- rep(groupN, length(ES))

ggplot(data = powerDat, aes(x=variable, y = groupN, fill=value>0.6)) +
  geom_tile(color="black") +
  theme_minimal() +
  labs(x = "Effect Size",
       y = "Sample Size",
       fill = "Power > 0.6") +
  theme(text = element_text(size = 18)) +
  scale_fill_manual(values = c("#C25120", "#03B9AA"))

# Now with informed prior

res2 <- array(NA, dim = c(1000, length(groupN), length(ES)))

for(i in seq_along(groupN)){
  for(j in seq_along(ES)){
    tmp <- BFDA_2sample_t(pop.ES = ES[j], alternative = "greater", group.n = groupN[i], prior.mu = 0.3, prior.var = 0.0225)
    res2[, i, j] <- tmp$BF
  }
  print(i)
}

power2 <- array(dim=c(length(groupN), length(ES)))

for(i in seq_along(groupN)){
  for(j in seq_along(ES)){
    power2[i, j] <- sum(res2[, i, j] > 10 | res2[, i, j] < 1/10)/1000
  }
}

rownames(power2) <- groupN
colnames(power2) <- ES

library(reshape2)
powerDat <- reshape2::melt(as.data.frame(power2))
powerDat$groupN <- rep(groupN, length(ES))

ggplot(data = powerDat, aes(x=variable, y = groupN, fill=value>0.8)) +
  geom_tile(color="black") +
  theme_minimal() +
  labs(x = "Effect Size",
       y = "Sample Size",
       fill = "Power > 0.8") +
  theme(text = element_text(size = 18)) +
  scale_fill_manual(values = c("#C25120", "#03B9AA"))

  
