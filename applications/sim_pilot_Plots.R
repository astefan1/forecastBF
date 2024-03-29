# ==============================================================================
# Plots for internal pilot study design
# ==============================================================================

rm(list = ls())
library(ggplot2)
library(ggpattern)
library(magick)

################## Model-averaged predictions, zero-centered prior #############

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilot_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilot_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

boundary <- c(1/10,10)
power <- 0.8

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilot_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0)
  res$conclusiveBF <- (res$BFFinal < boundary[1]) | (res$BFFinal > boundary[2])
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  simresults <- rbind(simresults,res)
  
}

simresults$SSconclusiveFinal <- simresults$SSFinal
simresults$SSconclusiveFinal[!simresults$conclusiveBF] <- simresults$SSconclusiveFinal[!simresults$conclusiveBF]*-1
simresults$SSconclusiveFinal <- as.factor(simresults$SSconclusiveFinal)
simresults$ES.pop <- as.factor(simresults$ES.pop)

countsMatrix <- simplify2array(tapply(simresults$SSconclusiveFinal, simresults$ES.pop, table))
countsDF <- as.data.frame(matrix(NA,nrow=length(SIM_names)*length(levels(simresults$SSconclusiveFinal))))
countsDF$ES.pop <- as.numeric(rep(levels(simresults$ES.pop), each = length(levels(simresults$SSconclusiveFinal))))
countsDF$SSconclusiveFinal <- as.numeric(rep(levels(simresults$SSconclusiveFinal), times =length(levels(simresults$ES.pop))))
countsDF$counts <- as.vector(countsMatrix)*sign(countsDF$SSconclusiveFinal)
countsDF$ES.pop <- as.factor(countsDF$ES.pop)
countsDF$SSconclusiveFinal <- factor(countsDF$SSconclusiveFinal, levels = c("-200", "-150", "-100", "-50", "-20", "20", "50", "100", "150", "200"), ordered = TRUE)

pat <- c(rep("right45", 5), rep("left30", 5))
mycolors <- c("#BEFFFF", "#A0E8EE", "#3FA0FF","#264FFF", "#292ADA","#F6EA89", "#F3D06E", "#FFAC72", "#F76E5E",  "#A61620")

theme_set(theme_classic())
update_geom_defaults("text", list(colour = "grey20", family = theme_get()$text$family, size = theme_get()$text$size))

p <- ggplot(countsDF, aes(x = ES.pop, y = counts/10)) +
  geom_bar_pattern(aes(pattern_type = SSconclusiveFinal, fill = SSconclusiveFinal), pattern = "magick", position = "stack", stat = "identity", size=0.1, pattern_key_scale_factor = 0.7,
                   colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = mycolors,
                    labels = c("N=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |", "N=20 inconclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE,
                    breaks = c( "20", "50", "100", "150", "200","-20", "-50", "-100","-150", "-200")
  ) +
  scale_pattern_type_manual(values = pat,
                    labels = c("N=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |", "N=20 inconclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE,
                    breaks = c( "20", "50", "100", "150", "200","-20", "-50", "-100","-150", "-200"),
  ) +
  theme_classic() +
  annotate("segment", x=0, y=-100, xend=0, yend=100,
           arrow=arrow(length=unit(0.3, "cm"), ends = "both")) +
  theme(legend.title = element_blank(),
        text = element_text(size = 18),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(margin = margin(r =25)),
        ) +
  labs(x = "Population Effect Size",
       y = "") +
  scale_y_continuous(limits = c(-100, 100), labels = c("100", "50", "0", "50", "100")) +
  scale_x_discrete(expand = c(-0.5, 1)) +
  coord_cartesian(clip = "off") 

p <- p +
  geom_text(aes(label = label), data.frame(label="% conclusive BF"), x = -1.3, y = 50, angle = 90, adj = 0.5, size = 7) +
  geom_text(aes(label = label), data.frame(label="% inconclusive BF"), x = -1.3, y = -50, angle = 90, adj = 0.5, size = 7) +
  annotate("segment", x = -0.1, y = 0, xend = 0.1, yend = 0)
  
p

# Produce the color bars as a legend because patterned legend doesn't work

ES.pop <- rep(seq(0, 1, by = 0.1), each=10)
SSconclusiveFinal <- rep(c(-200, -150, -100, -50, -20, 20, 50, 100, 150, 200), 11)
counts <- rep(1, 110)
countsDF2 <- data.frame(ES.pop = ES.pop,
                        SSconclusiveFinal = as.factor(SSconclusiveFinal),
                        counts=counts)
ggplot(countsDF2, aes(x = ES.pop, y = counts/10)) +
  geom_bar_pattern(aes(pattern_type = SSconclusiveFinal, fill = SSconclusiveFinal), pattern = "magick", position = "stack", stat = "identity", size=0.1, pattern_key_scale_factor = 0.7,
                   colour = 'black') +
  scale_fill_manual(values = mycolors) +
  scale_pattern_type_manual(values = pat) +
  theme_classic() 


################ Predictions under H1, zero-centered prior (Appendix) ##########

rm(list = ls())

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilotH1_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilotH1_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

boundary <- c(1/6, 6)
power <- 0.6
whichThreshold <- "both"

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilotH1_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0, whichThreshold = whichThreshold)
  res$conclusiveBF <- checkThreshold(res$BFFinal, boundary, which= whichThreshold)
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  simresults <- rbind(simresults,res)
  
}

simresults$SSconclusiveFinal <- simresults$SSFinal
simresults$SSconclusiveFinal[!simresults$conclusiveBF] <- simresults$SSconclusiveFinal[!simresults$conclusiveBF]*-1
simresults$SSconclusiveFinal <- as.factor(simresults$SSconclusiveFinal)
simresults$ES.pop <- as.factor(simresults$ES.pop)

countsMatrix <- simplify2array(tapply(simresults$SSconclusiveFinal, simresults$ES.pop, table))
countsDF <- as.data.frame(matrix(NA,nrow=length(SIM_names)*length(levels(simresults$SSconclusiveFinal))))
countsDF$ES.pop <- as.numeric(rep(levels(simresults$ES.pop), each = length(levels(simresults$SSconclusiveFinal))))
countsDF$SSconclusiveFinal <- as.numeric(rep(levels(simresults$SSconclusiveFinal), times =length(levels(simresults$ES.pop))))
countsDF$counts <- as.vector(countsMatrix)*sign(countsDF$SSconclusiveFinal)
countsDF$ES.pop <- as.factor(countsDF$ES.pop)
countsDF$SSconclusiveFinal <- factor(countsDF$SSconclusiveFinal, levels = c("-200", "-150", "-100", "-50", "-20", "20", "50", "100", "150", "200"))


ggplot(countsDF, aes(fill = SSconclusiveFinal, x = ES.pop, y = counts/10)) +
  geom_bar(position = "stack", stat = "identity", size=0.1, color="black") +
  scale_fill_manual(values = c("#A61620", "#F76E5E", "#FFAC72", "#F3D06E", "#F6EA89", "#BEFFFF", "#A0E8EE", "#3FA0FF","#264FFF", "#292ADA"),
                    labels = c("N=200 inconclusive", "N=150       |", "N=100       |", "N=50         |", "N=20         |\n", "\nN=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 18)) +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(limits = c(-100, 100), labels = c("100", "50", "0", "50", "100"))

################# Model-averaged predictions, Vohs prior (Appendix) ############

rm(list = ls())

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilotInformed_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilotInformed_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

boundary <- c(1/6, 6)
power <- 0.8
whichThreshold <- "both"

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilotInformed_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0, whichThreshold = whichThreshold)
  res$conclusiveBF <- checkThreshold(res$BFFinal, boundary, which= whichThreshold)
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  simresults <- rbind(simresults,res)
  
}

simresults$SSconclusiveFinal <- simresults$SSFinal
simresults$SSconclusiveFinal[!simresults$conclusiveBF] <- simresults$SSconclusiveFinal[!simresults$conclusiveBF]*-1
simresults$SSconclusiveFinal <- as.factor(simresults$SSconclusiveFinal)
simresults$ES.pop <- as.factor(simresults$ES.pop)

countsMatrix <- simplify2array(tapply(simresults$SSconclusiveFinal, simresults$ES.pop, table))
countsDF <- as.data.frame(matrix(NA,nrow=length(SIM_names)*length(levels(simresults$SSconclusiveFinal))))
countsDF$ES.pop <- as.numeric(rep(levels(simresults$ES.pop), each = length(levels(simresults$SSconclusiveFinal))))
countsDF$SSconclusiveFinal <- as.numeric(rep(levels(simresults$SSconclusiveFinal), times =length(levels(simresults$ES.pop))))
countsDF$counts <- as.vector(countsMatrix)*sign(countsDF$SSconclusiveFinal)
countsDF$ES.pop <- as.factor(countsDF$ES.pop)
countsDF$SSconclusiveFinal <- factor(countsDF$SSconclusiveFinal, levels = c("-200", "-150", "-100", "-50", "-20", "20", "50", "100", "150", "200"))

ggplot(countsDF, aes(fill = SSconclusiveFinal, x = ES.pop, y = counts/10)) +
  geom_bar(position = "stack", stat = "identity", size=0.1, color="black") +
  scale_fill_manual(values = c("#A61620", "#F76E5E", "#FFAC72", "#F3D06E", "#F6EA89", "#BEFFFF", "#A0E8EE", "#3FA0FF","#264FFF", "#292ADA"),
                    labels = c("N=200 inconclusive", "N=150       |", "N=100       |", "N=50         |", "N=20         |\n", "\nN=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 18)) +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(limits = c(-100, 100), labels = c("100", "50", "0", "50", "100")) 

############### Predictions under H0, zero-centered prior (Appendix) ###########

rm(list = ls())

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilotH0_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilotH0_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

boundary <- 1/6
power <- 0.8
whichThreshold <- "lower"

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilotH0_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0, whichThreshold = whichThreshold)
  res$conclusiveBF <- checkThreshold(res$BFFinal, boundary, which= whichThreshold)
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  simresults <- rbind(simresults,res)
  
}

simresults$SSconclusiveFinal <- simresults$SSFinal
simresults$SSconclusiveFinal[!simresults$conclusiveBF] <- simresults$SSconclusiveFinal[!simresults$conclusiveBF]*-1
simresults$SSconclusiveFinal <- as.factor(simresults$SSconclusiveFinal)
simresults$ES.pop <- as.factor(simresults$ES.pop)

countsMatrix <- simplify2array(tapply(simresults$SSconclusiveFinal, simresults$ES.pop, table))
countsDF <- as.data.frame(matrix(NA,nrow=length(SIM_names)*length(levels(simresults$SSconclusiveFinal))))
countsDF$ES.pop <- as.numeric(rep(levels(simresults$ES.pop), each = length(levels(simresults$SSconclusiveFinal))))
countsDF$SSconclusiveFinal <- as.numeric(rep(levels(simresults$SSconclusiveFinal), times =length(levels(simresults$ES.pop))))
countsDF$counts <- as.vector(countsMatrix)*sign(countsDF$SSconclusiveFinal)
countsDF$ES.pop <- as.factor(countsDF$ES.pop)
countsDF$SSconclusiveFinal <- factor(countsDF$SSconclusiveFinal, levels = c("-200", "-150", "-100", "-50", "-20", "20", "50", "100", "150", "200"))

ggplot(countsDF, aes(fill = SSconclusiveFinal, x = ES.pop, y = counts/10)) +
  geom_bar(position = "stack", stat = "identity", size=0.1, color="black") +
  scale_fill_manual(values = c("#A61620", "#F76E5E", "#FFAC72", "#F3D06E", "#F6EA89", "#BEFFFF", "#A0E8EE", "#3FA0FF","#264FFF", "#292ADA"),
                    labels = c("N=200 inconclusive", "N=150       |", "N=100       |", "N=50         |", "N=20         |\n", "\nN=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 18)) +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(limits = c(-100, 100), labels = c("100", "50", "0", "50", "100")) 

############## Model-averaged predictions, Vohs prior (Appendix) ###############

rm(list = ls())

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilotInformedH1_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilotInformedH1_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

boundary <- 10
power <- 0.8
whichThreshold <- "upper"

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilotInformedH1_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0, whichThreshold = whichThreshold)
  res$conclusiveBF <- checkThreshold(res$BFFinal, boundary, which= whichThreshold)
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  simresults <- rbind(simresults,res)
  
}

simresults$SSconclusiveFinal <- simresults$SSFinal
simresults$SSconclusiveFinal[!simresults$conclusiveBF] <- simresults$SSconclusiveFinal[!simresults$conclusiveBF]*-1
simresults$SSconclusiveFinal <- as.factor(simresults$SSconclusiveFinal)
simresults$ES.pop <- as.factor(simresults$ES.pop)

countsMatrix <- simplify2array(tapply(simresults$SSconclusiveFinal, simresults$ES.pop, table))
countsDF <- as.data.frame(matrix(NA,nrow=length(SIM_names)*length(levels(simresults$SSconclusiveFinal))))
countsDF$ES.pop <- as.numeric(rep(levels(simresults$ES.pop), each = length(levels(simresults$SSconclusiveFinal))))
countsDF$SSconclusiveFinal <- as.numeric(rep(levels(simresults$SSconclusiveFinal), times =length(levels(simresults$ES.pop))))
countsDF$counts <- as.vector(countsMatrix)*sign(countsDF$SSconclusiveFinal)
countsDF$ES.pop <- as.factor(countsDF$ES.pop)
countsDF$SSconclusiveFinal <- factor(countsDF$SSconclusiveFinal, levels = c("-200", "-150", "-100", "-50", "-20", "20", "50", "100", "150", "200"))

ggplot(countsDF, aes(fill = SSconclusiveFinal, x = ES.pop, y = counts/10)) +
  geom_bar(position = "stack", stat = "identity", size=0.1, color="black") +
  scale_fill_manual(values = c("#A61620", "#F76E5E", "#FFAC72", "#F3D06E", "#F6EA89", "#BEFFFF", "#A0E8EE", "#3FA0FF","#264FFF", "#292ADA"),
                    labels = c("N=200 inconclusive", "N=150       |", "N=100       |", "N=50         |", "N=20         |\n", "\nN=20 conclusive", "N=50         |", "N=100       |", "N=150       |", "N=200       |"),
                    drop = FALSE) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 18)) +
  labs(x = "Population Effect Size",
       y = "% Ending At") +
  scale_y_continuous(limits = c(-100, 100), labels = c("100", "50", "0", "50", "100")) 

########### Plots for internal pilot study design: Meta-analytic Bias ##########

#### Report everything ####

rm(list = ls())

source("applications/sim_pilot_Analyze.R")

SIM_names <- dir(".", pattern = "SIM_pilot_")
lapply(SIM_names, load, .GlobalEnv)

SIM_numbers <- gsub("SIM_pilot_", "", SIM_names)
SIM_numbers <- substr(SIM_numbers, 1, nchar(SIM_numbers)-6)

# Change the following parameters for different thresholds and desired power
boundary <- c(1/6,6)
power <- 0.6

simresults <- as.data.frame(matrix(NA, ncol=7, nrow=0))
meta <- as.data.frame(matrix(NA, ncol = 3, nrow = 0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilot_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0)
  res$conclusiveBF <- (res$BFFinal < boundary[1]) | (res$BFFinal > boundary[2])
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  resMeta <- metaPilot(res)
  
  simresults <- rbind(simresults,res)
  meta <- rbind(meta, cbind(coef(resMeta), resMeta$ci.lb, resMeta$ci.ub))
  
}

colnames(meta) <- c("Estimate", "CI_low", "CI_high")

ES.pop <- as.numeric(SIM_numbers)/10

#### Extract only the results from trajectories where sampling continues ####

meta2 <- as.data.frame(matrix(NA, ncol = 4, nrow = 0))

for(i in 1:(length(SIM_names))){
  
  res <- analyzePilot(get(paste0("SIM_pilot_", SIM_numbers[i])), boundary = boundary, power = power, abandon = TRUE, nmin = 0)
  res$conclusiveBF <- (res$BFFinal < boundary[1]) | (res$BFFinal > boundary[2])
  res$ES.pop <- as.numeric(SIM_numbers[i])/10
  
  resMeta <- metaPilot(res[(res$maxPower > power), ])
  meta2 <- rbind(meta2, cbind(coef(resMeta), resMeta$ci.lb, resMeta$ci.ub, resMeta$k.all))
}

colnames(meta2) <- c("Estimate", "CI_low", "CI_high", "N_Studies")

# Plot estimates for reporting everything
plot(ES.pop, meta$Estimate, pch = 19, bty="l", xlab = "Population Effect Size", 
     ylab = "Meta-Analytic Estimate", cex.lab = 1.5, cex.axis = 1.5, cex=1.5, ylim = c(-0.2,1.1), las=1)
abline(a=0, b=1, col = "grey")

# Add estimates for selective reporting
points(ES.pop, meta2$Estimate, pch = 19, bty="l", xlab = "Population Effect Size", ylab = "Meta-Analytic Estimate", cex.lab = 1.5, cex.axis = 1.5, cex=1.5, ylim = c(-0.2, 1), col = "red")

