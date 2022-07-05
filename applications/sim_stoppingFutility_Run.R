# Design analysis of the futility stopping design with different effect sizes

source("ttest_2sample_normalprior/mcmc_2sample_t.R")

SIM_stopfutil_00 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0,
                                             alternative = "greater")
save(SIM_stopfutil_00, file = "SIM_stopfutil_00.RData")

SIM_stopfutil_01 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.1,
                                             alternative = "greater")
save(SIM_stopfutil_01, file = "SIM_stopfutil_01.RData")

SIM_stopfutil_02 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.2,
                                             alternative = "greater")
save(SIM_stopfutil_02, file = "SIM_stopfutil_02.RData")

SIM_stopfutil_03 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.3,
                                             alternative = "greater")
save(SIM_stopfutil_03, file = "SIM_stopfutil_03.RData")

SIM_stopfutil_04 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.4,
                                             alternative = "greater")
save(SIM_stopfutil_04, file = "SIM_stopfutil_04.RData")

SIM_stopfutil_05 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.5,
                                             alternative = "greater")
save(SIM_stopfutil_05, file = "SIM_stopfutil_05.RData")

SIM_stopfutil_06 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.6,
                                             alternative = "greater")
save(SIM_stopfutil_06, file = "SIM_stopfutil_06.RData")

SIM_stopfutil_07 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.7,
                                             alternative = "greater")
save(SIM_stopfutil_07, file = "SIM_stopfutil_07.RData")

SIM_stopfutil_08 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.8,
                                             alternative = "greater")
save(SIM_stopfutil_08, file = "SIM_stopfutil_08.RData")

SIM_stopfutil_09 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 0.9,
                                             alternative = "greater")
save(SIM_stopfutil_09, file = "SIM_stopfutil_09.RData")

SIM_stopfutil_10 <- designAnalysis.stopfutil(n.min = 10, 
                                             n.max = 100, 
                                             stepsize = 10, 
                                             ES.pop = 1.0,
                                             alternative = "greater")
save(SIM_stopfutil_10, file = "SIM_stopfutil_10.RData")





