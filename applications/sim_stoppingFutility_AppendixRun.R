# ==============================================================================
# Stopping for futility: Additional analyses to explore
# ==============================================================================
source("applications/sim_stoppingFutility_Functions.R")

################################ Small N #######################################

SIM_stopfutil_00_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0,
                                                    alternative = "greater")
save(SIM_stopfutil_00_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_00_smallN.RData")

SIM_stopfutil_01_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.1,
                                                    alternative = "greater")
save(SIM_stopfutil_01_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_01_smallN.RData")

SIM_stopfutil_02_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.2,
                                                    alternative = "greater")
save(SIM_stopfutil_02_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_02_smallN.RData")

SIM_stopfutil_03_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.3,
                                                    alternative = "greater")
save(SIM_stopfutil_03_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_03_smallN.RData")

SIM_stopfutil_04_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.4,
                                                    alternative = "greater")
save(SIM_stopfutil_04_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_04_smallN.RData")

SIM_stopfutil_05_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.5,
                                                    alternative = "greater")
save(SIM_stopfutil_05_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_05_smallN.RData")

SIM_stopfutil_06_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.6,
                                                    alternative = "greater")
save(SIM_stopfutil_06_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_06_smallN.RData")

SIM_stopfutil_07_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.7,
                                                    alternative = "greater")
save(SIM_stopfutil_07_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_07_smallN.RData")

SIM_stopfutil_08_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.8,
                                                    alternative = "greater")
save(SIM_stopfutil_08_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_08_smallN.RData")

SIM_stopfutil_09_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 0.9,
                                                    alternative = "greater")
save(SIM_stopfutil_09_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_09_smallN.RData")

SIM_stopfutil_10_smallN <- designAnalysis.stopfutil(n.min = 2, 
                                                    n.max = 10, 
                                                    stepsize = 1, 
                                                    ES.pop = 1,
                                                    alternative = "greater")
save(SIM_stopfutil_10_smallN, file = "./SIM_stopfutil_smallN/SIM_stopfutil_10_smallN.RData")

#################### Different Prior ###########################################

SIM_stopfutil_00_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_00_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_00_differentPrior.RData")

SIM_stopfutil_01_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.1,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_01_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_01_differentPrior.RData")

SIM_stopfutil_02_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.2,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_02_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_02_differentPrior.RData")

SIM_stopfutil_03_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.3,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_03_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_03_differentPrior.RData")

SIM_stopfutil_04_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.4,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_04_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_04_differentPrior.RData")

SIM_stopfutil_05_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.5,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_05_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_05_differentPrior.RData")

SIM_stopfutil_06_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.6,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_06_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_06_differentPrior.RData")

SIM_stopfutil_07_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.7,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_07_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_07_differentPrior.RData")

SIM_stopfutil_08_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.8,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_08_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_08_differentPrior.RData")

SIM_stopfutil_09_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 0.9,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_09_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_09_differentPrior.RData")

SIM_stopfutil_10_differentPrior <- designAnalysis.stopfutil(n.min = 10, 
                                                            n.max = 100, 
                                                            stepsize = 10, 
                                                            ES.pop = 1,
                                                            alternative = "greater", 
                                                            prior.mu = 0.3,
                                                            prior.var = 0.01)

save(SIM_stopfutil_10_differentPrior, file = "./SIM_stopfutil_differentPrior/SIM_stopfutil_10_differentPrior.RData")

#################### Different Futility Threshold ##############################

SIM_stopfutil_00_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)
                                             
save(SIM_stopfutil_00_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_00_highthreshold.RData")

SIM_stopfutil_01_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.1,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_01_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_01_highthreshold.RData")

SIM_stopfutil_02_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.2,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_02_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_02_highthreshold.RData")

SIM_stopfutil_03_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.3,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_03_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_03_highthreshold.RData")

SIM_stopfutil_04_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.4,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_04_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_04_highthreshold.RData")

SIM_stopfutil_05_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.5,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_05_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_05_highthreshold.RData")

SIM_stopfutil_06_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.6,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_06_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_06_highthreshold.RData")

SIM_stopfutil_07_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.7,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_07_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_07_highthreshold.RData")

SIM_stopfutil_08_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.8,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_08_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_08_highthreshold.RData")

SIM_stopfutil_09_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 0.9,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_09_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_09_highthreshold.RData")

SIM_stopfutil_10_highthreshold <- designAnalysis.stopfutil(n.min = 10,
                                                           n.max = 100, 
                                                           stepsize = 10, 
                                                           ES.pop = 1,
                                                           alternative = "greater",
                                                           futilitythreshold = 0.1)

save(SIM_stopfutil_10_highthreshold, file = "./SIM_stopfutil_highthreshold/SIM_stopfutil_10_highthreshold.RData")
