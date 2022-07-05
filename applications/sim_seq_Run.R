# ==============================================================================
# Script to conduct design analyses for sequential design application (without
# futility stopping) using the BFDA package
# ==============================================================================

# Simulations need several hours to run. Results can be found on https://osf.io/xaqh4/
# To run the analysis script (sim_seq_Analyze.R), results need to
# be saved in the R-project folder.

library(BFDA)

SIM_seqMaxN_00 <- BFDA.sim(expected.ES = 0, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_00, file = "SIM_seqMaxN_00.RData")                                    

SIM_seqMaxN_01 <- BFDA.sim(expected.ES = 0.1, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_01, file = "SIM_seqMaxN_01.RData")                                    

SIM_seqMaxN_02 <- BFDA.sim(expected.ES = 0.2, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_02, file = "SIM_seqMaxN_02.RData")   

SIM_seqMaxN_03 <- BFDA.sim(expected.ES = 0.3, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_03, file = "SIM_seqMaxN_03.RData") 

SIM_seqMaxN_04 <- BFDA.sim(expected.ES = 0.4, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_04, file = "SIM_seqMaxN_04.RData") 

SIM_seqMaxN_05 <- BFDA.sim(expected.ES = 0.5, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_05, file = "SIM_seqMaxN_05.RData") 

SIM_seqMaxN_06 <- BFDA.sim(expected.ES = 0.6, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_06, file = "SIM_seqMaxN_06.RData") 

SIM_seqMaxN_07 <- BFDA.sim(expected.ES = 0.7, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_07, file = "SIM_seqMaxN_07.RData") 

SIM_seqMaxN_08 <- BFDA.sim(expected.ES = 0.8, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_08, file = "SIM_seqMaxN_08.RData") 

SIM_seqMaxN_09 <- BFDA.sim(expected.ES = 0.9, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_09, file = "SIM_seqMaxN_09.RData") 

SIM_seqMaxN_10 <- BFDA.sim(expected.ES = 1, type = "t.between", 
                           prior = list("normal",
                                        list(prior.mean = 0, prior.variance = 1)),
                           n.min = 10, n.max = 100, design = "sequential",
                           boundary = c(1/10,10), B = 1000, stepsize = 10, 
                           alternative = "greater", cores = 4, seed = 1234)
save(SIM_seqMaxN_10, file = "SIM_seqMaxN_10.RData") 
