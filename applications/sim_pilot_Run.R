# Design analysis of the pilot study design with different effect sizes

source("applications/sim_pilot_Functions.R")

SIM_pilot_00 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_00, file="SIM_pilot_00.RData")                          

SIM_pilot_01 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.1, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_01, file="SIM_pilot_01.RData")                          

SIM_pilot_02 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.2, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_02, file="SIM_pilot_02.RData")                          

SIM_pilot_03 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.3, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_03, file="SIM_pilot_03.RData")                          

SIM_pilot_04 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.4, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_04, file="SIM_pilot_04.RData")                          

SIM_pilot_05 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.5, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_05, file="SIM_pilot_05.RData")                          

SIM_pilot_06 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.6, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_06, file="SIM_pilot_06.RData")

SIM_pilot_07 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.7, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_07, file="SIM_pilot_07.RData")

SIM_pilot_08 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.8, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_08, file="SIM_pilot_08.RData")  

SIM_pilot_09 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 0.9, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_09, file="SIM_pilot_09.RData")   

SIM_pilot_10 <- designAnalysis.pilot(n.pilot = 20, 
                                     n.main = c(30, 80, 130, 180), 
                                     ES.pop = 1, 
                                     alternative = "greater", 
                                     prior.mu=0, 
                                     prior.var=1, 
                                     iter = 1000, 
                                     iter.DA = 1000,
                                     mcmc.setting = list(chains=2, burnin=0, iter=5000))
save(SIM_pilot_10, file="SIM_pilot_10.RData")                          

SIM_pilotH1_00 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))
                                     
save(SIM_pilotH1_00, file="SIM_pilotH1_00.RData")                         

SIM_pilotH1_01 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.1, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_01, file="SIM_pilotH1_01.RData") 

SIM_pilotH1_02 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.2, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_02, file="SIM_pilotH1_02.RData") 

SIM_pilotH1_03 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.3, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_03, file="SIM_pilotH1_03.RData") 

SIM_pilotH1_04 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.4, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_04, file="SIM_pilotH1_04.RData") 

SIM_pilotH1_05 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.5, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_05, file="SIM_pilotH1_05.RData") 

SIM_pilotH1_06 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.6, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_06, file="SIM_pilotH1_06.RData")    

SIM_pilotH1_07 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.7, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_07, file="SIM_pilotH1_07.RData")     

SIM_pilotH1_08 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.8, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_08, file="SIM_pilotH1_08.RData")  

SIM_pilotH1_09 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.9, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_09, file="SIM_pilotH1_09.RData")   

SIM_pilotH1_10 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 1, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H1",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH1_10, file="SIM_pilotH1_10.RData")                          

SIM_pilotInformed_00 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))
                                     
save(SIM_pilotInformed_00, file="SIM_pilotInformed_00.RData")                          

SIM_pilotInformed_01 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.1, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_01, file="SIM_pilotInformed_01.RData")                          

SIM_pilotInformed_02 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.2, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_02, file="SIM_pilotInformed_02.RData")

SIM_pilotInformed_03 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.3, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_03, file="SIM_pilotInformed_03.RData") 

SIM_pilotInformed_04 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.4, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_04, file="SIM_pilotInformed_04.RData")  

SIM_pilotInformed_05 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.5, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_05, file="SIM_pilotInformed_05.RData")   

SIM_pilotInformed_06 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.6, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_06, file="SIM_pilotInformed_06.RData")   

SIM_pilotInformed_07 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.7, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_07, file="SIM_pilotInformed_07.RData")  

SIM_pilotInformed_08 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.8, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_08, file="SIM_pilotInformed_08.RData")   

SIM_pilotInformed_09 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0.9, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_09, file="SIM_pilotInformed_09.RData")   

SIM_pilotInformed_10 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 1, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformed_10, file="SIM_pilotInformed_10.RData")                          

SIM_pilotH0_00 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_00, file="SIM_pilotH0_00.RData")   

SIM_pilotH0_01 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.1, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_01, file="SIM_pilotH0_01.RData")

SIM_pilotH0_02 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.2, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_02, file="SIM_pilotH0_02.RData")   

SIM_pilotH0_03 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.3, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_03, file="SIM_pilotH0_03.RData")   

SIM_pilotH0_04 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.4, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_04, file="SIM_pilotH0_04.RData")   

SIM_pilotH0_05 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.5, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_05, file="SIM_pilotH0_05.RData") 

SIM_pilotH0_06 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.6, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_06, file="SIM_pilotH0_06.RData") 

SIM_pilotH0_07 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.7, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_07, file="SIM_pilotH0_07.RData") 

SIM_pilotH0_08 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.8, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_08, file="SIM_pilotH0_08.RData")   

SIM_pilotH0_09 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 0.9, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_09, file="SIM_pilotH0_09.RData")   

SIM_pilotH0_10 <- designAnalysis.pilot(n.pilot = 20,
                                       n.main = c(30, 80, 130, 180), 
                                       ES.pop = 1, 
                                       alternative = "greater", 
                                       prior.mu=0, 
                                       prior.var=1, 
                                       iter = 1000, 
                                       iter.DA = 1000,
                                       forecastmodel = "H0",
                                       mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotH0_10, file="SIM_pilotH0_10.RData")   

SIM_pilotInformedH1_00 <- designAnalysis.pilot(n.pilot = 20,
                                             n.main = c(30, 80, 130, 180), 
                                             ES.pop = 0, 
                                             alternative = "greater", 
                                             prior.mu=0.3, 
                                             prior.var=0.0225, 
                                             iter = 1000, 
                                             iter.DA = 1000,
                                             forecastmodel = "H1",
                                             mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_00, file="SIM_pilotInformedH1_00.RData")                          

SIM_pilotInformedH1_01 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.1, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_01, file="SIM_pilotInformedH1_01.RData")                          

SIM_pilotInformedH1_02 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.2, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_02, file="SIM_pilotInformedH1_02.RData")     

SIM_pilotInformedH1_03 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.3, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_03, file="SIM_pilotInformedH1_03.RData")   

SIM_pilotInformedH1_04 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.4, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_04, file="SIM_pilotInformedH1_04.RData") 

SIM_pilotInformedH1_05 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.5, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_05, file="SIM_pilotInformedH1_05.RData")

SIM_pilotInformedH1_06 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.6, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_06, file="SIM_pilotInformedH1_06.RData")  

SIM_pilotInformedH1_07 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.7, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_07, file="SIM_pilotInformedH1_07.RData")    

SIM_pilotInformedH1_08 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.8, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_08, file="SIM_pilotInformedH1_08.RData")    

SIM_pilotInformedH1_09 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 0.9, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_09, file="SIM_pilotInformedH1_09.RData")   

SIM_pilotInformedH1_10 <- designAnalysis.pilot(n.pilot = 20,
                                               n.main = c(30, 80, 130, 180), 
                                               ES.pop = 1, 
                                               alternative = "greater", 
                                               prior.mu=0.3, 
                                               prior.var=0.0225, 
                                               iter = 1000, 
                                               iter.DA = 1000,
                                               forecastmodel = "H1",
                                               mcmc.setting = list(chains=2, burnin=0, iter=5000))

save(SIM_pilotInformedH1_10, file="SIM_pilotInformedH1_10.RData")                          

