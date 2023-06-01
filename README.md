# Interim Design Analysis Using Bayes Factor Forecasts

Repository containing code for the paper "Interim Design Analysis Using Bayes Factor Forecasts" by Angelika Stefan, Quentin Gronau, and Eric-Jan Wagenmakers.

## Code Structure

The code is organized in an R project. The project is organized as follows:

The code is separated in three main folders:


:file_folder: __ttest_2sample_normalprior__

This folder contains all necessary code to conduct a Bayes factor forecast for an independent-samples t-test with a normal prior on $\delta$ under $\mathcal{M}_1$. 

* __bf_2sample_t.R__ contains functions to compute Bayes factors
* __bfda_fixed_2sample_t.R__ contains functions to compute a fixed-N BFDA for the independent-samples t-test
* __forecast_2sample_t.R__ contains function to run the Bayes factor forecast
* __mcmc_2sample_t.R__ contains function to sample from joint posterior distribution
* __zzz_posterior_2sample_t.R__ contains function to sample from posterior distribution for $\delta$ 

:file_folder: __bf_forecast__

This folder contains all necessary code to reproduce the results and plots in the first section of the manuscript (introduction of interim BFDA).

* __bfda_plot.R__ contains the code to reproduce the visualization of the outcome of a fixed-N BFDA (density plot of Bayes factor distribution)
* __plot_forecast_layers.R__ contains the code to reproduce the visualization of the Bayes factor forecast (model-averaged and model-specific; fan-plot design)
* __zzz_forecast_fanplot.R__ contains a simpler version of the fan plot produced by `plot_forecast_layers.R`
* __zzz_forecast_noisyfanplot.R__ contains a function to plot the results of a Bayes factor forecast similar to `plot_forecast_layers.R`, but representing them as many trajectories instead of a fan
* __zzz_forecast_violinplot.R__ contains a function to plot the results of a model-averaged Bayes factor forecast as a violinplot, i.e., two vertical density plots attached to each other

:file_folder: __applications__

This folder contains all necessary code to reproduce the results and plots in the two sections presenting simulated application examples. All scripts with the prefix `sim_pilot` refer to the Bayesian internal pilot study design. Scripts with the prefix `sim_stoppingFutility` refer to the sequential Max-N design with futility stopping. In contrast, scripts with the prefix `sim_seq` refer to the sequential Max-N design (without futility stopping). `sim_both` refers to the direct comparison of the sequential design with and without futility stopping. 

Scripts should be executed according to the following order of postfixes: `Functions` $\rightarrow$ `Run` $\rightarrow$ `Analyze` $\rightarrow$ `Plots` .

*Bayesian Internal Pilot Studies*

* __sim_pilot_Functions.R__ contains functions to simulate a Bayesian internal pilot study design and conduct a design analysis
* __sim_pilot_Run.R__ is a script to run the design analyses for the Bayesian internal pilot study design
* __sim_pilot_Analyze.R__ contains functions to analyze the results of a design analysis for the Bayesian internal pilot study design. Allows setting of stopping thresholds, desired power, and stopping for efficiency / insufficient conditional power.
* __sim_pilot_Plots.R__ contains code to reproduce all plots in the manuscript section covering the Bayesian internal pilot study design (as well as in the appendix).
* __bfda_power.R__ reproduces the plot explaining the large proportion of early stopping in the Bayesian internal pilot study design using the results of fixed-N BFDAs.
* __zzz_sim_pilot_Plots_Scenarios.R__ contains code to make a plot similar to `sim_stoppingFutility_Plot-StoppingScenarios.R` that was not included in the manuscript
* __zzz_sim_pilot_pointEstimates.R__ contains code that bases the continuation decisions in pilot studies on point estimates. Not included in the manuscript.

*Futility Stopping in Sequential Bayesian Desings* 

* __sim_stoppingFutility_Functions.R__ contains functions to simulate a Bayesian sequential design with futility stopping and conduct a design analysis
* __sim_stoppingFutility_Run.R__ is a script to run the design analyses for the Bayesian sequential design with futility stopping
* __sim_stoppingFutility_Analyze.R__ reproduces the stacked barplot figures presented in the manuscript
* __sim_stoppingFutility_Plot-StoppingScenarios.R__ reproduces the panels in the figure depicting the stopping scenarios in the manuscript
* __tryout_sim_stoppingFutility_HighFutility.R__ contains an exploratory analysis where the futility stopping threshold is increased. Analysis was conducted mainly as a validity check for the results of `sim_stoppingFutility_Analyze.R`.
* __sim_stoppingFutility_AppendixRun.R__ contains script to run additional simulations of the futility stopping design presented in the online appendix

*Bayesian Sequential Max-N Design (no futility stopping)*
* __sim_seq_Run.R__ contains a script to run the design analyses for the Bayesian sequential Max-N design based on the [`BFDA` R-package](https://github.com/nicebread/BFDA)
* __sim_seq_Analyze.R__ contains a script to reproduce the stacked barplot figure presented in the manuscript for the sequential max-N design
* __sim_both_Analyze.R__ reproduces the comparison of expected sample sizes for the futility stopping and the regular sequential Max-N design

:file_folder: __appendix__

Scripts to reproduce the plots in the online appendix. The simulation results themselves can be reproduced using the script __sim_stoppingFutility_AppendixRun.R__ in the "applications" folder.

* __sim_stoppingFutility_differentPrior_plot.R__ Contains plotting functions for the simulation condition with the informed prior
* __sim_stoppingFutility_highthreshold_plot.R__ Contains plotting functions for the simulation condition with the higher futility threshold
* __sim_stoppingFutility_smallN_plot.R__ Contains plotting functions for the simulation condition with smaller max-N

## General Remarks

* Most figures in the paper were produced by post-processing the R-generated plots in Microsoft PowerPoint. PowerPoint slides and PDF versions of the figures can be found in the associated OSF repository [https://osf.io/xh2ep/](https://osf.io/xh2ep/)

* Simulation results overall take several days to run on a standard multi-core computer. For easier reproducibility, we saved simulation results in the associated OSF repository [https://osf.io/xaqh4/](https://osf.io/xaqh4/)

* Files labeled with the prefix `zzz` contain scripts and functions used in code development that were not included in the final manuscript. However, they may be useful in other circumstances, so we kept them. 

* Files labeled with the prefix `tryout` contain additional analyses usually conducted as validity checks.

```
> sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.6.8

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggpattern_1.0.1 ggplot2_3.4.2  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9           lattice_0.20-45      prettyunits_1.1.1    png_0.1-8            class_7.3-20        
 [6] ps_1.7.1             utf8_1.2.2           V8_4.2.1             plyr_1.8.7           R6_2.5.1            
[11] stats4_4.2.1         e1071_1.7-13         pillar_1.8.1         rlang_1.1.1          curl_4.3.2          
[16] pscl_1.5.5           rstudioapi_0.14      callr_3.7.2          magick_2.7.4         Matrix_1.5-1        
[21] labeling_0.4.2       stringr_1.4.1        loo_2.5.1            munsell_0.5.0        proxy_0.4-27        
[26] compiler_4.2.1       gridpattern_1.0.2    rstan_2.26.13        pkgconfig_2.0.3      pkgbuild_1.3.1      
[31] tidyselect_1.1.2     tibble_3.1.8         gridExtra_2.3        codetools_0.2-18     matrixStats_0.62.0  
[36] fansi_1.0.3          viridisLite_0.4.1    crayon_1.5.2         dplyr_1.0.10         withr_2.5.0         
[41] sf_1.0-13            MASS_7.3-57          grid_4.2.1           jsonlite_1.8.2       gtable_0.3.1        
[46] lifecycle_1.0.3      DBI_1.1.3            magrittr_2.0.3       units_0.8-2          StanHeaders_2.26.13 
[51] scales_1.2.1         RcppParallel_5.1.5   KernSmooth_2.23-20   stringi_1.7.8        cli_3.4.1           
[56] cachem_1.0.6         reshape2_1.4.4       farver_2.1.1         LaplacesDemon_16.1.6 generics_0.1.3      
[61] vctrs_0.6.2          tools_4.2.1          glue_1.6.2           purrr_0.3.5          processx_3.7.0      
[66] parallel_4.2.1       fastmap_1.1.0        inline_0.3.19        colorspace_2.0-3     classInt_0.4-9      
[71] memoise_2.0.1   
```