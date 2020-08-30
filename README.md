# flux_fanova
Bayesian functional analysis of variance for carbon flux estimates from multiple models

The processing pipeline for remote sensing data products involves several implementation choices that can influence the scientific conclusions drawn from the final products. Among the various levels of processing, Level 4 products often have the most scientific utility to the broadest community of users. At the same time, these products typically require the most complex processing algorithms and largest number of upstream dependencies. When possible, it is valuable to evaluate the robustness of the key features of the final product to different choices for these implementation factors. 

The software in this repository implements a statistical model known as functional analysis of variance (FANOVA) to estimate the common features and algorithm-specific anomalies in a collection of Level 4 products. The approach can be applied to global and regional carbon flux estimates produced from inversion systems that assimilate Level 2 and/or Level 3 carbon dioxide products. The [NASA MEaSUREs](https://climatesciences.jpl.nasa.gov/co2measures) project on Records of Fused and Assimilated Satellite Carbon Dioxide Observations and Fluxes from Multiple Instruments is developing these products.  

The functional ANOVA methodology is based on [Kaufman and Sain (2010)](https://doi.org/10.1214/10-BA505)

***

### Basic simulation

1. Simulate data fields and output with `fanova_sim_data.R`
2. Create burn-in sampling files with `fanova_sim_setup_burn.R`
3. Generate initial states with `fanova_sim_geninit.R`
4. Run burn-in for each chain with `fanova_sim_mcmc_burn.R`
    - Run within R console/RStudio:  
    `> config = 'fanova_example_config.csv'`  
    `> chain = 1`  
    `> system(paste("Rscript fanova_sim_mcmc_burn.R",config,chain))`
    - Run in batch (remote/cluster):  
    `Rscript fanova_sim_mcmc_burn.R fanova_example_config.csv 1`
    - Repeat for desired number of chains (4 in example)
5. Burn-in diagnostics with `fanova_sim_burn_diag.R`
6. Run post burn-in sampling with `fanova_sim_mcmc_post.R`  
Follow same convention as burn-in step
