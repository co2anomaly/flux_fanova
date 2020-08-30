# flux_fanova
Bayesian functional analysis of variance for carbon flux estimates from multiple models

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
