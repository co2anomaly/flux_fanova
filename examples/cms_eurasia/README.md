## CMS-Flux Data Source Example 

This implementation executes the functional ANOVA for the CMS-Flux inversion system for different time periods and data sources. The statistical model uses the [Vecchia approximation](https://doi.org/10.1214/19-STS755) for Gaussian processes for computational efficiency. Supporting MCMC routines are found in `func_anova_vecchia.R` and `func_anova_fns.R`

This example uses inversions using two different approaches to aggregation of atmospheric carbon dioxide data from OCO-2, including spatially aggregated retrievals available through a [NASA MEaSUREs project](https://doi.org/10.5067/582L7HTJ343N). Flux estimates are analyzed for the JJA season for 2015 and 2016 over Eurasia. 

**Compuational Note:** The MCMC sampling scripts (`mcmc_burn`, `mcmc_post`) can have long runtimes. These can be run from the R console/RStudio interactively or non-interactively in batch mode.

* Run within R console/RStudio:  
`> config = 'config/fanova_cms_eurasia_jja.csv'`  
`> chain = 1`  
`> system(paste("Rscript fanova_cms_mcmc_burn_rnd1.R",config,chain))`
* Run in batch (remote/cluster):  
`Rscript fanova_sim_mcmc_burn.R fanova_example_config.csv 1`
* Repeat for desired number of chains (4 in example)

The [posterior summaries and maps](#mcmc-stationary-summary) can be reproduced without running the full MCMC samplers by loading the archived posterior sampling data files.

### Data and configuration preparation

* The flux estimates are packaged in the NetCDF file `CMS_Eurasia_Flux_JJA_2Yr.nc`. This file should placed in the `examples/cms_eurasia` directory. The `Coast.csv` file should be placed in the same directory if maps will be plotted.
* Optional: generate maps of fluxes with `fanova_cms_map_flux.R` 
* The configuration files for this example are in the `config` directory
    - ANOVA config: `fanova_cms_eurasia_jja.csv`
    - MCMC and prior parameters: `CMS_Eurasia_JJA_Prior_Beta.json`
* Perform initial likelihood grid search for GP parameters with `fanova_cms_vecchia_like.R`  
Site-specific ANOVA estimates are computed and the grid search is performed for the intercept, main effects, and interaction.

### MCMC burn-in setup

* Setup output files for MCMC burn-in with `fanova_cms_setup_burn.R`
* Generate initial states for round 1 with `fanova_cms_geninit_rnd1.R`

### MCMC burn-in sampling and evaluation

The MCMC burn-in works best by running in three *rounds*, with different parameters fixed and sampled for each round.

* Run MCMC round 1 with `fanova_cms_mcmc_burn_rnd1.R`. For this round, only the ANOVA component fields are sampled and the GP parameters are fixed. Upon completion, diagnostic plots can be generated by running the notebook `CMS_Eurasia_Burn_JJA.Rmd`
* Generate initial states for round 2 with `fanova_cms_geninit_rnd2.R`
* Run MCMC round 2 with `fanova_cms_mcmc_burn_rnd2.R`. For this round, the ANOVA component fields are sampled along with the noise process (epsilon) GP parameters. Upon completion, diagnostic plots can be generated by running the notebook `CMS_Eurasia_Burn_JJA.Rmd`
* Generate initial states for round 3 with `fanova_cms_geninit_rnd3.R`
* Run MCMC round 3 with `fanova_cms_mcmc_burn_rnd3.R`. For this round, the ANOVA component fields are sampled along with all GP parameters. Upon completion, diagnostic plots can be generated by running the notebook `CMS_Eurasia_Burn_JJA.Rmd`

### MCMC stationary setup, sampling, and evaluation

* Setup output files for MCMC posterior with `fanova_cms_setup_post.R`
* Run stationary MCMC with `fanova_cms_mcmc_post.R`. Upon completion, diagnostic plots can be generated by running the notebook `CMS_Eurasia_Post_JJA.Rmd`

### MCMC stationary summary 

The posterior summary results can be produced with the scripts below and the external datasets with the posterior samples.

* The stationary posterior distributions are summarized with the final sections of the notebook `CMS_Eurasia_Post_JJA.Rmd`
* Posteior maps are produced with `fanova_postfield_cms.R` 
