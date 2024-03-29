## OCO-2 Flux Model Intercomparison Project

This implementation executes the functional ANOVA for model intercomparison projects (MIPs), such as the [OCO-2 flux MIP](https://gml.noaa.gov/ccgg/OCO2_v9mip/). The statistical model uses the [Vecchia approximation](https://doi.org/10.1214/19-STS755) for Gaussian processes for computational efficiency. Supporting MCMC routines are found in `func_anova_vecchia.R` and `func_anova_fns.R`

This example uses a subset of four models from the OCO-2 V9 MIP and examines the monthly flux estimates for June-July-August (JJA) of 2016 over Africa.

**Configuration Note:** All of the `fanova_mip` R scripts apply to the MIP ANOVA for any region. The ANOVA config CSV file may need to be changed in these scripts to point to the appropriate region's example.  
**Compuational Note:** The MCMC sampling scripts (`mcmc_burn`, `mcmc_post`) can have long runtimes. These can be run from the R console/RStudio interactively or non-interactively in batch mode. 

* Run within R console/RStudio:  
`> config = 'config/fanova_mip_africa_prdev.csv'`  
`> chain = 1`  
`> system(paste("Rscript fanova_mip_mcmc_burn_rnd1.R",config,chain))`
* Run in batch (remote/cluster):  
`Rscript fanova_sim_mcmc_burn.R fanova_mip_africa_prdev.csv 1`
* Repeat for desired number of chains (4 in example)

The [posterior summaries and maps](#mcmc-stationary-summary) can be reproduced without running the full MCMC samplers by loading the archived posterior sampling data files.

### Data and configuration preparation

* The MIP fluxes and priors should be requested and downloaded from the [OCO-2 V9 MIP website](https://gml.noaa.gov/ccgg/OCO2_v9mip/). In addition, the region mask file, linked at the same site under protocol details, should also be downloaded and stored in the same directory. For the analysis in this example, the following datasets should be downloaded from the [MIP downloads page](https://gml.noaa.gov/ccgg/OCO2_v9mip/download.php)
    - Gridded IS Fluxes
    - Gridded LNLG Fluxes
    - Gridded Prior Fluxes
* The configuration files for this example are in the `config` directory
    - ANOVA config: `fanova_mip_africa_prdev.csv`  
In addition to file and directory settings for the functional ANOVA, this file includes a setting `data_src_dir` that should identify the directory that stores the downloaded MIP results.
    - MCMC and prior parameters: `MIP_Africa_PrDev_Prior_Beta.json`
* The `Coast.csv` file should placed in the `examples/cms_eurasia` directory if maps will be plotted.
* The data are extracted from the MIP collection, mapped, and prepped for analysis with the `mip_data_africa_gdal.R` script. The spatial covariance functions in the GPvecchia package compute distances internally, so the data processing script converts location information to UTM coordinates.
* Perform initial likelihood grid search for GP parameters with `fanova_mip_vecchia_like.R`  
Site-specific ANOVA estimates are computed and the grid search is performed for the intercept, main effects, and interaction.

### MCMC burn-in setup

* Setup output files for MCMC burn-in with `fanova_mip_setup_burn.R`
* Generate initial states for round 1 with `fanova_mip_geninit_rnd1.R`

### MCMC burn-in sampling and evaluation

The MCMC burn-in works best by running in three *rounds*, with different parameters fixed and sampled for each round.

* Run MCMC round 1 with `fanova_mip_mcmc_burn_rnd1.R`. For this round, only the ANOVA component fields are sampled and the GP parameters are fixed. Upon completion, diagnostic plots can be generated by running the notebook `MIP_Africa_Burn_PrDev.Rmd`
* Generate initial states for round 2 with `fanova_mip_geninit_rnd2.R`
* Run MCMC round 2 with `fanova_mip_mcmc_burn_rnd2.R`. For this round, the ANOVA component fields are sampled along with the noise process (epsilon) GP parameters. Upon completion, diagnostic plots can be generated by running the notebook `MIP_Africa_Burn_PrDev.Rmd `
* Generate initial states for round 3 with `fanova_mip_geninit_rnd3.R`
* Run MCMC round 3 with `fanova_mip_mcmc_burn_rnd3.R`. For this round, the ANOVA component fields are sampled along with all GP parameters. Upon completion, diagnostic plots can be generated by running the notebook `MIP_Africa_Burn_PrDev.Rmd`

### MCMC stationary setup, sampling, and evaluation

* Setup output files for MCMC posterior with `fanova_mip_setup_post.R`
* Run stationary MCMC with `fanova_mip_mcmc_post.R`. Upon completion, diagnostic plots can be generated by running the notebook `MIP_Africa_Post_PrDev.Rmd`

### MCMC stationary summary 

The posterior summary results can be produced with the scripts below and the external datasets with the posterior samples.

* The stationary posterior distributions are summarized with the final sections of the notebook `MIP_Africa_Post_PrDev.Rmd`
* Posteior maps are produced with `fanova_postfield_mip_africa.R` 

