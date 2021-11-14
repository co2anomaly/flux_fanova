## OCO-2 Flux Model Intercomparison Project

This implementation executes the functional ANOVA for model intercomparison projects (MIPs), such as the [OCO-2 flux MIP](https://gml.noaa.gov/ccgg/OCO2_v9mip/). The statistical model uses the [Vecchia approximation](https://doi.org/10.1214/19-STS755) for Gaussian processes for computational efficiency. Supporting MCMC routines are found in `func_anova_vecchia.R` and `func_anova_fns.R`

This example uses a subset of four models from the MIP and examines the monthly flux estimates for June-July-August (JJA) of 2016.

### Development build of GPvecchia

The spatial covariance routines use the [GPvecchia R package](https://github.com/katzfuss-group/GPvecchia)

The internal matrix operations in the CRAN version of GPvecchia can sometimes cause the MCMC sampler to crash, but the development version of the package handles the errors without crashing. Therefore, it can be useful to build the package from the repository directly. Some general guidelines for the procedure can be found as part of the [R package tutorial](https://kbroman.org/pkg_primer/pages/build.html)

The procedure to build from the repository on a remote Linux system is

* Clone repository: `git clone https://github.com/katzfuss-group/GPvecchia.git`
* Add `markdown` as a dependency in DESCRIPTION
* Create/add the following line to `.Renviron` file
```
R_LIBS_USER=/home/user/Rlibs
```
The `R_LIBS_USER` directory is the directory where the library will be installed. It does not need to be the same location as the cloned repository.
* Build the package with `R CMD BUILD GPvecchia` or `R CMD build GPvecchia`
* Install the package in local directory with `R CMD INSTALL GPvecchia_0.1.3.tar.gz`

### Data and configuration preparation

* The data are extracted from the MIP collection with the `mip_data_north_amer_gdal.R` script in the `processing` directory. The spatial covariance functions in the GPvecchia package compute distances internally, so the data processing script converts location information to UTM coordinates.
* Configuration
* Perform initial likelihood grid search for GP parameters with `fanova_mip_mltlev_vec_like.R`  
Site-specific ANOVA estimates are computed and the grid search is performed for the intercept, main effects, and interaction.

### MCMC burn-in setup

### MCMC burn-in sampling and evaluation

### MCMC stationary (post burn-in) setup

### MCMC stationary sampling and evaluation
