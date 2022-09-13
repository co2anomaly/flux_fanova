# flux_fanova
Bayesian functional analysis of variance for carbon flux estimates from multiple models

The processing pipeline for remote sensing data products involves several implementation choices that can influence the scientific conclusions drawn from the final products. Among the various levels of processing, Level 4 products often have the most scientific utility to the broadest community of users. At the same time, these products typically require the most complex processing algorithms and largest number of upstream dependencies. When possible, it is valuable to evaluate the robustness of the key features of the final product to different choices for these implementation factors. 

The software in this repository implements a statistical model known as functional analysis of variance (FANOVA) to estimate the common features and algorithm-specific anomalies in a collection of Level 4 products. The approach can be applied to global and regional carbon flux estimates produced from inversion systems that assimilate Level 2 and/or Level 3 carbon dioxide products. The [NASA MEaSUREs](https://climatesciences.jpl.nasa.gov/co2measures) project on Records of Fused and Assimilated Satellite Carbon Dioxide Observations and Fluxes from Multiple Instruments is developing these products.  

The functional ANOVA methodology is based on [Kaufman and Sain (2010)](https://doi.org/10.1214/10-BA505)

***

### Examples

Several examples of the functional ANOVA approach for carbon flux inversions are found in the `examples` directory

* [CMS-Flux Eurasia](examples/cms_eurasia/README.md): CMS-Flux esimates over Eurasia for two years, using different aggregation methods
* [MIP North America](examples/mip_namer/README.md): OCO-2 V9 flux model intercomparison project estimates over North America
* [MIP Africa](examples/mip_africa/README.md): OCO-2 V9 flux model intercomparison project estimates over Africa

***

### Dependencies

The functional ANOVA implementations are in R. Additional R packages required include 

* `jsonlite`
* `reshape2`
* `plyr`
* `ggplot2`
* `colorspace`
* `fields`
* `ncdf4`
    - The system will also need NetCDF libraries installed
* `rgdal`
    - This library is only required if generating data files from the OCO-2 flux MIP, e.g. `mip_data_north_amer_gdal.R`
    - The system will also need libgdal installed
* `GPvecchia` (development build preferred)
    - For a development build, all of the packages GPvecchia depends on are needed: `Rcpp`, `methods`, `stats`, `sparseinv`, `fields`, `Matrix`, `parallel`, `GpGp`, `FNN`, `markdown`, `RcppArmadillo`, `BH`
    - Additional packages needed for building documentation `rmarkdown`, `knitr`, `mvtnorm`

***

### Development build of GPvecchia

The spatial covariance routines use the [GPvecchia R package](https://github.com/katzfuss-group/GPvecchia)

The internal matrix operations in the CRAN version of GPvecchia can sometimes cause the MCMC sampler to crash, but the development version of the package handles the errors without crashing. Therefore, it can be useful to build the package from the repository directly. Some general guidelines for the procedure can be found as part of the [R package tutorial](https://kbroman.org/pkg_primer/pages/build.html)

The procedure to build from the repository on a remote Linux system is

* Clone repository: `git clone https://github.com/katzfuss-group/GPvecchia.git`
* Add `markdown` to the Imports list in DESCRIPTION
* Create/add the following line to `.Renviron` file
```
R_LIBS_USER=/home/user/Rlibs
```
The `R_LIBS_USER` directory is the directory where the library will be installed. It does not need to be the same location as the cloned repository.
* Build the package with `R CMD BUILD GPvecchia` or `R CMD build GPvecchia`
* Install the package in local directory with `R CMD INSTALL GPvecchia_0.1.3.tar.gz`

