### Development build of GPvecchia

The spatial covariance routines use the [GPvecchia R package](https://github.com/katzfuss-group/GPvecchia)  
The `flux_fanova` MCMC routines generally work well with the recent v0.1.4 release of GPvecchia

**Older versions of GPvecchia**

For older versions of GPvecchia, the internal matrix operations in the CRAN version of GPvecchia can sometimes cause the MCMC sampler to crash, but the development version of the package handles the errors without crashing. Therefore, it can be useful to build the package from the repository directly. Some general guidelines for the procedure can be found as part of the [R package tutorial](https://kbroman.org/pkg_primer/pages/build.html)

The necessary R dependencies are
* For a development build, all of the packages GPvecchia depends on are needed: `Rcpp`, `methods`, `stats`, `sparseinv`, `fields`, `Matrix`, `parallel`, `GpGp`, `FNN`, `markdown`, `RcppArmadillo`, `BH`
* Additional packages needed for building documentation `rmarkdown`, `knitr`, `mvtnorm`

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

