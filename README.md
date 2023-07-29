# flux_fanova
Bayesian functional analysis of variance for carbon flux estimates from multiple models

The processing pipeline for remote sensing data products involves several implementation choices that can influence the scientific conclusions drawn from the final products. Among the various levels of processing, Level 4 products often have the most scientific utility to the broadest community of users. At the same time, these products typically require the most complex processing algorithms and largest number of upstream dependencies. When possible, it is valuable to evaluate the robustness of the key features of the final product to different choices for these implementation factors. 

The software in this repository implements a statistical model known as functional analysis of variance (FANOVA) to estimate the common features and algorithm-specific anomalies in a collection of Level 4 products. The approach can be applied to global and regional carbon flux estimates produced from inversion systems that assimilate Level 2 and/or Level 3 carbon dioxide products. The [NASA MEaSUREs](https://climatesciences.jpl.nasa.gov/co2measures) project on Records of Fused and Assimilated Satellite Carbon Dioxide Observations and Fluxes from Multiple Instruments is developing these products.  

The functional ANOVA methodology is based on [Kaufman and Sain (2010)](https://doi.org/10.1214/10-BA505)

***

### Examples

Several examples of the functional ANOVA approach for carbon flux inversions are found in the `examples` directory

* Examples with a spatio-temporal model for the functional ANOVA error term
    - [CMS-Flux Eurasia](examples/cms_eurasia_sptm/README.md): CMS-Flux esimates over Eurasia for two years, using different aggregation methods
    - [MIP North America](examples/mip_namer_sptm/README.md): OCO-2 V9 flux model intercomparison project estimates over North America
    - [MIP Africa](examples/mip_africa_sptm/README.md): OCO-2 V9 flux model intercomparison project estimates over Africa
* Examples for the spatial only functional ANOVA model
    - [CMS-Flux Eurasia](examples/cms_eurasia/README.md): CMS-Flux esimates over Eurasia for two years, using different aggregation methods
    - [MIP North America](examples/mip_namer/README.md): OCO-2 V9 flux model intercomparison project estimates over North America
    - [MIP Africa](examples/mip_africa/README.md): OCO-2 V9 flux model intercomparison project estimates over Africa

***

### Dependencies

The functional ANOVA implementations are in R. Additional R packages required include 

* `jsonlite`
* `tidyverse`
* `ggplot2`
* `colorspace`
* `fields`
* `ncdf4`
    - The system will also need NetCDF libraries installed
* `rgdal`
    - This library is only required if generating data files from the OCO-2 flux MIP, e.g. `mip_data_north_amer_tidy.R`
    - The system will also need libgdal installed
* `GPvecchia` (v0.1.4 or later preferred)
    - For earlier versions, a [development build](GPvecchia_devbuild.md) is needed
* `scales`
    - This library is only required for prior distribution visualizations

