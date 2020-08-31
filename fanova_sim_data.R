# Simulation for functional ANOVA
# Save dataset to NetCDF

library(fields)
library(ncdf4)
source("func_anova_fns.R")

exmpcfg = read.csv("fanova_example_config.csv")
set.seed(133204)

### dimensions and sizes
smdm = c(10,10)
lv2 = c(2,2)
nrp = 5
ic=jc=c(-1,1)
nab = (lv2[1]-1)*(lv2[2]-1)

### set true parameters
mu0=3
#sig.mu=4; sig.a=sig.b=3; sig.ab=2; sig.eps=1
#ra.mu=.4; ra.a=ra.b=.3; ra.ab=.2; ra.eps=.1

# covariance params
rnglst = list(mean = 4, mainA = 3, mainB = 3, interac = 2, noise = 1)
siglst = list(mean = 2, mainA = 1.5, mainB = 1.5, interac = 1, noise = 0.5)
simnu = 1.5

simout = spat_anova_sim_2fact(smdm,lv2,nrp, rnglst, siglst) 

# Save data (and generating parameters) in NetCDF
# Variables for contrasts?
dimlc = ncdim_def("location","",seq(1,nrow(simout$locfrm)),create_dimvar = FALSE)
dimlvA = ncdim_def("level_A","",seq(1,lv2[1]))
dimlvB = ncdim_def("level_B","",seq(1,lv2[2]))
dimeffA = ncdim_def("effect_A","",seq(1,lv2[1]-1),create_dimvar = FALSE)
dimeffB = ncdim_def("effect_B","",seq(1,lv2[2]-1),create_dimvar = FALSE)
dimeffAB = ncdim_def("effect_interact","",seq(1,nab))
dimrep = ncdim_def("replicate","",seq(1,nrp))

varx = ncvar_def("coord_easting","None",list(dimlc),-9999,
                 longname="West-East coordinate variable")
vary = ncvar_def("coord_northing","None",list(dimlc),-9999,
                 longname="South-North coordinate variable")
vardt = ncvar_def("data_field","None",list(dimlc,dimlvA,dimlvB,dimrep),-9999,
                  longname="Simulated data fields")
varmu = ncvar_def("true_mean","None",list(dimlc),-9999,
                  longname="True overall mean field")
varA = ncvar_def("true_mainA","None",list(dimlc,dimeffA),-9999,
                 longname="True A main effect field")
varB = ncvar_def("true_mainB","None",list(dimlc,dimeffB),-9999,
                 longname="True A main effect field")
varAB = ncvar_def("true_interact","None",list(dimlc,dimeffAB),-9999,
                 longname="True interaction effect field")
varctrstA = ncvar_def("contrast_A","None",list(dimlvA),-9999,
                      longname="Main effect A zero sum contrasts")
varctrstB = ncvar_def("contrast_B","None",list(dimlvB),-9999,
                      longname="Main effect B zero sum contrasts")

nc1 = nc_create("Examples/TwoFactor_SimData.nc",list(varx,vary,vardt,varmu,varA,varB,varAB,varctrstA,varctrstB))
ncvar_put(nc1,varx,simout$locfrm[,1])
ncatt_put(nc1,varx,"missing_value",-9999.0)
ncvar_put(nc1,vary,simout$locfrm[,2])
ncatt_put(nc1,vary,"missing_value",-9999.0)
ncvar_put(nc1,vardt,simout$datarr)
ncatt_put(nc1,vardt,"missing_value",-9999.0)
ncvar_put(nc1,varmu,as.vector(simout$mu))
ncatt_put(nc1,vardt,"missing_value",-9999.0)
ncvar_put(nc1,varA,as.vector(simout$alpha))
ncatt_put(nc1,varA,"missing_value",-9999.0)
ncvar_put(nc1,varB,as.vector(simout$beta))
ncatt_put(nc1,varB,"missing_value",-9999.0)
ncvar_put(nc1,varAB,as.vector(simout$ab))
ncatt_put(nc1,varAB,"missing_value",-9999.0)
ncvar_put(nc1,varctrstA,ic)
ncatt_put(nc1,varctrstA,"missing_value",-9999.0)
ncvar_put(nc1,varctrstB,ic)
ncatt_put(nc1,varctrstB,"missing_value",-9999.0)
nc_close(nc1)


