# Generate initial states for functional ANOVA
# This is a "second round" initialization, assuming that initial GP samples have been run with fixed params
# Here only the noise GP params are updated

library(fields)
library(ncdf4)
library(jsonlite)
library(plyr)
source("func_anova_fns.R")

cfgfile = "config/fanova_mip_namer_prdev.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

prinfo = fromJSON(txt = cfglst$json_file)

insd = as.integer(cfglst$init_seed)
set.seed(insd)
nchain = as.integer(cfglst$nchain)
nloc = as.integer(cfglst$nloc)

nmnA = as.integer(cfglst$nleva)
nmnB = as.integer(cfglst$nlevb)
nab = (nmnA-1) * (nmnB-1)
nrp = as.integer(cfglst$nrep)

dimlc = ncdim_def("location","",seq(1,nloc),create_dimvar = FALSE)
dimlvA = ncdim_def("level_A","",seq(1,nmnA))
dimlvB = ncdim_def("level_B","",seq(1,nmnB))
dimeffA = ncdim_def("effect_A","",seq(1,nmnA-1),create_dimvar = FALSE)
dimeffB = ncdim_def("effect_B","",seq(1,nmnB-1),create_dimvar = FALSE)
dimeffAB = ncdim_def("effect_interact","",seq(1,nab))
dimrep = ncdim_def("replicate","",seq(1,nrp))
dimscl = ncdim_def("scalar","",seq(1,1),create_dimvar = FALSE)

# Read in Dataset
nc1 = nc_open(cfglst$data_file)
dtarr = ncvar_get(nc1,cfglst$data_variable)
ctrstA = ncvar_get(nc1,"contrast_A")
ctrstB = ncvar_get(nc1,"contrast_B")
nc_close(nc1)

# Likelihood estimates
load(cfglst$like_est)

# groupings of GPs
vct = 0
varlst = list()
gpgrps = c("noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness")
gplng = c("replicate error")
for (g1 in seq(1,length(gpgrps))) {
  vct = vct + 1
  sgnm = paste0("gp_stddev_",gpgrps[g1])
  sgln = paste0("Gaussian process standard deviation for ",gpgrps[g1])
  varlst[[vct]] = sgnm

  vct = vct + 1
  rgnm = paste0("gp_range_",gpgrps[g1])
  rgln = paste0("Gaussian process range for ",gpgrps[g1])
  varlst[[vct]] = rgnm 
 
  vct = vct + 1
  nunm = paste0("gp_smoothness_",gpgrps[g1])
  nuln = paste0("Gaussian process smoothness for ",gpgrps[g1])
  varlst[[vct]] = nunm
}

# Read noise GP pars from final state file and store in initial 
for (j in seq(1,nchain)) {
    # Copy to burn/post final state files
    brnint = paste0(cfglst$burn_init_file,j,".nc")
    brnfnl = paste0(cfglst$burn_final_file,j,".nc")

    gpstor = rep(0,length(gpvnms)) 
    ncf = nc_open(brnfnl,write=FALSE)
    for (g1 in seq(1,length(gpgrps))) {
        for (k in seq(1,length(gpvnms))) {
            gvrnm = paste0(gpvnms[k],"_",gpgrps[g1])
            gpstor[k] = ncvar_get(ncf,gvrnm)
        }
    } 
    nc_close(ncf)

    # Generate GP options
    nc1 = nc_open(brnint,write=TRUE)
    for (g1 in seq(1,length(gpgrps))) {
        for (k in seq(1,length(gpvnms))) {
            gvrnm = paste0(gpvnms[k],"_",gpgrps[g1])
            print(gvrnm)
            ncvar_put(nc1,gvrnm,gpstor[k])
        }
    }
    nc_close(nc1)
}


