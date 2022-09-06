# Generate initial states for functional ANOVA

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

# Coef estimates to initialize
cfbsln = aaply(dtarr,1,.fun = fit_anova_mltlev_arr3d, .progress="text", contrstA=ctrstA, contrstB=ctrstB)
sdbsln = apply(dtarr,1:3,sd)
gpsd = median(sdbsln) * 0.5
cfsd = sd(cfbsln)

# groupings of GPs
vct = 0
varlst = list()
gpgrps = c("mean","mainA","mainB","interact","noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness")
gplng = c("overall mean","main effect factor A","main effect factor B","interaction","replicate error")
for (g1 in seq(1,length(gpgrps))) {
  vct = vct + 1
  sgnm = paste0("gp_stddev_",gpgrps[g1])
  sgln = paste0("Gaussian process standard deviation for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(sgnm,"None",list(dimscl),-9999,
                            longname=sgln)

  vct = vct + 1
  rgnm = paste0("gp_range_",gpgrps[g1])
  rgln = paste0("Gaussian process range for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(rgnm,"None",list(dimscl),-9999,
                            longname=rgln)
  
  vct = vct + 1
  nunm = paste0("gp_smoothness_",gpgrps[g1])
  nuln = paste0("Gaussian process smoothness for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(nunm,"None",list(dimscl),-9999,
                            longname=nuln)
}

# GP fields
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mean","None",list(dimlc),-9999,
                          longname="Gaussian process field for overall mean")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mainA","None",list(dimlc,dimeffA),-9999,
                          longname="Gaussian process field for main effect factor A")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mainB","None",list(dimlc,dimeffB),-9999,
                          longname="Gaussian process field for main effect factor B")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_interact","None",list(dimlc,dimeffAB),-9999,
                          longname="Gaussian process field for main effect factor B")


# Generate sampling files for each chain
for (j in seq(1,nchain)) {
    ncnm = paste0(cfglst$burn_init_file,j,".nc")
    nc1 = nc_create(ncnm,varlst)
    for (i in seq(1,length(varlst))) {
        ncatt_put(nc1,varlst[[i]],"missing_value",-9999.0)
    }
    
    # Initial values: from prior for GP pars
    for (g1 in seq(1,length(gpgrps))) {
        for (k in seq(1,length(gpvnms))) {
            gvrnm = paste0(gpvnms[k],"_",gpgrps[g1])
            prprs = prinfo$prior[[ gpgrps[g1] ]][[ gpvnms[k] ]]
            if (prprs$form == "halfnormal") {
                psmp = abs(rnorm(1,sd=prprs$stddev))
                if (gpvnms[k] == "gp_stddev") {
                    if (gpgrps[g1] == "noise") {
                        psmp = exp(log(psmp) * 0.4 + 0.6 * log(gpsd))
                    }
                    else {
                        psmp = exp(log(psmp) * 0.4 + 0.6 * log(gpsd))
                    }
                }
            }
            else if (prprs$form == "lognormal") {
                psmp = exp(rnorm(1, mean = prprs$mean, sd=prprs$stddev))
            }
            else if (prprs$form == "beta") {
                psmp = prprs$lim1 + (prprs$lim2 - prprs$lim1) * rbeta(1, shape1 = prprs$shape1, shape2 = prprs$shape2)
            }
            ncvar_put(nc1,gvrnm,psmp)
        }
    }
    # Priors for GP fields - need appropriate sequences
    cfct = 1
    cfinit = cfbsln[,cfct] + rnorm(nloc,sd = gpsd[1])
    ncvar_put(nc1,"gp_field_mean",cfinit)
    
    cfct = cfct + 1
    cffn = cfct + nmnA - 2
    cfinit = cfbsln[,cfct:cffn] + array(rnorm(nloc,sd = gpsd[1]),c(nloc,cffn-cfct+1))
    ncvar_put(nc1,"gp_field_mainA",cfinit)
    
    cfct = cffn + 1
    cffn = cfct + nmnB - 2
    cfinit = cfbsln[,cfct:cffn] + array(rnorm(nloc,sd = gpsd[1]),c(nloc,cffn-cfct+1))
    ncvar_put(nc1,"gp_field_mainB",cfinit)
    
    cfct = cffn + 1
    cffn = cfct + nab - 1
    cfinit = cfbsln[,cfct:cffn] + array(rnorm(nloc,sd = gpsd[1]),c(nloc,cffn-cfct+1))
    ncvar_put(nc1,"gp_field_interact",cfinit)
    
    nc_close(nc1)
      
    # Copy to burn/post final state files
    brnfnl = paste0(cfglst$burn_final_file,j,".nc")
    file.copy(ncnm,brnfnl,overwrite = TRUE)
    pstfnl = paste0(cfglst$post_final_file,j,".nc")
    file.copy(ncnm,pstfnl,overwrite = TRUE)
    
}

