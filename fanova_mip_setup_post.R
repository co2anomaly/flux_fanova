# Generate post burn-in sampling files for functional ANOVA

library(ncdf4)

cfgfile = "config/fanova_mip_namer_prdev.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

nchain = as.integer(cfglst$nchain)
niter = as.integer(cfglst$npost)
nthin = as.integer(cfglst$nthin)
nloc = as.integer(cfglst$nloc)
nsamp = floor(niter / nthin)
nadapt = as.integer(cfglst$post_mh_adapt)
nmhsv = floor(niter / nadapt)


nmnA = as.integer(cfglst$nleva)
nmnB = as.integer(cfglst$nlevb)
nab = (nmnA-1) * (nmnB-1)

dimlc = ncdim_def("location","",seq(1,nloc),create_dimvar = FALSE)
dimlvA = ncdim_def("level_A","",seq(1,nmnA))
dimlvB = ncdim_def("level_B","",seq(1,nmnB))
dimeffA = ncdim_def("effect_A","",seq(1,nmnA-1),create_dimvar = FALSE)
dimeffB = ncdim_def("effect_B","",seq(1,nmnB-1),create_dimvar = FALSE)
dimeffAB = ncdim_def("effect_interact","",seq(1,nab))
dimsmp = ncdim_def("iteration","",1:nsamp,create_dimvar = FALSE)
dimadp = ncdim_def("adapt","",1:nmhsv,create_dimvar = FALSE)

# groupings of GPs
vct = 0
varlst = list()
gpgrps = c("mean","mainA","mainB","interact","noise")
gplng = c("overall mean","main effect factor A","main effect factor B","interaction","replicate error")
for (g1 in seq(1,length(gpgrps))) {
  vct = vct + 1
  sgnm = paste0("gp_stddev_",gpgrps[g1])
  sgln = paste0("Gaussian process standard deviation for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(sgnm,"None",list(dimsmp),-9999,
                            longname=sgln)

  vct = vct + 1
  rgnm = paste0("gp_range_",gpgrps[g1])
  rgln = paste0("Gaussian process range for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(rgnm,"None",list(dimsmp),-9999,
                            longname=rgln)
  
  vct = vct + 1
  nunm = paste0("gp_smoothness_",gpgrps[g1])
  nuln = paste0("Gaussian process smoothness for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(nunm,"None",list(dimsmp),-9999,
                            longname=nuln)
}

# GP fields
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mean","None",list(dimsmp,dimlc),-9999,
                          longname="Gaussian process field for overall mean")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mainA","None",list(dimsmp,dimlc,dimeffA),-9999,
                          longname="Gaussian process field for main effect factor A")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_mainB","None",list(dimsmp,dimlc,dimeffB),-9999,
                          longname="Gaussian process field for main effect factor B")
vct = vct + 1
varlst[[vct]] = ncvar_def("gp_field_interact","None",list(dimsmp,dimlc,dimeffAB),-9999,
                          longname="Gaussian process field for main effect factor B")

# MH variables
for (g1 in seq(1,length(gpgrps))) {
  vct = vct + 1
  sgnm = paste0("mh_stddev_",gpgrps[g1])
  sgln = paste0("MH jumping standard deviation for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(sgnm,"None",list(dimadp),-9999,
                            longname=sgln)
  
  vct = vct + 1
  rgnm = paste0("mh_acc_rate_",gpgrps[g1])
  rgln = paste0("MH acceptance rate for ",gpgrps[g1])
  varlst[[vct]] = ncvar_def(rgnm,"None",list(dimadp),-9999,
                            longname=rgln)
}


# Generate sampling files for each chain
for (j in seq(1,nchain)) {
    ncnm = paste0(cfglst$post_samp_file,j,".nc")
    nc1 = nc_create(ncnm,varlst)
    for (i in seq(1,length(varlst))) {
        ncatt_put(nc1,varlst[[i]],"missing_value",-9999.0)
    }
    nc_close(nc1)
}

