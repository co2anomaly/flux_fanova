# Generate initial states for functional ANOVA
# This is a "second round" initialization, assuming that initial GP samples have been run with fixed params

library(ncdf4)
library(jsonlite)

cfgfile = "config/fanova_mip_namer_prdev_sptm.csv"
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

# groupings of GPs
vct = 0
varlst = list()
gpgrps = c("mean","mainA","mainB","interact","noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness")
gplng = c("overall mean","main effect factor A","main effect factor B","interaction","replicate error")

# Copy original init files to storage, final to init 
for (j in seq(1,nchain)) {
    # Copy to burn/post final state files
    brnint = paste0(cfglst$burn_init_file,j,".nc")
    brntmp = paste0(cfglst$burn_init_file,"_Orig",j,".nc")
    brnfnl = paste0(cfglst$burn_final_file,j,".nc")
    file.copy(brnint,brntmp,overwrite = TRUE)
    file.copy(brnfnl,brnint,overwrite = TRUE)

    # Perturb inits
    nc1 = nc_open(brnint,write=TRUE)
    for (g1 in seq(1,length(gpgrps))) {
        for (k in seq(1,length(gpvnms))) {
            gvrnm = sprintf("%s_%s",gpvnms[k],gpgrps[g1])
            crvl = ncvar_get(nc1,gvrnm)
            if (gpvnms[k] == "gp_stddev") {
                crlg = 0.5 * log(crvl)
            }
            else {
                crlg = log(crvl)
            }
            psmp = exp(rnorm(1,mean = crlg, sd=0.1))
            ncvar_put(nc1,gvrnm,psmp)
        }
        if (gpgrps[g1] == "noise") {
            ncnm = sprintf("ar_correlation_%s",gpgrps[g1])
            crvl = ncvar_get(nc1,ncnm)
            psmp = runif(1,min=crvl-0.15,max=crvl+0.15)
            ncvar_put(nc1,ncnm,psmp)
        }
    }
    nc_close(nc1)
}


