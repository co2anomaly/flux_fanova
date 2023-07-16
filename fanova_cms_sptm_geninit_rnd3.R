# Generate initial states for functional ANOVA
# This is a "second round" initialization, assuming that initial GP samples have been run with fixed params
# Here only the noise GP params are updated

library(ncdf4)
library(jsonlite)

exmpcfg = read.csv("config/fanova_cms_eurasia_jja_sptm.csv",header = TRUE,stringsAsFactors = FALSE)
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
gpgrps = c("noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness","ar_correlation")
gplng = c("replicate error")

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


