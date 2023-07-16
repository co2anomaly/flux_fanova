# Generate initial states for functional ANOVA
# This is a "second round" initialization, assuming that initial GP samples have been run with fixed params

library(ncdf4)
library(jsonlite)

cfgfile = "config/fanova_cms_eurasia_jja_sptm.csv"
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

}


