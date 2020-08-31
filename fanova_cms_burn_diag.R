# Generate burn-in sampling files for functional ANOVA

library(fields)
library(ncdf4)
library(ggplot2)
library(plyr)
source("func_anova_fns.R")

theme_mat = theme_bw() 
theme_mat$strip.text$size = 12
theme_mat$axis.title.x$size = 12
theme_mat$axis.text.x$size = 11
theme_mat$axis.title.y$size = 12
theme_mat$axis.text.y$size = 11
theme_mat$plot.title$size = 12
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 10
theme_mat$legend.title$size = 10


exmpcfg = read.csv("fanova_cms_namer.csv",header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

# GP parameters
gpgrps = c("mean","mainA","mainB","interact","noise")
gpbnms = paste0("Example_CMS_NAmer/sp_pars_",gpgrps,"_burn.png")
gprfrm = data.frame(vrnm=gpgrps,ofile=gpbnms,nrow=rep(2,5),width=rep(1200,5),height=rep(900,5),stringsAsFactors = FALSE)
d_ply(gprfrm,c("vrnm"),.fun=gp_partrace_burn,.progress="text",cflst=cfglst,pltthm=theme_mat)

# MH Acceptance
mhfrm = gprfrm
mhfrm$width=rep(750,5)
mhfrm$height=rep(600,5)
mhfrm$ofile = paste0("Example_CMS_NAmer/mh_sp_pars_",gpgrps,"_burn.png")
d_ply(mhfrm,c("vrnm"),.fun=ascatter,.progress="text",cflst=cfglst,pltthm=theme_mat)

# ANOVA fields at selected locations
fldfrm = gprfrm[1:4,]
fldfrm$nctrst = c(0,1,1,1)
fldfrm$locst = rep(2,4)
fldfrm$locfn = rep(98,4)
fldfrm$locinc = rep(12,4)
fldfrm$nrow = rep(3,4)
fldfrm$ndigit = rep(1,4)
fldfrm$ofile = paste0("Example_CMS_NAmer/field_locs_",gpgrps[1:4],"_burn.png")
d_ply(fldfrm,c("vrnm"),.fun=gp_fldtrace_burn,.progress="text",cflst=cfglst,pltthm=theme_mat)


