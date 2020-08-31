# Plot CMS posterior fields

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 14
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 14
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14

smryfrm = function(dfrm,vlvr = "value") {
    # Summarize dataset with quantiles
    dtsb = dfrm[,vlvr]
    dtsb = dtsb[!is.na(dtsb)]
    qs = quantile(dtsb,probs = c(0.05,0.5,0.95))
    frmout = data.frame(Q05 = qs[1], Q50 = qs[2], Q95 = qs[3])
}

# Read CMS Locs
nc1 = nc_open("Example_CMS_NAmer/CMS_NAmer_Flux.nc")
latcms = as.vector(ncvar_get(nc1,"latitude"))
loncms = as.vector(ncvar_get(nc1,"longitude"))
nc_close(nc1)

nloc = length(latcms)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = loncms,Latitude = latcms)

# Read posterior results
fctnms = c("mean","mainA","mainB","interact","noise")
fctlbs = c("Mean","Season","DataSrc","Interact","Noise")
nfct = length(fctnms)
fctord = seq(1,nfct)
fctfrm = data.frame(ANOVACmp=fctord,ANOVALbl=fctlbs)

niter = 5000
nchn =  4

rngarr = array(0,c(niter,nfct,nchn))
sdvarr =  array(0,c(niter,nfct,nchn))
nuarr = array(0,c(niter,nfct,nchn))
#gp_stddev_

for (j in seq(1,nchn)) {
    
    pstnm = paste0("Example_CMS_NAmer/NAmer_PostSamp",j,".nc")
    ncpst = nc_open(pstnm)
    for (k in seq(1,nfct)) {
        gpnm = paste0("gp_stddev_",fctnms[k])
        sdvarr[,k,j] = ncvar_get(ncpst,gpnm)
        gpnm = paste0("gp_smoothness_",fctnms[k])
        nuarr[,k,j] = ncvar_get(ncpst,gpnm)
        gpnm = paste0("gp_range_",fctnms[k])
        rngarr[,k,j] = ncvar_get(ncpst,gpnm)
    }
    nc_close(ncpst)
  
}

sdmlt = melt(sdvarr,varnames = c("Iteration","ANOVACmp","Chain"),value.name = "GPStdDev")
sdmlt$GPStdDev_gYr = sdmlt$GPStdDev * 1314.9
sdsmry = ddply(sdmlt,c("ANOVACmp"),.fun=smryfrm, vlvr="GPStdDev_gYr")
sdsmry$ANOVACmp = factor(sdsmry$ANOVACmp,labels = fctfrm$ANOVALbl)

rngmlt = melt(rngarr,varnames = c("Iteration","ANOVACmp","Chain"),value.name = "GPRange")
rngsmry = ddply(rngmlt,c("ANOVACmp"),.fun=smryfrm, vlvr="GPRange")
rngsmry$ANOVACmp = factor(rngsmry$ANOVACmp,labels = fctfrm$ANOVALbl)

gsig = ggplot(sdsmry,aes(x=ANOVACmp,y=Q50)) + geom_errorbar(aes(ymin=Q05,ymax=Q95),width=0.4,size=1.2) + 
       geom_point(size=5) + scale_y_log10("GP Std Deviation",limits=c(0.01,15)) + xlab("ANOVA Component") + 
       ggtitle("GP Standard Deviation Posterior") + theme_mat 
pnm = "CMS_GPStdDev_gYr_Posterior_NAmer_Log10.pdf"
pdf(pnm,width = 8,height=6)
print(gsig)
dev.off()

gsig = ggplot(sdsmry,aes(x=ANOVACmp,y=Q50)) + geom_errorbar(aes(ymin=Q05,ymax=Q95),width=0.4,size=1.2) + 
  geom_point(size=5) + scale_y_continuous("GP Std Deviation",limits=c(0.01,15)) + xlab("ANOVA Component") + 
  ggtitle("GP Standard Deviation Posterior") + theme_mat 
pnm = "CMS_GPStdDev_gYr_Posterior_NAmer.pdf"
pdf(pnm,width = 8,height=6)
print(gsig)
dev.off()

sdmlt$ANOVACmp = factor(sdmlt$ANOVACmp,labels = fctfrm$ANOVALbl)
gsig2 = ggplot(sdmlt,aes(x=ANOVACmp,y=GPStdDev_gYr)) + geom_violin() + 
  scale_y_continuous("GP Std Deviation",limits=c(0.01,15)) + xlab("ANOVA Component") + 
  ggtitle("GP Standard Deviation Posterior") + theme_mat 
pnm = "CMS_GPStdDev_gYr_Violin_NAmer.pdf"
pdf(pnm,width = 8,height=6)
print(gsig2)
dev.off()

grng = ggplot(rngsmry,aes(x=ANOVACmp,y=Q50)) + geom_errorbar(aes(ymin=Q05,ymax=Q95),width=0.4,size=1.2) + 
  geom_point(size=5) + scale_y_log10("GP Range [km]",limits=c(20,2000)) + xlab("ANOVA Component") + 
  ggtitle("GP Range Posterior") + theme_mat 
pnm = "CMS_GPRange_Posterior_NAmer_Log10.pdf"
pdf(pnm,width = 8,height=6)
print(grng)
dev.off()

grng = ggplot(rngsmry,aes(x=ANOVACmp,y=Q50)) + geom_errorbar(aes(ymin=Q05,ymax=Q95),width=0.4,size=1.2) + 
  geom_point(size=5) + scale_y_continuous("GP Range [km]",limits=c(0,1600)) + xlab("ANOVA Component") + 
  ggtitle("GP Range Posterior") + theme_mat 
pnm = "CMS_GPRange_Posterior_NAmer.pdf"
pdf(pnm,width = 8,height=6)
print(grng)
dev.off()

