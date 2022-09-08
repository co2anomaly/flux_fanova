# Plot CMS posterior fields

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("mapfunctions.R")

cfgfile = "config/fanova_cms_eurasia_jja.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

wrldfl = paste0(cfglst$output_dir,"/Coast.csv")
wrld = centermap2(wrldfl,center=90.0,intrarg = FALSE)
lnlb = lonlbs(seq(0,180,by=60),center=90)
ltlb = latlbs(seq(0,80,by=20))

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14

# Read CMS Locs
nc1 = nc_open(cfglst$data_file)
latcms = as.vector(ncvar_get(nc1,"latitude"))
loncms = as.vector(ncvar_get(nc1,"longitude"))
ctrstA = ncvar_get(nc1,"contrast_A")
ctrstB = ncvar_get(nc1,"contrast_B")
yrtxt = ncvar_get(nc1,"year")
nc_close(nc1)

nloc = length(latcms)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = loncms,Latitude = latcms)

# Read posterior results
fctnms = c("mean","mainA","mainB","interact")
gpnms = c("stddev","range","smoothness")
fctlbs = c("Mean","Year","Aggregation","Interact")
fctord = seq(1,4)
fctfrm = data.frame(ANOVACmp=fctord,ANOVALbl=fctlbs)

niter = 5000
nchn =  4
pstflds = array(0,c(niter,nloc,4,4))
amns = array(0,c(niter,nloc,2))
gparr = array(0,c(niter,3,4,5))

for (j in seq(1,nchn)) {
    
    pstnm = paste0(cfglst$post_samp_file,j,".nc")
    ncpst = nc_open(pstnm)
    for (k in seq(1,4)) {
        gpnm = paste0("gp_field_",fctnms[k])
        if (k == 1) {
            pstflds[,,j,k] = ncvar_get(ncpst,gpnm,start=c(1,1),count=c(niter,nloc))
        }
        else {
            pstflds[,,j,k] = ncvar_get(ncpst,gpnm,start=c(1,1,1),count=c(niter,nloc,1))
        }
        for (t in seq(1,3)) {
            prnm = paste("gp",gpnms[t],fctnms[k],sep="_")
            gparr[,t,j,k] = ncvar_get(ncpst,prnm,start=c(1),count=c(niter))
        }
    }
    for (t in seq(1,3)) {
      prnm = paste("gp",gpnms[t],"noise",sep="_")
      gparr[,t,j,5] = ncvar_get(ncpst,prnm,start=c(1),count=c(niter))
    }
    nc_close(ncpst)
  
}

pstfldmn = apply(pstflds,c(2,4),FUN = mean)
pstmlt = melt(pstfldmn,varnames = c("LocID","ANOVACmp"))

pstmrg = merge(pstmlt,fctfrm)
pstmrg = merge(pstmrg,lcfrm)
pstmrg = fklon(pstmrg,lonvar = "Longitude",center=90)
pstmrg$ANOVACmp = factor(pstmrg$ANOVACmp,labels = fctfrm$ANOVALbl)
pstmrg$FluxPlt = pstmrg$value
pstmrg$FluxPlt[pstmrg$FluxPlt < -400] = -399.0
pstmrg$FluxPlt[pstmrg$FluxPlt >  400] = 399.0

#pstmrg$Flux360 = pstmrg$value * 1314.9
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)

gmn = ggplot(pstmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_wrap(~ ANOVACmp) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-400,400)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Flux ANOVA") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_ANOVA_Posterior_Eurasia_JJA.pdf")
pdf(pnm,width = 12,height=7)
print(gmn)
dev.off()

# Plot Effects Individually
lmsq = c(400,5,2,2)
lwlm = -1.0 * lmsq
for (k in seq(1,length(fctlbs)) ) {
  pstmsb = pstmrg[(pstmrg$ANOVACmp == fctlbs[k]),]
  pstmsb$ANOVACmp = factor(pstmsb$ANOVACmp)
  print(summary(pstmsb$FluxPlt))
  gbi = ggplot(pstmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
    geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
    facet_wrap(~ ANOVACmp,nrow=1) + 
    scale_fill_gradientn("Flux",colors=r7,limits=c(lwlm[k],lmsq[k])) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
    theme_mat + coord_equal()
  pnm = paste0(cfglst$output_dir,"/CMS_ANOVA_Posterior_",fctlbs[k],"_Eurasia_JJA.pdf")
  pdf(pnm,width = 7.5,height=3)
  print(gbi)
  dev.off()
}

# Plot only Aggregation and Interaction Effects
pstmsb = pstmrg[(pstmrg$ANOVACmp == "Aggregation") | (pstmrg$ANOVACmp == "Interact"),]
pstmsb$ANOVACmp = factor(pstmsb$ANOVACmp)
gbi = ggplot(pstmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_wrap(~ ANOVACmp,nrow=2) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2,2)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Flux ANOVA Estimates") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"CMS_ANOVA_Posterior_Agg_Intr_Eurasia_JJA.pdf")
pdf(pnm,width = 8,height=6)
print(gbi)
dev.off()

# Plot Year Effect 
pstmsb = pstmrg[(pstmrg$ANOVACmp == "Aggregation") | (pstmrg$ANOVACmp == "Interact"),]
pstmsb$ANOVACmp = factor(pstmsb$ANOVACmp)
gbi = ggplot(pstmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_wrap(~ ANOVACmp,nrow=2) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2,2)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Flux ANOVA Estimates") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"CMS_ANOVA_Posterior_Agg_Intr_Eurasia_JJA.pdf")
pdf(pnm,width = 8,height=6)
print(gbi)
dev.off()


# Season Means
mncr = array(0,c(nloc,2,2))
for(i in 1:2) { 
  for(j in 1:2) {
    ictrst = ctrstA[i] * ctrstB[j]
    tst1 = pstfldmn[,1] + pstfldmn[,2] * as.vector(ctrstA[i]) + pstfldmn[,3] * as.vector(ctrstB[j]) + pstfldmn[,4] * as.vector(ictrst)
    mncr[,i,j] = as.vector(tst1)
  }
}
mnA = apply(mncr,1:2,FUN=mean)
dimnames(mnA)[2] = list(yrtxt)
mnAfrm = melt(mnA,varnames = c("LocID","Year"))
amrg = merge(mnAfrm,lcfrm)
amrg = fklon(amrg,lonvar = "Longitude",center=90)
amrg$FluxPlt = amrg$value
amrg$FluxPlt[amrg$FluxPlt < -400] = -399.0
amrg$FluxPlt[amrg$FluxPlt >  400] = 399.0
gseas = ggplot(amrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_wrap(~ Year,nrow=2) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-400,400)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Flux ANOVA Estimates") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_ANOVA_Posterior_YearMean_Eurasia_JJA.pdf")
pdf(pnm,width = 8,height=6)
print(gseas)
dev.off()


## Posterior probabilities
pstindmn = abs(pstflds[,,,1]) - abs(pstflds[,,,3])
pstindmnb = array(0,dim(pstindmn))
pstindmnb[pstindmn > 0] = 1
pstmnprb = apply(pstindmnb,2,mean)
pstindyr = abs(pstflds[,,,2]) - abs(pstflds[,,,3])
pstindyrb = array(0,dim(pstindyr))
pstindyrb[pstindyr > 0] = 1
pstyrprb = apply(pstindyrb,2,mean)


prbmlt = data.frame(LocID=seq(1,length(pstmnprb)), PstMnPrb=pstmnprb, PstYrPrb=pstyrprb )
prbmrg = merge(prbmlt,lcfrm)
prbmrg = fklon(prbmrg,lonvar = "Longitude",center=90)

gbmn = ggplot(prbmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PstMnPrb)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_continuous("Prob", type="viridis", limits=c(0.9,1.0)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("Posterior Probability (Mean > Aggregation)") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_PostProb_MeanAgg","_Eurasia_JJA.pdf")
pdf(pnm,width = 7.5,height=3)
print(gbmn)
dev.off()

gbyr = ggplot(prbmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PstYrPrb)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_continuous("Prob", type="viridis", limits=c(0.5,0.6)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("Posterior Probability (Year > Aggregation)") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_PostProb_YearAgg","_Eurasia_JJA.pdf")
pdf(pnm,width = 7.5,height=3)
print(gbyr)
dev.off()
