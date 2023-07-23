# Plot CMS posterior fields

library(ggplot2)
library(ncdf4)
library(colorspace)
library(dplyr)
library(tidyr)

source("mapfunctions_tidy.R")

cfgfile = "config/fanova_cms_eurasia_jja_sptm.csv"
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
fctlbs = c("Mean","Year","Aggregation","Interaction")
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
dimnames(pstfldmn) = list(LocID = seq(1,nloc),
                       CompID = sprintf("Comp%d",seq(1,4)) )
fldfrm = as_tibble(pstfldmn, rownames="LocID") 
fldfrm = fldfrm %>% pivot_longer(cols = starts_with("Comp"), names_to = "ANOVACmp",
                             names_prefix="Comp", values_to = "PostMean")
fldfrm = fldfrm %>% mutate_at(c('LocID','ANOVACmp'), as.numeric)

pstmrg = fldfrm %>% left_join(fctfrm, by=c("ANOVACmp"))
pstmrg = fldfrm %>% left_join(lcfrm, by="LocID")
pstmrg = fklon_tbl(pstmrg,lonvar = "Longitude",center=90)


pstmrg$ANOVACmp = factor(pstmrg$ANOVACmp,labels = fctfrm$ANOVALbl)
pstmrg$FluxPlt = pstmrg$PostMean
pstmrg$FluxPlt[pstmrg$FluxPlt < -400] = -399.0
pstmrg$FluxPlt[pstmrg$FluxPlt >  400] = 399.0

r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
flxlb2 = bquote(atop("Flux Increment", "gC" ~ "m" ^ -2 ~ "yr" ^-1))

gmn = ggplot(pstmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  facet_wrap(~ ANOVACmp) + 
  scale_fill_gradientn(flxlb2,colors=r7,limits=c(-400,400)) + 
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
    geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
    facet_wrap(~ ANOVACmp,nrow=1) + 
    scale_fill_gradientn(flxlb2,colors=r7,limits=c(lwlm[k],lmsq[k])) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
    theme_mat + coord_equal()
  pnm = paste0(cfglst$output_dir,"/CMS_ANOVA_Posterior_",fctlbs[k],"_Eurasia_JJA.pdf")
  pdf(pnm,width = 7.5,height=3)
  print(gbi)
  dev.off()
}


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
dimnames(mnA) = list(LocID = seq(1,nloc),
                          YearID = sprintf("Year%d",yrtxt) )
afrm = as_tibble(mnA, rownames="LocID") 
afrm = afrm %>% pivot_longer(cols = starts_with("Year"), names_to = "Year",
                                 names_prefix="Year", values_to = "PostMean")
afrm = afrm %>% mutate_at(c('LocID'), as.numeric)

amrg = afrm %>% left_join(lcfrm, by="LocID")
amrg = fklon_tbl(amrg,lonvar = "Longitude",center=90)
amrg$FluxPlt = amrg$PostMean
amrg$FluxPlt[amrg$FluxPlt < -400] = -399.0
amrg$FluxPlt[amrg$FluxPlt >  400] = 399.0

gseas = ggplot(amrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  facet_wrap(~ Year,nrow=2) + 
  scale_fill_gradientn(flxlb2,colors=r7,limits=c(-400,400)) + 
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


prbmlt = tibble(LocID=seq(1,length(pstmnprb)), PstMnPrb=pstmnprb, PstYrPrb=pstyrprb )
prbmrg = prbmlt %>% left_join(lcfrm, by="LocID")
prbmrg = fklon_tbl(prbmrg,lonvar = "Longitude",center=90)

gbmn = ggplot(prbmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PstMnPrb)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  scale_fill_continuous("Prob", type="viridis", limits=c(0.9,1.0)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("Posterior Probability (Mean > Aggregation)") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_PostProb_MeanAgg","_Eurasia_JJA.pdf")
pdf(pnm,width = 7.5,height=3)
print(gbmn)
dev.off()

gbyr = ggplot(prbmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PstYrPrb)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  scale_fill_continuous("Prob", type="viridis", limits=c(0.55,0.65)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("Posterior Probability (Year > Aggregation)") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_PostProb_YearAgg","_Eurasia_JJA.pdf")
pdf(pnm,width = 7.5,height=3)
print(gbyr)
dev.off()
