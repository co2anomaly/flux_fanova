# Map CMS-Flux datasets, JJA consecutive years

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("mapfunctions.R")
source("func_anova_fns.R")

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


# Read data file
nc1 = nc_open(cfglst$data_file)
lon = ncvar_get(nc1,"longitude")
lat = ncvar_get(nc1,"latitude")
flux = ncvar_get(nc1,"CO2_Flux")
fluxpr = ncvar_get(nc1,"CO2_Flux_Prior")
fluxdv = ncvar_get(nc1,"CO2_Flux_Dev_Prior")
nc_close(nc1)

locsq = seq(1,length(lon))
locfrm = data.frame(LocIdx=locsq,Longitude=lon,Latitude=lat)

dsrc = c("Control","Fused")
yrs = c(2015,2016)


# Plot original fluxes
dimnames(flux) = list(LocIdx=locsq,Year=yrs,DataSource=dsrc,Month=1:3)
flxall = melt(flux,value.name = "Flux")
flxmrg = merge(flxmrg,locfrm)

flxmrg = fklon(flxmrg,lonvar = "Longitude",center=90)
flxmrg$Season = "JJA"
flxmrg$Trt = paste(flxmrg$Year,flxmrg$DataSource,sep="_")
flxmrg$MFct = factor(flxmrg$Month,levels = 1:3,labels=c("Jun","Jul","Aug"))

flxmrg$FluxPlt = flxmrg$Flux
flxmrg$FluxPlt[flxmrg$Flux < -2000] = -1999.0
flxmrg$FluxPlt[flxmrg$Flux >  2000] = 1999.0
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
gmn = ggplot(flxmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
      geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.6, color="#444444") +
      facet_grid(Trt ~ MFct) + 
      scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) +  
      scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
      scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
      ggtitle("CMS Land Flux") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_Eurasia_JJA_2Yr.pdf")
pdf(pnm,width = 12,height=8 )
print(gmn)
dev.off()

# Flux prior
dimnames(fluxpr) = list(LocIdx=locsq,Year=yrs,DataSource=dsrc,Month=1:3)
prall = melt(fluxpr,value.name = "PrAnn")
prmrg = merge(prall,locfrm)

prmrg = fklon(prmrg,lonvar = "Longitude",center=90)
prmrg$Season = "JJA"
prmrg$Trt = paste(prmrg$Year,prmrg$DataSource,sep="_")
prmrg$MFct = factor(prmrg$Month,levels = 1:3,labels=c("Jun","Jul","Aug"))

prmrg$PrPlt = prmrg$PrAnn
prmrg$PrPlt[prmrg$PrAnn < -2000] = -1999.0
prmrg$PrPlt[prmrg$PrAnn >  2000] = 1999.0
prplt = prmrg[prmrg$DataSource == dsrc[1],]
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
gmn = ggplot(prmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PrPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.6, color="#444444") +
  facet_grid(Year ~ MFct) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) +  
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Land Flux") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_Eurasia_Prior_JJA_2Yr.pdf")
pdf(pnm,width = 12,height=4 )
print(gmn)
dev.off()

# Map deviations
dimnames(fluxdv) = list(LocIdx=locsq,Year=yrs,DataSource=dsrc,Month=1:3)
mrgall = melt(fluxdv,value.name = "PrDev")
dvmrg = merge(mrgall,locfrm)

dvmrg = fklon(dvmrg,lonvar = "Longitude",center=90)
dvmrg$Season = "JJA"
dvmrg$Trt = paste(dvmrg$Year,dvmrg$DataSource,sep="_")
dvmrg$MFct = factor(dvmrg$Month,levels = 1:3,labels=c("Jun","Jul","Aug"))

dvmrg$FluxPlt = dvmrg$PrDev
dvmrg$FluxPlt[dvmrg$PrDev < -2000] = -1999.0
dvmrg$FluxPlt[dvmrg$PrDev >  2000] = 1999.0
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
gmn = ggplot(dvmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.6, color="#444444") +
  facet_grid(Trt ~ MFct) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) +  
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Land Flux") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_Eurasia_PrDev_JJA_2Yr.pdf")
pdf(pnm,width = 12,height=8 )
print(gmn)
dev.off()

