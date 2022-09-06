# Process MIP datasets for functional ANOVA, map monthly values
# North America, JJA 2016
# Save locations as UTM coordinates

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)
library(rgdal)

source("mapfunctions.R")
source("func_anova_fns.R")

cfgfile = "config/fanova_mip_africa_prdev.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

dtdr = cfglst$data_src_dir 
wrldfl = paste0(cfglst$output_dir,"/Coast.csv")

wrld = centermap2(wrldfl,center=15.0,intrarg = FALSE)
lnlb = lonlbs(seq(-20,40,by=20),center=15)
ltlb = latlbs(seq(-40,40,by=20))

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14

read_mip_subset = function(dfrm, mtchfrm, nlon = 360,nlat = 180) {
    # Read flux MIP results and subset to selected regions
    # Input data frame contains
    #     FileName: Name of flux file
    #     LonVar: Longitude variable name
    #     LatVar: Latitude variable name
    #     LandFlux: Land flux variable name
    #     OceanFlux: 
    #     TimeStart
    #     TimeCount
    # mtchfrm is the data frame with region assignments
    # Total flux is sum of land and ocean, although regions generally delineate these anyway
    # Flux variables should be 3D (lat,lon,time)
  
    nc1 = nc_open(dfrm$FileName)
    lon = ncvar_get(nc1,dfrm$LonVar)
    lat = ncvar_get(nc1,dfrm$LatVar)
    flxocn = ncvar_get(nc1,dfrm$OceanFlux,start=c(1,1,dfrm$TimeStart),count=c(nlon,nlat,dfrm$TimeCount))
    flxlnd = ncvar_get(nc1,dfrm$LandFlux,start=c(1,1,dfrm$TimeStart),count=c(nlon,nlat,dfrm$TimeCount))
    nc_close(nc1)
    
    flxnet = flxocn + flxlnd
    #flxnet = flxlnd
    dimnames(flxnet)[1] = list(lon)
    dimnames(flxnet)[2] = list(lat)
    flxmlt = melt(flxnet,varnames=c("Longitude","Latitude","Time"),value.name = "Flux")
    
    flxmrg = merge(flxmlt,mtchfrm)
    return(flxmrg)
    
}

mdnms = c('Ames','Baker','CMS-Flux','OU')
lndvrs = c('land','land','land','land')
ocnvrs = c('ocean','ocean','ocean','ocean')
expt = c('LNLG','IS')

mipfrm = data.frame(Model=rep(mdnms,2), LandFlux=rep(lndvrs,2), OceanFlux=rep(ocnvrs,2), Expt=rep(expt,each=4))
mipfrm$FileName = sprintf("%s/%s_gridded_fluxes_%s.nc4",dtdr,mipfrm$Model,mipfrm$Expt)
mipfrm$TimeStart = 18
mipfrm$TimeCount = 3
mipfrm$LonVar = "longitude"
mipfrm$LatVar = "latitude"


# Region Masks
mskfl = sprintf("%s/oco2_regions_l4mip_v7.nc",dtdr) 
nc1 = nc_open(mskfl)
mlat = ncvar_get(nc1,"latitude")
mlon = ncvar_get(nc1,"longitude")
rgn = ncvar_get(nc1,"mip_regions")
mipnm = ncvar_get(nc1,"mip_names")
tcm = ncvar_get(nc1,"transcom_regions")
tcmnm = ncvar_get(nc1,"transcom_names")
nc_close(nc1)
dimnames(rgn)[1] = list(mlon)
dimnames(rgn)[2] = list(mlat)
miprgns = melt(rgn,varnames = c("Longitude","Latitude"),value.name = "Region")

# Subset
miprgns = miprgns[(miprgns$Region == 6) | (miprgns$Region == 7) | (miprgns$Region == 8) | (miprgns$Region == 9),]
miprgns = miprgns[order(miprgns$Latitude,miprgns$Longitude),]

mipfnl = ddply(mipfrm,c("Model","Expt"),.fun = read_mip_subset,.progress = "text",mtchfrm = miprgns)
mipfnl = fklon(mipfnl,lonvar = "Longitude",center=15)
tmlbls = c("6_Jun","7_Jul","8_Aug")
mipfnl$Time = factor(mipfnl$Time,levels=c(1,2,3),labels = tmlbls)

# Plots, impose high end cutoff
mipfnl$FluxPlt = mipfnl$Flux
mipfnl$FluxPlt[mipfnl$Flux < -2000] = -1999.0
mipfnl$FluxPlt[mipfnl$Flux >  2000] = 1999.0
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
for (i in seq(1,length(expt))) {
  frmsb = mipfnl[mipfnl$Expt == expt[i],]  
  ttlstr = paste0('Africa JJA 2016: ',expt[i])
  gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
    geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
    facet_grid(Time ~ Model) + 
    scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-35,35)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(-40,40)) + 
    ggtitle(ttlstr) + theme_mat + coord_equal()
  paste0(cfglst$output_dir,"/Africa_2016_JJA_",expt[i],".pdf")
  pdf(pnm,width = 12,height=7)
  print(gmn)
  dev.off()
}

# Prior fields
prfrm = data.frame(Model=rep(mdnms,2), LandFlux=rep(lndvrs,2), OceanFlux=rep(ocnvrs,2), Expt=rep(expt,each=4))
prfrm$FileName = sprintf("%s/%s_gridded_fluxes_Prior.nc4",dtdr,prfrm$Model)
prfrm$TimeStart = 18
prfrm$TimeCount = 3
prfrm$LonVar = "longitude"
prfrm$LatVar = "latitude"

prfnl = ddply(prfrm,c("Model","Expt"),.fun = read_mip_subset,.progress = "text",mtchfrm = miprgns)
prfnl = fklon(prfnl,lonvar = "Longitude",center=15)
tmlbls = c("6_Jun","7_Jul","8_Aug")
prfnl$Time = factor(prfnl$Time,levels=c(1,2,3),labels = tmlbls)

prfnl$PrAnn = prfnl$Flux
prfnl$FluxPlt = prfnl$Flux
prfnl$FluxPlt[prfnl$Flux < -2000] = -1999.0
prfnl$FluxPlt[prfnl$Flux >  2000] = 1999.0
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
frmsb = prfnl[prfnl$Expt == expt[1],]  
ttlstr = 'Africa JJA 2016: Prior'
gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_grid(Time ~ Model) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-35,35)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(-40,40)) + 
  ggtitle(ttlstr) + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/Africa_2016_JJA_Prior.pdf")
pdf(pnm,width = 12,height=7)
print(gmn)
dev.off()

# Merge and map deviations
mrgall = merge(mipfnl[,c("Longitude","Latitude","Time","Region","Model","fk360","Expt","Flux")],
               prfnl[,c("Longitude","Latitude","Time","Region","Model","fk360","Expt","PrAnn")])
mrgall$PrDev = mrgall$Flux - mrgall$PrAnn

mrgall$FluxPlt = mrgall$PrDev
mrgall$FluxPlt[mrgall$PrDev < -2000] = -1999.0
mrgall$FluxPlt[mrgall$PrDev >  2000] = 1999.0
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)

for (i in seq(1,length(expt))) {
  frmsb = mrgall[mrgall$Expt == expt[i],]  
  ttlstr = paste0('Africa JJA 2016 Post-Prior: ',expt[i])
  gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
    geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
    facet_grid(Time ~ Model) + 
    scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000)) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-35,35)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(-40,40)) + 
    ggtitle(ttlstr) + theme_mat + coord_equal()
  pnm = paste0(cfglst$output_dir,"/Africa_2016_JJA_PrDev_",expt[i],".pdf")
  pdf(pnm,width = 12,height=7 )
  print(gmn)
  dev.off()
}

# Prepare Output
nloc = nrow(miprgns)
nrep = 3
nlva = length(mdnms)
nlvb = length(expt)

y=array(dim=c(nloc,nlva,nlvb,nrep))
ypr=array(dim=c(nloc,nlva,nlvb,nrep))
for(i in seq(1,nlva)){
  for(j in seq(1,nlvb)){
    for(k in 1:nrep){
      frmsb = mipfnl[(mipfnl$Model == mdnms[i]) & (mipfnl$Expt == expt[j]) & (mipfnl$Time == tmlbls[k]),]
      frmsb = frmsb[order(frmsb$Latitude,frmsb$Longitude),]
      y[,i,j,k] = frmsb$Flux
      prsb = mrgall[(mrgall$Model == mdnms[i]) & (mrgall$Expt == expt[j]) & (mrgall$Time == tmlbls[k]),]
      prsb = prsb[order(prsb$Latitude,prsb$Longitude),]
      ypr[,i,j,k] = prsb$PrAnn
    }
  }
}
yprdv = y - ypr

# Plot summary mean
flxmn = apply(y,1,FUN=mean)
ymnfrm = frmsb
ymnfrm$Mu=as.vector(flxmn)
ymnfrm$Mu[ymnfrm$Mu < -2000] = -1999.0
ymnfrm$Mu[ymnfrm$Mu >  2000] = 1999.0

gmn = ggplot(ymnfrm,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Mu)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-2000,2000))  + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-35,35)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(-40,40)) + 
  ggtitle("MIP Flux Summary Statistic: Mean") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Africa_Summary_Mean_JJA_2016.pdf")
pdf(pnm,width = 12,height=8)
print(gmn)
dev.off()

# Get UTM coordinates
cord.dec <- SpatialPoints(cbind(frmsb$Longitude, frmsb$Latitude), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +lon_0=15e +zone=33 +ellps=GRS80"))


nab = (nlva-1) * (nlvb-1)
# Contrasts
hmtb = gp_scale_helmert(nlvb)
hmta = gp_scale_helmert(nlva)
hlst = list(hmta,hmtb)

dimlc = ncdim_def("location","",seq(1,nloc),create_dimvar = FALSE)
dimlvA = ncdim_def("model_A","",seq(1,nlva),create_dimvar = FALSE)
dimlvB = ncdim_def("data_source_B","",seq(1,nlvb),create_dimvar = FALSE)
dimeffA = ncdim_def("effect_A","",seq(1,nlva-1),create_dimvar = FALSE)
dimeffB = ncdim_def("effect_B","",seq(1,nlvb-1),create_dimvar = FALSE)
dimeffAB = ncdim_def("effect_interact","",seq(1,nab))
dimrep = ncdim_def("replicate","",seq(1,nrep))
dimchr = ncdim_def("char","",1:40,create_dimvar = FALSE)

varx = ncvar_def("longitude","degrees_east",list(dimlc),-9999,
                 longname="Longitude")
vary = ncvar_def("latitude","degrees_north",list(dimlc),-9999,
                 longname="Latitude")
varux = ncvar_def("utm_easting","km",list(dimlc),-9.999e6, 
                  longname="UTM Easting Coordinate")
varuy = ncvar_def("utm_northing","km",list(dimlc),-9.999e6, 
                  longname="UTM Northing Coordinate")
vardt = ncvar_def("CO2_Flux","None",list(dimlc,dimlvA,dimlvB,dimrep),-9999,
                  longname="Flux data fields")
varpr = ncvar_def("CO2_Flux_Prior","None",list(dimlc,dimlvA,dimlvB,dimrep),-9999,
                  longname="Prior flux")
vardv = ncvar_def("CO2_Flux_Dev_Prior","None",list(dimlc,dimlvA,dimlvB,dimrep),-9999,
                  longname="Posterior-Prior flux")
varctrstA = ncvar_def("contrast_A","None",list(dimlvA,dimeffA),-9999,
                      longname="Main effect A zero sum contrasts")
varctrstB = ncvar_def("contrast_B","None",list(dimlvB,dimeffB),-9999,
                      longname="Main effect B zero sum contrasts")
varA = ncvar_def("model","",list(dimchr,dimlvA),prec="char")
varB = ncvar_def("data_source","",list(dimchr,dimlvB),prec="char")

nc1 = nc_create(cfglst$data_file,list(varx,vary,varux,varuy,
                                      vardt,varpr,vardv,varctrstA,varctrstB,varA,varB))
ncvar_put(nc1,varx,frmsb$Longitude)
ncatt_put(nc1,varx,"missing_value",-9999.0)
ncvar_put(nc1,vary,frmsb$Latitude)
ncatt_put(nc1,vary,"missing_value",-9999.0)
ncvar_put(nc1,varux,cord.UTM@coords[,1]*0.001)
ncatt_put(nc1,varux,"missing_value",-9.999e6)
ncvar_put(nc1,varuy,cord.UTM@coords[,2]*0.001)
ncatt_put(nc1,varuy,"missing_value",-9.999e6)

ncvar_put(nc1,vardt,y)
ncatt_put(nc1,vardt,"missing_value",-9999.0)
ncatt_put(nc1,vardt,"units","g C/m^2/sec")
ncvar_put(nc1,varpr,ypr)
ncatt_put(nc1,varpr,"missing_value",-9999.0)
ncatt_put(nc1,varpr,"units","g C/m^2/sec")
ncvar_put(nc1,vardv,yprdv)
ncatt_put(nc1,vardv,"missing_value",-9999.0)
ncatt_put(nc1,vardv,"units","g C/m^2/sec")
ncvar_put(nc1,varctrstA,hmta)
ncatt_put(nc1,varctrstA,"missing_value",-9999.0)
ncvar_put(nc1,varctrstB,hmtb)
ncatt_put(nc1,varctrstB,"missing_value",-9999.0)
ncvar_put(nc1,varA,mdnms)
ncvar_put(nc1,varB,expt)
nc_close(nc1)


