# Process MIP datasets for functional ANOVA, map monthly values
# North America, JJA 2016
# Save locations as UTM coordinates

library(ggplot2)
library(ncdf4)
library(colorspace)
library(dplyr)
library(tidyr)
library(sp)

source("mapfunctions_tidy.R")
source("func_anova_tidy.R")

cfgfile = "config/fanova_mip_namer_prdev_sptm.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

dtdr = cfglst$data_src_dir 
wrldfl = paste0(cfglst$output_dir,"/Coast.csv")
wrld = centermap2(wrldfl,center=-90.0,intrarg = FALSE)
lnlb = lonlbs(seq(-150,-30,by=30),center=-90)
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

dimnames(rgn) = list(Longitude = mlon,Latitude = sprintf("Lat%.1f",mlat))
rgn_tbl = as_tibble(rgn,rownames = "Longitude")
rgn_tbl = rgn_tbl %>% pivot_longer(cols = starts_with("Lat"), names_to = "Latitude", 
                                  names_prefix = "Lat", values_to = "Region")
rgn_tbl = rgn_tbl %>% mutate_at(c('Longitude','Latitude'), as.numeric)
miprgns = rgn_tbl %>% filter(Region %in% c(1,2))


read_mip_subset = function(dfrm, mtchfrm, nlon = 360,nlat = 180,netoptn = FALSE) {
    # Read flux MIP results and subset to selected regions
    # Input data frame contains
    #     FileName: Name of flux file
    #     LonVar: Longitude variable name
    #     LatVar: Latitude variable name
    #     LandFlux: Land flux variable name
    #     OceanFlux: Ocean flux variable name
    #     TimeStart: Time index start
    #     TimeCount: Time index count
    # mtchfrm is the data frame with region assignments
    # Total flux is sum of land and ocean, although regions generally delineate these anyway
    # Flux variables should be 3D (lat,lon,time)
  
    nc1 = nc_open(dfrm$FileName)
    lon = ncvar_get(nc1,dfrm$LonVar)
    lat = ncvar_get(nc1,dfrm$LatVar)
    if (netoptn) {
        flxnet = ncvar_get(nc1,"net",start=c(1,1,dfrm$TimeStart),count=c(nlon,nlat,dfrm$TimeCount))
    } else {
        flxocn = ncvar_get(nc1,dfrm$OceanFlux,start=c(1,1,dfrm$TimeStart),count=c(nlon,nlat,dfrm$TimeCount))
        flxlnd = ncvar_get(nc1,dfrm$LandFlux,start=c(1,1,dfrm$TimeStart),count=c(nlon,nlat,dfrm$TimeCount))
        flxnet = flxocn + flxlnd
    }
    nc_close(nc1)
    
    dimnames(flxnet) = list(Longitude = mlon,Latitude = sprintf("Lat%.1f",mlat),Time = sprintf("Time%02d",seq(1,dfrm$TimeCount)))
    flx_tbl = as_tibble(flxnet,rownames = "Longitude")
    flx_tbl = flx_tbl %>% pivot_longer(cols = starts_with("Lat"), names_to = "LatTime", 
                                       names_prefix = "Lat", values_to = "Flux")
    flx_tbl = flx_tbl %>% separate(LatTime, c("Latitude","Time"), sep = "\\.Time", remove=FALSE)
    flx_tbl = flx_tbl %>% mutate_at(c('Longitude','Latitude'), as.numeric) %>% select(-c("LatTime"))

    flxmrg = flx_tbl %>% inner_join(mtchfrm,by=c("Longitude","Latitude"))
    return(flxmrg)
}

mdnms = c('Ames','Baker','CMS-Flux','OU')
lndvrs = c('land','land','land','land')
ocnvrs = c('ocean','ocean','ocean','ocean')
expt = c('LNLG','IS')

# Prior fields
prfrm = data.frame(Model=rep(mdnms,2), LandFlux=rep(lndvrs,2), OceanFlux=rep(ocnvrs,2), Expt=rep(expt,each=4))
prfrm$FileName = sprintf("%s/Prior/%s_gridded_fluxes_Prior.nc4",dtdr,prfrm$Model)
prfrm$TimeStart = 18
prfrm$TimeCount = 3
prfrm$LonVar = "longitude"
prfrm$LatVar = "latitude"

prfnl = prfrm %>% group_by(Model,Expt) %>% group_modify(~ read_mip_subset(.x,mtchfrm=miprgns)) %>% ungroup
prfnl = fklon_tbl(prfnl,lonvar = "Longitude",center=-90)
tmlbls = c("6_Jun","7_Jul","8_Aug")
prfnl$Time = factor(prfnl$Time,levels=c("01","02","03"),labels = tmlbls)

prfnl = prfnl %>% mutate(FluxPlt = case_when(Flux < -2000 ~ -1999.0,
                                             Flux > 2000 ~ 1999.0,
                                             .default = Flux),
                         PrAnn = Flux)
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
frmsb = prfnl %>% filter(Expt == expt[1]) 

flxlb2 = bquote(atop("Flux", "gC" ~ "m" ^ -2 ~ "yr" ^-1))
ttlstr = 'North America JJA 2016: Prior'
gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  facet_grid(Time ~ Model) + 
  scale_fill_gradientn(flxlb2,colors=r7,limits=c(-2000,2000)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle(ttlstr) + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_NorthAmerica_2016_JJA_Prior.pdf")
pdf(pnm,width = 12,height=7)
print(gmn)
dev.off()

# MIP Experiment Results
mipfrm = data.frame(Model=rep(mdnms,2), LandFlux=rep(lndvrs,2), OceanFlux=rep(ocnvrs,2), Expt=rep(expt,each=4))
mipfrm$FileName = sprintf("%s/%s_gridded_fluxes_%s.nc4",dtdr,mipfrm$Model,mipfrm$Expt)
mipfrm$TimeStart = 18
mipfrm$TimeCount = 3
mipfrm$LonVar = "longitude"
mipfrm$LatVar = "latitude"

mipfnl = mipfrm %>% group_by(Model,Expt) %>% group_modify(~ read_mip_subset(.x,mtchfrm=miprgns)) %>% ungroup
mipfnl = fklon_tbl(mipfnl,lonvar = "Longitude",center=-90)
tmlbls = c("6_Jun","7_Jul","8_Aug")
mipfnl$Time = factor(mipfnl$Time,levels=c("01","02","03"),labels = tmlbls)

mipfnl = mipfnl %>% mutate(FluxPlt = case_when(Flux < -2000 ~ -1999.0,
                                               Flux > 2000 ~ 1999.0,
                                               .default = Flux))
for (i in seq(1,length(expt))) {
  frmsb = mipfnl %>% filter(Expt == expt[i]) 
  ttlstr = paste0('North America JJA 2016: ',expt[i])
  gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
    geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
    facet_grid(Time ~ Model) + 
    scale_fill_gradientn(flxlb2,colors=r7,limits=c(-2000,2000)) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
    ggtitle(ttlstr) + theme_mat + coord_equal()
  pnm = paste0(cfglst$output_dir,"/MIP_NorthAmerica_2016_JJA_",expt[i],".pdf")
  pdf(pnm,width = 12,height=7)
  print(gmn)
  dev.off()
}

# Merge and map deviations
mipfnl = mipfnl %>% select(c("Longitude","Latitude","Time","Region","Model","fk360","Expt","Flux"))
prfnl = prfnl %>% select(c("Longitude","Latitude","Time","Region","Model","fk360","Expt","PrAnn"))
mrgall = mipfnl %>% inner_join(prfnl,by=c("Longitude","Latitude","Time","Region","Model","fk360","Expt"))
mrgall = mrgall %>% mutate(PrDev = Flux - PrAnn)
mrgall = mrgall %>% mutate(FluxPlt = case_when(PrDev < -2000 ~ -1999.0,
                                               PrDev > 2000 ~ 1999.0,
                                               .default = PrDev))

for (i in seq(1,length(expt))) {
  frmsb = mrgall[mrgall$Expt == expt[i],]  
  ttlstr = paste0('North America JJA 2016 Post-Prior: ',expt[i])
  gmn = ggplot(frmsb,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
    geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
    facet_grid(Time ~ Model) + 
    scale_fill_gradientn(flxlb2,colors=r7,limits=c(-2000,2000)) + 
    scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
    scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
    ggtitle(ttlstr) + theme_mat + coord_equal()
  pnm = paste0(cfglst$output_dir,"/MIP_NorthAmerica_2016_JJA_PrDev_",expt[i],".pdf")
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

# Get UTM coordinates
cord.dec <- SpatialPoints(cbind(frmsb$Longitude, frmsb$Latitude), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +lon_0=93w +zone=15 +ellps=GRS80"))


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


