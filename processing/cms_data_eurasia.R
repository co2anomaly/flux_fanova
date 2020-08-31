# Map CMS-Flux datasets, DJF and JJA

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("mapfunctions.R")

wrld = centermap2("Coast.csv",center=90.0,intrarg = FALSE)
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

read_cms = function(dfrm,flxvr = "CO2_Flux") {
    # Read CMS Flux results
    flstr = sprintf("%s/%04d/%02d.nc",dfrm$CMSDir,dfrm$Year,dfrm$Month)
    nc1 = nc_open(flstr)
    flx = ncvar_get(nc1,flxvr)
    lon = ncvar_get(nc1,"lon")
    lat = ncvar_get(nc1,"lat")
    nc_close(nc1)
    
    dimnames(flx)[1] = list(lon)
    dimnames(flx)[2] = list(lat)
    flxmlt = melt(flx,varnames=c("Longitude","Latitude"),value.name = "Flux")
    return(flxmlt)
}

# Construct a data frame with structure
cmsdir = c("CMS/MEASURES-GOSAT+OCO2-control-land","CMS/MEASURES-GOSAT+OCO2-fused-land")
dsrc = c("Control","Fused")

d2 = data.frame(Year=rep (rep(c(2014,rep(2015,5)) ),2), Month=rep(c(12,1,2,6,7,8),2),
                CMSDir = rep(cmsdir,each=6), DataSource = rep(dsrc,each=6))
flxall = ddply(d2,c("Year","Month","DataSource"),.fun = read_cms)

flxall = fklon(flxall,lonvar = "Longitude",center=90)
flxall$Season = "DJF"
flxall$Season[(flxall$Month <= 8) & (flxall$Month >= 6)] = "JJA"
flxall$Trt = paste(flxall$Season,flxall$DataSource,sep="_")
flxall$MSeq = 1
flxall$MSeq[(flxall$Season == "DJF") & (flxall$Month < 10)] = flxall$Month[(flxall$Season == "DJF") & (flxall$Month < 10)] + 1
flxall$MSeq[(flxall$Season == "JJA")] = flxall$Month[(flxall$Season == "JJA")] - 5
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)

# Region Masks
rgns = read.csv("CMS/CMS_Region_Match.csv")
flxmrg = merge(flxall,rgns)
flxmrg = flxmrg[(flxmrg$Region == 10) | (flxmrg$Region == 11) | (flxmrg$Region == 16),]

gmn = ggplot(flxmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Flux)) + 
      geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.6, color="#444444") +
      facet_grid(Trt ~ MSeq) + 
      scale_fill_gradientn("Flux",colors=r7,limits=c(-0.2,0.2)) + 
      scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
      scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
      ggtitle("CMS Land Flux") + theme_mat + coord_equal()
pnm = "CMS_Eurasia.pdf"
pdf(pnm,width = 12,height=8 )
print(gmn)
dev.off()

# Prepare Output
dsrc = c("Control","Fused")
ssns = c("DJF","JJA")
frmsb = flxmrg[(flxmrg$DataSource == dsrc[1]) & (flxmrg$Season == ssns[1]) & (flxmrg$MSeq == 1),]
nloc = nrow(frmsb)
nrep = 3

y=array(dim=c(nloc,length(ssns),length(dsrc),nrep))
for(i in seq(1,length(ssns))){
  for(j in seq(1,length(dsrc))){
    for(k in 1:nrep){
      frmsb = flxmrg[(flxmrg$DataSource == dsrc[j]) & (flxmrg$Season == ssns[i]) & (flxmrg$MSeq == k),]
      frmsb = frmsb[order(frmsb$Latitude,frmsb$Longitude),]
      y[,i,j,k] = frmsb$Flux
    }
  }
}

nab = 1
ic=jc=c(-1,1)

dimlc = ncdim_def("location","",seq(1,nrow(frmsb)),create_dimvar = FALSE)
dimlvA = ncdim_def("season_A","",seq(1,length(ssns)),create_dimvar = FALSE)
dimlvB = ncdim_def("data_source_B","",seq(1,length(dsrc)),create_dimvar = FALSE)
dimeffA = ncdim_def("effect_A","",seq(1,1),create_dimvar = FALSE)
dimeffB = ncdim_def("effect_B","",seq(1,1),create_dimvar = FALSE)
dimeffAB = ncdim_def("effect_interact","",seq(1,nab))
dimrep = ncdim_def("replicate","",seq(1,nrep))
dimchr = ncdim_def("char","",1:40,create_dimvar = FALSE)

varx = ncvar_def("longitude","degrees_east",list(dimlc),-9999,
                 longname="Longitude")
vary = ncvar_def("latitude","degrees_north",list(dimlc),-9999,
                 longname="Latitude")
vardt = ncvar_def("CO2_Flux","None",list(dimlc,dimlvA,dimlvB,dimrep),-9999,
                  longname="Simulated data fields")
varctrstA = ncvar_def("contrast_A","None",list(dimlvA),-9999,
                      longname="Main effect A zero sum contrasts")
varctrstB = ncvar_def("contrast_B","None",list(dimlvB),-9999,
                      longname="Main effect B zero sum contrasts")
varA = ncvar_def("season","",list(dimchr,dimlvA),prec="char")
varB = ncvar_def("data_source","",list(dimchr,dimlvB),prec="char")

nc1 = nc_create("../FunctionalANOVA/Example_CMS_Eurasia/CMS_Eurasia_Flux.nc",list(varx,vary,vardt,varctrstA,varctrstB,varA,varB))
ncvar_put(nc1,varx,frmsb$Longitude)
ncatt_put(nc1,varx,"missing_value",-9999.0)
ncvar_put(nc1,vary,frmsb$Latitude)
ncatt_put(nc1,vary,"missing_value",-9999.0)
ncvar_put(nc1,vardt,y)
ncatt_put(nc1,vardt,"missing_value",-9999.0)
ncvar_put(nc1,varctrstA,ic)
ncatt_put(nc1,varctrstA,"missing_value",-9999.0)
ncvar_put(nc1,varctrstB,ic)
ncatt_put(nc1,varctrstB,"missing_value",-9999.0)
ncvar_put(nc1,varA,ssns)
ncvar_put(nc1,varB,dsrc)
nc_close(nc1)


