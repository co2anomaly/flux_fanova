# Map CMS-Flux datasets, DJF and JJA

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("mapfunctions.R")

wrld = centermap2("Coast.csv",center=0.0,intrarg = FALSE)
lnlb = lonlbs(seq(-180,120,by=60),center=0)
ltlb = latlbs(seq(-90,90,by=30))

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
cmsdir = c("CMS/Ocean/MEASURES-GOSAT+OCO2-control","CMS/Ocean/MEASURES-GOSAT+OCO2-fused")
dsrc = c("Control","Fused")

d2 = data.frame(Year=rep (rep(c(2014,rep(2015,5)) ),2), Month=rep(c(12,1,2,6,7,8),2),
                CMSDir = rep(cmsdir,each=6), DataSource = rep(dsrc,each=6))
flxall = ddply(d2,c("Year","Month","DataSource"),.fun = read_cms)

flxall = fklon(flxall,lonvar = "Longitude",center=0)
flxall$Season = "DJF"
flxall$Season[(flxall$Month <= 8) & (flxall$Month >= 6)] = "JJA"
flxall$Trt = paste(flxall$Season,flxall$DataSource,sep="_")
flxall$MSeq = 1
flxall$MSeq[(flxall$Season == "DJF") & (flxall$Month < 10)] = flxall$Month[(flxall$Season == "DJF") & (flxall$Month < 10)] + 1
flxall$MSeq[(flxall$Season == "JJA")] = flxall$Month[(flxall$Season == "JJA")] - 5
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.5)

gmn = ggplot(flxall,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Flux)) + 
      geom_path(aes(x=fk360,y=Y,group=group2), data=wrld) +
      facet_grid(Trt ~ MSeq) + 
      scale_fill_gradientn("Flux",colors=r7,limits=c(-0.01,0.01)) + 
      scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-180,180)) + 
      scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(-90,90)) + 
      ggtitle("CMS Ocean Flux") + theme_mat + coord_equal()
pnm = "CMS_Ocean_Season.pdf"
pdf(pnm,width = 12,height=10)
print(gmn)
dev.off()


