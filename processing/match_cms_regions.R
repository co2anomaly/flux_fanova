# Setup flux region definitions for CMS grid

library(ncdf4)
library(fields)
library(reshape2)

rdist_ref_min = function(tstloc,rffrm,rfvar = "Region") {
    # Find minimum distance
    mtrf = as.matrix(tstloc[,c("Longitude","Latitude")])
    mt1 = as.matrix(rffrm[,c("Longitude","Latitude")])
    dsts = rdist.earth(mtrf,mt1,miles = FALSE)
    
    tfrm = rffrm
    tfrm$Dist = as.vector(dsts)
    
    sbfrm = tfrm[ (tfrm$Dist == min(dsts)), c(rfvar,"Dist") ]
    if (nrow(sbfrm) > 0) {
        sbfrm = sbfrm[1,]
    }
    return(sbfrm)
}

cmslnd = "/Users/jhobbs/Documents/MEASURES_CO2/Data/CMS/MEASURES-GOSAT+OCO2-control-land/2015/03.nc"
nc1 = nc_open(cmslnd)
flxlnd = ncvar_get(nc1,"CO2_Flux")
lon = ncvar_get(nc1,"lon")
lat = ncvar_get(nc1,"lat")
nc_close(nc1)
dimnames(flxlnd)[1] = list(lon)
dimnames(flxlnd)[2] = list(lat)
lndmlt = melt(flxlnd, varnames=c("Longitude","Latitude"),value.name = "LandFlux")

cmsocn = "/Users/jhobbs/Documents/MEASURES_CO2/Data/CMS/Ocean/MEASURES-GOSAT+OCO2-control/2015/03.nc"
nc1 = nc_open(cmsocn)
flxocn = ncvar_get(nc1,"CO2_Flux")
lon = ncvar_get(nc1,"lon")
lat = ncvar_get(nc1,"lat")
nc_close(nc1)
dimnames(flxocn)[1] = list(lon)
dimnames(flxocn)[2] = list(lat)
ocnmlt = melt(flxocn, varnames=c("Longitude","Latitude"),value.name = "OceanFlux")

flxmrg = merge(lndmlt,ocnmlt)

flxbth = flxmrg[(flxmrg$OceanFlux != 0.0) & (flxmrg$LandFlux != 0.0),]

# Get flux regions
nc1 = nc_open("oco2_regions_l4mip_v7.nc")
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

cmsmtch = ddply(flxmrg, c("Longitude","Latitude"), .fun = rdist_ref_min, .progress = "text", rffrm = miprgns)

# Output as CSV
write.csv(cmsmtch,file="CMS/CMS_Region_Match.csv",quote = FALSE,row.names = FALSE)
