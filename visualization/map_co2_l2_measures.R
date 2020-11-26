# Map CO2

library(ncdf4)
library(ggplot2)
library(maps)
library(reshape2)
library(colorspace)
library(viridisLite)

source("mapfunctions.R")

# Orbit 4974

usst = centermap2("States.csv",center=-115,intrarg = FALSE)
lnlb = lonlbs(c(-120,-110,-90),center=-115.0)
ltlb = latlbs(seq(30,50,by=5))

nc1 = nc_open("oco2_LtCO2_150608_B9003r_180928050743s.nc4")
orb1 = ncvar_get(nc1,"Sounding/orbit")
lat1 = ncvar_get(nc1,"latitude")
lon1 = ncvar_get(nc1,"longitude")
sdg1 = ncvar_get(nc1,"sounding_id")
time = ncvar_get(nc1,"time")
vrtlt = ncvar_get(nc1,"vertex_latitude")
vrtln = ncvar_get(nc1,"vertex_longitude")
qflg = ncvar_get(nc1,"xco2_quality_flag")
#print(nc1)
xco2 = ncvar_get(nc1,"xco2")
nc_close(nc1)

sq1 = seq(1,length(orb1))

sborb = sq1[orb1 == 4974]
sborb2 = sq1[(lat1 > 38) & (lat1 < 40) & (lon1 > -130) & (lon1 < -100) & (qflg == 0)]

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 12
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 12
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 14
theme_mat$legend.position = "right"
theme_mat$plot.title$hjust = 0.5
theme_mat$strip.text$size = 12
theme_mat$legend.text$size = 11
theme_mat$legend.title$size = 12
theme_mat$panel.grid.major = element_blank()
theme_mat$panel.grid.minor = element_blank()

dt1 = data.frame(Longitude=lon1[sborb],Latitude=lat1[sborb])

# MEASURES
nc2 = nc_open("OCO2GriddedXCO2_20150608_v1_1588978942.nc")
lonms = ncvar_get(nc2,"longitude")
latms = ncvar_get(nc2,"latitude")
xco2ms = ncvar_get(nc2,"xco2")
nc_close(nc2)

sq2 = seq(1,length(lonms))
sbmeas = sq2[latms > 30 & latms < 60 & lonms > -130 & lonms < -100]

dt3$XCO2[(dt3$XCO2 < 390) & !(is.na(dt3$XCO2))] = 390.01
dt3 = data.frame(Longitude=as.vector(lonms[sbmeas]),Latitude=as.vector(latms[sbmeas]),
                 XCO2=as.vector(xco2[sbmeas]))
dt3 = fklon(dt3,lonvar = "Longitude",center=-115)

r9 = viridis(9)
g2 = ggplot(dt3,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=XCO2)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=usst) +
  scale_fill_gradientn("XCO2",colors=r9,limits=c(390,405)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-14,18)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(32,52)) + 
  ggtitle("Level 3 XCO2 2015-06-08") + theme_mat + coord_equal()

pdfms = "XCO2Map_MEASURES.pdf"
pdf(pdfms,width=8,height=6)
print(g2)
dev.off()


# Back to Level 2
sdstr = sprintf("%.0f",sdg1)
vltsb = vrtlt[,sborb2]
vlnsb = vrtln[,sborb2]

dimnames(vltsb)[2] = list(sdstr[sborb2])
dimnames(vlnsb)[2] = list(sdstr[sborb2])
ltmlt = melt(vltsb,varnames=c("Vertex","SoundStr"),value.name="VertLat")
lnmlt = melt(vlnsb,varnames=c("Vertex","SoundStr"),value.name="VertLon")
vrtmrg = merge(lnmlt,ltmlt)
vrtmrg$SoundStr = sprintf("%.0f",vrtmrg$SoundStr)
vrtmrg = fklon(vrtmrg,lonvar = "VertLon",center=-115)

lnlbl2 = lonlbs(c(-113),center=-115.0)
ltlbl2 = latlbs(seq(38,40,by=1))

dt2 = data.frame(Longitude=as.vector(lon1[sborb2]),Latitude=as.vector(lat1[sborb2]),
                 XCO2=as.vector(xco2[sborb2]),SoundStr=sdstr[sborb2])

dtmrg = merge(dt2,vrtmrg)
dtmrg = dtmrg[order(dtmrg$SoundStr,dtmrg$Vertex),]


glv2 = ggplot(dtmrg,aes(x=fk360,y=VertLat,group = SoundStr,fill=XCO2)) + 
       geom_polygon(aes(fill=XCO2)) +
       scale_fill_gradientn("XCO2",colors=r9,limits=c(390,405)) + 
       scale_x_continuous("",breaks=lnlbl2$fk360,labels=parse(text=lnlbl2$labxpr),limits=c(1.4,2.6)) + 
       scale_y_continuous("",breaks=ltlbl2$origlat,labels=parse(text=ltlbl2$labxpr),limits=c(37.9,40.1)) + 
       coord_equal() + theme_mat + ggtitle("Level 2 XCO2 2015-06-08")
pdfout = "XCO2Map_LiteL2.pdf"
pdf(pdfout,width=5.5,height=6.5)
print(glv2)
dev.off()

