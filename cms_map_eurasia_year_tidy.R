# Process MIP datasets for functional ANOVA, map monthly values
# CMS-Flux data aggregation example
# Save locations as UTM coordinates

library(ggplot2)
library(ncdf4)
library(colorspace)
library(dplyr)
library(tidyr)
library(sp)

source("mapfunctions_tidy.R")
source("func_anova_tidy.R")

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
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)


# Read CMS data directly
nc1 = nc_open(cfglst$data_file)
flx = ncvar_get(nc1,"CO2_Flux")
flxpr = ncvar_get(nc1,"CO2_Flux_Prior")
lon = ncvar_get(nc1,"longitude")
lat = ncvar_get(nc1,"latitude")
yrA = ncvar_get(nc1,"year")
dtsrcB = ncvar_get(nc1,"data_source")

nc_close(nc1)


nloc = length(lat)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = lon,Latitude = lat)


# Prior fields, collapse data source dimension
prarr = flxpr[,,1,]
dimnames(prarr) = list(LocID = seq(1,nloc),YearID = sprintf("Year%d",yrA),Time = sprintf("Time%02d",seq(1,3)))
pr_tbl = as_tibble(prarr,rownames = "LocID")
pr_tbl = pr_tbl %>% pivot_longer(cols = starts_with("Year"), names_to = "YearTime", 
                                   names_prefix = "Year", values_to = "Flux")
pr_tbl = pr_tbl %>% separate(YearTime, c("Year","Time"), sep = "\\.Time", remove=FALSE)
pr_tbl = pr_tbl %>% mutate_at(c('LocID','Time'), as.integer) %>% select(-c("YearTime"))
pr_tbl$MFct = factor(pr_tbl$Time,levels = 1:3,labels=c("Jun","Jul","Aug"))

prfnl = pr_tbl %>% left_join(lcfrm, by="LocID")
prfnl = prfnl %>% mutate(FluxPlt = case_when(Flux < -2000 ~ -1999.0,
                                             Flux > 2000 ~ 1999.0,
                                             .default = Flux),
                         PrAnn = Flux)
prfnl = fklon_tbl(prfnl,lonvar = "Longitude",center=90)

flxlb2 = bquote(atop("Flux", "gC" ~ "m" ^ -2 ~ "yr" ^-1))
ttlstr = 'CMS Land Prior Flux'
gmn = ggplot(prfnl,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  facet_grid(Year ~ MFct) + 
  scale_fill_gradientn(flxlb2,colors=r7,limits=c(-2000,2000)) +  
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle(ttlstr) + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_Eurasia_Prior_JJA_2Yr.pdf")
pdf(pnm,width = 12,height=4 )
print(gmn)
dev.off()

# CMS Experiment Results
dimnames(flx) = list(LocID = seq(1,nloc),DatSrc = sprintf("Src%s",dtsrcB),
                     YearID = sprintf("Year%d",yrA),Time = sprintf("Time%02d",seq(1,3)))
cms_tbl = as_tibble(flx,rownames = "LocID")
cms_tbl = cms_tbl %>% pivot_longer(cols = starts_with("Src"), names_to = "SrcYearTime", 
                                   names_prefix = "Src", values_to = "Flux")
cms_tbl = cms_tbl %>% separate(SrcYearTime, c("DatSrc","YearTime"), sep = "\\.Year", remove=FALSE)
cms_tbl = cms_tbl %>% separate(YearTime, c("Year","Time"), sep = "\\.Time", remove=FALSE)
cms_tbl = cms_tbl %>% mutate_at(c('LocID','Time'), as.integer) %>% select(-c("YearTime","SrcYearTime"))
cms_tbl$MFct = factor(cms_tbl$Time,levels = 1:3,labels=c("Jun","Jul","Aug"))

flxmrg = cms_tbl %>% left_join(lcfrm, by="LocID")
flxmrg = flxmrg %>% mutate(FluxPlt = case_when(Flux < -2000 ~ -1999.0,
                                             Flux > 2000 ~ 1999.0,
                                             .default = Flux))
flxmrg = fklon_tbl(flxmrg,lonvar = "Longitude",center=90)
flxmrg$Trt = paste(flxmrg$Year,flxmrg$DatSrc,sep="_")

gmn = ggplot(flxmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=FluxPlt)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, linewidth=0.8, color="#444444") +
  facet_grid(Trt ~ MFct) + 
  scale_fill_gradientn(flxlb2,colors=r7,limits=c(-2000,2000)) +  
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-110,110)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(0,84)) + 
  ggtitle("CMS Land Flux") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/CMS_Eurasia_JJA_2Yr.pdf")
pdf(pnm,width = 12,height=8 )
print(gmn)
dev.off()


