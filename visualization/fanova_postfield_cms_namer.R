# Plot CMS posterior fields

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("../Data/mapfunctions.R")

wrld = centermap2("../Data/Coast.csv",center=-90.0,intrarg = FALSE)
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

# Read CMS Locs
nc1 = nc_open("Example_CMS_NAmer/CMS_NAmer_Flux.nc")
latcms = as.vector(ncvar_get(nc1,"latitude"))
loncms = as.vector(ncvar_get(nc1,"longitude"))
nc_close(nc1)

nloc = length(latcms)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = loncms,Latitude = latcms)

# Read posterior results
fctnms = c("mean","mainA","mainB","interact")
fctlbs = c("Mean","Season","DataSrc","Interact")
fctord = seq(1,4)
fctfrm = data.frame(ANOVACmp=fctord,ANOVALbl=fctlbs)

niter = 5000
nchn =  4
pstflds = array(0,c(niter,nloc,4,4))

for (j in seq(1,nchn)) {
    
    pstnm = paste0("Example_CMS_NAmer/NAmer_PostSamp",j,".nc")
    ncpst = nc_open(pstnm)
    for (k in seq(1,4)) {
        gpnm = paste0("gp_field_",fctnms[k])
        pstflds[,,j,k] = ncvar_get(ncpst,gpnm)
    }
    nc_close(ncpst)
  
}

pstfldmn = apply(pstflds,c(2,4),FUN = mean)
pstmlt = melt(pstfldmn,varnames = c("LocID","ANOVACmp"))

pstmrg = merge(pstmlt,fctfrm)
pstmrg = merge(pstmrg,lcfrm)
pstmrg = fklon(pstmrg,lonvar = "Longitude",center=-90)
pstmrg$ANOVACmp = factor(pstmrg$ANOVACmp,labels = fctfrm$ANOVALbl)
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)

gmn = ggplot(pstmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  facet_wrap(~ ANOVACmp) + 
  scale_fill_gradientn("Flux",colors=r7,limits=c(-0.05,0.05)) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("CMS Flux ANOVA") + theme_mat + coord_equal()
pnm = "CMS_ANOVA_Posterior_NAmer.pdf"
pdf(pnm,width = 12,height=8)
print(gmn)
dev.off()

