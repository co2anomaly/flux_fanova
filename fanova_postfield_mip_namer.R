# Plot MIP posterior fields (or chunks)

library(ggplot2)
library(reshape2)
library(ncdf4)
library(colorspace)
library(plyr)

source("mapfunctions.R")

cfgfile = "config/fanova_mip_namer_prdev.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

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

# Read MIP Locs, metadata
nc1 = nc_open(cfglst$data_file)
latmip = as.vector(ncvar_get(nc1,"latitude"))
lonmip = as.vector(ncvar_get(nc1,"longitude"))
mdlA = ncvar_get(nc1,"model")
dtsrcB = ncvar_get(nc1,"data_source")
ctrA = ncvar_get(nc1,"contrast_A")
ctrB = ncvar_get(nc1,"contrast_B")
nc_close(nc1)

if (is.na(ncol(ctrA))) {
  ctrA = matrix(ctrA,ncol=1)
}
if (is.na(ncol(ctrB))) {
  ctrB = matrix(ctrB,ncol=1)
}



nloc = length(latmip)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = lonmip,Latitude = latmip)

# Read posterior results
fctnms = c("mean","mainA","mainB","interact")
fctlbs = c("Mean","Season","DataSrc","Interact")
fctord = seq(1,4)
fctfrm = data.frame(ANOVACmp=fctord,ANOVALbl=fctlbs)

itst = 1
niter = 5000
nchn =  4
neffA = ncol(ctrA)
neffB = ncol(ctrB)
pstmu = array(0,c(niter,nloc,nchn))
pstA = array(0,c(niter,nloc,neffA,nchn))
pstB = array(0,c(niter,nloc,neffB,nchn))

for (j in seq(1,nchn)) {
    
    pstnm = paste0(cfglst$post_samp_file,j,".nc")
    ncpst = nc_open(pstnm)
    pstmu[,,j] = ncvar_get(ncpst,"gp_field_mean",start=c(itst,1),count=c(niter,nloc))
    pstA[,,,j] = ncvar_get(ncpst,"gp_field_mainA",start=c(itst,1,1),count=c(niter,nloc,neffA))
    pstB[,,,j] = ncvar_get(ncpst,"gp_field_mainB",start=c(itst,1,1),count=c(niter,nloc,neffB))
    nc_close(ncpst)
  
}

pstfldmu = apply(pstmu,2,FUN = mean)
pstfldsd = apply(pstmu,2,FUN = sd)
mufrm = data.frame(LocID = seq(1,nloc), Mu=as.vector(pstfldmu), StdDev=as.vector(pstfldsd))

mumrg = merge(mufrm,lcfrm)
mumrg = fklon(mumrg,lonvar = "Longitude",center=-90)
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
std9 = sequential_hcl(9,h=340,c.=c(0,80),l=c(100,50)) 

gmn = ggplot(mumrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Mu)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-800,800))  + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_Mu.pdf")
pdf(pnm,width = 12,height=8)
print(gmn)
dev.off()

gstd = ggplot(mumrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=StdDev)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Std Dev",colors=std9,limits=c(65,75))  + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean Posterior Std Dev") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_Mu_StdDev.pdf")
pdf(pnm,width = 12,height=8)
print(gstd)
dev.off()



# A Effect
pstfldA = apply(pstA,2:3,FUN = mean)
pstAsm = pstfldA %*% t(ctrA)
dimnames(pstAsm)[2] = list(as.vector(mdlA))
amlt = melt(pstAsm,varnames = c("LocID","FactorA"))
amrg = merge(amlt,lcfrm)
amrg = fklon(amrg,lonvar = "Longitude",center=-90)
geffA = ggplot(amrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-60,60))  + 
  facet_wrap( ~ FactorA) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Model Effect") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_AEff.pdf")
pdf(pnm,width = 12,height=8)
print(geffA)
dev.off()

# B Effect
pstfldB = apply(pstB,2:3,FUN = mean)
pstBsm = pstfldB %*% t(ctrB)
dimnames(pstBsm)[2] = list(as.vector(dtsrcB))
bmlt = melt(pstBsm,varnames = c("LocID","FactorB"))
bmrg = merge(bmlt,lcfrm)
bmrg = fklon(bmrg,lonvar = "Longitude",center=-90)
geffB = ggplot(bmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-60,60))  + 
  facet_wrap( ~ FactorB,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Data Source Effect") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_BEff.pdf")
pdf(pnm,width = 12,height=6)
print(geffB)
dev.off()


# Credible intervals - Overall mean
pstqmu = apply(pstmu,2,FUN = quantile,probs=c(0.025,0.975))
muqmlt = melt(pstqmu,varnames = c("Pctile","LocID"))
muqmrg = merge(muqmlt,lcfrm)
muqmrg = fklon(muqmrg,lonvar = "Longitude",center=-90)
gmuq = ggplot(muqmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-1000,1000))  + 
  facet_wrap( ~ Pctile,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean Interval") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_Mu_Cred.pdf")
pdf(pnm,width = 12,height=6)
print(gmuq)
dev.off()


# Transport contrast
#pstfldA = apply(pstA,2:3,FUN = mean)
#pstAsm = pstfldA %*% t(ctrA)
ctr1_4 = c(0.5, -0.5, 0.5, -0.5)

#ctr1_4 = c(1, 0, -1, 0)
ctr1_3 = t(ctr1_4) %*% ctrA

pstctr = ctr1_3[1] * pstA[,,1,] + ctr1_3[2] * pstA[,,2,] + ctr1_3[3] * pstA[,,3,] 
pstqctr = apply(pstctr,2,FUN = quantile,probs=c(0.025,0.975))
ctrqmlt = melt(pstqctr,varnames = c("Pctile","LocID"))
ctrqmrg = merge(ctrqmlt,lcfrm)
ctrqmrg = fklon(ctrqmrg,lonvar = "Longitude",center=-90)
gctrq = ggplot(ctrqmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-60,60))  + 
  facet_wrap( ~ Pctile,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Model Contrast Interval") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_AEff_Contrast.pdf")
pdf(pnm,width = 12,height=6)
print(gctrq)
dev.off()


# 
#ctr2_2 = c(-1.0, 1.0)
#ctr2_1 = t(ctr2_2) %*% ctrB
pstqctrb = apply(pstB[,,1,],2,FUN = quantile,probs=c(0.025,0.975))
ctrqmlt = melt(pstqctrb,varnames = c("Pctile","LocID"))
ctrqmrg = merge(ctrqmlt,lcfrm)
ctrqmrg = fklon(ctrqmrg,lonvar = "Longitude",center=-90)
gctrqb = ggplot(ctrqmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=value)) + 
  geom_path(aes(x=fk360,y=Y,group=group2), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-100,100))  + 
  facet_wrap( ~ Pctile,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Data Source Interval") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_PrDev_BEff_Contrast.pdf")
pdf(pnm,width = 12,height=6)
print(gctrqb)
dev.off()

