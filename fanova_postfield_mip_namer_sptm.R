# Plot MIP posterior fields (or chunks)

library(ggplot2)
library(ncdf4)
library(colorspace)
library(dplyr)
library(tidyr)

source("mapfunctions_tidy.R")

cfgfile = "config/fanova_mip_namer_prdev_sptm.csv"
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
mufrm = tibble(LocID = seq(1,nloc), Mu=as.vector(pstfldmu), StdDev=as.vector(pstfldsd))

mumrg = mufrm %>% left_join(lcfrm, by="LocID")
mumrg = fklon_tbl(mumrg,lonvar = "Longitude",center=-90)
r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=0.3)
std9 = sequential_hcl(9,h=340,c.=c(0,80),l=c(100,50)) 

gmn = ggplot(mumrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Mu)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-800,800))  + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_SpTm_Mu.pdf")
pdf(pnm,width = 12,height=8)
print(gmn)
dev.off()

gstd = ggplot(mumrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=StdDev)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Std Dev",colors=std9,limits=c(65,75))  + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean Posterior Std Dev") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_SpTm_Mu_StdDev.pdf")
pdf(pnm,width = 12,height=8)
print(gstd)
dev.off()



# A Effect
pstfldA = apply(pstA,2:3,FUN = mean)
pstAsm = pstfldA %*% t(ctrA)

dimnames(pstAsm) = list(LocID = seq(1,nloc),
                          ModelID = sprintf("Model%s",mdlA) )
afrm = as_tibble(pstAsm, rownames="LocID") 
afrm = afrm %>% pivot_longer(cols = starts_with("Model"), names_to = "Model",
                                 names_prefix="Model", values_to = "PostMean")
afrm = afrm %>% mutate_at(c('LocID'), as.numeric)

amrg = afrm %>% left_join(lcfrm, by="LocID")
amrg = fklon_tbl(amrg,lonvar = "Longitude",center=-90)
geffA = ggplot(amrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PostMean)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-60,60))  + 
  facet_wrap( ~ Model) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Model Effect") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_SpTm_AEff.pdf")
pdf(pnm,width = 12,height=8)
print(geffA)
dev.off()

# B Effect
pstfldB = apply(pstB,2:3,FUN = mean)
pstBsm = pstfldB %*% t(ctrB)

dimnames(pstBsm) = list(LocID = seq(1,nloc),
                        ModelID = sprintf("Src%s",dtsrcB) )
bfrm = as_tibble(pstBsm, rownames="LocID") 
bfrm = bfrm %>% pivot_longer(cols = starts_with("Src"), names_to = "DatSrc",
                             names_prefix="Src", values_to = "PostMean")
bfrm = bfrm %>% mutate_at(c('LocID'), as.numeric)

bmrg = bfrm %>% left_join(lcfrm, by="LocID")
bmrg = fklon_tbl(bmrg,lonvar = "Longitude",center=-90)
geffB = ggplot(bmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=PostMean)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-60,60))  + 
  facet_wrap( ~ DatSrc,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Data Source Effect") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_SpTm_BEff.pdf")
pdf(pnm,width = 12,height=6)
print(geffB)
dev.off()


# Credible intervals - Overall mean
pstqmu = apply(pstmu,2,FUN = quantile,probs=c(0.025,0.975))
qchr = dimnames(pstqmu)[[1]]
dimnames(pstqmu) = list(Pctile = sprintf("Q%s",qchr), LocID = seq(1,nloc) )

muqfrm = as_tibble(t(pstqmu), rownames="LocID") 
muqfrm = muqfrm %>% pivot_longer(cols = starts_with("Q"), names_to = "Pctile",
                                 names_prefix="Q", values_to = "Flux")
muqfrm = muqfrm %>% mutate_at(c('LocID'), as.numeric)

muqmrg = muqfrm %>% left_join(lcfrm, by="LocID")
muqmrg = fklon_tbl(muqmrg,lonvar = "Longitude",center=-90)
gmuq = ggplot(muqmrg,aes(x=fk360,y=Latitude)) + geom_tile(aes(fill=Flux)) + 
  geom_path(aes(x=fk360,y=Y,group=grp2chr), data=wrld, size=0.8, color="#444444") +
  scale_fill_gradientn("Flux",colors=r7,limits=c(-1000,1000))  + 
  facet_wrap( ~ Pctile,nrow=1) + 
  scale_x_continuous("",breaks=lnlb$fk360,labels=parse(text=lnlb$labxpr),limits=c(-80,40)) + 
  scale_y_continuous("",breaks=ltlb$origlat,labels=parse(text=ltlb$labxpr),limits=c(10,84)) + 
  ggtitle("MIP Flux ANOVA Overall Mean Interval") + theme_mat + coord_equal()
pnm = paste0(cfglst$output_dir,"/MIP_ANOVA_Posterior_NAmer_SpTm_Mu_Cred.pdf")
pdf(pnm,width = 12,height=6)
print(gmuq)
dev.off()

