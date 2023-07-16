library(ggplot2)
library(maps)
library(purrr)

valmn = function(dtfrm) {
  meanout = mean(dtfrm$value,na.rm=TRUE)
  medout = median(dtfrm$value,na.rm=TRUE)
  frameout = data.frame(mean=meanout,median=medout)
  return(frameout)
}

chckgrp2 = function(mpfrm) {
  mpfrm = mpfrm[order(mpfrm$plyord),]
  mpfrm$dif1 = c(NA,diff(mpfrm$fk360))
  return(mpfrm)
}

grpct = function(dfrm) {
	frmout = data.frame(ctgrp = length(unique(dfrm$group2)))
	return(frmout)
}

spltna2 = function(dfrm) {
  # Split at NA's in the difference variable
  sq1 = seq(1,length(dfrm$dif1))	
  tsq = sq1[is.na(dfrm$dif1)]
  if (length(tsq) > 1) {
    tls = c(tsq[seq(2,length(tsq))] - 1, length(sq1))
  }
  else {
    tls = length(sq1)
  }
  dfrm$group2 = 1
  for (i in seq(1,length(tsq))) {
    dfrm$group2[tsq[i]:tls[i]] = i		
  }
  return(dfrm)
}


centermap2 = function(flnm,center=180,intrarg=TRUE) {
  # Make a map centered on a chosen longitude, using map_data
  frm1 = read.csv(flnm,header = TRUE)
  frm1$fk360 = (frm1$X - center) %% 360.0
  frm1$fk360[frm1$fk360 > 180] = frm1$fk360[frm1$fk360 > 180] - 360
  
  frm2 = frm1 %>% group_by(plyidx) %>% group_modify(~ chckgrp2(.x))
  frm2$dif1[frm2$dif1 < -100 & !is.na(frm2$dif1)] = NA
  frm2$dif1[frm2$dif1 > 100 & !is.na(frm2$dif1)] = NA
  frm3 = frm2 %>% group_modify(~ spltna2(.x)) %>% ungroup()
  frm3$grp2chr = paste(frm3$plyidx,frm3$group2,sep="_")
  return(frm3)
}

lonlbs = function(lonbrks,center=180) {
	# Return nicely formatted labels, with arbitrary map center
	# Assume lonbrks are (-180,180) or (0,360)
	frm1 = data.frame(origlon=lonbrks)
	frm1$fk360 = (frm1$origlon - center) %% 360.0
	frm1$fk360[frm1$fk360 > 180] = frm1$fk360[frm1$fk360 > 180] - 360
    frm1$labnum = frm1$origlon
    frm1$labnum[frm1$origlon < 0] = abs(frm1$origlon[frm1$origlon < 0])
    frm1$labnum[frm1$origlon > 180] = frm1$origlon[frm1$origlon > 180] - 360
    frm1$labdir = "* W"
    frm1$labdir[frm1$origlon == 0 | frm1$origlon == 180.0 | frm1$origlon == -180.0] = ""
    frm1$labdir[frm1$origlon > 0 & frm1$origlon < 180.0] = "* E"
    frm1$lablab = paste(sprintf("%d",frm1$labnum),"* degree",frm1$labdir)
    frm1$labxpr = as.expression(frm1$lablab)
    frm1 = frm1[order(frm1$fk360),]
    return(frm1)
    
}

latlbs = function(latbrks) {
	# Return nicely formatted labels
	# Assume latbrks are (-90,90)
	frm1 = data.frame(origlat=latbrks)
    frm1$labnum = frm1$origlat
    frm1$labnum[frm1$origlat < 0] = abs(frm1$origlat[frm1$origlat < 0])
    frm1$labdir = "* N"
    frm1$labdir[frm1$origlat == 0] = ""
    frm1$labdir[frm1$origlat < 0] = "* S"
    frm1$lablab = paste(sprintf("%d",frm1$labnum),"* degree",frm1$labdir)
    frm1$labxpr = as.expression(frm1$lablab)
    frm1 = frm1[order(frm1$origlat),]
    return(frm1)
}

fklon_tbl = function(dtbl,lonvar,center=180) {
    tblout = dtbl %>% mutate(tmp360 = (eval(parse(text=lonvar)) - center) %% 360.0  )
    tblout = tblout %>% mutate(fk360 = case_when(tmp360 > 180 ~ tmp360 - 360,
                                                 tmp360 <= 180 ~ tmp360))
    tblout = tblout %>% select(-c(tmp360))                            
    return(tblout)
}

fklon = function(dfrm,lonvar,center=180) {
	# Convert a data frame with a longitude coordinate to a system centered at center
	dfrm$fk360 = (dfrm[,lonvar] - center) %% 360.0
	dfrm$fk360[dfrm$fk360 > 180] = dfrm$fk360[dfrm$fk360 > 180] - 360
    return(dfrm)
}
