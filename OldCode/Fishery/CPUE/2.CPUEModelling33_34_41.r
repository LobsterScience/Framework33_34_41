###########################
#start with an overview of all fishing year by year
require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(devtools)
require(sf)
require(ggplot2)
require(mgcv)
require(ggeffects)
require(ggforce)
fd = file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUE')
setwd(fd)

gr = readRDS('C:/Users/Cooka/Documents/git/bio.lobster.data/mapping_data/GridPolyLand.rds')
st_crs(gr) <- 4326
gr = st_transform(gr,32620) 
st_geometry(gr) <- st_geometry(st_as_sf(gr$geometry/1000)) 
st_crs(gr) <- 32620


###beginning models

aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2024 & LFA %in% c(33,34))

a4 = lobster.db('process.logs.41')

aa = split(aT,f=list(aT$LFA,aT$SYEAR))
oo = list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  first.day<-min(tmp$DATE_FISHED)
  tmp$time<-julian(tmp$DATE_FISHED,origin=first.day-1)
  oo[[i]] = tmp
}

aT = do.call(rbind,oo)

aT$fYear = as.factor(aT$SYEAR)
aT$leffort = log(aT$NUM_OF_TRAPS)


l33 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(aT,LFA==33),family = Gamma(link='log'),method='REML')
l33a = gam(WEIGHT_KG~s(time)+fYear+offset(leffort),data=subset(aT,LFA==33),family = Gamma(link='log'),method='REML')
l33b = gam(WEIGHT_KG~s(time)+s(time,by=fYear)+fYear+offset(leffort),data=subset(aT,LFA==33),family = Gamma(link='log'),method='REML')

saveRDS(list(l33,l33a,l33b),'first3CPUEmodels33.rds')
w = readRDS('first3CPUEmodels33.rds')
l35 = w[[1]]
l35a = w[[2]]
l35b = w[[3]]

aT$rLID = as.factor(aT$LICENCE_ID)
l33c = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=subset(aT,LFA==33),family = Gamma(link='log'),method='REML') #random intercept 
saveRDS(l33c,'CPUEmodel_LicenceRE33.rds')
l33c = readRDS('CPUEmodel_LicenceRE33.rds')


aT$rV = as.factor(aT$VR_NUMBER)
l33d = gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=subset(aT,LFA==33),family = Gamma(link='log'),method='REML') # random intercept
saveRDS(l33d,'CPUEmodel_vesselRE33.rds')
#l33d = readRDS('CPUEmodel_vesselRE33.rds')

#################trip weighted marginal means

ind = aggregate(SD_LOG_ID~time+SYEAR,data=subset(aT,LFA ==33),FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
names(ind1)[2] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fYear=as.factor(ind$SYEAR)


b = emmeans::emmeans(l33,~fYear,offset=0,type='response',data=subset(aT,LFA==33))
bmm = as.data.frame(summary(b))[c('response', 'SE')]
bmm$Model = 'Base'
bmm$SYEAR = 2006:2023
names(bmm)[1:2]=c('wemm','wse')

d = emmeans::emmeans(l33a,~fYear+time,offset=0,type='response',data=subset(aT,LFA==33),at=list(time=c(1:193)))

dmm = as.data.frame(summary(d))
dmm = merge(dmm,ind)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BD'
dmme=dom


f = emmeans::emmeans(l33b,~fYear+time,offset=0,type='response',data=subset(aT,LFA==33),at=list(time=c(1:193)))
dmm = as.data.frame(summary(f))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxD'
fmm=dom


g = emmeans::emmeans(l33c,~fYear+time,offset=0,type='response',data=subset(aT,LFA==33),at=list(time=c(1:193)))
dmm = as.data.frame(summary(g))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxDRL'
gmm=dom

h = emmeans::emmeans(l35d,~fYear+time,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
dmm = as.data.frame(summary(h))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fYear))
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,wse=NA,wemm=NA)
for(i in 1:length(ii)){
  z = subset(dmm,Year==ii[i])
  usm = sum(z$SD_LOG_ID)
  dom[i,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[i,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[i,'SYEAR'] = ii[i]
}
dom$Model = 'BxDRV'
hmm=dom

vb = readRDS(file='unBIASED_CPUE.rds')
unb = subset(vb[[1]],LFA==35,select=c(unBCPUE, unBVar,SYEAR))
unb$unBVar = sqrt(unb$unBVar)
names(unb)=c('wemm','wse','SYEAR')
unb$Model='unbiased'
oo = do.call(rbind,list(unb,bmm,dmme,fmm,gmm,hmm))
write.csv(oo,'LFA35MultipleModelsCPUE.csv')
#marginal mean, internally consistent

ggplot(oo, aes(x = SYEAR, y = wemm ,colour=Model,fill=Model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = wemm-wse, ymax = wemm+wse), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
##################################
# with temperature

#glorys from git/bio.lobster.glorys/inst/Analysis/GeneratingWeeklyClimatologiesLobGrids.r
cpy = lobster.db('process.logs')
cpy$woy = lubridate::week(cpy$DATE_FISHED)
cpy$yr = lubridate::year(cpy$DATE_FISHED)
cpy$GRID_NO = cpy$GRID_NUM
reag = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+GRID_NO+WOS+SYEAR+yr+woy,data=cpy,FUN=sum)
grb = readRDS(file.path(project.datadirectory('bio.lobster.glorys'),'WeeklyTemp_by_grid_05-24_withAnom.rds'))
grbu=grb
grbu$geometry <- NULL
grbu = grbu[!duplicated(grbu),]
rea = merge(reag,grbu,by=c('LFA','woy','yr','GRID_NO'))

reag = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+GRID_NO+Anomaly+bottomT[,3]+WOS+SYEAR,data=rea,FUN=sum)
