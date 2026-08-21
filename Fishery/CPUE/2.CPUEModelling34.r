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

require(sdmTMB)
require(purrr)
la()

fd = file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUE')
setwd(fd)

m = readRDS(file='data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
ca = m[[3]]

#aT = lobster.db('process.logs')
ca = subset(ca,SYEAR>2005 & SYEAR<2026 & LFA %in% c(33,34))
ca$mn = lubridate::month(ca$DATE_FISHED)
ca$SYEAR = ifelse(ca$mn %in% c(11,12), ca$yr+1, ca$yr)
ca$fYear = as.factor(ca$SYEAR)
ca$leffort = log(ca$NUM_OF_TRAPS)


l34 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')
l34a = gam(WEIGHT_KG~s(DOS)+fYear+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')
l34b = gam(WEIGHT_KG~s(DOS)+s(DOS,by=fYear)+fYear+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')
l34c = gam(WEIGHT_KG~s(DOS)+s(DOS,by=fYear)+fYear+s(bcT)+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')

saveRDS(list(l34,l34a,l34b,l34c),'first4CPUEmodels34.rds')
w = readRDS('first4CPUEmodels34.rds')
l34 = w[[1]]
l34a = w[[2]]
l34b = w[[3]]
l34c = w[[4]]


#################trip weighted marginal means
aT=ca
ind = aggregate(SD_LOG_ID~DOS+SYEAR,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
names(ind1)[2] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fYear=as.factor(ind$SYEAR)

b = emmeans::emmeans(l34,~fYear,offset=0,type='response',data=subset(as.data.frame(aT),LFA==34))
bmm = as.data.frame(summary(b))[c('response', 'SE')]
bmm$Model = 'Base'
bmm$SYEAR = 2006:2025
names(bmm)[1:2]=c('wemm','wse')

d = emmeans::emmeans(l34a,~fYear+DOS,offset=0,type='response',data=subset(aT,LFA==34),at=list(DOS=c(1:193)))

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


f = emmeans::emmeans(l34b,~fYear+DOS,offset=0,type='response',data=subset(aT,LFA==34),at=list(DOS=c(1:193)))
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


vb = readRDS(file='unBIASED_CPUE.rds')
unb = subset(vb[[1]],LFA==34,select=c(unBCPUE, unBVar,SYEAR))
unb$unBVar = sqrt(unb$unBVar)
names(unb)=c('wemm','wse','SYEAR')
unb$Model='unbiased'
oo = do.call(rbind,list(unb,bmm,dmme,fmm))


##l34c

tmps = round(as.numeric(with(subset(aT,LFA==34), quantile(bcT,na.rm=T,c(0.025,seq(0.1,0.9,by=.2),0.975)))),1)
f = emmeans::emmeans(l34c,~fYear+DOS+bcT,offset=0,type='response',data=subset(aT,LFA==34),at=list(DOS=seq(1,193,by=3),bcT=tmps))
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
dom$Model = 'BxDxt'
gmm=dom

oo = do.call(rbind,list(gmm,oo))


write.csv(oo,'LFA34MultipleModelsCPUE.csv')
#marginal mean, internally consistent

ggplot(oo, aes(x = SYEAR, y = wemm ,colour=Model,fill=Model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = wemm-wse, ymax = wemm+wse), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)

