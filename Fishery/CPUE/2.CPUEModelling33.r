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
ca = subset(ca,SYEAR>2005 & SYEAR<2025 & LFA %in% c(33,34))
ca$mn = lubridate::month(ca$DATE_FISHED)
ca$SYEAR = ifelse(ca$mn %in% c(11,12), ca$yr+1, ca$yr)
ca$fYear = as.factor(ca$SYEAR)
ca$leffort = log(ca$NUM_OF_TRAPS)


l33 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca,LFA==33),family = Gamma(link='log'),method='REML')
l33a = gam(WEIGHT_KG~s(time)+fYear+offset(leffort),data=subset(ca,LFA==33),family = Gamma(link='log'),method='REML')
l33b = gam(WEIGHT_KG~s(time)+s(time,by=fYear)+fYear+offset(leffort),data=subset(ca,LFA==33),family = Gamma(link='log'),method='REML')

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


# #plot index
# #these are conditional
# pr1 = ggpredict(l35,terms='fYear',condition = c(leffort=0))
# pr2 = ggpredict(l35a,terms=c('fYear'),condition = c(leffort=0, time=10))
# pr3 = ggpredict(l35b,terms='fYear',condition = c(leffort=0,time=10))
# pr4 = ggpredict(l35c,terms='fYear',condition = c(leffort=0,time=10))
# pr5 = ggpredict(l35d,terms='fYear',condition = c(leffort=0,time=10))
# 
# pr1t = data.frame(x=pr1$x,predicted=pr1$predicted,Model='Base',conf.low=pr1$conf.low,conf.high=pr1$conf.high,se = pr1$std.error)
# pr2t = data.frame(x=pr2$x,predicted=pr2$predicted,Model='BD',conf.low=pr2$conf.low,conf.high=pr2$conf.high,se = pr2$std.error)
# pr3t = data.frame(x=pr3$x,predicted=pr3$predicted,Model='BxD',conf.low=pr3$conf.low,conf.high=pr3$conf.high,se = pr3$std.error)
# pr4t = data.frame(x=pr4$x,predicted=pr4$predicted,Model='BxDRL',conf.low=pr4$conf.low,conf.high=pr4$conf.high,se = pr4$std.error)
# pr5t = data.frame(x=pr5$x,predicted=pr5$predicted,Model='BxDRV',conf.low=pr5$conf.low,conf.high=pr5$conf.high,se = pr5$std.error)
# 
# pp = rbind(pr5t, pr4t,pr3t,pr2t,pr1t)
# pp$x = as.numeric(as.character(pp$x))
# 
# ggplot(pp, aes(x = x, y = predicted,colour=Model,fill=Model)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Standardized CPUE')+
#   theme_test(base_size = 14)
# 
# 
# # ##these are marginal means and have equal weighting across all days of the season... which is why they are lower
# # require(ggeffects) #new as of 
# b = emmeans::emmeans(l35,~fYear,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
# bmm = as.data.frame(summary(b))[c('response', 'lower.CL','upper.CL')]
# bmm$Model = 'Base'
# bmm$Year = 2006:2023
# d = emmeans::emmeans(l35a,~fYear,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
# dmm = as.data.frame(summary(d))[c('response', 'lower.CL','upper.CL')]
# dmm$Model = 'BD'
# dmm$Year = 2006:2023
# f = emmeans::emmeans(l35b,~fYear,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
# fmm = as.data.frame(summary(f))[c('response', 'lower.CL','upper.CL')]
# fmm$Model = 'BxD'
# fmm$Year = 2006:2023
# g = emmeans::emmeans(l35c,~fYear,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
# gmm = as.data.frame(summary(g))[c('response', 'lower.CL','upper.CL')]
# gmm$Model = 'BxDRL'
# gmm$Year = 2006:2023
# 
# h = emmeans::emmeans(l35d,~fYear,offset=0,type='response',data=subset(aT,LFA==35),at=list(time=c(1:60,160:292)))
# hmm = as.data.frame(summary(h))[c('response', 'lower.CL','upper.CL')]
# hmm$Model = 'BxDRV'
# hmm$Year = 2006:2023
# 
# 
# 
# 
# oo = do.call(rbind,list(bmm,dmm,fmm,gmm,hmm))
# write.csv(oo,'LFA35MarginalMeanCPUE.csv')
# 
# #marginal mean, internally consistenat
# 
# ggplot(oo, aes(x = Year, y = response ,colour=Model,fill=Model)) +
#   geom_line(linewidth=1) +
#   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width=0.2)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Marginal Mean CPUE')+
#   theme_test(base_size = 14)
# 
# # 
# #
# 

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


# # ##
# # packageurl = "https://cran.r-project.org/src/contrib/Archive/ggeffects/ggeffects_1.3.4.tar.gz"
# # 
# # install.packages(packageurl, repos=NULL, type="source")
# # 
# # #turning conditional into marginal weighting by prop of total catch
#  pr2 = ggpredict(l35a,terms=c('fYear','time[1:235]'),condition = c(leffort=0), back_transform = F)
#  pr2ta = data.frame(fYear=pr2$x,predicted=pr2$predicted,Model='BD',se=pr2$std.error,time=pr2$group)
# # 
# # 
# # ##marginal weighted by effort
# 
# ind = aggregate(NUM_OF_TRAPS~time+SYEAR,data=subset(aT,LFA ==35),FUN=sum)
# ind1 = aggregate(NUM_OF_TRAPS~SYEAR,data=subset(aT,LFA ==35),FUN=sum)
# names(ind1)[2] = 'SumTraps'
# ind = merge(ind,ind1)
# ind$prop = ind$NUM_OF_TRAPS/ind$SumTraps
# ind$fYear=as.factor(ind$SYEAR)
# ppr = merge(pr2ta,ind[,c('fYear','time','prop')])
# ppr$wpre = ppr$predicted*ppr$prop
# ppr$wvar = (ppr$se)^2*ppr$prop
# p2 = aggregate(cbind(wpre,wvar)~fYear+Model,data=ppr,FUN=sum)
# p2$wse = sqrt(p2$wvar)
# 
# pr3 = ggpredict(l35b,terms=c('fYear','time[c(1:235]'),condition = c(leffort=0))
# pr3ta = data.frame(fYear=pr3$x,predicted=pr3$predicted,Model='BxD',se=pr3$std.error,time=pr3$group)
# pp3r = merge(pr3ta,ind[,c('fYear','time','prop')])
# pp3r$wpre = pp3r$predicted*pp3r$prop
# pp3r$wvar = (pp3r$se)^2*pp3r$prop
# p3 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp3r,FUN=sum)
# p3$wse = sqrt(p3$wvar)
# 
# pr4 = ggpredict(l35c,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
# pr4ta = data.frame(fYear=pr4$x,predicted=pr4$predicted,Model='BxDRL',se=pr4$std.error,time=pr4$group)
# pp4r = merge(pr4ta,ind[,c('fYear','time','prop')])
# pp4r$wpre = pp4r$predicted*pp4r$prop
# pp4r$wvar = (pp4r$se)^2*pp4r$prop
# p4 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp4r,FUN=sum)
# p4$wse = sqrt(p4$wvar)
# 
# pr5 = ggpredict(l35d,terms=c('fYear','time[1:235]'),condition = c(leffort=0))
# pr5ta = data.frame(fYear=pr5$x,predicted=pr5$predicted,Model='BxDRV',se=pr5$std.error,time=pr5$group)
# pp5r = merge(pr5ta,ind[,c('fYear','time','prop')])
# pp5r$wpre = pp5r$predicted*pp5r$prop
# pp5r$wvar = (pp5r$se)^2*pp5r$prop
# p5 = aggregate(cbind(wpre,wvar)~fYear+Model,data=pp5r,FUN=sum)
# p5$wse = sqrt(p5$wvar)
# 
# 
# 
# p1 = pr1t
# names(p1)[c(1,2,6)] =c('fYear','wpre','wse')
# 
# pS = rbind(p1[,c(1,2,3,6)],p2[,c(1,2,3,5)],p3[,c(1,2,3,5)],p4[,c(1,2,3,5)],p5[,c(1,2,3,5)])
# pS$Year = as.numeric(as.character(pS$fYear))
# 
# ggplot(pS, aes(x = Year, y = wpre ,colour=Model,fill=Model)) +
#   geom_line(linewidth=1) +
#   geom_errorbar(aes(ymin = wpre-wse, ymax = wpre+wse), width=0.2)+
#   #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
#   xlab('Year')+
#   ylab('Standardized CPUE')+
#   theme_test(base_size = 14)
# 
# 
# 
# ####cross valiation
# library(tidyverse)
# library(modelr)
# library(mgcv)
# 
# #100 samples
# cv_df = crossv_mc(subset(aT,LFA==35), 50) 
# 
# #names samples
# cv_df =
#   cv_df |> 
#   mutate(
#     train = map(train, as_tibble),
#     test = map(test, as_tibble))
# 
# # run
# cv_res = 
#   cv_df |> 
#   mutate(
#     sim_mod  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+offset(leffort),data=df,family = Gamma(link='log'))),
#     mod2  = map(train, \(df) mgcv::gam(WEIGHT_KG~fYear+s(time)+offset(leffort),data=df,family = Gamma(link='log'))), 
#     mod3  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+offset(leffort),data=df,family = Gamma(link='log'))),
#     mod4  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rLID,bs='re')+offset(leffort),data=df,family = Gamma(link='log'))), 
#     mod5  = map(train, \(df) mgcv:: gam(WEIGHT_KG~fYear+s(time,by=fYear)+s(rV,bs='re')+offset(leffort),data=df,family = Gamma(link='log')))) |>
#    mutate(
#     rmse_M1 = map2_dbl(sim_mod, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod2, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod3, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod4, test, \(mod, df) rmse(model = mod, data = df)),
#     rmse_M2 = map2_dbl(mod5, test, \(mod, df) rmse(model = mod, data = df)))
# 
# #plot cv
# cv_res |> 
#   select(starts_with("rmse")) |> 
#   pivot_longer(
#     everything(),
#     names_to = "model", 
#     values_to = "rmse",
#     names_prefix = "rmse_") |> 
#   mutate(model = fct_inorder(model)) |> 
#   ggplot(aes(x = model, y = rmse)) + geom_violin()
# 
# ############################
# #spatial
# 
# aa = subset(aT,SYEAR>2005 & SYEAR<2023)
# xx = split(aa,f=aa$SYEAR)
# m=0
# ou = list()
# for(i in 1:length(xx)){
#      w=   xx[[i]]
#      y = unique(w$SYEAR)
#      ww = unique(w$GRID_NUM)
#     for(j in 1:length(ww)){
#       m=m+1
#       www = subset(w,GRID_NUM==ww[j])
#       g = unique(www$GRID_NUM)
#       uL = length(unique(www$LICENCE_ID))
#       uLFA = length(unique(www$LFA))
#       #grArea = as.numeric(sum(st_area(subset(gr,grid==ww[j]))))#these are wrong do no deal with holes right
#       wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
#       wa$CUME = cumsum(wa$NUM_OF_TRAPS)
#       wa$CUML = cumsum(wa$WEIGHT_KG)
#       wa$rT = round(wa$temp,2)
#       wa = merge(wa,preT,by.x='rT',by.y='temp')
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$SCPUE = wa$pred*wa$CPUE
#       wa$cpue=wa$SCPUE
#       wa$effort = wa$NUM_OF_TRAPS
#       wa$catch = wa$WEIGHT_KG
#       wa$weight = wa$NUM_OF_TRAPS
#       le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
#       ou[[m]] = c(y,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
#     }
#   
# }
# out = as.data.frame(do.call(rbind,ou))
# names(out) = c('Year','Grid','NLics','NLFAs','q','No','Land','Effort','Trips')
# out$CPUE = out$Land/out$Effort
# 
# ################ grouping grids
# ou = list()
# m=0
# LFAgrid <- read.csv(file.path( project.datadirectory("bio.lobster"), "data","maps","GridPolys.csv"))
# grL = (LFAgrid)
# grL$area=grL$PID
# attr(grL,'projection') <- 'LL'
# 
# #grs is the grids
# #making the list of connections to grids
# require(sp)
# require(spdep)
# 
# 
#   aa = subset(aT,SYEAR>2005 & SYEAR<2023)
#   xx = split(aa,f=aa$SYEAR)
#   m=0
#   ou = list()
#   for(i in 1:length(xx)){
#     w=   xx[[i]]
#     y = unique(w$SYEAR)
#     ww = unique(w$GRID_NUM)
#     grL1 = subset(grL,SID %in% ww)
#     g1 = split(grL1,grL1$SID)
#     nm = c()
#     gp = list()
#     for(k in 1:length(g1)){
#       gp[[k]] = Polygons(list(Polygon(g1[[k]][,c('X','Y')])),unique(g1[[k]]$SID))
#     }
#     gpp = SpatialPolygons(gp,proj4string=CRS("+proj=longlat +datum=WGS84"))
#     gpnb = poly2nb(gpp,row.names=names(gpp))
#     names(gpnb)=names(gpp)
#     
#     for(j in 1:length(gpnb)){
#       gr = names(gpnb)[gpnb[[j]]]
#       ss = ifelse(length(gr)>=3,3,length(gr))
#       gr = sample(gr,ss)
#       m=m+1
#       www = subset(w,GRID_NUM==gr)
#       g = length(unique(www$GRID_NUM))
#       gn = paste(gr,collapse="_")
#       uL = length(unique(www$LICENCE_ID))
#       uLFA = length(unique(www$LFA))
#       bgf = subset(grL,SID %in%gr)
#       attr(bgf,'projection') <-"LL"
#       grArea = sum(calcArea(convUL(bgf))$area)
#       wa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~GRID_NUM+DATE_FISHED+temp,data=www,FUN=sum)
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$CPUE = wa$WEIGHT_KG / wa$NUM_OF_TRAPS
#       wa$CUME = cumsum(wa$NUM_OF_TRAPS)
#       wa$CUML = cumsum(wa$WEIGHT_KG)
#       wa$rT = round(wa$temp,2)
#       wa = merge(wa,preT,by.x='rT',by.y='temp')
#       wa = wa[order(wa$DATE_FISHED),]
#       wa$SCPUE = wa$pred*wa$CPUE
#       wa$cpue=wa$SCPUE
#       wa$effort = wa$NUM_OF_TRAPS
#       wa$catch = wa$WEIGHT_KG
#       wa$weight = wa$NUM_OF_TRAPS
#       le = deluryLeslie(y=wa,estimate='leslie',method='robust',plot=F,weight = wa$weight)
#       ou[[m]] = c(y,gn,g,uL,uLFA,le$q,le$N0,max(wa$CUML),max(wa$CUME),nrow(www))
#     }
#     
#   }
#   outa = as.data.frame(do.call(rbind,ou))
#   names(outa) = c('Year','Grid','NGrids','NLics','NLFAs','q','No','Land','Effort','Trips')
#   outa$CPUE = outa$Land/outa$Effort
#   outa = toNums(outa,c(1,3:10))
#   outas = subset(outa,No>0)
#   outas$Rexp = outas$Land/outas$No
#   ggplot(outas,aes(x=Land,y=No))+geom_point()+facet_wrap(~Year,scales='free')
#   ggplot(outas,aes(x=Year,y=Rexp))
#   
#   
# head(outa)  
