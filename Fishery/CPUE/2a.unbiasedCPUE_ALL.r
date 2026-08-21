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
la()
gr = readRDS('C:/Users/Cooka/Documents/git/bio.lobster.data/mapping_data/GridPolyLand.rds')
st_crs(gr) <- 4326
gr = st_transform(gr,32620) 
st_geometry(gr) <- st_geometry(st_as_sf(gr$geometry/1000)) 
st_crs(gr) <- 32620

#need the temperature - catch relationship using IP/catchability/TempCatch.r output
#preT = readRDS(file=file.path(project.datadirectory('bio.lobster'),'analysis','ClimateModelling','tempCatchability.rds'))
#preT$temp = round(preT$Temperature,2)
#preT = aggregate(pred~temp,data=preT,FUN=mean)

aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2026 & LFA %in% c(33,34))

#start with an overall CPUE trend

aa = split(aT,f=list(aT$LFA,aT$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-round(julian(tmp$date,origin=first.day-1))
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)))
  cpue.ann[[i]]=g
}

ca =as.data.frame(do.call(rbind,cpue.ann))
ggplot(subset(ca,LFA %in% c(33,34)),aes(x=SYEAR,y=unBCPUE))+geom_point()+geom_line()+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.3)+
  facet_wrap(~LFA)+xlab('Fishing Year (ending)')+ylab('unbiased CPUE')+theme_test()

cc =as.data.frame(do.call(rbind,cpue.lst))
cc$dyear = cc$yr+cc$t/365
#cc = merge(cc,preT)

ggplot(subset(cc,lfa %in% c(33) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()
ggsave(file.path(fig_dir,'33_daily_CPUE.png'))


ggplot(subset(cc,lfa %in% c(34) ),aes(x=t,y=unBCPUE))+geom_point(size=.1)+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.1)+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('unbiased CPUE')+ylim(c(0,8))+theme_test()
ggsave(file.path(fig_dir,'34_daily_CPUE.png'))

saveRDS(list(ca,cc),file='unBIASED_CPUE.rds')


####logbook errors

a = lobster.db('process.logs')
b=lobster.db('process.logs.unfiltered')

aa = aggregate(SD_LOG_ID~SYEAR+LFA,data=a,FUN=function(x) length(unique(x)))
bb = aggregate(SD_LOG_ID~SYEAR+LFA,data=b,FUN=function(x) length(unique(x)))

names(bb)[3]='SD_LOG_UNF'

ab = merge(aa,bb)
ab$prop = 1-ab$SD_LOG_ID/ab$SD_LOG_UNF

ab=subset(ab,LFA %in% c(33,34) & SYEAR>2005 & SYEAR<2026)
agg = aggregate(prop~LFA,data=ab,FUN=mean)
ggplot(ab,aes(x=SYEAR,y=prop))+geom_line()+facet_wrap(~LFA)+theme_test(base_size = 14)+xlab('Fishing Season')+ylab('Proportion Logs Containing Errors')+geom_hline(data=agg,aes(yintercept=prop),color='red')
ggsave(file.path(fig_dir,'ProportionOfLogsContainingErrors.png'))




######within season diffs
cd = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~DOS+LFA,data=aT,FUN=sum)
cd$oCPUE = cd$WEIGHT_KG/cd$NUM_OF_TRAPS
cd = subset(cd, select=c(DOS,oCPUE,LFA))

cc = merge(cc,cd,by.x=c('t','lfa'),by.y=c('DOS','LFA'))
cc$ANO = cc$CPUE - cc$oCPUE
cc$cols=  ifelse(cc$ANO>0,'red','blue')
ggplot(subset(cc,lfa %in% c(33) ))+geom_segment(aes(x=t,xend=t,y=0,yend=ANO,colour=cols))+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('Daily CPUE Anomaly')+theme_test()
ggsave(file.path(fig_dir,'33_DailyCPUEAnomalies_Errors.png'))


ggplot(subset(cc,lfa %in% c(34) ))+geom_segment(aes(x=t,xend=t,y=0,yend=ANO,colour=cols))+
  facet_wrap(~yr)+xlab('Day of Fishing Year')+ylab('Daily CPUE Anomaly')+theme_test()
ggsave(file.path(fig_dir,'34_DailyCPUEAnomalies_Errors.png'))

########################################################################################
###temp too



aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2026 & LFA %in% c(33,34))


aa = split(aT,f=list(aT$LFA,aT$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-round(julian(tmp$date,origin=first.day-1)/7)+1
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)))
  cpue.ann[[i]]=g
}

ca =as.data.frame(do.call(rbind,cpue.ann))
ggplot(subset(ca,LFA %in% c(33,34)),aes(x=SYEAR,y=unBCPUE))+geom_point()+geom_line()+geom_errorbar(aes(ymin=l95,ymax=u95),width=0,alpha=.3)+
  facet_wrap(~LFA)+xlab('Fishing Year (ending)')+ylab('unbiased CPUE')+theme_test()

cc =as.data.frame(do.call(rbind,cpue.lst))
cc$dyear = cc$yr+cc$t/365


m = readRDS(file='data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
aT = m[[3]]
aT$WOS = floor((aT$DOS-1)/7)+1
aT$catchwttemp = aT$WEIGHT_KG*aT$bcT


cd = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~WOS+LFA,data=aT,FUN=sum)
cd$oCPUE = cd$WEIGHT_KG/cd$NUM_OF_TRAPS
cd = subset(cd, select=c(WOS,oCPUE,LFA))
cc = merge(cc,cd,by.x=c('t','lfa'),by.y=c('WOS','LFA'))
cc$ANO = cc$CPUE - cc$oCPUE
cc$cols=  ifelse(cc$ANO>0,'red','blue')

td = aggregate(cbind(WEIGHT_KG,catchwttemp)~WOS+LFA,data=aT,FUN=sum)
td$ot = td$catchwttemp/td$WEIGHT_KG
td = subset(td, select=c(WOS,ot,LFA))

ttd = aggregate(cbind(WEIGHT_KG,catchwttemp)~WOS+LFA+SYEAR,data=aT,FUN=sum)
ttd$at = ttd$catchwttemp/ttd$WEIGHT_KG
ttd = subset(ttd, select=c(WOS,at,LFA,SYEAR))

ta = merge(ttd,td,by=c('LFA','WOS'))
ta$ANO_T = ta$at - ta$ot

x = merge(cc,ta,by.x=c('t','lfa','yr'),by.y=c('WOS','LFA','SYEAR'))

ggplot(subset(x,lfa %in% c(33) ))+geom_segment(aes(x=t,xend=t,y=0,yend=ANO,colour=cols),linewidth=1)+geom_line(aes(x=t,y=ANO_T))+
  facet_wrap(~yr)+xlab('Week of Fishing Year')+ylab('Anomaly (bars = CPUE, line=temperature)')+theme_test(base_size =14)+theme(
    legend.position = "none"  )


ggsave(file.path(fig_dir,'33_DailyCPUEAnomalies_Errors_temp.png'))



ggplot(subset(x,lfa %in% c(34) ))+geom_segment(aes(x=t,xend=t,y=0,yend=ANO,colour=cols),linewidth=1)+geom_line(aes(x=t,y=ANO_T))+
  facet_wrap(~yr)+xlab('Week of Fishing Year')+ylab('Anomaly (bars = CPUE, line=temperature)')+theme_test(base_size =14)+theme(
    legend.position = "none"  )
ggsave(file.path(fig_dir,'34_DailyCPUEAnomalies_Errors.png'))


library(dplyr)
xall = x
x <- subset(xall, lfa==34) %>%
  mutate(
    quadrant = case_when(
      ANO_T >= 0 & ANO >= 0 ~ "Q1",
      ANO_T <  0 & ANO >= 0 ~ "Q2",
      ANO_T <  0 & ANO <  0 ~ "Q3",
      ANO_T >= 0 & ANO <  0 ~ "Q4"
    )
  )
quad_props <- x %>%
  count(quadrant) %>%
  mutate(
    prop = n / sum(n),
    label = paste0(
      scales::percent(prop, accuracy = 0.1),
      "\n(n=", n, ")"
    )
  )

xr <- range(x$ANO_T, na.rm = TRUE)
yr <- range(x$ANO, na.rm = TRUE)

quad_labels <- data.frame(
  quadrant = c("Q1", "Q2", "Q3", "Q4"),
  x = c(0.75*xr[2], 0.75*xr[1], 0.75*xr[1], 0.75*xr[2]),
  y = c(0.75*yr[2], 0.75*yr[2], 0.75*yr[1], 0.75*yr[1])
) %>%
  left_join(quad_props, by = "quadrant")


require(ggdensity)
ggplot(x, aes(ANO_T, ANO)) +
  geom_point(shape = 21, fill = "white", color = "black", size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_text(
    data = quad_labels,
    aes(x, y, label = label),
    fontface = "bold",
    size = 5
  ) +
  geom_hdr_lines(probs = c(0.5, 0.75, 0.9)) +
  labs(
    x = "Temperature Anomaly",
    y = "CPUE Anomaly"
  ) +
  theme_test(base_size = 14)


library(ggplot2)

ggplot(df, aes(temp_anom, cpue_anom)) +
  geom_point(shape = 21, fill = "white", size = 2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    breaks = c(0.50, 0.75, 0.90),
    alpha = 0.2
  ) +
  stat_density_2d(
    breaks = c(0.50, 0.75, 0.90),
    linewidth = 0.8
  ) +
  theme_bw()