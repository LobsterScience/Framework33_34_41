## temperature optimal interpolation
#using ILTS and RV for 4VWX
require(bio.lobster)
require(devtools)
require(gstat)
require(sf)
la()
require(dplyr)
require(ggplot2)
require(tidyr)
lobster.db('survey')

isf = subset(surveyCatch,select=c(TRIP_ID,SET_NO,SET_LAT,SET_LONG)) %>%
        distinct()

#summer ILTS
y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(70,300),aggregate=T,species=2550,biomass = F)
y = subset(y, !is.na(temp) & month(SET_DATE) %in% c(6,7),select=c(YEAR,SET_LONG,SET_LAT,temp))
y$X = y$SET_LONG
y$Y = y$SET_LAT
y$year = y$YEAR
#summer RV

##prediction grid from bathy
ba = readRDS('~/git/bio.lobster.data/mapping_data/bathymetrySF.rds')
ba = subset(ba,z>15)
pa = readRDS('~/git/bio.lobster.data/mapping_data/Summer_RV_strata.rds')
pa  = st_make_valid(subset(pa,PID>=470))
pau = st_union(pa)
pau = st_transform(pau,crs=32620)
ba_i = ba[st_within(ba,st_as_sf(pau), sparse=FALSE),]



d = groundfish.db('gsinf.odbc')
d = subset(d,strat %in% unique(pa$PID) & month(sdate) %in% c(6,7,8) & !is.na(bottom_temperature))
d$year = year(d$sdate)
d = subset(d,select=c(year,slat,slong,bottom_temperature))
d$X = convert.dd.dddd(d$slong)* -1
d$Y = convert.dd.dddd(d$slat)
d$temp = d$bottom_temperature

yd = rbind(subset(d,select=c(year,X,Y,temp)),subset(y,select=c(year,X,Y,temp)))

yd = st_as_sf(yd,coords=c('X','Y'),crs=4326)
ydu = st_transform(yd,crs=32620)
sf_use_s2(FALSE)

#do this is you want external drift and change the gstat formula to be temp ~ z

#bathy_model <- gstat(id = "z", formula = z ~ 1, locations = ba, nmax = 7, set = list(idp = 2))
#ydu$z = predict(bathy_model,newdata = ydu)$z.pred

yrs=unique(ydu$year)
for(i in 1:length(yrs)){
  v = subset(ydu,year==yrs[i])
  im = gstat(id = "temp", formula = temp ~ 1, locations = v, nmax = 7, set = list(idp = 2))
  pp = predict(im, newdata=ba_i)
  st_geometry(pp) <- NULL
  ba_i$new = pp$temp.pred
  ba_i = rename.df(ba_i,'new',paste0('X',yrs[i]))
}

##climatology 1991-2020

ba_i$clim <- ba_i %>%
        st_set_geometry(NULL) %>%
        dplyr::select(dplyr::matches("^X(199[1-9]|200[0-9]|201[0-9]|2020)")) %>%
        rowMeans( na.rm = TRUE)


ggplot(ba_i,aes(colour=clim,fill=clim))+geom_sf()+
  scale_color_viridis_c()+
  scale_fill_viridis_c()

ba_i$id=1:nrow(ba_i)

long_df <- ba_i %>%
  st_set_geometry(NULL) %>%
  gather(key = "year", value = "value", matches("^X[0-9]{4}")) %>%
  mutate(year = as.numeric(sub("X", "", year)))


ld = st_as_sf(merge(long_df,ba_i[,c('id')],all=T))

ld$Anomaly = ld$value - ld$clim


#slows it down
ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs

ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -63, ymax = 45.5))))

ns_coast <- st_transform(ns_coast, 32620)

ggplot()+geom_sf(data=ns_coast)+
 geom_sf(data=subset(ld,year%in%c(2012,2016,2020,2023,2024,2025)),aes(fill=Anomaly,colour=Anomaly),size=1)+
 scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
 scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
   theme_test_adam()+
 facet_wrap(~year)+
 coord_sf()


ag = aggregate(Anomaly~year,data=ld,FUN=mean)
agm <- ag %>%
  arrange(year) %>%
  mutate(running_mean = zoo::rollmean(Anomaly, k = 5, fill = NA, align = "right"))


require(rshift)
# Run regime shift analysis (forward)
result <- Rodionov(ag,'Anomaly','year',l=5,merge=T)
result = result %>%
  mutate(index=cumsum(RSI>0))
result$index = result$index+1
rf = aggregate(Anomaly~index,data=result,FUN=mean)
rfyn = aggregate(year~index,data=result,FUN=min)
rfyx = aggregate(year~index,data=result,FUN=max)
names(rfyx)[2]='maxy'
rf = merge(merge(rf,rfyn),rfyx)


agr = ag[rev(rownames(ag)),]
#reverse
result_r <- Rodionov(agr,'Anomaly','year',l=5,merge=T)
result_r = result_r %>%
            mutate(index=cumsum(RSI>0))

result_r$index = result_r$index+1
rr = aggregate(Anomaly~index,data=result_r,FUN=mean)
rryn = aggregate(year~index,data=result_r,FUN=min)
rryx = aggregate(year~index,data=result_r,FUN=max)
names(rryx)[2]='maxy'
rr = merge(merge(rr,rryn),rryx)

rff = aggregate(Anomaly~index,data=result,FUN=mean)

mean(ld$clim)

tm = ggplot(agm,aes(x=year,y=Anomaly))+geom_line()+geom_point(shape=1,size=2.5)+
  geom_line(aes(y=running_mean),colour='black',size=1.4)+
  geom_segment(data=rr,aes(x=year,xend=maxy,y=Anomaly,yend=Anomaly),colour='red',size=1.3)+
  geom_segment(data=rf,aes(x=year,xend=maxy,y=Anomaly,yend=Anomaly),colour='blue',size=1.3)+
  geom_hline(yintercept=0,colour='grey',size=1.4)+
  geom_hline(yintercept=sd(ld$clim)*.5,colour='grey',size=1.4,linetype='dashed')+
  geom_hline(yintercept=sd(ld$clim)*-.5,colour='grey',size=1.4,linetype='dashed')+
  annotate("text", x = 2019, y = -1.7, label = 'Mean = 7.63\u00B0C ' ,
            size = 5, color = "black")+
  theme_test(base_size = 14)+labs(x='Year',y='Temperature Anomaly')

saveRDS(list(agm,rr,rf,sd(ld$clim)),file='BTAnom4x.rds')
