require(bio.lobster)
require(devtools)
require(sf)
require(ggplot2)
la()

td = file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'Glorys1994_2025wBiasCorrColumn_doy_sept24.rds')
td = readRDS(td)
td$gloryTemp = td$bottomT + td$pred
td = subset(td,select=c(Date,gloryTemp, X,Y,z,YR))

tds = st_as_sf(td,coords = c('X','Y'),crs=32620)
po = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFA41OffAreas.rds')))
po = st_transform(po,crs=32620)
st_geometry(po) = st_geometry(po)/1000
st_crs(po) <- 32620

tp = st_join(tds,po)
tps = subset(tp,!is.na(OFFAREA ))


ggplot()+geom_sf(data=po)+geom_sf(data=subset(tps,Date=='1994-12-31'))

tps$Month = lubridate::month(tps$Date)
aggTemp = aggregate(gloryTemp~Month+YR+OFFAREA,data=tps,FUN=mean)

saveRDS(tps,'OffshoreGlorysTemperatures_withBC.rds')
ggplot(subset(aggTemp,Month %in% c(3,7)),aes(YR, gloryTemp,colour=OFFAREA))+geom_point()+geom_smooth()+facet_wrap(~Month)


cw = lobster.db('cw.temperature.data')
cw = st_as_sf(subset(cw,!is.na(LAT_DD)),coords=c('LON_DD','LAT_DD'),crs=4326)
cw = st_transform(cw,crs=32620)
st_geometry(cw) = st_geometry(cw)/1000
st_crs(cw) <- 32620

cw = st_join(cw,po)
cw = subset(cw,!is.na(OFFAREA))
cw$Month = lubridate::month(cw$DateSet)
cw$YR = lubridate::year(cw$DateSet)
cwa = aggregate(Temp~YR+Month+OFFAREA,data=cw,FUN=mean)
