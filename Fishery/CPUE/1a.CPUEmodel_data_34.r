require(mgcv)
require(bio.lobster)
require(spdep)
require(sf)
require(devtools)
require(dplyr)
la()

setwd(file.path(project.datadirectory('Framework_LFA33_34_41')))

##temperature data

d = readRDS(file=file.path(project.datadirectory('bio.lobster.glorys'),'Glorys1994_2025wBiasCorrColumn_doy_grid_agg_sept2.rds'))
d$Date = as.Date(d$Date)
d$bcT = d$bcT[,3]
d$z = d$z[,1]
d = subset(d,select=c(LFA,GRID_NO,yr,doy,bcT,z))



g = lobster.db('process.logs')
g = subset(g, SYEAR<=2026 & LFA==34)
po = lobster.db('port_location')
#bring in voluntary log data to populate <2005
v = lobster.db('process.vlog')
v = subset(v,LFA %in% c('34' ))

#assign port code to most common grid
vPC = unique(v$PORT_CODE)
gPC = aggregate(SD_LOG_ID~GRID_NUM+LFA+COMMUNITY_CODE, data=subset(g, LFA==34 & SYEAR<=2020 & COMMUNITY_CODE %in% vPC),FUN=length)
library(dplyr)

com_grid = gPC %>%
  group_by(COMMUNITY_CODE, GRID_NUM) %>%
  summarise(total_SD_LOG_ID = sum(SD_LOG_ID), .groups = "drop") %>%
  group_by(COMMUNITY_CODE) %>%
  slice_max(total_SD_LOG_ID, n = 1, with_ties = FALSE)

v = merge(v,com_grid[,c('COMMUNITY_CODE','GRID_NUM')],by.x=c('PORT_CODE'),by.y='COMMUNITY_CODE')
v = subset(v,select=c(LFA,DATE_FISHED,SYEAR,UNIQUE_FISHER, WOS, DOS, GRID_NUM,PORT_CODE, N_TRP,W_KG ))
v$SD_LOG_ID = 1:nrow(v)
v$SOURCE = 'vlog'

g = subset(g,select=c(LFA,DATE_FISHED,SYEAR,LICENCE_ID, WOS,DOS, GRID_NUM,COMMUNITY_CODE, NUM_OF_TRAPS, WEIGHT_KG,SD_LOG_ID ))
g$SOURCE = 'marfis'
names(v) = names(g)

aT = bind_rows(v,g)
names(aT)[c(4)] = 'UNIQUE_HARVESTER'

aT$GRID_NO = aT$GRID_NUM
aT$DOY = lubridate::yday(aT$DATE_FISHED)
aT$yr = lubridate::year(aT$DATE_FISHED)

aT = merge(aT,d,by.x=c('LFA','DOY','GRID_NUM','yr'),by.y=c('LFA','doy','GRID_NO','yr'),all.x=T)

#####
saveRDS(aT, 'ExtendedTS_LFA34_logs.rds')

#####




gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
coa = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
coa= st_make_valid(coa)
coa = subset(coa,PROVINCE=='Nova Scotia')
coa = subset(coa,st_area(coa)==max(st_area(coa))) # remove islands and cape breton

gr$GRID_NO = as.numeric(gr$GRID_NO)
#remove the islands or multipart polygons and keep only the biggest ones
gr <- gr %>%
  mutate(area = st_area(geometry)) %>%
  group_by(GRID_NO, LFA) %>%
  slice_max(order_by = area, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-area)
gall = gtot = gr
gtot$centroid = st_centroid(gtot)
cen_coords = st_coordinates(gtot$centroid)
gtot$X = cen_coords[,1] 
gtot$Y = cen_coords[,2] 
gtot$centroid <- NULL
gtot$geometry <- NULL
gtots = st_as_sf(gtot,coords=c('X','Y'),crs=4326)
gtot = subset(gtots,LFA %in% c(33))
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
gtot$geometry <- NULL

po3 = subset(po, PORT_CODE %in% unique(aT$COMMUNITY_CODE) & LFA ==33)
correctPorts_aT = subset(aT, COMMUNITY_CODE %in% po3$PORT_CODE) 

cwg = merge(correctPorts_aT,gtot)
cwg = merge(cwg,po3[,c('PORT_CODE','CENTLAT','CENTLON')],by.x='COMMUNITY_CODE',by.y='PORT_CODE')
require(geosphere)

cwg$distance <- distHaversine(
  p1 = cwg[, c("X", "Y")],
  p2 = cwg[, c("CENTLON", "CENTLAT")]
)/1000

##center of gravity for both catch and effort
cc = split(subset(cwg,!is.na(NUM_OF_TRAPS)),f=subset(cwg,!is.na(NUM_OF_TRAPS))$SYEAR)
o = list()
for(i in 1:length(cc)){
  xcw = sum(cc[[i]]$X*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  xce = sum(cc[[i]]$X*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  ycw = sum(cc[[i]]$Y*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  yce = sum(cc[[i]]$Y*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  d2s = sum(cc[[i]]$distance*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  o[[i]] = c(unique(cc[[i]]$SYEAR),xce,xcw,yce,ycw,d2s)      
}

w = as.data.frame(do.call(rbind,o))
names(w) = c('SYEAR','Xef','Xca','Yef','Yca','Dist2Shore')

g = ggLobsterMap('33',return.object = T)
gca = g+geom_point(data=w,aes(x=Xca,y=Yca,colour='red'))+geom_path(data=w,aes(x=Xca,y=Yca,colour='red'))
gef = g+geom_point(data=w,aes(x=Xef,y=Yef,colour='red'))+geom_path(data=w,aes(x=Xef,y=Yef,colour='red'))

#total distance and bearing for catch and effort
require(geosphere)

dist.effort = (distVincentySphere(w[-1,c('Xef','Yef')],w[-nrow(w),c('Xef','Yef')]))
bear.effort = (bearing(w[-nrow(w),c('Xef','Yef')],w[-1,c('Xef','Yef')]))

bear.effort.weight = sum(dist.effort*bear.effort) / sum(dist.effort)
ov_dist.ef = distVincentySphere(w[1,c('Xef','Yef')],w[nrow(w),c('Xef','Yef')])
ov_bear.ef = bearing(w[1,c('Xef','Yef')],w[nrow(w),c('Xef','Yef')])
dest.ef <- destPoint(p=c(w$Xef[1],w$Yef[1]),  b=ov_bear.ef, d= ov_dist.ef)
def = data.frame(sx = w$Xef[1],sy = w$Yef[1], ex = dest.ef[1],ey = dest.ef[2])

dist.cat = distVincentySphere(w[-1,c('Xca','Yca')],w[-nrow(w),c('Xca','Yca')])
bear.cat = bearing(w[-1,c('Xca','Yca')],w[-nrow(w),c('Xca','Yca')])
bear.cat.weight = sum(dist.cat*bear.cat) / sum(dist.cat)

ov_dist.ca = distVincentySphere(w[1,c('Xca','Yca')],w[nrow(w),c('Xca','Yca')])
ov_bear.ca = bearing(w[1,c('Xca','Yca')],w[nrow(w),c('Xca','Yca')])
dest.ca <- destPoint(p=c(w$Xca[1],w$Yca[1]),  b=ov_bear.ca, d= ov_dist.ca)
dca = data.frame(sx = w$Xca[1],sy = w$Yca[1], ex = dest.ca[1],ey = dest.ca[2])

plot(w$Xca,w$Yca,type='l')
text(w$Xca,w$Yca,labels=w$SYEAR)
arrows(x0=dca$sx,x1=dca$ex,y0=dca$sy,y1=dca$ey,col='red',lwd=2,length=.1)
plot(w$Xef,w$Yef,type='l')
text(w$Xef,w$Yef,labels=w$SYEAR)



gca+geom_segment(data=dca,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')
gef+geom_segment(data=def,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')

ggplot(w,aes(SYEAR, Dist2Shore))+geom_line()


###temperature modelling
d = readRDS(file=file.path(project.datadirectory('bio.lobster.glorys'),'Glorys2000_2025wBiasCorrColumn_doy_grid_agg_july29.rds'))
d$Date = as.Date(d$Date)
d$bcT = d$bcT[,3]
d$z = d$z[,1]
d = subset(d,select=c(LFA,GRID_NO,yr,doy,bcT,z))
#grids missing temperature adn will use the adjacent grid
temgrids = as.data.frame(cbind(missing=c(42,43,44,69,77,107,141,159,196,214,1278,1344,1346,1409,1591),fills=c(53,54,55,81,78,108,140,158,195,213,1279,1343,1347,1408,167)))
v=subset(d,GRID_NO %in% temgrids$fills& LFA %in% c(33,34,41))
v$GRID_NO <- temgrids$missing[match(v$GRID_NO, temgrids$fills)]
d = rbind(d,v)
cwga$yr = lubridate::year(cwga$DATE_FISHED)
cda = merge(cwga,d,by.x=c('LFA','GRID_NO','yr','DOY'), by.y=c('LFA','GRID_NO','yr','doy'))

#common start data at the end of november early december (41 is all year, but need to start somewhere)

ca <- cda %>%
  mutate(         month = month(DATE_FISHED)) %>%
  group_by(SYEAR) %>%
  mutate(start = min(DATE_FISHED[month %in% 11 : 12 ], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DOS = as.integer(DATE_FISHED - start)/24/60/60+1)

ca$fyr = as.factor(ca$SYEAR)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca$fGRID_NO = as.factor(ca$GRID_NO)


require(sdmTMB)
require(splines)
crs_utm20 <- 32620
gs <- st_transform(goo, crs_utm20)
st_geometry(gs) = st_geometry(gs)/1000
st_crs(gs) <- crs_utm20
gs$X = st_coordinates(gs)[,1]
gs$Y = st_coordinates(gs)[,2]
gto = as_tibble(gs)
ns_coast =readRDS(file.path( bio.directory, "bio.lobster.data","mapping_data","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ylim=c(41.1,48); 		xlim=c(-67.8,-57.8)

sf_use_s2(FALSE) #needed for cropping
ns_coast = suppressWarnings(suppressMessages(st_crop(ns_coast,xmin=xlim[1],ymin=ylim[1],xmax=xlim[2],ymax=ylim[2])))

#remove island
nsc <- ns_coast %>%
  mutate(area = st_area(geometry)) %>%
  slice_max(order_by = area, n = 5, with_ties = FALSE) %>%
  select(-area)
nsc <- st_transform(nsc, crs_utm20)

ns_coast = st_simplify(nsc,dTolerance = 10000)

st_geometry(ns_coast) = st_geometry(ns_coast)/1000
st_crs(ns_coast) <- crs_utm20

ca = subset(ca, X< -60.5)
cas = st_as_sf(ca,coords = c('X','Y'),crs=4326)
cas = st_transform(cas,crs=crs_utm20)
st_geometry(cas) = st_geometry(cas)/1000
st_crs(cas) <- crs_utm20
cas$X = st_coordinates(cas)[,1]
cas$Y = st_coordinates(cas)[,2]
ca = as_tibble(cas)

saveRDS(list(gto,ns_coast,ca),file=file.path('CPUE','data_33_34_41_sdmtmb.rds'))
