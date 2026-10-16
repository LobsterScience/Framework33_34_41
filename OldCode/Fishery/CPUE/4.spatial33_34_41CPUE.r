require(mgcv)
require(bio.lobster)
require(spdep)
require(sf)
require(devtools)
require(dplyr)
la()


################ grouping grids
aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>2005 & SYEAR<2025 & LFA %in% c(33,34))
aT$GRID_NO = aT$GRID_NUM
aT$DOY = lubridate::yday(aT$DATE_FISHED)
a4 = lobster.db('process.logs41')
a4 = subset(a4,!is.na(GRID_NO) & yr>2005 & yr<2025)
a4$SYEAR = a4$yr
a4$DATE_FISHED = a4$FV_FISHED_DATETIME
a4$DOY = lubridate::yday(a4$FV_FISHED_DATETIME)
a4$WEIGHT_KG = a4$ADJCATCH_KG
u4 = unique(a4$GRID_NO)

gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
gr41 = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFA41_grid_polys.rds')))
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
gr41$LFA = as.character(gr41$LFA)
gr41 = subset(gr41,GRID_NO %in% u4)
gtot = bind_rows(gr,gr41)
gtot$centroid = st_centroid(gtot)
cen_coords = st_coordinates(gtot$centroid)
gtot$X = cen_coords[,1] 
gtot$Y = cen_coords[,2] 
gtot$centroid <- NULL
gtot$geometry <- NULL
gtots = st_as_sf(gtot,coords=c('X','Y'),crs=4326)
gtot = subset(gtots,LFA %in% c(33,34,41))
gtot$dist_to_shore = as.numeric(st_distance(gtot,coa))/1000
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
gtot = subset(gtot,X< -60.5)
goo = gtot
gtot$geometry <- NULL

logs = subset(aT,select=c(LFA, GRID_NO, SYEAR, DATE_FISHED,WEIGHT_KG,NUM_OF_TRAPS,DOY))
l41 = subset(a4,select=c(LFA, GRID_NO, SYEAR, DATE_FISHED,WEIGHT_KG,NUM_OF_TRAPS,DOY))
st_geometry(l41) = NULL

combined = bind_rows(logs,l41)
cwga = cwg = merge(combined,gtot)
cwg$rDist = round(cwg$dist_to_shore/50) * 50
##center of gravity for both catch and effort
cc = split(cwg,f=cwg$SYEAR)
o = list()
for(i in 1:length(cc)){
      xcw = sum(cc[[i]]$X*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
      xce = sum(cc[[i]]$X*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
      ycw = sum(cc[[i]]$Y*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
      yce = sum(cc[[i]]$Y*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
      d2s = sum(cc[[i]]$dist*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
    o[[i]] = c(unique(cc[[i]]$SYEAR),xce,xcw,yce,ycw,d2s)      
  }

w = as.data.frame(do.call(rbind,o))
names(w) = c('SYEAR','Xef','Xca','Yef','Yca','Dist2Shore')

g = ggLobsterMap('SWN',return.object = T)
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

gca+geom_segment(data=dca,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')
gef+geom_segment(data=def,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')

ggplot(w,aes(SYEAR, Dist2Shore))+geom_point()


# without LFA 41

combined = bind_rows(logs)
cwg = merge(combined,gtot)
cwg$rDist = round(cwg$dist_to_shore/50) * 50
##center of gravity for both catch and effort
cwg = subset(cwg,LFA %ni% 41)
cc = split(cwg,f=cwg$SYEAR)
o = list()
for(i in 1:length(cc)){
  xcw = sum(cc[[i]]$X*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  xce = sum(cc[[i]]$X*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  ycw = sum(cc[[i]]$Y*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  yce = sum(cc[[i]]$Y*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  d2s = sum(cc[[i]]$dist*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  o[[i]] = c(unique(cc[[i]]$SYEAR),xce,xcw,yce,ycw,d2s)      
}

w = as.data.frame(do.call(rbind,o))
names(w) = c('SYEAR','Xef','Xca','Yef','Yca','Dist2Shore')

g = ggLobsterMap('SWN',return.object = T)
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

gca+geom_segment(data=dca,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')
gef+geom_segment(data=def,aes(x=sx,y=sy,xend=ex,yend=ey),arrow=arrow(length=unit(0.1,'cm')),colour='blue')

ggplot(w,aes(SYEAR, Dist2Shore))+geom_point()


####
####modelling
#bring in temperature

d = readRDS(file.path(project.datadirectory('bio.lobster.glorys'),'DailyTemp_by_grid_05-24_withAnom.rds'))
d = subset(d,select=c(LFA,GRID_NO,yr,doy,bottomT,Anomaly))
d$bT = d$bottomT[,3]
d$bottomT <- NULL
cwga$yr = lubridate::year(cwga$DATE_FISHED)
cda = merge(cwga,d,by.x=c('LFA','GRID_NO','yr','DOY'), by.y=c('LFA','GRID_NO','yr','doy'))

require(mgcv)
summary(cda)
cda$fyr = as.factor(cda$yr)
cda$logWt = log(cda$WEIGHT_KG)
cdaa = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~bT+DOY+yr+GRID_NO+LFA+X+Y+dist_to_shore,data=cda,FUN=sum)
cdaa$logWt = log(cdaa$WEIGHT_KG)
cdaa$fyr = as.factor(cdaa$yr)

m1 = bam(WEIGHT_KG~s(bT)+s(DOY,by=fyr)+fyr,offset = log(NUM_OF_TRAPS),data=cda,correlation=corExp(form= ~X +Y),family=scat())

m2 = bam(logWt~s(bT)+s(DOY,by=fyr)+fyr+GRID_NO,offset = log(NUM_OF_TRAPS),data=cdaa,correlation=corExp(form= ~X +Y),family=scat())

###focus on LFA 33 and 34 and use day of season
cdas = subset(cda,LFA %in% 33:34)

ca <- cdas %>%
  mutate(         month = month(DATE_FISHED)) %>%
  group_by(SYEAR) %>%
  mutate(start = min(DATE_FISHED[month %in% 11 : 12 ], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DOS = as.integer(DATE_FISHED - start)/24/60/60+1)
ca$fyr = as.factor(ca$SYEAR)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca$fGRID_NO = as.factor(ca$GRID_NO)
m3 = gam(WEIGHT_KG~fyr+s(DOS,by=fyr)+s(bT)+LFA+offset(leffort),data=ca,family = 'nb')

m4 = gam(WEIGHT_KG~fyr+s(DOS,by=fyr)+s(DOS,by=GRID_NO)+GRID_NO+s(bT)+LFA+offset(leffort),data=ca,family = 'nb')


saveRDS(list(ca,m3),file=file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUEModelM3.rds'))
m3 = readRDS(file=file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUEModelM3.rds'))


#weighted marginal means
ind = aggregate(SD_LOG_ID~DOS+SYEAR+LFA,data=aT,FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR+LFA,data=ind,FUN=sum)
names(ind1)[3] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fyr=as.factor(ind$SYEAR)

da = da[order(round(da))]
da = unique(da)

h = emmeans::emmeans(m3,~fyr+DOS+LFA,offset=0,type='response',data=ca,at=list(DOS=da))
dmm = as.data.frame(summary(h))
dmm = merge(dmm,ind)
usm = sum(dmm$SD_LOG_ID)
dmm$Year = as.numeric(as.character(dmm$fyr))
jj=33:34
ii = unique(dmm$Year)
dom = data.frame(SYEAR=NA,LFA=NA, wse=NA,wemm=NA)
m=0
for(i in 1:length(ii)){
  for(j in 1:length(jj)){
    m=m+1
  z = subset(dmm,Year==ii[i] & LFA==jj[j])
  usm = sum(z$SD_LOG_ID)
  dom[m,'wse'] = sqrt(sum((z$SD_LOG_ID/usm)^2  * (z$SE^2)))
  dom[m,'wemm'] = sum(z$SD_LOG_ID/usm*z$response) 
  dom[m,'SYEAR'] = ii[i]
  dom[m,'LFA'] = jj[j]
  }
}
dom$Model = 'm3'
ggplot(dom,aes(x=SYEAR,y=wemm,colour=as.character(LFA)))+geom_point()+geom_smooth()


###########sdmtmb

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
ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)/1000
st_crs(ns_coast) <- crs_utm20

ca=m3[[1]]
ca = subset(ca, X< -60.5)
cas = st_as_sf(ca,coords = c('X','Y'),crs=4326)
cas = st_transform(cas,crs=crs_utm20)
st_geometry(cas) = st_geometry(cas)/1000
st_crs(cas) <- crs_utm20
cas$X = st_coordinates(cas)[,1]
cas$Y = st_coordinates(cas)[,2]
ca = as_tibble(cas)

gtos = subset(gto, GRID_NO %in% unique(ca$GRID_NO))

mes = sdmTMB::make_mesh(ca,xy_cols = c('X','Y'),n_knots=nrow(gtos)-1)


# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  mes, ns_coast, range_fraction = .1,
  proj_scaling = 1, plot = TRUE
)

m4 = sdmTMB(WEIGHT_KG~ 0+
              s(bT)+s(DOS)+fyr,
            offset = 'leffort',
            data=ca,
            family = nbinom2(link='log') ,
#            time_varying = ~ 1 + bs(DOS, degree=6, intercept=T),
#            time_varying_type = "rw0",
            mesh = bspde,
            spatial='on',
#            time='SYEAR',
 #           spatiotemporal='AR1'
)
