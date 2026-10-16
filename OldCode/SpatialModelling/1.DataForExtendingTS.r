##ILTS and RV survey overviews

require(bio.lobster)
require(devtools)
require(bio.utilities)
require(sf)
require(ggplot2)
require(tidyr)
theme_set(theme_test(base_size = 14))

setwd(file.path(project.datadirectory('Framework_LFA33_34_41')))
la()


#recruit numb proportion

y = ILTS_ITQ_All_Data(redo_base_data = F,size=c(70,82),aggregate=T,species=2550,biomass = F,extend_ts = F)
y = subset(y, SA_CORRECTED>0 & month(y$SET_DATE)<8 & LFA %in% c('L34','L35','L36','L37','L38'),select=c(TRIP_ID,SET_NO,SET_LONG,SET_LAT,SA_CORRECTED_PRORATED_N,tot))

x = RV_sets()
x = subset(x,month(x$DATE) %in% c(6,7,8))

xle = x %>% pivot_longer(starts_with('P'))
xle$Length = as.numeric(substr(xle$name,3,8))
#xle= na.zero(xle,cols='value')
xxa = aggregate(value~mission+setno+Lobster+LATITUDE+LONGITUDE,data=subset(xle,Length %in% 70:82),FUN=sum)

xs = st_as_sf(xxa,coords = c('LONGITUDE','LATITUDE'),crs=4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326

xs = st_join(xs,rL,join=st_within)
xs$X=st_coordinates(xs)[,1]
xs$Y=st_coordinates(xs)[,2]
st_geometry(xs) <- NULL
xs = subset(xs,LFA %in% c(34,35,36,37,38,41,40,33) & X< -64,select=c(mission,setno,X,Y,value,Lobster))
names(xs)[5] = 'Recruits'
names(y)=names(xs)
y$mission = as.character(y$mission)
y$setno = as.character(y$setno)

xy = bind_rows(y,xs)
saveRDS(xy,'proportions70_82.rds')
####################

#comparing to Scallop and NEFSC older sets

ou = readRDS(file.path( bio.directory,'bio.lobster.data', 'survey_corrections','modelled_recruit_Proportions_34-38.rds'))
#if(all(size==c(82,300))){ou = readRDS(file.path( bio.directory,'bio.lobster.data','survey_corrections' ,'modelled_commercial_Proportions_34-38.rds'));u=1}
#if(all(size==c(70,300))){ou = readRDS(file.path( bio.directory,'bio.lobster.data','survey_corrections' ,'modelled_commercial_and_recruit_Proportions_34-38.rds'));u=1}
#if(all(size==c(82,300))&biomass){ou = readRDS(file.path( bio.directory,'bio.lobster.data','survey_corrections' ,'modelled_commercial_Proportions_wt_34-38.rds'));u=1; biomass=F}

n = NEFSC_sets()
n = st_as_sf(n,coords=c('LONGITUDE','LATITUDE'),crs = 4326)
rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
ns = st_join(n,rL,join=st_within)
ns = subset(ns,!is.na(LFA))


ss = st_nearest_feature(ns,ou)
ns$di = st_distance(ns,ou[ss,],by_element=T)
st_geometry(ou) = NULL
ns$prop = ou$Modelled_Proportion[ss]
ns$obsProp = ns$Recruit/ns$Lobster

nn = subset(ns,Lobster>0)
ggplot(nn,aes(x=obsProp,y=prop))+geom_point()
##no good
