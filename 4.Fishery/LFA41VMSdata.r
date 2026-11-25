#VMS_plots
require(sf)
require(devtools)
require(ggplot2)
require(bio.lobster)
la()

lobster.db('logs41')
vms.data = lobster.db('lfa41.vms')
logs41$YEAR = year(logs41$FV_FISHED_DATETIME)

logs41r = subset(logs41, !is.na(DDLON) & !is.na(DDLAT) & !is.na(FV_FISHED_DATETIME))
logs41r$DDLON = logs41r$DDLON * -1

l = st_as_sf(logs41r,coords = c('DDLON','DDLAT'),crs=4326)
v = st_as_sf(vms.data,coords = c('LON','LAT'),crs=4326)
  
l$DATE  =  as.Date(format(strftime(l$FV_FISHED_DATETIME,format="%Y-%m-%d"), tz="America/Halifax",usetz=TRUE))
v$DATE = as.Date(v$DATE)

l$id = paste(l$VR_NUMBER,l$DATE,sep="-")
v$id = paste(v$VRN,v$DATE,sep="-")

da = unique(l$id)
dav = unique(v$id)
ids = intersect(dav,da)

#something wrong, keeps filtering out good vms data
figs = F
if(figs){
pdf(file=file.path(project.figuredirectory('bio.lobster'),'lfa41LogsVMS.pdf'), onefile = T)
for(d in 1: length(ids)) {
 ll = subset(l,id == ids[d])
 vv = subset(v,id == ids[d])
gall = ggplot()+
  geom_sf(data=ll,colour='red')+
  geom_sf(data=vv,colour='blue')+
  labs(title=paste(ids[d],'all'))

vvv = subset(vv,SPEED_KNOTS<5)
dist_matrix <- st_distance(vvv,ll)
points_to_remove <- apply(dist_matrix, 1, function(row) all(row < 20000))
vvv <- vvv[points_to_remove, ]

gs = ggplot()+
  geom_sf(data=ll,colour='red')+
  geom_sf(data=vvv,colour='blue')+
  labs(title=paste(ids[d],'speed and distance filtered'))

print(gridExtra::grid.arrange(gall,gs,ncol=2))
}
dev.off()
}

####
lobster.db('logs41')
vms.data = lobster.db('lfa41.vms')
logs41$YEAR = year(logs41$FV_FISHED_DATETIME)

logs41r = subset(logs41, !is.na(DDLON) & !is.na(DDLAT) & !is.na(FV_FISHED_DATETIME))
logs41r$DDLON = logs41r$DDLON * -1

l = st_as_sf(logs41r,coords = c('DDLON','DDLAT'),crs=4326)
v = st_as_sf(vms.data,coords = c('LON','LAT'),crs=4326)
v$X = st_coordinates(v)[,1]
v$Y = st_coordinates(v)[,2]
v = subset(v,X> -67 & X < -55 & Y>40 & Y<47)

l$DATE  =  as.Date(format(strftime(l$FV_FISHED_DATETIME,format="%Y-%m-%d"), tz="America/Halifax",usetz=TRUE))
v$DATE = as.Date(v$DATE)

lat <- l |>
  as_tibble() |>
  group_by(TRIP_ID,VESSEL_NAME) |>
  summarise(
    MinDate = min(DATE, na.rm = TRUE),
    MaxDate = max(DATE, na.rm = TRUE)
  ) |>
  ungroup()

lat = data.table(lat)
library(data.table)

setDT(lat)
setDT(v)


v[lat,
  TRIP_ID := i.TRIP_ID,
  on = .(
    VESSEL_NAME,
    DATE >= MinDate, 
    DATE <= MaxDate)
]


v_trip = subset(v, !is.na(TRIP_ID))
v_trip = st_as_sf(v_trip,coords=c('X','Y'),crs=4326)
ids = unique(v_trip$TRIP_ID)
da=c()
for(d in 1: length(ids)) {
  pts = subset(v_trip,TRIP_ID == ids[d])
  pts = pts[order(pts$VMSDATElocal),]
  if(nrow(pts)>4){
    x <- st_distance(
      pts[-nrow(pts), ],
      pts[-1, ],
      by_element = TRUE
    )
    x=sum(x)
    da=  rbind(da,c(ids[d],x))
  }
}

dad = as.data.frame(da)
dad$V2 = as.numeric(dad$V2)
names(dad) = c('TRIP_ID','Distance')
dak = merge(dad,lat)
dak = subset(dak,Distance < 1e7)
dak$YR = lubridate::year(dak$MinDate)
dak$MN = lubridate::month(dak$MinDate)
dak$DAYS = dak$MaxDate - dak$MinDate
dag = aggregate(Distance~YR+MN,data=subset(dak),FUN=median)

ggplot(dag,aes(x=MN,y=Distance/1000))+geom_point()+geom_line()+facet_wrap(~YR)
daga = aggregate(Distance~YR,data=subset(dak),FUN=median)
ggplot(daga,aes(x=YR,y=Distance/1000))+geom_point()+geom_line()


dak$DAYS = dak$DAYS+1
dlen = aggregate(DAYS~YR,data=subset(dak),FUN=median)


ggplot(dlen,aes(x=YR,y=DAYS))+geom_point()+geom_line()+labs(y='Length of Trip')



dak$kmpd = dak$Distance/as.numeric(dak$DAYS)

dkm = aggregate(kmpd~YR,data=subset(dak),FUN=median)

ggplot(dkm,aes(x=YR,y=kmpd/1000))+geom_point()+geom_line()+labs(y='KM per Day')



xy = read.csv(file=file.path(project.datadirectory('bio.lobster'),'analysis','indicators','commercialCatchrates.csv'))
v = merge(xy,dkm,by.x='fishingYear',by.y='YR')

library(ggrepel)
yrs <- sort(unique(v$fishingYear))

lab_yrs <- yrs[
  ((yrs - min(yrs)) %% 5 == 0) |
    yrs == max(yrs)
]

labs <- subset(v, fishingYear %in% lab_yrs)
ggplot(v, aes(kmpd, lL)) +
  geom_path() +
  geom_point() +
  geom_text_repel(
    data = labs,
    aes(label = fishingYear)
  )+
  labs(x='km per day',y='CPUE')+
  theme_test()
