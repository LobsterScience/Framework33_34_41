require(bio.lobster)
require(devtools)
require(sf)
require(ggplot2)

la()

outdir = file.path(project.datadirectory('Framework_LFA33_34_41'),'Connectivity')


#refresh logbooks before running this
lobster.db('logs41')


#polygons
po = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFAPolysSF.rds')))
pou = st_transform(po,crs=32620)

logs41$DDLON = logs41$DDLON*-1
lo = st_as_sf(subset(logs41,!is.na(DDLON)),coords = c('DDLON','DDLAT'),crs=4326)
lou = st_transform(lo,crs=32620)

#prune out errors so logs are only in LFA 41

lous = st_join(lou,subset(pou,LFA==41))
lous = subset(lous,!is.na(LFA))

#distance to LFA 33, 34 and 40 for each set in LFA 41
pous = subset(pou,LFA %in% c(33,34,40))
pousl = st_cast(pous,'MULTILINESTRING')
dm = st_distance(lous,pousl)

lous$closest_LFA <- apply(dm, 1, function(x) pousl$LFA[which.min(x)])
lous$min_distance <- apply(dm, 1, min)
lous$yr = lubridate::year(lous$FV_FISHED_DATETIME)

require(dplyr)
lousF = lous %>%
  group_by(closest_LFA) %>%
  filter( min_distance <= quantile(min_distance, 0.95))

xx = aggregate(min_distance~closest_LFA+yr,data=subset(lousF,OFFAREA %ni% c('GBANK','GBASIN')),FUN=median)
xxn = aggregate(min_distance~closest_LFA+yr,data=lousF,FUN=length)

ggplot(subset(xx,yr<2025),aes(x=yr,y=min_distance/1000,colour=closest_LFA))+
  geom_point()+geom_line()+
  labs(x='Year',y='Median Distance to adjacent LFA Boundary (km)')+
  theme_test(base_size = 14)


ggplot(subset(lousF,yr<2025 & closest_LFA==33),aes(x=as.character(yr),y=min_distance/1000))+
  geom_boxplot()+
  labs(x='Year',y='Median Distance to adjacent LFA Boundary (km)')+
  theme_test(base_size = 14)+
  theme(axis.text.x  = element_text(angle=90))+
  facet_wrap(~closest_LFA)

ggplot(subset(lousF,yr<2025 & closest_LFA==34),aes(x=as.character(yr),y=min_distance/1000))+
  geom_boxplot()+
  labs(x='Year',y='Median Distance to adjacent LFA Boundary (km)')+
  theme_test(base_size = 14)+
  theme(axis.text.x  = element_text(angle=90))+
  facet_wrap(~closest_LFA)

ggplot(subset(lousF,yr<2025 & closest_LFA==40),aes(x=as.character(yr),y=min_distance/1000))+
  geom_boxplot()+
  labs(x='Year',y='Median Distance to adjacent LFA Boundary (km)')+
  theme_test(base_size = 14)+
  theme(axis.text.x  = element_text(angle=90))+
  facet_wrap(~closest_LFA)



