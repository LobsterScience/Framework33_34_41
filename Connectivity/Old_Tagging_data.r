#Tagging

require(bio.lobster)
require(devtools)
require(sf)
require(ggplot2)
require(dplyr)
require(tidyr)
require(stars)
require(raster)
la()
sf_use_s2(FALSE) #needed for cropping


###################################################
#mapping objects

#other mapping objects
po = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFAPolysSF.rds')))
co = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
us = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','US_FishingAreas.rds')))
us = subset(us,Stock %ni% 'CAN')
us = st_transform(us,crs=4326)

co = suppressWarnings(suppressMessages(st_crop(co,xmin=-74,ymin=40,xmax=-57,ymax=48)))

po = suppressWarnings(suppressMessages(st_crop(po,xmin=-74,ymin=40,xmax=-60,ymax=48)))

#bathymetry
#working with the bathymetry data
ras = terra::rast(file.path(git.repo,'bio.lobster.data','mapping_data','bathymetryRaster.tif'))

ras1 = ras
ras1[ras1<3 | ras1>600] <- NA

pou = st_transform(po,crs=32620)
cou = st_transform(co,crs=32620)
usu = st_transform(us,crs=32620)
zras1 = crop(ras1,extent(-153000,550000,4300000,5300000))
ff = as.data.frame(zras1,xy=TRUE)

cf = coord_sf(xlim=c(-149000,500000),ylim=c(4500000,5100000))

bmap = ggplot()+
  geom_raster(data=ff,aes(x=x,y=y,fill=bathymetryRaster))+
  scale_fill_viridis_c()+
  geom_sf(data=pou, fill=NA,colour='orange',linewidth=1.2)+
  geom_sf(data=cou, fill='wheat')+
  geom_sf(data=subset(usu), fill=NA,colour='red', linewidth=1.2)+
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
  cf

########################################################

#digitized data
aab =read.csv( file.path(project.datadirectory('bio.lobster'),'data','inputs','Tagging','Campbell1984_tags.csv'))
aab$Source='Campbell84'
aac = read.csv( file.path(project.datadirectory('bio.lobster'),'data','inputs','Tagging','Campbell1985_tags.csv'))
aac$Source='Campbell85'

aa = bind_rows(aac,aab)
a1 = subset(aa,select=c(TAGNO,Release_X,Release_Y))
a2 = subset(aa,select=c(TAGNO,Recapture_X,Recapture_Y))
names(a1)[1:3] = names(a2)[1:3] <- c('TAGNUM','LON','LAT')
a1$RelCap = 'Rel'
a2$RelCap = 'Cap'
a = bind_rows(a1,a2)
a$DATE = NA

x = read.csv(file.path(project.datadirectory('bio.lobster'),'data','inputs','Tagging','crtags.csv'))
x1 = subset(x,select=c('TAGNUM','TAGTYPE','RELDATE','RLAT','RLON','RCLMMS','RSEX','REGG'))
x2 =subset(x,select=c('TAGNUM','TAGTYPE','CAPDATE','CLAT','CLON','CAPCLMMS','CAPSEX','CAPEGG'))
x1$RelCap = 'Rel'
x2$RelCap = 'Cap'

names(x1) = names(x2) = c('TAGNUM','TAGTYPE','DATE','LAT','LON','CL','SEX','EGG','RelCap')

xa = bind_rows(x1,x2)
xa$TAGNUM = as.character(xa$TAGNUM)
xa$Source = 'CrTags'
xa = bind_rows(xa,a)

xa$LON = ifelse(xa$LON>0,  xa$LON*-1,xa$LON)
xa = subset(xa,!is.na(LAT) & LON< -60 & LAT>40)

xas = st_as_sf(xa,coords=c('LON','LAT'),crs=4326)

#remove on land
cc = st_join(xas,co)

xas = subset(cc,is.na(COUNTRY))
xas$DATE = as.Date(xas$DATE)
xasu=st_transform(xas,crs=32620)

#ours and SWLSS tagging
fb = "R:\\Science\\Population Ecology Division\\Shared\\!PED_Unit17_Lobster\\Lobster Unit Shared\\Projects and Programs\\Tagging\\Master_data"
v = readxl::read_xlsx(file.path(fb,"SWLSS_releases20240229.xlsx" ))
w = readxl::read_xlsx(file.path(fb,"XY_releases_2021.2022.2023.2024.2025_master.xlsx"  ))
wv = readxl::read_xlsx(file.path(fb,"LBT_RECAPTURES.xlsx"  ))

w$LAT_DD = w$LAT_DEGREES+w$LAT_MINUTES/60
w$LON_DD = (w$LON_DEGREES-(w$LON_MINUTES/60) )

v$LAT_DD = v$LAT_DEGREES+v$LAT_MINUTES/60
v$LON_DD = (v$LON_DEGREES-(v$LON_MINUTES/60)) * -1

#only recaps

wv = subset(wv,select=c(TAG_PREFIX,TAG_NUMBER, DAY, MONTH, YEAR, LAT_DD,LON_DD,CAPTURE_LENGTH,SEX))
wv$Relcap = 'Cap'
w = subset(w,select=c(TAG_PREFIX,TAG_NUM, DAY, MONTH, YEAR, LAT_DD,LON_DD,CARAPACE_LENGTH,SEX))
w$Relcap = 'Rel'

v = subset(v,select=c(TAG_PREFIX,TAG_NUM, CARAPACE_LENGTH, SEX,DAY,MONTH,YEAR,LON_DD,LAT_DD))
v$Relcap = 'Rel'

names(wv)[c(2,8)] = c('TAG_NUM','CARAPACE_LENGTH')

wv$id = paste(wv$TAG_PREFIX,wv$TAG_NUM,sep="_")
v$id = paste(v$TAG_PREFIX,v$TAG_NUM,sep="_")
w$id = paste(w$TAG_PREFIX,w$TAG_NUM,sep="_")

wv = wv %>% mutate_at(2:9, as.numeric)
w = w %>% mutate_at(2:9, as.numeric)
v = v %>% mutate_at(2:9, as.numeric)

vp = subset(v,id %in% unique(wv$id))
wp = subset(w,id %in% unique(wv$id))

cp = bind_rows(list(wv,wp,vp))
cp$DATE = as.Date(paste(cp$YEAR,cp$MONTH,cp$DAY,sep="-"))
pd <- cp %>%
      group_by(id) %>%
        filter(n()>1) %>%
      arrange(id,DATE)

pd$LON_DD = ifelse(pd$LON_DD<0,pd$LON_DD,pd$LON_DD*-1)

pds = st_as_sf(pd,coords=c('LON_DD','LAT_DD'),crs=4326)

pdd = st_transform(pds,crs=32620)

bmap+
  geom_sf(data=subset(pdd,Relcap=='Rel' ),colour='pink',size=0.8)+
 geom_sf(data=subset(pdd,Relcap=='Cap'),colour='red',size=0.8)+
  cf



##summary data

library(geosphere)
require(gdistance)

# Pair release and recapture events
paired_data <- xasu %>%
  group_by(TAGNUM) %>%
  filter(n() == 2) %>%  # Ensure each animal has both events
  arrange(DATE) %>%     # Sort by date
  summarise(
    release = first(geometry),
    recapture = last(geometry),
  #  release_date = first(DATE),
  #  recapture_date = last(DATE),
    distance_km = st_distance(first(geometry), last(geometry), by_element = TRUE) / 1000,  # Convert to km
  #  release_distance_to_shore = min(st_distance(first(geometry),cou)),
  # recapture_distance_to_shore = min(st_distance(last(geometry),cou)),
  #  duration_days = as.numeric(last(DATE) - first(DATE)),  # Time difference in days
#    direction = bearing(st_coordinates(first(geometry)), st_coordinates(last(geometry)))  # Calculate bearing
  )


#allocate mark to LFAs
po$Id = po$LFA
po$Stock = 'CAN'
po = subset(po,select=c(Id,Stock))
pd = st_join(paired_data,po)
pd = subset(pd,!is.na(Id))
pd = subset(pd,select=c(TAGNUM,Id,Stock))
st_geometry(pd) <- NULL
names(pd) = c('TAGNUM','Mark_ID','Mark_Stock')

#allocated recapture to Fishing areas
pr = paired_data
st_geometry(pr) = NULL
st_geometry(pr) = pr$recapture
prd = st_join(pr,po)
prd1 = subset(prd,is.na(Id),select=c(-Stock, -Id))
prd2 = subset(prd,!is.na(Id))
prd2 = subset(prd2,select=c(TAGNUM,Id,Stock))
st_geometry(prd2) <- NULL
names(prd2) = c('TAGNUM','Recap_ID','Recap_Stock')

prd1 = st_join(prd1,us)
prd1 = subset(prd1,select=c(TAGNUM,Id,Stock))
st_geometry(prd1) <- NULL
names(prd1) = c('TAGNUM','Recap_ID','Recap_Stock')
prd1$Recap_ID = as.character(prd1$Recap_ID)

prd = bind_rows(prd1,prd2)

pd = merge(pd,prd)

paired_data1=left_join(paired_data,as.data.frame(pd))


sPd = aggregate(TAGNUM~Mark_ID+Mark_Stock+Recap_ID+Recap_Stock,data=paired_data1,FUN=function(x) length(unique(x)))
names(sPd)[5] = 'N_Recaptures'

####Figures 
###########
#summary figures
##########
#locations of releases
bmap+
  geom_sf(data=subset(xasu,RelCap=='Rel'),colour='pink',size=0.8)+
  geom_sf(data=subset(xasu,RelCap=='Rel' & TAGNUM %in% unique(paired_data$TAGNUM)),colour='red',size=0.8)+
  cf

bmap+
  geom_sf(data=subset(xasu,RelCap=='Cap'),colour='pink',size=0.8)+
  cf



ggplot(sPd,aes(x=Mark_ID,y=Recap_ID,fill=N_Recaptures))+
  scale_fill_viridis_c(trans='log', labels = scales::label_number(accuracy = 1))+
  geom_raster()+
  theme_test(base_size = 14)+
labs(x='Release LFA', y='Recapture Fishing Area')


pos = subset(po, Id>32)
# Compute centroids

cents = readRDS( file.path(bio.directory,'bio.lobster.data','mapping_data',"LFALabelsSF.rds"))
cents = subset(cents,PID %in% c(33,34,35,36,37,38,40,41))
sf_cent_us <- st_centroid(us)

# Create the plot
ggplot() +
  geom_sf(data = us) +  # Adjust fill as needed
  geom_text(data = sf_cent_us, aes(x = st_coordinates(geometry)[,1], 
                                     y = st_coordinates(geometry)[,2], 
                                     label = Id), size = 2, color = "black") +
  geom_sf(data = pos) +  # Adjust fill as needed
  geom_text(data = cents, aes(x = st_coordinates(geometry)[,1], 
                              y = st_coordinates(geometry)[,2], 
                              label = label), size = 2, color = "black") +
  
  theme_test_adam()


###LFA 38B
ggplot() +
  geom_sf(data = subset(us,Id==511),fill='blue',alpha=.5) +  # Adjust fill as needed
  geom_text(data = subset(sf_cent_us,Id==511), aes(x = st_coordinates(geometry)[,1], 
                                   y = st_coordinates(geometry)[,2], 
                                   label = Id), size = 5, color = "black") +
  geom_sf(data = subset(pos,Id==38),fill='red',alpha=.5) +  # Adjust fill as needed
  geom_text(data = subset(cents,label==38), aes(x = st_coordinates(geometry)[,1], 
                              y = st_coordinates(geometry)[,2], 
                              label = label), size = 5, color = "black") +
  
  theme_test_adam()


#########################################################################################################################
#least cost distance

paired_data$releaseUTM = st_transform(paired_data$release,crs=32620)
paired_data$recaptureUTM = st_transform(paired_data$recapture,crs=32620)


 transition_function <- function(x){
            ifelse(x>350 | x<3,0.0001,1)
 }

findMoments(lo=3,up=400,dist = 'weibull')

#depth dependent transitions
 transition_function_depth <- function(x){
            a = dweibull(mean(x),1.018,110.98) / 0.007 #scalar 
            if(mean(x)<3) a = 0.000001 #less than 3m no go 
            return(a)
         }
bathy_trans <- transition(ras, transition_function_depth, directions = 16, symm = TRUE) 
bt = geoCorrection(bathy_trans,type='r')
#paths
outpaths = list()
paired_data$path = NA
for(i in 1:nrow(paired_data)){
        start_coords <- as.matrix(st_coordinates(paired_data$releaseUTM[i]))
        names(start_coords) = c('x','y')
        end_coords <- as.matrix(st_coordinates(paired_data$recaptureUTM[i]))
        names(end_coords) = c('x','y')
        path <- try(shortestPath(bt, start_coords, end_coords, output = "SpatialLines"),silent=T)
        if (inherits(path, "try-error")) {
          next  # Skip to the next iteration
        }
        outpaths[[i]] <- st_as_sf(path)
        paired_data$path[i] = st_length(outpaths[[i]])
}

oo = bind_rows(outpaths)

# Plot Results

ggplot()+
  geom_raster(data=ff,aes(x=x,y=y,fill=layer))+
            scale_fill_viridis_c()+
  geom_sf(data=oo,colour='red')

