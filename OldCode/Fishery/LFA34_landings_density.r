
require(bio.lobster)
require(devtools)
require(sf)
require(ggplot2)
require(dplyr)
require(tidyr)
require(stars)
require(raster)
require(terra)
la()
sf_use_s2(FALSE) #needed for cropping

###historic landings time series

a = lobster.db('annual.landings.redo')
b = lobster.db('seasonal.landings.redo')
a = lobster.db('annual.landings')
b = lobster.db('seasonal.landings')


d = lobster.db('historic.landings')
d = subset(d, !is.na(LFA))
d$SYEAR = d$YEAR
d = aggregate(LANDINGS_MT~SYEAR+LFA,data=d,FUN=sum)
d = subset(d, SYEAR<1947)
names(d) = c('YR','LFA','LAND')

b$YR = substr(b$SYEAR,6,9)
a = subset(a,YR<1976)
b = subset(b,YR>1975 )
ff = "LFA34Update2026"
fpf1 = file.path(project.datadirectory('bio.lobster'),"analysis",ff)

LFA = c('LFA34')
for(i in LFA){
  aa = a[,c('YR',i)]
  bb = b[,c('YR',i)]
  dd = subset(d,LFA==34)
  dd$LFA <- NULL
  names(dd)[2] <- i
  aa = rbind(rbind(aa,bb),dd)
  aa = aa[order(aa$YR),]
  plot(aa$YR,aa[,i],type='h',lwd=4,col='black',xlab='Year',ylab='Landings (t)')
  lines(aa$YR,runmed(aa[,i],3),col='salmon',lwd=3)
}

##aa is landing time series





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

pou = st_transform(po,crs=32620)
cou = st_transform(co,crs=32620)
usu = st_transform(us,crs=32620)

v = readRDS('C:\\Users\\Cooka\\Documents\\bio_data\\bio.lobster\\data\\OFI_lobster_covariates/ofi_pointData500m.rds')
vs = st_as_sf(v,crs=32620)
vs$z = vs$depth *-1
vs = subset(vs,z< 500)
l = readRDS(file.path(git.repo,'bio.lobster.data/mapping_data/LFAPolysSF.rds'))
lp = st_transform(l,crs=32620)
pip <- st_join(vs, lp)
pip = subset(pip,!is.na(LFA))

ggplot()+
  geom_sf(data=subset(pip,LFA %in% c(33,34,40,41)),aes(fill=(z),colour=(z)))+
  geom_sf(data=subset(pou,LFA %in% c(33,34,40,41)),colour='orange',fill=NA)

#### make the 55m contour

v <- vect(pip)
r <- rast(ext(v),resolution=500)

rz <- rasterize(v, r, field = "z",resolution=500)
c55 <- as.contour(rz, levels = 55)

plot(rz)
lines(c55, col = "red", lwd = 2)

#just the line close to shore
cons = subset(cou,PROVINCE=='Nova Scotia')
library(sf)

coast <- st_union(cons)
shore <- st_boundary(coast)
shore_lines <- st_cast(shore, "MULTILINESTRING") |>
  st_cast("LINESTRING")
shore_lines$len <- st_length(shore_lines)
#longest line to get just shore
main_shore <- shore_lines[which.max(shore_lines$len), ]

c55_sf <- st_as_sf(c55,crs=32620)
ggplot()+geom_sf(data=main_shore)+geom_sf(data=c55_sf,colour='red')
##################################################################################################
#area between 3 and 55m for LFA 34 Duggan and pezzack id'ed in 

v34 <- vect(subset(pip,LFA==34))
r34 <- rast(ext(v34),resolution=500)

rz34 <- rasterize(v34, r34, field = "z",resolution=500)
c55_34 <- as.contour(rz34, levels = 55)

plot(rz34)
lines(c55_34, col = "red", lwd = 2)

#just the line close to shore
cons = subset(cou,PROVINCE=='Nova Scotia')
library(sf)

coast <- st_union(cons)
shore <- st_boundary(coast)
shore_lines <- st_cast(shore, "MULTILINESTRING") |>
  st_cast("LINESTRING")
shore_lines$len <- st_length(shore_lines)
#longest line to get just shore
main_shore <- shore_lines[which.max(shore_lines$len), ]

c55_sf <- st_as_sf(c55_34,crs=32620)
ggplot()+geom_sf(data=main_shore)+geom_sf(data=c55_sf,colour='red')

#area between  55m for LFA 34
shallow_34 <- rz34 <= 55 
sh34 <- st_as_sf(shallow_34,crs=32620)
sh34 <- as.polygons(shallow_34, dissolve = TRUE)
sh34 <- st_as_sf(sh34,crs=32620)
sh34 = subset(sh34,last==1)

# area (m2)
cell_area_34 <- cellSize(shallow_34, unit = "m")

area_m2_34 <- global(cell_area_34 * shallow_34, "sum", na.rm = TRUE)[1,1]

area_km2 <- area_m2_34 / 1e6

aa$Density = aa$LFA34 / area_km2
aa$Density[which(aa$YR>1972)]<-NA


##1995 onward
fd = file.path(project.datadirectory('Framework_LFA33_34_41'))
setwd(fd)
m = readRDS(file='ExtendedTS_LFA34_logs.rds')
mm = aggregate(WEIGHT_KG~GRID_NO+SYEAR,data=m,FUN=sum)

#which grids are the sme area

ggLobsterMap('34',addGrids = T,addGridLabels = T)+geom_sf(data=sh34,colour='red',fill=NA)
library(sf)
library(dplyr)
grids = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','Grids_updated_2026.rds'))
# make geometries valid
grids <- st_transform(st_make_valid(subset(grids,LFA==34)),crs=32620)

red   <- st_make_valid(sh34)

grids$grid_area <- st_area(grids)
int <- st_intersection(
  grids %>% select(GRID_NO, grid_area),
  red
)

cov <- int %>%
  mutate(red_area = st_area(geometry)) %>%
  st_drop_geometry() %>%
  group_by(GRID_NO) %>%
  summarise(red_area = sum(red_area))

grids2 <- grids %>%
  left_join(cov, by = "GRID_NO") %>%
  mutate(
    red_area = ifelse(is.na(red_area), units::set_units(0, m^2), red_area),
    prop_red = as.numeric(red_area / grid_area)
  )

majority_red <- grids2 %>%
  filter(prop_red > 0.35)

sl = aggregate(WEIGHT_KG~GRID_NO+SYEAR,data=subset(mm,SYEAR>1998 & GRID_NO %in% majority_red$GRID_NO),FUN=sum)

slm = merge(sl,majority_red)
slm$grid_area = as.numeric(slm$grid_area)
slm$red_area = as.numeric(slm$red_area)
dt = aggregate(cbind(WEIGHT_KG,grid_area,red_area)~SYEAR,data=slm,FUN=sum)
dt$grid_area = dt$grid_area/1e6
dt$red_area = dt$red_area/1e6

dt$Densityred = dt$WEIGHT_KG/dt$red_area/1000
dt$Densitygrid = dt$WEIGHT_KG/dt$grid_area/1000

ii <- is.na(aa$Density)

aa$Density[ii] <- dt$Densitygrid[
  match(aa$YR[ii], dt$SYEAR)
]

aa$YR <- as.numeric(aa$YR)

ggplot(aa, aes(YR, Density)) +
  geom_point() +
  scale_x_continuous(
    breaks = seq(1840, 2020, by = 20)
  ) +
  theme_bw()


