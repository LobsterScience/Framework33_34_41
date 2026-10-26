require(bio.lobster)
require(devtools)
require(sf)
require(bio.utilities)
require(dplyr)
require(ggplot2)

la()

p = list()
p$yrs = 2026

#survey data
g = compileAbundPresAbs_vessel_corr(redo = F,size=F) #turn to redo=T for a new run
g = subset(g,!is.na(LATITUDE) | !is.na(LONGITUDE))

gs = st_as_sf(g,coords=c('LONGITUDE','LATITUDE'),crs=4326)
po = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFAPolysSF.rds')))
co = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
us = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','US_FishingAreas.rds')))
us = subset(us,Stock %ni% 'CAN')
us = st_transform(us,crs=4326)
sf_use_s2(FALSE) #needed for cropping

co = suppressWarnings(suppressMessages(st_crop(co,xmin=-74,ymin=40,xmax=-57,ymax=48)))
po = suppressWarnings(suppressMessages(st_crop(po,xmin=-74,ymin=40,xmax=-56,ymax=48)))
go = suppressWarnings(suppressMessages(st_crop(gs,xmin=-70,ymin=41,xmax=-62,ymax=46)))

#full data set
g = subset(g,SOURCE %in% c('ILTS','DFO_RV','NEFSC_RV','Snow Crab survey','Scallop Survey','MNR'))


#assign GLORYS temperatures
#ois = assignGlorys(x=g, temp=file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'Glorys1994-2025wBiasCorrColumn_doy_Sept24.rds'))
#survey data where I have glorys to match 

ois = readRDS(file=file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'lobsterData_withGlorys_oct19.rds'))



