# set up data from modelling
require(bio.lobster)
require(devtools)
require(sf)
require(PBSmapping)
require(bio.utilities)
require(dplyr)
require(ggplot2)

la()


lobster.db('logs41.redo')
lobster.db('greyzone_logs.redo')
lobster.db('fsrs.redo')
lobster.db('atSea.redo')
groundfish.db('odbc.redo',datayrs = 2023:2024)
groundfish.db('gs_trawl_conversions_redo')
lobster.db("survey.redo")
lobster.db('scallop.redo')
nefsc.db('odbc.dump.redo')
nefsc.db('clean.redo')
ILTS_ITQ_All_Data(redo_set_data = T, redo_base_data = T, biomass = F)  

#need to udpate the MNR sets via https://mainedmr.shinyapps.io/MaineDMR_Trawl_Survey_Portal/'



g  =compileAbundPresAbs_vessel_corr(size=F)
gs = st_as_sf(g,coords=c('LONGITUDE','LATITUDE'),crs=4326)
po = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFAPolysSF.rds')))
co = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
us = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','US_FishingAreas.rds')))
us = subset(us,Stock %ni% 'CAN')
us = st_transform(us,crs=4326)

co = suppressWarnings(suppressMessages(st_crop(co,xmin=-74,ymin=40,xmax=-57,ymax=48)))

po = suppressWarnings(suppressMessages(st_crop(po,xmin=-74,ymin=40,xmax=-56,ymax=48)))
ggplot()+geom_sf(data=co) +geom_sf(data=subset(gs,Gear %in% c('Commercial','RecruitmentTrap')), aes(colour=as.factor(Empty)),size=.3)

