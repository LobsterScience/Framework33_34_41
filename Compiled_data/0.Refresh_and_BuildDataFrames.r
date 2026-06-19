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
groundfish.db('odbc.redo',datayrs = 2023:2025)
groundfish.db('gs_trawl_conversions_redo')
lobster.db("survey.redo")
lobster.db('scallop.redo')
nefsc.db('odbc.dump.redo')
nefsc.db('clean.redo')
ILTS_ITQ_All_Data(redo_set_data = T, redo_base_data = T, biomass = F)  

#need to udpate the MNR sets via https://mainedmr.shinyapps.io/MaineDMR_Trawl_Survey_Portal/'



g  =compileAbundPresAbs_vessel_corr(redo = F,size=F)
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
ggplot()+geom_sf(data=co) +geom_sf(data=subset(gs,Gear %in% c('Commercial','RecruitmentTrap')), aes(colour=as.factor(Empty)),size=.3)

#assign GLORYS temperatures

##assignGlorys(x=g, temp=file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'Glorys2000-2025wBiasCorrColumn_doy_june15.rds'))


#data for raph

g = subset(g,SOURCE %in% c('ILTS', 
                           'DFO_RV', 
                           'NEFSC_RV',
                            "Scallop Survey"
                           ,'MNR'
                           ))


saveRDS(g,'AllSurveyData.rds')


#####length data

sc = scallop_sets(length.group=1)
sc$EMPTY = ifelse(sc$Lobster==0,1,0)
sc$SOURCE = 'Scallop Survey'
sc$Gear = 'Dredge'
sc$OFFSET_METRIC="TowedDist x wing spread m2"
sc$DATE = sc$TOW_DATE
sc$LONGITUDE = sc$X
sc$LATITUDE = sc$Y
sc$id = paste('Scall',sc$TOW_SEQ,sep="_")

sc = sc %>%
  select(id, Lobster, Legal, Legal_wt, Berried, Recruit, Recruit_wt,YEAR, DATE,Juv, EMPTY, starts_with("P."), LONGITUDE, LATITUDE,SOURCE, OFFSET, OFFSET_METRIC, Gear)
ilts = ILTS_ITQ_All_Data(biomass = F,aggregate=F,redo_base_data = F)  
ilts$N = ilts$SA_CORRECTED_PRORATED_N * ilts$sweptArea ###covertt back to raw #s
#sc1=seq(13,253,by=1)
ilts$SZ = ilts$FISH_LENGTH
ilts$Berried= ilts$Recruit = ilts$Legal=  ilts$Juv =0
ilts$Berried = ifelse(ilts$SEX==3,ilts$N,ilts$Berried)
ilts$Recruit = ifelse(ilts$FISH_LENGTH %in% 70:81,ilts$N,ilts$Recruit)
ilts$Recruit = ifelse(ilts$FISH_LENGTH ==82,ilts$N/2,ilts$Recruit)
ilts$Juv = ifelse(ilts$FISH_LENGTH <=60,ilts$N,ilts$Juv)

ilts$Legal = ifelse(ilts$FISH_LENGTH >82,ilts$N,ilts$Legal)
ilts$Legal = ifelse(ilts$FISH_LENGTH ==82,ilts$N/2,ilts$Legal)
ilts$Legal_wt = (lobLW(CL=ilts$FISH_LENGTH,sex=ilts$SEX) * ilts$Legal)/1000

ilts$Recruit_wt = (lobLW(CL=ilts$FISH_LENGTH,sex=ilts$SEX) * ilts$Recruit)/1000

ilts$ID = paste(ilts$TRIP_ID,ilts$SET_NO,sep="_")
dA = aggregate(N~SZ+ID,data=ilts,FUN=sum)

dS = aggregate(cbind(Berried,Legal,N,Legal_wt,Recruit,Juv,Recruit_wt)~TRIP_ID+SET_NO+ID,data=ilts,FUN=sum)
dS$Lobster = dS$N
dS$N = NULL
dA$P = dA$N
aa = aggregate(P~ID+SZ,data=dA,FUN=sum)
bb = reshape(aa[,c('ID','SZ','P')],idvar='ID',timevar='SZ', direction='wide')
bb = na.zero(bb)

ca = merge(dS,bb)
set = ilts %>%
  distinct(TRIP_ID,SET_NO,SET_DATE,SET_LONG,SET_LAT,sweptArea,YEAR,temp,Length_comps)
ilt = merge(set,ca,all.x=T)
ilt$EMPTY = ifelse(ilt$Lobster>0,0,1)
ilt$id = paste(ilt$TRIP_ID,ilt$SET_NO,sep="_")
ilt$OFFSET = ilt$sweptArea *1e6
ilt$OFFSET_METRIC = "TowedDist x wing spread m2"
ilt$LONGITUDE = ilt$SET_LONG
ilt$LATITUDE = ilt$SET_LAT
ilt$DATE = ilt$SET_DATE
ilt$SOURCE = 'ILTS'
ilt$Gear = 'NEST'

p_cols <- grep("^P\\.", names(ilt), value = TRUE)
extra_cols <- c("Berried", "Legal", "Legal_wt", "Recruit", "Juv","Recruit_wt")

# Combine all target columns
target_cols <- c(p_cols, extra_cols)
ilt[ilt$Length_comps == 0, target_cols] <- NA   

ilt = ilt %>%
  #  filter(YEAR>1998) %>%
  select(id, Lobster, Legal, Legal_wt, Berried, Recruit,Recruit_wt, YEAR, DATE,Juv, EMPTY, starts_with("P."), LONGITUDE, LATITUDE,SOURCE, OFFSET, OFFSET_METRIC, Gear)
#RV Survey  ###covertt back to raw #s
rv = RV_sets(length_group = 1)
rv$TEMP = rv$bottom_temperature
rv = subset(rv,YEAR>1998)
rv$id = paste(rv$mission,rv$setno,sep='_')
rv$OFFSET = rv$OFFSET*1e6
rv$OFFSET_METRIC = "TowedDist x wing spread m2"

rv = rv %>%
  select(id,Lobster, Legal, Legal_wt,Berried, Recruit, Recruit_wt, Juv,YEAR,DATE,EMPTY,starts_with("P."), LONGITUDE, LATITUDE, SOURCE, OFFSET, OFFSET_METRIC,Gear )


xx = bind_rows(list(ilt,rv,sc))
xx = na.zero(xx)


saveRDS(xx,'AllSurveyData_all_lengths_may26.rds')


#with glorys bias correction surface june2026
ois = readRDS(file=file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'lobsterData_withGlorys.rds'))
 


