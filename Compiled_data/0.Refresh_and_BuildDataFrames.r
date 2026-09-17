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

####sept 9 need to work through this to get an id for each station

g  =compileAbundPresAbs_vessel_corr(redo = F,size=F)
g = subset(g,!is.na(LATITUDE) | !is.na(LONGITUDE))
g = subset(g,SOURCE %in% c("ILTS", "DFO_RV", "NEFSC_RV","Scallop Survey")) 

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

#assign GLORYS temperatures


ois = assignGlorys(x=go, temp=file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),'Glorys2000-2025wBiasCorrColumn_doy_june15.rds'))


#data for raph 
saveRDS(ois,file.path(project.datadirectory('Framework_LFA33_34_41'),'Outputs','SURVEYS','SurveyOnlyData.rds'))

v = readRDS(file.path(project.datadirectory('Framework_LFA33_34_41'),'Outputs','SURVEYS','SurveyOnlyData.rds'))
v = subset(v,SOURCE %in% c('ILTS','DFO_RV') & lubridate::month(DATE) %in% c(6,7,8))
#adding in shell hardness
lobster.db('survey')
il = subset(surveyMeasurements,SPECCD_ID==2550 & !is.na(SHELL) ,select=c(TRIP_ID,SET_NO,FISHSET_ID,SET_DATE,SET_LAT,SET_LONG,SHELL, FISH_LENGTH,SEX,CALC_WT_G))
il$CALC_WT_G = lobLW(CL=il$FISH_LENGTH,sex = il$SEX)
library(dplyr)
library(tidyr)
il$id = paste(il$TRIP_ID,il$SET_NO,sep="_")
ilts_shell <- il %>%                 ``
  group_by(
id,
    SHELL
  ) %>%
  summarise(
    shell_wt = sum(CALC_WT_G, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
  id) %>%
  mutate(
    prop_wt = shell_wt / sum(shell_wt)
  ) %>%
  ungroup() %>%
  select(
    id,SHELL, prop_wt
  ) %>%
  pivot_wider(
    names_from = SHELL,
    values_from = prop_wt,
    names_prefix = "shell_",
    values_fill = 0
  )


#vx = merge(v,shell_props_wide,by.x=c('Date','Y'),by.y=c('SET_DATE','SET_LAT'),all.x=T)

x = groundfish.db('special.lobster.sampling')
xx = groundfish.db('gs_trawl_conversions')
set = xx$gsinf
cas = xx$gscat
de = xx$gsdet
de = subset(de,spec==2550)

me = subset(x,select=c(flen,mission,setno,fsex,specimen_id,molt_stage))
mei = subset(me, !is.na(molt_stage) & molt_stage %ni% c(0,50))
mei$calwt = lobLW(CL=mei$flen,sex=mei$fsex)
mei$id = paste(mei$mission,mei$setno,sep='_')

library(dplyr)
library(tidyr)

molt_props <- mei %>%
  group_by(
id,
    molt_stage
  ) %>%
  summarise(
    molt_wt = sum(calwt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
  id) %>%
  mutate(
    prop_wt = molt_wt / sum(molt_wt)
  ) %>%
  ungroup() %>%
  select(
    id,
    molt_stage, prop_wt
  ) %>%
  pivot_wider(
    names_from = molt_stage,
    values_from = prop_wt,
    names_prefix = "shell_",
    values_fill = 0
  )

cxx = rbind(molt_props,ilts_shell)

v1 = merge(v,cxx,by='id',all.x=T)
saveRDS(v1,file = file.path(project.datadirectory('Framework_LFA33_34_41'),'Outputs','SURVEYS','SurveyOnlyData_shellhard.rds')) #sent to raph sept 9 2026

#does hardness affect catch, doing it local
v1$SOFT = v1$shell_1+v1$shell_2+v1$shell_7
v1$HARD = v1$shell_3+v1$shell_4+v1$shell_5+v1$shell_6

#cluster to 40km areas
v1 = st_as_sf(v1)
d <- st_distance(v1)
hc <- hclust(as.dist(d))
v1$cluster40 <- cutree(hc, h = 40)

library(lme4)

mod <- lmer(
  Lobster ~ SOFT + (1|cluster40),
  data = v1
)
#clusters
library(lwgeom)

hulls2 <- v1 %>%
  group_by(cluster40) %>%
  summarise() %>%
  st_convex_hull()

cc <- v1 %>%
group_by(cluster40) %>%
summarise(
  mean_soft = mean(SOFT, na.rm = TRUE),
  mean_lob = mean(Lobster, na.rm = TRUE),
  n = n()
)
ggplot() +
  geom_sf(data = hulls2,
          fill ='red',
          colour = "black",
          linewidth = 0.3) +
  geom_sf(data = subset(hulls2,cluster40 %in% subset(cc,mean_soft>.1)$cluster40),
        fill ='blue',
        colour = "red",
        linewidth = 0.3) +
geom_sf(data = v1,
        aes(),
        size = .2)
  

require(ggeffects)
library(ggeffects)
library(ggplot2)

p <- ggpredict(mod, terms = "SOFT")

ggplot(p, aes(x, predicted)) +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1.2,
            colour = "blue") +
  theme_bw() +
  labs(
    x = "Proportion Soft-Shelled",
    y = "Predicted Lobster Density",
    title = "Mixed Model Estimated Effect of Shell Softness"
  )

#######################################################################################################
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

