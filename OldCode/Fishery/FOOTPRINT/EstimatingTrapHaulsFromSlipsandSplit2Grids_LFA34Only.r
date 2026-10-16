require(mgcv)
require(bio.lobster)
require(spdep)
require(sf)
require(devtools)
require(dplyr)
la()
require(bio.utilities)
load_all('C:/Users/Cooka/Documents/git/bio.utilities')

setwd(file.path(project.datadirectory('Framework_LFA33_34_41')))

aT = lobster.db('process.logs')
aT = subset(aT,SYEAR>=1999 & SYEAR<2026 & LFA==34)
aT$GRID_NO = aT$GRID_NUM
aT$DOY = lubridate::yday(aT$DATE_FISHED)
gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
coa = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
coa= st_make_valid(coa)
coa = subset(coa,PROVINCE=='Nova Scotia')
coa = subset(coa,st_area(coa)==max(st_area(coa))) # remove islands and cape breton

gr$GRID_NO = as.numeric(gr$GRID_NO)
#remove the islands or multipart polygons and keep only the biggest ones
gr1 <- gr %>%
  mutate(area = abs(st_area(geometry))) %>%
  group_by(GRID_NO, LFA) %>%
  slice_max(order_by = area, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-area)
gall = gtot = gr
gtot$centroid = st_centroid(gtot)
cen_coords = st_coordinates(gtot$centroid)
gtot$X = cen_coords[,1] 
gtot$Y = cen_coords[,2] 
gtot$centroid <- NULL
gtot$geometry <- NULL
gtots = st_as_sf(gtot,coords=c('X','Y'),crs=4326)
gtot = subset(gtots,LFA %in% c(34))
gtot$dist_to_shore = as.numeric(st_distance(gtot,coa))/1000
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
gtot = subset(gtot,X< -60.5)
goo = gtot
gtot$geometry <- NULL

logs = subset(aT,select=c(LFA, GRID_NO, SYEAR, DATE_FISHED,WEIGHT_KG,NUM_OF_TRAPS,DOY))


layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")
r<-readRDS(file.path( layerDir,"GridPolysSF.rds"))
r = st_as_sf(r)

b = lobster.db('seasonal.landings')
b = subset(b,!is.na(SYEAR))
b$SYEAR = 1976:2026
b$LFA38B <- NULL
b = subset(b,SYEAR>=1999 & SYEAR<2026, select=c(SYEAR,LFA34))
b$LFA=34
names(b)[1:2]=c('YR','SlipLand')


bup = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~SYEAR+LFA,data=logs,FUN=sum)
bup$CPUE = bup$WEIGHT_KG/bup$NUM_OF_TRAPS
bAll = merge(b,bup,by.x=c('YR','LFA'),by.y=c('SYEAR','LFA'))

sL= split(logs,f=list(logs$LFA, logs$SYEAR))
sL = rm.from.list(sL)
cpue.lst<-list()
cpue.ann = list()

  for(i in 1:length(sL)){
    tmp<-sL[[i]]
    tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
    names(tmp)<-c('time','catch','effort')
    tmp$date<-as.Date(tmp$time)
    first.day<-min(tmp$date)
    tmp$time<-julian(tmp$date,origin=first.day-1)
    g<-biasCorrCPUE(tmp,by.time = F)
    cpue.lst[[i]] <- c(lfa=unique(sL[[i]]$LFA),yr = unique(sL[[i]]$SYEAR),g)
  }
  
  cc =as.data.frame(do.call(rbind,cpue.lst))
  
cAll = merge(bAll,cc,by.x=c('LFA','YR'),by.y=c('lfa','yr'))

cAll$NTRAPs = cAll$SlipLand*1000/as.numeric(cAll$unBCPUE)
cAll$NTRAPSU = cAll$SlipLand*1000/as.numeric(cAll$l95)
cAll$NTRAPSL = cAll$SlipLand*1000/as.numeric(cAll$u95)


###########################################
#part the effort and landings to grids

partEffort = list()

for(i in 1:length(sL)){
  tmp = sL[[i]]
  tTH = aggregate(NUM_OF_TRAPS~LFA,data=tmp,FUN=sum)
  tC = subset(cAll, LFA==unique(tmp$LFA) & YR == unique(tmp$SYEAR)) 
  pTH = aggregate(NUM_OF_TRAPS~GRID_NO+LFA+SYEAR,data=tmp,FUN=sum)
  pTH$BTTH = pTH$NUM_OF_TRAPS / tTH$NUM_OF_TRAPS * tC$NTRAPs
  pTH$BlTH = pTH$NUM_OF_TRAPS / tTH$NUM_OF_TRAPS * tC$NTRAPSL
  pTH$BuTH = pTH$NUM_OF_TRAPS / tTH$NUM_OF_TRAPS * tC$NTRAPSU
  tTH = aggregate(WEIGHT_KG~LFA,data=tmp,FUN=sum)
  tC = subset(cAll, LFA==unique(tmp$LFA) & YR == unique(tmp$SYEAR)) 
  pTHa = aggregate(WEIGHT_KG~GRID_NO+LFA+SYEAR,data=tmp,FUN=sum)
  pTH$BL = pTHa$WEIGHT_KG / (tTH$WEIGHT_KG )* (tC$SlipLand*1000)
  
  partEffort[[i]] = pTH
}

partEffort = do.call(rbind, partEffort)

names(partEffort)[c(5,8)]=c('TrapHauls','Landings')
lan = subset(partEffort,select=c(LFA,GRID_NO,SYEAR,TrapHauls,Landings))


saveRDS(lan,'Landings_byGrid_year_34.rds')

v = ggLobsterMap('34',addGrids = T)
gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
gr = subset(gr,LFA==34)
grc = st_centroid(gr)
lg = st_as_sf(merge(lan,gr))

v+geom_sf(data=subset(lg,SYEAR %in% 1999:2025),aes(fill=Landings))+    scale_fill_distiller(trans='sqrt',palette='Spectral') +
  facet_wrap(~SYEAR)+theme_test_adam()

lg$cpue = lg$Landings/lg$TrapHauls
v+geom_sf(data=subset(lg,SYEAR %in% 1999:2025 & cpue<4),aes(fill=cpue))+    scale_fill_distiller(trans='sqrt',palette='Spectral') +
  facet_wrap(~SYEAR)+theme_test_adam()


grlc = st_as_sf(merge(lan,grc))

grlc$lon = st_coordinates(grlc)[,1]
grlc$lat = st_coordinates(grlc)[,2]

rlc=grlc

cog_year <- rlc %>%
  filter(!is.na(Landings), Landings > 0) %>%
  group_by(LFA, SYEAR) %>%
  summarise(
    cog_lon = sum(lon * Landings) / sum(Landings),
    cog_lat = sum(lat * Landings) / sum(Landings),
    total_landings = sum(Landings),
    n_cells = n_distinct(GRID_NO),
    .groups = "drop"
  )
