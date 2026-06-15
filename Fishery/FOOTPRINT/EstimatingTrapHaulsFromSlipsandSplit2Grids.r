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
aT = subset(aT,SYEAR>2005 & SYEAR<2026)
aT$GRID_NO = aT$GRID_NUM
aT$DOY = lubridate::yday(aT$DATE_FISHED)
a4 = lobster.db('process.logs41')
a4 = subset(a4,!is.na(ID) & yr>2005 & yr<2026)
a4$SYEAR = a4$yr
a4$DATE_FISHED = a4$FV_FISHED_DATETIME
a4$DOY = lubridate::yday(a4$FV_FISHED_DATETIME)
a4$WEIGHT_KG = a4$ADJCATCH_KG
u4 = unique(a4$ID)

gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
gr41 = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFA41_grid_polys.rds')))
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
gr41$LFA = as.character(gr41$LFA)
gr41 = subset(gr41,ID %in% u4 & ID %ni% c(860,1199))
names(gr41)[2] = c('GRID_NO')
gall = gtot = bind_rows(gr,gr41)
gtot$centroid = st_centroid(gtot)
cen_coords = st_coordinates(gtot$centroid)
gtot$X = cen_coords[,1] 
gtot$Y = cen_coords[,2] 
gtot$centroid <- NULL
gtot$geometry <- NULL
gtots = st_as_sf(gtot,coords=c('X','Y'),crs=4326)
gtot = subset(gtots,LFA %in% c(33,34,41))
gtot$dist_to_shore = as.numeric(st_distance(gtot,coa))/1000
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
gtot = subset(gtot,X< -60.5)
goo = gtot
gtot$geometry <- NULL

logs = subset(aT,select=c(LFA, GRID_NO, SYEAR, DATE_FISHED,WEIGHT_KG,NUM_OF_TRAPS,DOY))
l41 = subset(a4,select=c(LFA, ID, SYEAR, DATE_FISHED,WEIGHT_KG,NUM_OF_TRAPS,DOY))
st_geometry(l41) = NULL


layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")
r<-readRDS(file.path( layerDir,"GridPolysSF.rds"))
r = st_as_sf(r)

b = lobster.db('seasonal.landings')
b = subset(b,!is.na(SYEAR))
b$SYEAR = 1976:2026
b$LFA38B <- NULL
b = subset(b,SYEAR>2004 & SYEAR<2026)
b = reshape(b,idvar='SYEAR', varying=list(2:6),direction='long')
b$LFA=rep(c(33,34,35,36,38),each=21)
b$time <- NULL
names(b)[1:2]=c('YR','SlipLand')


d = lobster.db('annual.landings')
d = subset(d,YR>2004 & YR<2026, select=c(YR,LFA27,LFA28,LFA29,LFA30,LFA31A,LFA31B,LFA32))
d = reshape(d,idvar='YR', varying=list(2:8),direction='long')
d$LFA=rep(c(27,28,29,30,'31A','31B',32),each=21)
d$time <- NULL
names(d)[1:2]=c('YR','SlipLand')
bd = rbind(d,b)

bup = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~SYEAR+LFA,data=logs,FUN=sum)
bup$CPUE = bup$WEIGHT_KG/bup$NUM_OF_TRAPS
bAll = merge(bd,bup,by.x=c('YR','LFA'),by.y=c('SYEAR','LFA'))

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

l4 = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+ID+SYEAR,data=l41,FUN=sum)
names(l4)[c(2,4,5)]=c('GRID_NO','Landings','TrapHauls')
l4$LFA = as.character(l4$LFA)
l4 = dplyr::bind_rows(lan,l4)


l4g = merge(l4,gall)

saveRDS(l4g,'Landings_byGrid_year.rds')


