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
figpath = 'C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/FootPrints'
dir.create(figpath,recursive = T)

aT = lobster.db('process.logs')
a34 = subset(aT,LFA==34 &SYEAR>1998 & SYEAR<2026)
a33 = subset(aT,LFA==33 &SYEAR>2005 & SYEAR<2026)
arest = subset(aT,LFA %ni% c(33,34) &SYEAR>2005 & SYEAR<2026)
aT = bind_rows(a34,a33,arest)
aT$GRID_NO = aT$GRID_NUM
aT$DOY = lubridate::yday(aT$DATE_FISHED)
a4 = lobster.db('process.logs41')
a4 = subset(a4,!is.na(ID) & yr>2005 & yr<2026)
a4$SYEAR = a4$yr
a4$DATE_FISHED = a4$FV_FISHED_DATETIME
a4$DOY = lubridate::yday(a4$FV_FISHED_DATETIME)
a4$WEIGHT_KG = a4$ADJCATCH_KG
u4 = unique(a4$ID)

layerDir=file.path(project.datadirectory("bio.lobster"), "data","maps")
r<-readRDS(file.path( layerDir,"GridPolysSF.rds"))
gr = st_make_valid(st_as_sf(r))
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
b1 = b = subset(b, SYEAR<2026)
b = reshape(b,idvar='SYEAR', varying=list(2:6),direction='long')
b$LFA=rep(c(33,34,35,36,38),each=nrow(b1))
b$time <- NULL
names(b)[1:2]=c('YR','SlipLand')
sliplan = b

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

saveRDS(l4g,'Landings_byGrid_year_allLFAs.rds')

v = ggLobsterMap('34',addGrids = T)
gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
gr = subset(gr,LFA==34)
grc = st_centroid(gr)
lg = st_as_sf(merge(lan,gr))

v+geom_sf(data=subset(lg,SYEAR %in% 1999:2025),aes(fill=Landings))+    scale_fill_distiller(trans='sqrt',palette='Spectral') +
  facet_wrap(~SYEAR)+theme_test_adam()+  theme_test_adam() +
  theme(
    strip.background = element_rect(fill = "white", colour = "grey50" ),
      strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave(file.path(figpath,'LFA34Footprints.png'))




#which are the core grids
core <- lg %>%
  group_by(SYEAR, GRID_NO) %>%
  summarise(L = sum(Landings), .groups = "drop") %>%
  arrange(SYEAR, desc(L)) %>%
  group_by(SYEAR) %>%
  mutate(
    prop = L / sum(L),
    cumprop = cumsum(prop),
    core = cumsum(prop) - prop < 0.8
  )
core <- sf::st_make_valid(core)

persistence <- core %>%
  group_by(GRID_NO) %>%
  summarise(
    n_core = sum(core),
    prop_core = mean(core)
  )

ggplot(persistence,aes( fill=prop_core))+ geom_sf()+ scale_fill_distiller(trans='sqrt',palette='Spectral',name='Prop Core') +
  theme_test_adam()
ggsave(file.path(figpath,'ProportionYRSCore34.png'))

overlap <- core %>%
  filter(core) %>%
  select(SYEAR, GRID_NO) %>%
  group_by(SYEAR) %>%
  summarise(grids = list(GRID_NO))

#how many grids by year
overlapx = overlap %>%
  mutate(n_grids = lengths(grids)) %>%
  select(SYEAR, n_grids)
overlapx$LFA=34

#grids and total landings
overlapx = merge(overlapx,sliplan,by.x=c('SYEAR','LFA'),by.y=c('YR','LFA'))
yrs <- sort(unique(overlapx$SYEAR))

lab_yrs <- yrs[
  ((yrs - min(yrs)) %% 5 == 0) |
    yrs == max(yrs)
]

labs <- subset(overlapx, SYEAR %in% lab_yrs)
require(ggrepel)
ggplot(overlapx,aes(x=n_grids,y=SlipLand))+geom_point()+geom_path()+
  labs(x='Core Grids',y = 'Landings')+
  geom_text_repel(data = labs,aes(label = SYEAR))+
  theme_test()
ggsave(file.path(figpath,'NCore_2_landings34.png'))

########tracking landings in same grids since 1999
i = core$GRID_NO[which(core$SYEAR==1999 & core$core)]
core_sub = subset(core,GRID_NO %in% i)
core_land = aggregate(L~SYEAR,data=core_sub,FUN=sum)
tot_land = aggregate(L~SYEAR,data=subset(core, GRID_NO %ni% i),FUN=sum)
ggplot()+
  geom_path(data = tot_land,aes(x=SYEAR,y=L/1000), linewidth=1.1)+
  geom_path(data = core_land,aes(x=SYEAR,y=L/1000),colour='red',linewidth = 1.1)+
  labs(x="Year",y='Landings (t)')+
  theme_test(base_size = 14)
ggsave(file.path(figpath,'Core2NonCorelandings34.png'))

ggplot(subset(core,SYEAR==1999),aes( fill=core))+ geom_sf()+ scale_fill_distiller(trans='sqrt',palette='Spectral',name='Core') +
  theme_test_adam()+ theme(legend.position = "none")
ggsave(file.path(figpath,'CoreGrids1999_34.png'))

###cpue
v = aggregate(cbind(TrapHauls,Landings)~SYEAR,data=subset(l4,GRID_NO %in% i & LFA ==34),FUN=sum)
v$cpue = v$Landings/v$TrapHauls
w = aggregate(cbind(TrapHauls,Landings)~SYEAR,data=subset(l4, GRID_NO %ni% i &LFA ==34),FUN=sum)
w$cpue = w$Landings/w$TrapHauls
ggplot()+
  geom_path(data = w,aes(x=SYEAR,y=cpue),linewidth = 1.1)+
  geom_path(data = v,aes(x=SYEAR,y=cpue),colour='red',linewidth = 1.1)+
  labs(x="Year",y='CPUE (kg/th)')+
  theme_test(base_size = 14)
ggsave(file.path(figpath,'Core2noncorecpue34.png'))



####grids that are lost

jaccard <- function(a,b){
  length(intersect(a,b)) /
    length(union(a,b))
}

js = (data.frame(
  YEAR = overlap$SYEAR[-1],
  Similarity =
    mapply(
      jaccard,
      overlap$grids[-nrow(overlap)],
      overlap$grids[-1]
    )
))
jm = mean(js$Similarity)

ggplot(js,aes(x=YEAR,y=Similarity))+geom_point()+geom_path()+geom_hline(yintercept = jm,linetype='dashed')+labs(x='Year', y= 'Jaccard Similarity')+theme_test()
ggsave(file.path(figpath,'Jaccard34.png'))

library(dplyr)
library(purrr)

yrs <- sort(unique(core$SYEAR))
yrs <- yrs[yrs >= 2020]
library(dplyr)
library(purrr)

changes <- map_dfr(sort(unique(core$SYEAR[core$SYEAR > 2020])), function(y) {
  core <- core %>%
    ungroup()
  prev <- core %>%
    filter(SYEAR == y - 1, core) %>%
    pull(GRID_NO) %>%
    unique()
  
  curr <- core %>%
    filter(SYEAR == y, core) %>%
    pull(GRID_NO) %>%
    unique()
  
  tibble(
    SYEAR = y,
    GRID_NO = c(setdiff(prev, curr),
                setdiff(curr, prev)),
    Status = c(
      rep("Lost", length(setdiff(prev, curr))),
      rep("Gained", length(setdiff(curr, prev)))
    )
  )
})

library(dplyr)

turnover <- core %>%
  arrange(GRID_NO, SYEAR) %>%
  group_by(GRID_NO) %>%
  mutate(
    prev_core = lag(core),
    status = case_when(
      is.na(prev_core) ~ NA_character_,
      prev_core & core ~ "Stayed core",
      prev_core & !core ~ "Lost",
      !prev_core & core ~ "Gained",
      TRUE ~ "Not core"
    )
  ) %>%
  ungroup()

ggplot(filter(turnover, SYEAR %in% c(2001:2025))) +
  geom_sf(aes(fill = status), color = "grey40", linewidth = 0.1) +
  scale_fill_manual(
    values = c(
      "Lost" = "red",
      "Gained" = "dodgerblue",
      "Stayed core" = "grey70",
      "Not core" = "white"
    )
  ) +
  facet_wrap(~SYEAR)+
  theme_test_adam()+
  theme(
    strip.background = element_rect(fill = "white", colour = "grey50" ),
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave(file.path(figpath,'Gain_loss_core_34.png'))

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

ggplot(cog_year,aes(cog_lon, cog_lat))+geom_point()+geom_path()
library(ggrepel)
yrs <- sort(unique(cog_year$SYEAR))

lab_yrs <- yrs[
  ((yrs - min(yrs)) %% 5 == 0) |
    yrs == max(yrs)
]

labs <- subset(cog_year, SYEAR %in% lab_yrs)
ggplot(cog_year, aes(cog_lon, cog_lat)) +
  geom_path() +
  geom_point() +
  geom_text_repel(
    data = labs,
    aes(label = SYEAR)
  )+
  labs(x='Longitude',y = 'Latitude')+
  theme_test()
ggsave(file.path(figpath,'CenterofGravity_34.png'))


v+geom_path(data=cog_year,aes(x=cog_lon,y=cog_lat),col='red')

##########################################################

gr$GRID_NO = as.numeric(gr$GRID_NO)
#remove the islands or multipart polygons and keep only the biggest ones
gr <- gr %>%
  mutate(area = st_area(geometry)) %>%
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
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
gtot$geometry <- NULL
po = lobster.db('port_location')
po3 = subset(po, PORT_CODE %in% unique(aT$COMMUNITY_CODE) & LFA ==34)
correctPorts_aT = subset(aT, COMMUNITY_CODE %in% po3$PORT_CODE) 

cwg = merge(correctPorts_aT,gtot)
cwg = merge(cwg,po3[,c('PORT_CODE','CENTLAT','CENTLON')],by.x='COMMUNITY_CODE',by.y='PORT_CODE')
require(geosphere)

cwg$distance <- distHaversine(
  p1 = cwg[, c("X", "Y")],
  p2 = cwg[, c("CENTLON", "CENTLAT")]
)/1000

##center of gravity for both catch and effort
cc = split(subset(cwg,!is.na(NUM_OF_TRAPS)),f=subset(cwg,!is.na(NUM_OF_TRAPS))$SYEAR)
o = list()
for(i in 1:length(cc)){
  xcw = sum(cc[[i]]$X*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  xce = sum(cc[[i]]$X*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  ycw = sum(cc[[i]]$Y*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  yce = sum(cc[[i]]$Y*cc[[i]]$NUM_OF_TRAPS)/sum(cc[[i]]$NUM_OF_TRAPS)
  d2s = sum(cc[[i]]$distance*cc[[i]]$WEIGHT_KG)/sum(cc[[i]]$WEIGHT_KG)
  o[[i]] = c(unique(cc[[i]]$SYEAR),xce,xcw,yce,ycw,d2s)      
}

w = as.data.frame(do.call(rbind,o))
names(w) = c('SYEAR','Xef','Xca','Yef','Yca','Dist2Shore')

ggplot(w,aes(SYEAR, Dist2Shore))+geom_point()+
  labs(x='Year',y='Catch weighted Distance to Port')+
  theme_test(base_size = 14)
ggsave(file.path(figpath,'CatchWeightgedDist2Port_34.png'))


###########################################
####compared to historic
lf = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data',"LFAPolys37Split.rds"))
grids <- st_make_valid(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data',"Grids_updated_2026.rds")))
grids = subset(grids,LFA==34)
co = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','CoastSF.rds')))
co = suppressWarnings(suppressMessages(st_crop(co,xmin=-74,ymin=40,xmax=-57,ymax=48)))
ras = terra::rast(file.path(git.repo,'bio.lobster.data','mapping_data','bathymetryRaster.tif'))
grids = st_transform(grids,crs=32620)
cou = st_transform(co,crs=32620)
lf = st_transform(lf,crs=32620)
ff2 = as.data.frame(ras,xy=TRUE)
ff2 = subset(ff2,bathymetryRaster>3 &bathymetryRaster<600)
ff = st_as_sf(as.data.frame(ras,xy=TRUE),coords=c('x','y'),crs=32620)

f50 = subset(ff,bathymetryRaster< 55)
f34in34 = st_join(ff,lf)
f34in34 = subset(f34in34,PID==34)
cent <- st_centroid(grids)

d <- st_distance(cent, cou)

grids$dist_shore_km <- apply(d, 1, min) / 1000

grids$within20km <- grids$dist_shore_km <= 20

int <- st_intersection(
  grids |> select(GRID_NO),
  f34in34
)

dep <- int %>%
  st_drop_geometry() %>%
  group_by(GRID_NO) %>%
  summarise(
    mean_depth = mean(bathymetryRaster, na.rm = TRUE),
    median_depth = median(bathymetryRaster, na.rm = TRUE),
    prop_lt55 = mean(bathymetryRaster < 55, na.rm = TRUE),
    n = n()
  )

gridsless20 = subset(grids,within20km)
grids2 = subset(gridsless20, GRID_NO %in% dep$GRID_NO[which(dep$prop_lt55>0)])
ba=subset(f34in34, bathymetryRaster<55)
da =st_distance(ba,cou)
dist_shore_km <- apply(da, 1, min) / 1000

ba = ba[which(dist_shore_km<20),]
ggplot(ba)+geom_sf()


#1km resolution
dim(ba)[1]
####3982 km2 is this habitat for 90% of landings around 1985

lam = lobster.db('seasonal.landings')
lam$LFA34[which(lam$SYEAR=='1984-1985')]
#removals per km2
(5937-(5937*.11))/3982


grids$within20km <- grids$dist_shore_km <= 20

cf = coord_sf(xlim=c(144824.6,317770.3),ylim=c( 4727283.7 , 4979714.4 ))

ggplot()+
  geom_raster(data=f50,aes(x=x,y=y,fill=bathymetryRaster))+
  scale_fill_viridis_c()+
  geom_sf(data=grids, fill=NA,colour='orange',linewidth=1.2)+
  geom_sf(data=cou, fill='wheat')+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  cf

ggplot()+
  geom_raster(data=subset(ff2, bathymetryRaster<55),aes(x=x,y=y,fill=bathymetryRaster))+
  scale_fill_viridis_c()+
  geom_sf(data=grids, fill=NA,colour='orange',linewidth=1)+
  geom_sf(data=cou, fill='wheat')+
  geom_sf(data=grids2,fill=NA,colour='black', linewidth=1.01 )+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  cf

grids3 = merge(grids2,dep)

#merge gridded landings
lan4 = l4g |>
        as_tibble()  |>
        mutate(geometry = NULL)
lag2 = merge(grids3,lan4)
ugr = as.numeric(unique(subset(lag2,prop_lt55>.5)$GRID_NO))

lan4$Inshore = ifelse(lan4$GRID_NO%in% ugr, 1,0)

v =aggregate(cbind(Landings*Inshore,Landings)~SYEAR+LFA,data=subset(lan4,LFA==34),FUN=sum)
v$PropMidshore = 1-v$V1/v$Landings
vs = subset(v,select=c(SYEAR,PropMidshore))
vs = rbind(data.frame(SYEAR=c(1985,1995),PropMidshore=c(0.11,0.25)),vs)

ggplot(vs,aes(x=SYEAR,y=PropMidshore))+geom_point()+
  geom_point(data=subset(vs,SYEAR<1999),aes(x=SYEAR,y=PropMidshore),shape=8,size=3,col='red')+
  labs(x='Year',y='Proportion of Landings from Midshore areas')+
  theme_test( base_size = 14)
ggsave(file.path(figpath,'ProportionOfLandingsMidshore.png'))

ggplot()+
  geom_raster(data=subset(ff2),aes(x=x,y=y,fill=bathymetryRaster))+
  scale_fill_viridis_c()+
  geom_sf(data=grids, fill=NA,colour='orange',linewidth=1)+
  geom_sf(data=cou, fill='wheat')+
  geom_sf(data=subset(grids, GRID_NO %in% ugr),fill=NA,colour='black', linewidth=1.01 )+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  cf
