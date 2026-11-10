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
gr = st_as_sf(r)
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
  )
ggsave(file.path(figpath,'CenterofGravity_34.png'))


v+geom_path(data=cog_year,aes(x=cog_lon,y=cog_lat),col='red')
###################
#lfa 33
v = ggLobsterMap('33',addGrids = T)
gr = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','GridPolys_DepthPruned_37Split.rds'))
gr = subset(gr,LFA==33)
grc = st_centroid(gr)
lg = st_as_sf(merge(lan,gr))

v+geom_sf(data=subset(lg,SYEAR %in% 2006:2025),aes(fill=Landings))+    scale_fill_distiller(trans='sqrt',palette='Spectral') +
  facet_wrap(~SYEAR)+theme_test_adam()+  theme_test_adam() +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "grey50"
    ),
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave(file.path(figpath,'LFA33Footprints.png'))

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
ggsave(file.path(figpath,'ProportionYRSCore33.png'))

overlap <- core %>%
  filter(core) %>%
  select(SYEAR, GRID_NO) %>%
  group_by(SYEAR) %>%
  summarise(grids = list(GRID_NO))



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
ggsave(file.path(figpath,'Jaccard33.png'))


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

ggplot(filter(turnover, SYEAR %in% c(2006:2025) & !is.na(status))) +
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
ggsave(file.path(figpath,'Gain_loss_core_33.png'))





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
  labs(x='Longitude',y='Latitude')+
  theme_test()
ggsave(file.path(figpath,'CenterofGravity_33.png'))


v+geom_path(data=cog_year,aes(x=cog_lon,y=cog_lat),col='red')



#####41


a4 = lobster.db('process.logs41')
a4$SYEAR = a4$yr
a4$DATE_FISHED = a4$FV_FISHED_DATETIME
a4$DOY = lubridate::yday(a4$FV_FISHED_DATETIME)
a4$WEIGHT_KG = a4$ADJCATCH_KG
u4 = unique(a4$ID)

gr41 = st_as_sf(readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFA41_grid_polys.rds')))
gr41$LFA = as.character(gr41$LFA)
gr41 = subset(gr41,ID %in% u4 & ID %ni% c(860,1199))
names(gr41)[2] = c('GRID_NO')
east <- gr41 %>%
  mutate(lon = st_coordinates(st_centroid(.))[,1]) %>%
  arrange(desc(lon)) %>%
  slice(1:6)
gr41 = subset(gr41,GRID_NO %ni% east$GRID_NO)
gtot = gr41


gtot$centroid = st_centroid(gtot)
cen_coords = st_coordinates(gtot$centroid)
gtot$X = cen_coords[,1] 
gtot$Y = cen_coords[,2] 
gtot$centroid <- NULL
gtot$geometry <- NULL
gtots = st_as_sf(gtot,coords=c('X','Y'),crs=4326)
gtot = subset(gtots,LFA %in% c(41))
gtot$dist_to_shore = as.numeric(st_distance(gtot,coa))/1000
gtot$X <- st_coordinates(gtot)[,1]
gtot$Y <- st_coordinates(gtot)[,2]
goo = gtot
gtot$geometry <- NULL

b = aggregate(WEIGHT_KG~SYEAR+ID,data=a4,FUN=sum)
bb=st_as_sf(merge(b,gr41,by.x='ID',by.y='GRID_NO'))
v= ggLobsterMap('41')
v+geom_sf(data=subset(bb,SYEAR %in% 1999:2025),aes(fill=WEIGHT_KG/1000))+    scale_fill_distiller(trans='sqrt',palette='Spectral') +
  facet_wrap(~SYEAR)+theme_test_adam()+  theme_test_adam() +
  theme(
    strip.background = element_rect(fill = "white", colour = "grey50" ),
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave(file.path(figpath,'LFA41Footprints.png'))

#how many grids for 80% for all LFAs

