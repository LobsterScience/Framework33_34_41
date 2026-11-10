require(mgcv)
require(bio.lobster)
require(sf)
require(devtools)
require(dplyr)
la()
require(bio.utilities)
load_all('C:/Users/Cooka/Documents/git/bio.utilities')

setwd(file.path(project.datadirectory('Framework_LFA33_34_41')))
figpath = 'C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/FootPrints'
dir.create(figpath,recursive = T)

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

ann_land = lobster.db('annual.landings')
ann_land = subset(ann_land,!is.na(LFA41),select=c(YR,LFA41))
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


core <- bb %>%
  group_by(SYEAR, ID) %>%
  summarise(L = sum(WEIGHT_KG), .groups = "drop") %>%
  arrange(SYEAR, desc(L)) %>%
  group_by(SYEAR) %>%
  mutate(
    prop = L / sum(L),
    cumprop = cumsum(prop),
    core = cumsum(prop) - prop < 0.8
  )
core <- sf::st_make_valid(core)

persistence <- core %>%
  group_by(ID) %>%
  summarise(
    n_core = sum(core),
    prop_core = mean(core)
  )

ggplot(persistence,aes( fill=prop_core))+ geom_sf()+ scale_fill_distiller(trans='sqrt',palette='Spectral',name='Prop Core') +
  theme_test_adam()
ggsave(file.path(figpath,'ProportionYRSCore41.png'))

overlap <- core %>%
  filter(core) %>%
  select(SYEAR, ID) %>%
  group_by(SYEAR) %>%
  summarise(grids = list(ID))



#how many grids by year
overlapx = overlap %>%
  mutate(n_grids = lengths(grids)) %>%
  select(SYEAR, n_grids)
overlapx$LFA=41

#grids and total landings
overlapx = merge(overlapx,ann_land,by.x=c('SYEAR'),by.y=c('YR'))
yrs <- sort(unique(overlapx$SYEAR))

lab_yrs <- yrs[
  ((yrs - min(yrs)) %% 5 == 0) |
    yrs == max(yrs)
]

labs <- subset(overlapx, SYEAR %in% lab_yrs)

ggplot(overlapx,aes(x=n_grids,y=LFA41))+geom_point()+geom_path()+
  labs(x='Core Grids',y = 'Landings')+
  geom_text_repel(data = labs,aes(label = SYEAR))+
  theme_test()
ggsave(file.path(figpath,'NCore_2_landings41.png'))

js = (data.frame(
  YEAR = overlap$SYEAR[-1],
  Similarity =
    mapply(
      jaccard,
      overlap$grids[-nrow(overlap)],
      overlap$grids[-1]
    )
))
jm1 = mean(js$Similarity[which(js$YEAR<=2012)])

jm2 = mean(js$Similarity[which(js$YEAR>2012)])

jm = mean(js$Similarity)

ggplot(subset(js,YEAR<2026),aes(x=YEAR,y=Similarity))+geom_point()+geom_path()+
  geom_hline(yintercept = jm,linetype='dashed',linewidth=1)+
  geom_segment(x=2002,xend=2012,y=jm1,yend=jm1,linetype='dashed',col='red',linewidth=1)+
  geom_segment(x=2013,xend=2025,y=jm2,yend=jm2,linetype='dashed',col='blue',,linewidth=1)+
    labs(x='Year', y= 'Jaccard Similarity')+theme_test()
ggsave(file.path(figpath,'Jaccard41.png'))


yrs <- sort(unique(core$SYEAR))
yrs <- yrs[yrs >= 2020]
library(dplyr)
library(purrr)

changes <- map_dfr(sort(unique(core$SYEAR[core$SYEAR > 2020])), function(y) {
  core <- core %>%
    ungroup()
  prev <- core %>%
    filter(SYEAR == y - 1, core) %>%
    pull(ID) %>%
    unique()
  
  curr <- core %>%
    filter(SYEAR == y, core) %>%
    pull(ID) %>%
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
  arrange(ID, SYEAR) %>%
  group_by(ID) %>%
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
ggsave(file.path(figpath,'Gain_loss_core_41.png'))





grlc = st_centroid(bb)

grlc$lon = st_coordinates(grlc)[,1]
grlc$lat = st_coordinates(grlc)[,2]

rlc=grlc

cog_year <- rlc %>%
  mutate(Landings = WEIGHT_KG) %>%
  filter(!is.na(WEIGHT_KG), WEIGHT_KG > 0) %>%
  group_by(LFA, SYEAR) %>%
  summarise(
    cog_lon = sum(lon * Landings) / sum(Landings),
    cog_lat = sum(lat * Landings) / sum(Landings),
    total_landings = sum(Landings),
    n_cells = n_distinct(ID),
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
ggsave(file.path(figpath,'CenterofGravity_41.png'))


v+geom_path(data=cog_year,aes(x=cog_lon,y=cog_lat),col='red')



