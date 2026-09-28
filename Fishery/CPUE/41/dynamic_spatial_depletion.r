
require(bio.survey)
require(bio.lobster)
require(bio.utilities)
require(devtools)
la()

library(sf)
library(dplyr)
library(purrr)
library(broom)
library(dbscan)
library(slider)
library(lubridate)


lobster.db('logs41') #make sure to do a database recapture through logs41.redo before moving on

logs41 = rename.df(logs41,c('FV_FISHED_DATETIME'),c('DATE_FISHED'))


off41 = subset(off41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJ_LOB_LBS))
ziff41 = subset(ziff41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJCATCH))
logs41 = subset(logs41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJCATCH))
names(off41) = names(ziff41) = names(logs41) 
#oct16-oct15 fishing year until 2005 switch to Jan 1 to Dec 31



a41 = rbind(off41,ziff41,logs41)
a41 = subset(a41,!is.na(ADJCATCH))
a41$DDLON = a41$DDLON*-1
a41$fishingYear = sapply(a41$DATE_FISHED,offFishingYear)
a41 = makePBS(a41,polygon=F)
a41$ADJCATCH = a41$ADJCATCH *0.45359237



dat2 <- subset(a41, !is.na(X) & fishingYear>2004) %>%
    mutate( t = as.numeric(DATE_FISHED -
                     min(DATE_FISHED))/(60*60*24))


cell_size <- 10000      
gap_days  <- 21       

sfdat <- dat2 %>%
  filter(X< -64) %>%
  mutate(
    DATE_FISHED = as.Date(DATE_FISHED)
  ) %>%
  st_as_sf(
    coords = c("X","Y"),
    crs = 4326,
    remove = FALSE
  ) %>%
  st_transform(32620)

ln = readRDS(file.path(git.repo,'bio.lobster.data','mapping_data','LFAPolysSF.rds'))
lnp <- subset(ln,LFA==41) %>%
          st_transform(32620)

sfdat = sfdat %>%
      st_join(lnp) %>%
      filter(!is.na(LFA))



grid <- st_make_grid(
  sfdat,
  cellsize = cell_size
)

grid <- st_sf(
  cell = seq_along(grid),
  geometry = grid
)

sfdat <- st_join(  sfdat,  grid,  join = st_intersects)

#map
ggplot()+geom_sf(data=grid)+geom_sf(data=subset(grid,cell %in% sfdat$cell),fill='red')

episodes <- sfdat %>%
  st_drop_geometry() %>%
  arrange(cell, DATE_FISHED) %>%
  group_by(cell) %>%
  mutate(
    gap = c(0, diff(DATE_FISHED)),
    episode = cumsum(gap > gap_days)
  ) %>%
  ungroup()

#--------------------------------------------------
# DAILY CPUE
#--------------------------------------------------

daily_dat <- episodes %>%
  group_by(cell, episode, DATE_FISHED) %>%
  summarise(
    catch = sum(ADJCATCH, na.rm = TRUE),
    effort = sum(NUM_OF_TRAPS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cpue = catch / effort
  )

#--------------------------------------------------
# FIT ROBUST CPUE TREND WITHIN EACH EPISODE
#--------------------------------------------------

ep_split <- split(
  daily_dat,
  list(daily_dat$cell, daily_dat$episode),
  drop = TRUE
)

ep_trends <- lapply(ep_split, function(x){
  
  x <- x[order(x$DATE_FISHED), ]
  
  if(nrow(x) < 5) return(NULL)
  
  x$day <- as.numeric(
    x$DATE_FISHED - min(x$DATE_FISHED)
  )
  
  if(length(unique(x$cpue)) < 4) return(NULL)
  
  mod <- try(
    MASS::rlm(
      cpue ~ day,
      data = x
    ),
    silent = TRUE
  )
  
  if(inherits(mod, "try-error")) return(NULL)
  
  b <- coef(mod)
  
  start_cpue <- predict(
    mod,
    newdata = data.frame(day = 0)
  )
  
  end_cpue <- predict(
    mod,
    newdata = data.frame(day = max(x$day))
  )
  
  pct_change <- 100 *
    (end_cpue - start_cpue) /
    start_cpue
  
  data.frame(
    cell = x$cell[1],
    episode = x$episode[1],
    start_date = min(x$DATE_FISHED),
    end_date = max(x$DATE_FISHED),
    duration_days =
      as.numeric(max(x$DATE_FISHED) -
                   min(x$DATE_FISHED)) + 1,
    n_days = nrow(x),
    slope_cpue_day = b["day"],
    start_cpue = start_cpue,
    end_cpue = end_cpue,
    pct_change_cpue = pct_change
  )
})

ep_trends <- rbindlist(
  ep_trends,
  fill = TRUE
)

ep_trends[
  ,
  yr := year(start_date)
]

#--------------------------------------------------
# YEARLY SUMMARIES
#--------------------------------------------------

annual <- ep_trends[
  ,
  .(
    median_duration =
      median(duration_days, na.rm = TRUE),
    
    median_slope =
      median(slope_cpue_day, na.rm = TRUE),
    
    median_pct_change =
      median(pct_change_cpue, na.rm = TRUE),
    
    n_episodes = .N
  ),
  by = yr
]

ggplot(annual,aes(x=yr,y=median_pct_change))+geom_point()+geom_path()+labs(x='Year',y= "Median % CPUE Change Before Leaving")

cells_year <- sfdat %>%
  st_drop_geometry() %>%
  distinct(
    fishingYear,
    cell
  ) %>%
  count(fishingYear)

plot(
  n ~ fishingYear,
  data = cells_year,
  type = "b",
  xlab = "Year",
  ylab = "Unique grids fished"
)

search_area <- sfdat %>%
  st_drop_geometry() %>%
  distinct(
    fishingYear,
    cell
  ) %>%
  group_by(fishingYear) %>%
  summarise(
    n_cells = n()
  )

search_area$area_km2 <-
  search_area$n_cells *
  (cell_size/1000)^2

occ <- leave_cpue %>%
  group_by(yr) %>%
  summarise(
    med_days =
      median(duration_days),
    mean_days =
      mean(duration_days)
  )
plot(
  med_days ~ yr,
  data = occ,
  type = "b"
)


episodes_per_year <- leave_cpue %>%
  group_by(yr) %>%
  summarise(
    n_episodes = n()
  )

plot(
  n_episodes ~ yr,
  data = episodes_per_year,
  type = "b"
)

search_index <- merge(
  search_area,
  episodes_per_year,
  by.x="fishingYear",
  by.y="yr"
)

search_index <- merge(
  search_index,
  occ,
  by.x="fishingYear",
  by.y="yr"
)

search_index$searching_index <-
  search_index$n_episodes /
  search_index$med_days
