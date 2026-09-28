library(sf)
library(dplyr)
library(data.table)
library(dbscan)
library(lubridate)
library(adehabitatHR)
library(sp)

#=========================================================
# PREP DATA
#=========================================================

lobster.db('logs41') #make sure to do a database recapture through logs41.redo before moving on

logs41 = rename.df(logs41,c('FV_FISHED_DATETIME'),c('DATE_FISHED'))


off41 = subset(off41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJ_LOB_LBS))
ziff41 = subset(ziff41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJCATCH))
logs41 = subset(logs41,select=c(MON_DOC_ID,VR_NUMBER,DATE_FISHED,DDLAT,DDLON,NUM_OF_TRAPS,ADJCATCH))
names(off41) = names(ziff41) = names(logs41) 
#oct16-oct15 fishing year until 2005 switch to Jan 1 to Dec 31



a41 = do.call(rbind,(list(off41,ziff41,logs41)))
a41 = subset(a41,!is.na(ADJCATCH))
a41$DDLON = a41$DDLON*-1
a41$fishingYear = sapply(a41$DATE_FISHED,offFishingYear)
a41 = makePBS(a41,polygon=F)
a41$ADJCATCH = a41$ADJCATCH *0.45359237



dat2 <- subset(a41, !is.na(X) ) %>%
  mutate( t = as.numeric(DATE_FISHED -
                           min(DATE_FISHED))/(60*60*24))


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

sfdat$SOAK_DAYS = as.numeric(sfdat$SOAK_DAYS)


setDT(sfdat)
setorder(sfdat, VR_NUMBER, DATE_FISHED,TRIP_ID)

sfdat[, tripdur := cumsum(
  c(TRUE, diff(DATE_FISHED) > 1)
), by = VR_NUMBER]

trip_summary <- sfdat[
  ,
  {
    
    list(
      start_date = min(DATE_FISHED),
      end_date   = max(DATE_FISHED),
      trip_days  = as.numeric(max(DATE_FISHED) - min(DATE_FISHED)) + 1,
      n_records  = .N,
      traps      = sum(NUM_OF_TRAPS, na.rm = TRUE),
      catch_kg   = sum(ADJCATCH, na.rm = TRUE)
    )
  },
  by = .(VR_NUMBER, TRIP_ID)
]

ggplot(sfdat,aes(SOAK_DAYS))+geom_histogram()+facet_wrap(~fishingYear)


LFA41 = readRDS(file.path( git.repo,'bio.lobster.data','mapping_data',"LFA41Offareas.rds"))
lf = st_transform(LFA41,crs=32620)

sda = st_join(sfdat,lf)
