#Create Mesh Survey Data only

require(tidyr)
require(sdmTMB)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(PBSmapping)
require(SpatialHub)
require(sf)
la()
fd=file.path(project.datadirectory('Framework_LFA33_34_41'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping

###data in

crp = c(xmin = -71, ymin = 40.75, xmax = -57, ymax = 47.5)

#survey data
g  =compileAbundPresAbs_vessel_corr(size=F)
    addTemp2CompileAbun(g,temp.source='GLORYS')
g = subset(g, SOURCE %in% c("ILTS","DFO_RV", "NEFSC_RV","Snow crab survey","Scallop Survey","MNR") & YEAR>1999)
gs = st_as_sf(g,coords=c('LONGITUDE','LATITUDE'),crs=4326)
gs <- suppressWarnings(suppressMessages(
  st_crop(gs,
          crp)))
gs$X = st_coordinates(gs)[,1]
gs$Y = st_coordinates(gs)[,2]
i = which(gs$X< -70 & gs$Y<41.5)
j = which(gs$X< -69.4 & gs$Y <41)
ij=unique(c(i,j))
gs = gs[-ij,]
gs <- st_transform(gs, crs_utm20)
st_geometry(gs) = st_geometry(gs)/1000
st_crs(gs) <- crs_utm20


ns_coast =readRDS(file.path( bio.directory, "bio.lobster.data","mapping_data","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          crp)))
ns_coast <- st_transform(ns_coast, crs_utm20)
st_geometry(ns_coast) = st_geometry(ns_coast)/1000
st_crs(ns_coast) <- crs_utm20

ba = readRDS(file.path(bio.directory, "bio.lobster.data","mapping_data",'bathymetrySF.rds'))
ba = ba %>% st_as_sf() 
st_geometry(ba) = st_geometry(ba)/1000
st_crs(ba) = crs_utm20

rL = readRDS(file.path( bio.directory, "bio.lobster.data","mapping_data","LFAPolys37Split.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
rL = st_transform(rL,crs_utm20) 
st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
st_crs(rL) <- crs_utm20
rL = subset(rL,PID %in% c(30,311,312,32,33,34,35,36,38,40,41))

gs1 = st_join(gs, ns_coast)
gs = subset(gs1,is.na(COUNTRY), select=c(-COUNTRY, -PROVINCE))
bathy= ba
LFApolys = rL

ss = st_nearest_feature(gs,ba)
ds = st_distance(gs,ba[ss,],by_element=T)
st_geometry(ba) = NULL
gs$z = ba$z[ss]
gs$z_dist = as.numeric(ds)

gs$

  #lobster year Sept 1 as start of new year
decimal_year_sep1 <- function(date) {
    year <- year(date)
    mn = month(date)
    if(mn>=9){
    sep1 <- as.Date(paste0(year, "-09-01")) # Define September 1 of the same year
    next_sep1 <- as.Date(paste0(year + 1, "-09-01")) # Next year's September 1
    } else {
      sep1 <- as.Date(paste0(year-1, "-09-01")) # Define September 1 of the same year
      next_sep1 <- as.Date(paste0(year, "-09-01")) # Next year's September 1
      
    }
    # Compute fraction of the "year" (Sept 1 to next Sept 1)
    adjusted_year <- as.numeric(as.Date(date) - sep1) / as.numeric(next_sep1 - sep1)
    return(adjusted_year)
  }

lobster_year_sep1 <- function(date) {
  year <- year(date)
  mn = month(date)
  if(mn>=9){
    yr = year
    } else {
    yr = year-1
  }
  # Compute fraction of the "year" (Sept 1 to next Sept 1)
  return(yr)
}


# Apply function to vector
gs$DYEAR <- sapply(gs$DATE, decimal_year_sep1)
gs$LYEAR <- sapply(gs$DATE, lobster_year_sep1)
hist(gs$DYEAR)
  
  
saveRDS(list(gs,rL,ns_coast),'SurveyOnlyData.rds')

