re(tidyr)
require(sdmTMB)
require(sdmTMBextra)
require(splines)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(PBSmapping)
require(sf)
require(blockCV)
la()
fd=file.path(project.datadirectory('Framework_LFA33_34_41'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping

a = readRDS('Survey_data_package.rds')
aT = a[[1]]
bB = a[[2]]
cC = a[[3]]

#aT= subset(aT,SOURCE %ni% c('MNR','NEFSC_RV') & lubridate::month(DATE)%in% c(6,7,8) & !is.na(Legal) & z>0 & !is.na(Lobster))
aT= subset(aT,SOURCE %ni% c('MNR') & !is.na(Legal) & z>0 & !is.na(Lobster) & z<500)
#accounting for rounding with the coverserions

sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]
aT <- st_crop( aT, xmin = 124, xmax = 665, ymin = 4540, ymax = 5094)

survey = as_tibble(aT)
survey$of = log(as.numeric(survey$OFFSET))
ii = which(survey$Legal >= survey$Lobster)
survey$Lobster[ii] = survey$Legal[ii]


spde <- make_mesh(survey, xy_cols = c("X1000", "Y1000"),
                  cutoff = 12)
#plot(spde)
#n vertices
spde$mesh$n

# Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, cC, range_fraction = .1,
  proj_scaling = 1, plot = F
)

survey$lz = log(survey$z)
survey$IDS = "I"
or1 = as_tibble(survey)

path=file.path('Model_outputs/models_Nov2_seasons')
dir.create(path,recursive = T)


mod.label <- "m6"
  m <- sdmTMB(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE+s(Glor)+s(lz)+s(DYEAR,bs='cc'),
    offset = or1$of,
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1"
  )


saveRDS(m,file=file.path(path,paste0('Proplegal_',mod.label,'.rds')))
  
