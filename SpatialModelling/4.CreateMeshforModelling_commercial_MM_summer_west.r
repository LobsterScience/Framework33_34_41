#Create Mesh Survey Data only

require(tidyr)
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
la()
fd=file.path(project.datadirectory('Framework_LFA33_34_41'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping


a = readRDS('SurveyOnlyData.rds')
aT = a[[1]]
bB = a[[2]]
cC = a[[3]]

#crop

aT <- suppressWarnings(suppressMessages(
  st_crop(aT,
          c(xmin = 100, ymin = 4598, xmax = 430, ymax = 5100))))
ggplot(aT)+geom_sf()+geom_sf(data=aTa,colour='red')
aT = subset(aT,SOURCE %ni% 'MNR')
aT = subset(aT, lubridate::month(DATE) %in% c(5,6,7,8))
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]

survey = as_tibble(subset(aT,z>0 & !is.na(Legal)))
survey$of = log(as.numeric(survey$OFFSET))

spde <- make_mesh(survey, xy_cols = c("X1000", "Y1000"),
                  cutoff = 20)
#plot(spde)
#n vertices
spde$mesh$n

# Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, cC, range_fraction = .1,
  proj_scaling = 1, plot = F
)



i = which(survey$SOURCE=='Snow crab survey')
survey$Legal[i] = survey$Lobster[i]
survey$st = (survey$Glor-mean(survey$Glor)) / sd(survey$Glor)
survey$lz = log(survey$z)


survey$IDS = "I"
or1 = as_tibble(survey)
or1 = cv_SpaceTimeFolds(or1,idCol = 'IDS', nfolds=5)
path=file.path('Model_outputs/models_June23_26')
dir.create(path,recursive = T)
source(('C:/Users/cooka/Documents/git/Framework_LFA33_34_41/SpatialModelling/setupMultimodelTable.r'))
#source(file.path('~/git/Framework_LFA33_34_41/SpatialModelling/setupMultimodelTable.r'))

models = c('m1','m3','m5','m6')
################################################################################################################################
if('m1' %in% models){
  mod.label <- "m1" 
  m <- sdmTMB(
    data = or1,
    formula = Legal ~ SOURCE+s(lz), 
    offset = survey$of,
    mesh = bspde,
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = or1,
    formula = Legal ~ SOURCE+s(lz), 
    offset = survey$of,
    mesh = bspde,
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "YEAR",
    spatiotemporal = "ar1",
    fold_ids='fold_id',
    k_folds = 5
  )
  m1 = m
  ca <-mod.select.fn()
  mod.select <- rbind(mod.select, ca)
  saveRDS(m,file=file.path(path,paste0('commercialN_',mod.label,'.rds')))
  saveRDS(mod.select,file=file.path(path,'model_selection.rds'))
}



if('m2' %in% models){
  mod.label <- "m2" 
  
m2 <- sdmTMB(
  data = or1,
  formula = Legal ~ SOURCE+s(lz)+s(Glor), 
  offset = survey$of,
  mesh = bspde,
  spatial = "on",
  family =  delta_gamma(link1='logit',link2 = 'log'),
  time = "YEAR",
  spatiotemporal = "ar1"
    )

m_cv <- sdmTMB(
  data = or1,
  formula = Legal ~ SOURCE+s(lz)+s(Glor), 
  offset = survey$of,
  mesh = bspde,
  spatial = "on",
  family =  delta_gamma(link1='logit',link2 = 'log'),
  time = "YEAR",
  spatiotemporal = "ar1",
  fold_ids='fold_id',
  k_folds = 5
)
m2 = m
ca <-mod.select.fn()
mod.select <- rbind(mod.select, ca)
saveRDS(m,file=file.path(path,paste0('commercialN_',mod.label,'.rds')))
saveRDS(mod.select,file=file.path(path,'model_selection.rds'))
}
if('m3' %in% models){
  mod.label <- "m3" 
  m3 <- sdmTMB(
    data = or1,
    formula = Legal ~ SOURCE+s(lz), 
    offset = survey$of,
    mesh = bspde,
    spatial = "on",
    spatial_varying = ~(Glor),
    family =  delta_gamma(link1='logit',link2 = 'log'),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = or1,
    formula = Legal ~ SOURCE+s(lz), 
    offset = survey$of,
    mesh = bspde,
    spatial = "on",
    spatial_varying = ~(Glor),
    family =  delta_gamma(link1='logit',link2 = 'log'),
    time = "YEAR",
    spatiotemporal = "ar1",
    fold_ids='fold_id',
    k_folds = 5
  )
  m3 = m
  ca <-mod.select.fn()
  mod.select <- rbind(mod.select, ca)
  saveRDS(m,file=file.path(path,paste0('commercialN_',mod.label,'.rds')))
  saveRDS(mod.select,file=file.path(path,'model_selection.rds'))
}
