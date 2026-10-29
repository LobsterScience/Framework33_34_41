
# source('/home/cooka@ENT/git/Framework33_34_41/3.STModel/PropLegal/multimodelPropLegal.r')


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

aT= subset(aT,SOURCE %ni% c('MNR','NEFSC_RV') & lubridate::month(DATE)%in% c(6,7,8) & !is.na(Legal) & z>0 & !is.na(Lobster))
#accounting for rounding with the coverserions

sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]

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
#or1 = cv_SpaceTimeFolds(or1,idCol = 'IDS', nfolds=5)
#moving to blockCV for folds
surv_sf <- st_as_sf(survey, coords = c("X1000", "Y1000"), crs = 32620)
set.seed(3)
blo <- cv_spatial(surv_sf, k = 5, size = 50)
cv_plot(blo, surv_sf)

or1$fold_id = blo$folds_ids

path=file.path('Model_outputs/models_Oct26')
dir.create(path,recursive = T)
#source(('C:/Users/cooka/Documents/git/Framework_LFA33_34_41/3.STModel/setupMultimodelTable.r'))
# source('/home/cooka@ENT/git/Framework33_34_41/3.STModel/setupMultimodelTable.r')
models = c('m1','m3','m5')
################################################################################################################################
if('m1' %in% models){
  mod.label <- "m1" 
  m <- sdmTMB(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE, 
    offset = or1$of,
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE, 
    offset = 'of',
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1",
    fold_ids='fold_id' ,
    k_folds = 5
  )
  m1 = m
  ca <-sdmTMB_cv_stats(response_var = 'Legal',denom_var = 'Lobster',offset_var = 'of')
  mod.select <- ca
  saveRDS(m,file=file.path(path,paste0('Proplegal_',mod.label,'.rds')))
  saveRDS(mod.select,file=file.path(path,'Proplegal_model_selection.rds'))
}


if('m3' %in% models){
  mod.label <- "m3" 
  m <- sdmTMB(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE+s(Glor), 
    offset = or1$of,
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE+s(Glor), 
    offset = 'of',
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1",
    fold_ids='fold_id',
    k_folds = 5
  )
  m3 = m
  ca <-sdmTMB_cv_stats(response_var = 'Legal',denom_var = 'Lobster',offset_var = 'of')
  mod.select <- bind_rows(mod.select, ca)
  saveRDS(m,file=file.path(path,paste0('Proplegal_',mod.label,'.rds')))
  saveRDS(mod.select,file=file.path(path,'Proplegal_model_selection.rds'))
}

if('m5' %in% models){
  mod.label <- "m5" 
  m <- sdmTMB(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE+s(Glor)+s(lz), 
    offset = or1$of,
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
  
  m_cv <- sdmTMB_cv(
    data = or1,
    formula = cbind(Legal,Lobster-Legal) ~ 0+SOURCE+s(Glor)+s(lz), 
    offset = 'of',
    mesh = bspde,
    spatial = "on",
    family =  binomial(link='logit'),
    time = "YEAR",
    spatiotemporal = "ar1",
    fold_ids='fold_id',
    k_folds = 5
  )
  m5 = m
  ca <-sdmTMB_cv_stats(response_var = 'Legal',denom_var = 'Lobster',offset_var = 'of')
  mod.select <- bind_rows(mod.select, ca)
  saveRDS(m,file=file.path(path,paste0('Proplegal_',mod.label,'.rds')))
  saveRDS(mod.select,file=file.path(path,'Proplegal_model_selection.rds'))
}
