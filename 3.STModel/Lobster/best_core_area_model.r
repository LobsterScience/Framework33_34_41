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
require(tidyr)
la()
fd=file.path(project.datadirectory('Framework_LFA33_34_41'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping

a = readRDS('Survey_data_package.rds')
aT = a[[1]]
bB = a[[2]]
cC = a[[3]]

aT = subset(aT,SOURCE %ni% c('MNR','NEFSC_RV') & lubridate::month(DATE)%in% c(6,7,8))

sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]
aT <- st_crop( aT, xmin = 124, xmax = 665, ymin = 4540, ymax = 5094)

survey = as_tibble(subset(aT,z>0 & !is.na(Lobster) & z<500))
survey$of = log(as.numeric(survey$OFFSET))

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



i = which(survey$SOURCE=='Snow crab survey')
survey$lz = log(survey$z)
path=file.path('Model_outputs/models_Nov2')
dir.create(path,recursive = T)

mod.label <- "m5"
  m <- sdmTMB(
    data = survey,
    formula = Lobster ~ 0+SOURCE+s(Glor)+s(lz),
    offset = survey$of,
    mesh = bspde,
    spatial = "on",
    family =  nbinom2(),
    time = "YEAR",
    spatiotemporal = "ar1"
  )
saveRDS(m,file=file.path(path,paste0('totlobster_',mod.label,'.rds')))
 
m = readRDS(file=file.path(path,paste0('totlobster_',mod.label,'.rds')))

re = residuals(m)
qqnorm(re)
qqline(re)


ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(survey$z) & z>5)

temp=readRDS(file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),"Glorys1994_2025wBiasCorrColumn_doy_oct19.rds"))
library(data.table)
setDT(temp)
#temp -glory = est
#glory+est =temp
# July + August
out <- temp[DOY >= 182 & DOY <= 243,.(
              temp = mean(EST+GLORYS_TEMP, na.rm = TRUE),
              depth = first(depth)
              ),by = .(X, Y)]

out = as.data.frame(out)
names(out) = c('X1000','Y1000','Glor','z')

out = st_as_sf(out,coords=c('X1000','Y1000'),crs=32620)

#downscaling
library(gstat)

# ensure same CRS
ff <- st_transform(ff, st_crs(out))

# IDW interpolation of Glor onto ff locations
idw.pred <- gstat::idw(
  formula = Glor ~ 1,
  locations = out,
  newdata = ff,
  idp = 2
)

ff$Glor <- idw.pred$var1.pred


ff <- ff |>
  tidyr::crossing(YEAR = unique(survey$YEAR)) |>
  st_as_sf()

sg = st_coordinates(ff)
ff$X1000 = sg[,1]
ff$Y1000 = sg[,2]
ff$SOURCE='ILTS'
f = as_tibble(ff)
f$lz = log(f$z)

f34 = subset(f,LFA==34 &z>=5)
g34 = predict(m,newdata=f34,return_tmb_object = T)
cog34 = get_cog(g34,bias_correct = T)
wcog34 <- cog34 |>
  pivot_wider(
    id_cols = c(YEAR, type),
    names_from = coord,
    values_from = c(est, lwr, upr, se),
    names_glue = "{coord}_{.value}"
  )

ind34g50 = get_index(g34,bias_correct = T)

a = lm(est_y~est_x,data=cog34)

b = lm(Y~X,data=coy)


 ggplot(subset(ind34g50),aes(x=YEAR,y=est,ymin=lwr,ymax=upr)+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')

f41 = subset(f,LFA==41 & z>=5)
g41 = predict(m,newdata=f41,return_tmb_object = T)
ind41g5 = get_index(g41,bias_correct = T)
ggplot(subset(ind41g5),aes(x=YEAR,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


f3538 = subset(f,LFA %in% 38 & z>=50)
g3538 = predict(model,newdata=f3538,return_tmb_object = T)
ind38g50 = get_index(g3538,bias_correct = T)
g8 = ggplo(subset(ind38g50,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


saveRDS(list(ind35,ind36,ind38),'IndicesFromFullComboModelSept26.rds')



f4 = subset(f,LFA %in% 34)
g3538 = predict(model,newdata=f4,return_tmb_object = T)
ind34 = get_index(g3538,bias_correct = T)
g8 = ggplot(subset(ind34,year<2024),aes(x=year,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')



###everywhere
g = predict(model,newdata=f,return_tmb_object = T)
h = g$data
h$pred = model$family$linkinv(h$est)
gsf = st_as_sf(h)


