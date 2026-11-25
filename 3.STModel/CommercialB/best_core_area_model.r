
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

aT = subset(aT,SOURCE %ni% c('MNR', 'Scallop Survey'))

sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]
aT <- st_crop( aT, xmin = 124, xmax = 665, ymin = 4540, ymax = 5094)

survey = as_tibble(subset(aT,z>3 & !is.na(Legal_wt) & z<500))
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



survey$lz = log(survey$z)

path=file.path('Model_outputs/models_Nov13')
dir.create(path,recursive = T)
mod.label <- "m7"

m <- sdmTMB(
  data = survey,
  formula = Legal_wt ~ 0+SOURCE+s(Glor)+s(lz)+s(DYEAR,bs='cc'), 
  offset = survey$of,
  mesh = bspde,
  spatial = "on",
  family =  nbinom2(),
  time = "LYEAR",
  spatiotemporal = "ar1"
)
saveRDS(m,file=file.path(path,paste0('legallobster_',mod.label,'.rds')))
 
m = readRDS(file=file.path(path,paste0('legallobster_',mod.label,'.rds')))
survey$resid = residuals(m)

par(mfrow = c(2,2))

hist(survey$resid, breaks = 30)

qqnorm(survey$resid)
qqline(survey$resid)

plot(fitted(m), residuals(m))
lines(lowess(fitted(m), residuals(m)), col="red", lwd=2)

plot(survey$LYEAR, survey$resid)
abline(h = 0, lty = 2)

plot(survey$Glor, survey$resid)
plot(survey$z, survey$resid)
plot(survey$DYEAR, survey$resid)

ggplot(survey,
       aes(X1000, Y1000, colour = resid)) +
  geom_point() +
  scale_colour_gradient2()


survey$fit <- exp(predict(m)$est)
ggplot(survey,
       aes(fit, Legal_wt)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "red")




ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(survey$z) & z>5)

temp=readRDS(file.path(bio.lobster::project.datadirectory('bio.lobster.glorys'),"Glorys1994_2025wBiasCorrColumn_doy_oct19.rds"))
library(data.table)
setDT(temp)
if(SummerMonths){
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

}

##if you need to define dyear
if(DYEAR){
  temp[, DATE := as.Date(DOY - 1,
                         origin = paste0(YR, "-01-01"))]
    temp[, sept1 := as.Date(
    ifelse(
      format(DATE, "%m-%d") >= "09-01",
      paste0(YR, "-09-01"),
      paste0(YR - 1, "-09-01")
    )
  )]
  temp[, DYEAR := as.numeric(DATE - sept1) / 365.25]
  temp[, DYEAR_BIN := round(DYEAR*10) / 10]
  temp[, LYEAR := ifelse(month(DATE) >= 9,
                         year(DATE),
                         year(DATE) - 1)]
  temp[, Glory := EST+GLORYS_TEMP]
  
  temp = as.data.frame(temp)
  temp = subset(temp, select=c(X,Y,DYEAR_BIN,Glory,depth))
  names(temp) = c('X1000','Y1000','DYEAR','Glor','z')
  
  out = st_as_sf(temp,coords=c('X1000','Y1000'),crs=32620)
  #downscaling
  library(gstat)
  
  # ensure same CRS
  ff <- ff |>
    tidyr::crossing(
      DYEAR = seq(0, 1, by = 0.1)
    )
  ff <- st_transform(st_as_sf(ff), st_crs(out))
  
  idw.pred <- list()
    d = sort(unique(ff$DYEAR))
    for(i in 1:length(d)){
      x = subset(out, floor(DYEAR*10) == floor(d[i]*10))
      y = subset(ff, floor(DYEAR*10) == floor(d[i]*10))
      y$Glor  = gstat::idw(    Glor ~ 1,    locations =x ,  newdata   = y,  idp = 2)$var1.pred
      idw.pred[[i]] = y
     }
  
  dw = bind_rows(idw.pred)
  }
saveRDS(dw,'downscaledDYEAR_glorys4preds.rds')
dw = readRDS('downscaledDYEAR_glorys4preds.rds')

ff = dw

ff <- ff |>
  tidyr::crossing(LYEAR = unique(survey$LYEAR)) |>
  st_as_sf()

sg = st_coordinates(ff)
ff$X1000 = sg[,1]
ff$Y1000 = sg[,2]
ff$SOURCE='ILTS'
f = as_tibble(ff)
f$lz = log(f$z)

f34 = subset(f,LFA==34 &z>=5 & DYEAR == .8)
g34 = predict(m,newdata=f34,return_tmb_object = T)
g34I = get_index(g34,bias_correct = T)
cog34 = get_cog(g34,format = 'wide')
wcog34 <- cog34 |>
  pivot_wider(
    id_cols = c(YEAR, type),
    names_from = coord,
    values_from = c(est, lwr, upr, se),
    names_glue = "{coord}_{.value}"
  )
yrs = sort(unique(cog34$LYEAR))
lab_yrs <- yrs[
  ((yrs - min(yrs)) %% 5 == 0) |
    yrs == max(yrs)
]

labs <- subset(cog34, LYEAR %in% lab_yrs)
ggplot(cog34, aes(est_x, est_y)) +
  geom_path() +
  geom_point() +
  geom_text_repel(
    data = labs,
    aes(label = LYEAR)
  )+
  labs(x='Longitude',y = 'Latitude')+
  theme_test()




ind34g502 = get_index(g34,bias_correct = F)

a = lm(est_y~est_x,data=cog34)

b = lm(Y~X,data=coy)


 ggplot(subset(ind34g501),aes(x=LYEAR,y=est,ymin=lwr,ymax=upr))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000) ')

f41 = subset(f,LFA==41 & z>=5)
g41 = predict(m,newdata=f41,return_tmb_object = T)
ind41g5 = get_index(g41,bias_correct = F)
ggplot(subset(ind41g5),aes(x=LYEAR,y=est/1000,ymin=lwr/1000,ymax=upr/1000))+geom_point()+geom_line()+geom_ribbon(alpha=.25)+theme_test(base_size = 14)+labs(x='Year',y='Commercial Abundance (x000)')


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


