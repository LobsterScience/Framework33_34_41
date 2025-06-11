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
require(SpatialHub)
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


sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]

survey = as_tibble(subset(aT,z>0))
survey$of = log(survey$OFFSET)

spde <- make_mesh(survey, xy_cols = c("X1000", "Y1000"),
                  cutoff = 18)
plot(spde)
#n vertices
spde$mesh$n

# Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, cC, range_fraction = .1,
  proj_scaling = 1, plot = TRUE
)


mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]
ggplot() +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green")


k <- 3 # 4 basis splines
k3 <- Hmisc::wtd.quantile(survey$z, weights=survey$Legal, probs = seq(0, 1, len = k)[-c(1,k)]) ##quantiles of depth weighted by catch, to have splines at depth where catch located

k <- 4 # 5 basis splines
k4 <- Hmisc::wtd.quantile(survey$z, weights=survey$Legal, probs = seq(0, 1, len = k)[-c(1,k)]) ##quantiles of depth weighted by catch, to have splines at depth where catch located

m <- sdmTMB(
  data = survey,
  formula = Legal ~ 0+SOURCE, 
  offset = survey$off,
  mesh = bspde,
  time_varying = ~ 1 + bs(z, knots=k3, degree=3, intercept=FALSE),
  time_varying_type = "rw0",
  spatial = "on",
  family =  tweedie(link = "log"),
  time = "YEAR",
  spatiotemporal = "ar1"
  
)


survey$LO = log(survey$OFFSET)
fit = sdmTMB(Legal~
               s(lZ,k=5) + s(DYEAR,k=3),
             data=as_tibble(survey),
             offset = 'LO',
             time=YEAR, 
             mesh=bspde,
             extra_time = 17,
             family=tweedie(link='log'),
             spatial='on',
             spatiotemporal='ar1')

go =predict(fit) 
go$pred = fit$family$linkinv(go$est)



rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
st_crs(rL) <- 4326
crs_utm20 <- 32620
#rL = rL[-which(!(st_is_valid(rL))),]
rL <- suppressWarnings(suppressMessages(
  st_crop(rL,
          c(xmin = -68, ymin = 41, xmax = -56.5, ymax = 47.5))))
rL <- st_transform(rL, crs_utm20)
rL = st_union(rL)

baT <- baXY %>% st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
  st_crop(c(xmin = -68, ymin = 42, xmax = -57.5, ymax = 47)) %>%		
  st_transform(crs_utm20)				

baT = subset(baT,z>0 & z<400)
  #x
b = st_coordinates(baT)
baT$X1000 = b[,1]/1000
baT$Y1000 = b[,2]/1000
baT$X = b[,1]
baT$Y = b[,2]
baT$Depth = baT$z

ba = baT[,c('X','Y','Depth','X1000','Y1000')]
ba$geometry <- NULL
be = as.data.frame(sapply(ba,rep.int,27))
be$W = rep(0:26,each=dim(ba)[1])
be$geometry = NULL
be$lZ = log(be$Depth)
g = predict(fit,newdata=be)

g1 = fit$family$linkinv(g)

be$pred = apply(g1,1,median)
be$sd = apply(g1,1,sd)
be$lQ = apply(g1,1,quantile,0.25)
be$uQ = apply(g1,1,quantile,0.75)



gsf = st_as_sf(be,coords = c("X","Y"),crs=32620,remove=F)



#Maps
png('Figures/ModelOutput/lobstersdmTMBwk1-12.png', width = 10, height = 12,units='in',pointsize=12, res=300,type='cairo')
mm = c(0.001,max(gsf$pred))
ggplot(subset(gsf,W %in% 1)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c(trans='log',limits=mm) +
  scale_color_viridis_c(trans='log',limits=mm) +
  facet_wrap(~W) +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()
dev.off()



saveRDS(list(data=aT,grid=bspde,preds=be),file='results/dataForLFA33-35.rds')

