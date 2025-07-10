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
                  cutoff = 12)
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
i = which(survey$SOURCE=='Snow crab survey')
survey$Legal[i] = survey$Lobster[i]
m <- sdmTMB(
  data = survey,
  formula = Legal ~ 0+SOURCE, 
  offset = survey$of,
  mesh = bspde,
  time_varying = ~ 1 + bs(z, knots=k3, degree=3, intercept=FALSE),
  time_varying_type = "rw0",
  spatial = "on",
  family =  tweedie(link = "log"),
  time = "YEAR",
  spatiotemporal = "ar1"
  
)

s_m = simulate(m,nsim=500,type='mle-mvn')
r_m = sdmTMB::dharma_residuals(s_m,m,return_DHARMa=T)
plot(r_m)
DHARMa::testResiduals(r_m)
DHARMa::testSpatialAutocorrelation(r_m,x=survey$X1000,y=survey$Y1000)
DHARMa::testZeroInflation(r_m)


survey$resids <- residuals(m) # randomized quantile residuals
qqnorm(survey$resids)
qqline(survey$resids)



ff = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","bathy_LFAPolysSF.rds"))
ff = subset(ff,z<max(survey$z) & z>5)
ff$X1000 = st_coordinates(ff)[,1]
ff$Y1000 = st_coordinates(ff)[,2]

yy = unique(survey$YEAR)
o = list()
for(i in 1:length(yy)){
  f = ff
  f$YEAR=yy[i]
  o[[i]]=f
}
ff = do.call(rbind,o)

ff$SOURCE='ILTS'
f = as_tibble(ff)

g = predict(m,newdata=f,offset=rep(log(1000000),times=nrow(f)))

g1 = m$family$linkinv(g$est)

f$pred = g1
gsf = st_as_sf(f,coords = c("X1000","Y1000"),crs=32620,remove=F)



#Maps
png('Figures/ModelOutput/lobstersdmTMBwk1-12.png', width = 10, height = 12,units='in',pointsize=12, res=300,type='cairo')
mm = c(0.00001,quantile(gsf$pred,0.999))
ggplot(subset(gsf,YEAR==2016)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c(limits=mm) +
  scale_color_viridis_c(limits=mm) +
#  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()
dev.off()


survey$LegalPA = ifelse(survey$Legal>0,1,0)
m1 = sdmTMB(
         data = survey,
	   formula = LegalPA ~ 0+SOURCE,
	   offset = survey$of,
	     mesh = bspde,
	     time_varying = ~ 1 + bs(z, knots=k3, degree=3, intercept=FALSE),
	       time_varying_type = "rw0",
	       spatial = "on",
	         family =  binomial(link = "logit"),
	         time = "YEAR",
		   spatiotemporal = "ar1"
		 
		 )

saveRDS(list(data=aT,grid=bspde,preds=be),file='results/dataForLFA33-35.rds')

