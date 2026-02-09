
require(mgcv)
require(bio.lobster)
require(spdep)
require(sf)
require(devtools)
require(dplyr)
require(sdmTMB)
require(purrr)
la()

setwd(file.path(project.datadirectory('Framework_LFA33_34_41')))



m = readRDS(file='CPUE/data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
ca = m[[3]]


#one year for testing
ca = subset(ca,SYEAR%in% 2014:2024 & LFA %ni% 41 & WEIGHT_KG>0 & !is.na(bcT) & NUM_OF_TRAPS>10)
gtos = subset(gto, GRID_NO %in% unique(ca$GRID_NO))

mes = sdmTMB::make_mesh(ca,xy_cols = c('X','Y'),n_knots=nrow(gtos)-1)


# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  mes, ns_coast,range_fraction = .2,
  proj_scaling = 1, plot = TRUE
)


m4 = sdmTMB(WEIGHT_KG~ s(bcT)+s(DOS),
            offset = 'leffort',
            data=ca,
            family = nbinom2(link='log') ,
            mesh = bspde,
            spatial='on',
            time='SYEAR',
           spatiotemporal='IID'
)

#include a trap saturation effect

m5 = sdmTMB(WEIGHT_KG~ s(bcT)+s(DOS)+I(NUM_OF_TRAPS/(1+NUM_OF_TRAPS)),
            offset = 'leffort',
            data=ca,
            family = nbinom2(link='log') ,
            mesh = bspde,
            spatial='on',
            time='SYEAR',
            spatiotemporal='IID'
)

cAIC(m4)
cAIC(m5)

# b = visreg::visreg(m5,xvar='bcT',scale='response')
# t

 ##predictions
 gtos = subset(gto, LFA %ni% 41, select=c(LFA,GRID_NO,X,Y))
 temps = seq(quantile(ca$bcT,0.025),quantile(ca$bcT,0.975),length.out=10)
 dos = seq(1,max(ca$DOS),length.out=15)
 yr = seq(min(ca$SYEAR),max(ca$SYEAR))
 NUM_OF_TRAPS=100
 leffort=log(100)
 t1 = expand.grid(DOS=dos,bcT=temps,SYEAR=yr,NUM_OF_TRAPS=NUM_OF_TRAPS,leffort=leffort)
 pre = merge(gtos,t1)
 require(purrr)
base_subsets <- map(1:800, function(i) {
			pre %>%
			group_by(SYEAR) %>%
			slice_sample(n = 1, replace = TRUE) %>%
			ungroup()
																							                      })
sampled_ids <- bind_rows(base_subsets) %>% distinct()
remaining_df <- anti_join(pre, sampled_ids)

# Step 3: Randomly distribute remaining rows across the 200 subsets
remaining_split <- split(remaining_df, rep(1:800, length.out = nrow(remaining_df)))

# Step 4: Combine base samples with remaining rows
final_subsets <- map2(base_subsets, remaining_split, bind_rows)

years = unique(ca$SYEAR)
for(i in 1:length(final_subsets)) {
	        fs = final_subsets[[i]]
          fs = subset(fs, SYEAR %in% years)
      	  g = predict(m5,newdata=fs,se_fit=F,offset=fs$leffort)
	        fs$pred = m5$family$linkinv(g$est)
		        final_subsets[[i]] = fs
		        saveRDS(fs, file=paste0('cpue_predictions',i,'.rds'))
			        rm(fs,g)
}

fin = bind_rows(final_subsets)

saveRDS(fin,'compiled_cpue_predictions.rds')

fin = readRDS('compiled_cpue_predictions.rds')

fin$temp = round(fin$bcT*2)/2
fi = aggregate(pred~temp,data=fin,FUN=function(x) quantile(x,c(0.25,0.5,0.75)))
ggplot(fi,aes(x=temp,y=pred[,2]/100,ymin=pred[,1]/100,ymax=pred[,3]/100))+
  geom_point(color='steelblue4')+
  geom_line(color='steelblue4')+
  geom_ribbon(fill='steelblue',alpha=.3)+
  labs(x='Temperature',y='Marginal CPUE')


