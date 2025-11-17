
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



m = readRDS(file='data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
ca = m[[3]]


#one year for testing
ca = subset(ca,SYEAR%in% 2014:2024 & LFA %ni% 41 & WEIGHT_KG>0)
gtos = subset(gto, GRID_NO %in% unique(ca$GRID_NO))

mes = sdmTMB::make_mesh(ca,xy_cols = c('X','Y'),n_knots=nrow(gtos)-1)


# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  mes, ns_coast,range_fraction = .2,
  proj_scaling = 1, plot = TRUE
)


require(sdmTMB)
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



# b = visreg::visreg(m4,xvar='DOS',scale='response')
# g = visreg::visreg(m4,xvar='GlorBC_mean',scale='response')

 ##predictions
 gtos = subset(gto, LFA %ni% 41, select=c(LFA,GRID_NO,X,Y))
 temps = seq(quantile(ca$GlorBC_mean,0.025),quantile(ca$GlorBC_mean,0.975),length.out=10)
 dos = seq(1,max(ca$DOS),length.out=15)
 yr = seq(min(ca$SYEAR),max(ca$SYEAR))
 
 t1 = expand.grid(DOS=dos,GlorBC_mean=temps,SYEAR=yr)
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
      	  g = predict(m4,newdata=fs,se_fit=F)
	        fs$pred = m4$family$linkinv(g$est)
		        final_subsets[[i]] = fs
		        saveRDS(fs, file=paste0('cpue_predictions',i,'.rds'))
			        rm(fs,g)
}

 mod = predict(m4, newdata=pre, return_tmb_object=T)
 
 
