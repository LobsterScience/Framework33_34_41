
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

ca = subset(ca,SYEAR%in% 2003:2025 & LFA %ni% 41 & WEIGHT_KG>0 & !is.na(bcT) & NUM_OF_TRAPS>10)

gtos = subset(gto, GRID_NO %in% unique(ca$GRID_NO) & LFA %ni% 41) 
ca = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+DOS+X+Y+bcT+GRID_NO+z+dist_to_shore,data=ca,FUN=sum)
mes = sdmTMB::make_mesh(ca,xy_cols = c('X','Y'),n_knots=nrow(gtos)-1)


# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  mes, ns_coast,range_fraction = .2,
  proj_scaling = 1, plot = TRUE
)
cat = as_tibble(ca)
cat$leffort = log(cat$NUM_OF_TRAPS)
m3 = sdmTMB(WEIGHT_KG~ s(bcT)+s(DOS)+leffort
            ,
            data=cat,
            family = nbinom2(link='log') ,
            mesh = bspde,
            spatial='on',
            time='SYEAR',
           spatiotemporal='IID'
)

m4 = sdmTMB(WEIGHT_KG~ s(bcT)+s(DOS),
            offset = 'leffort',
            data=cat,
            family = nbinom2(link='log') ,
            mesh = bspde,
            spatial='on',
            time='SYEAR',
           spatiotemporal='IID'
)

#spatially varying smooth
X <- as.data.frame(smoothCon(	s(DOS, k = 5),	data = cat)[[1]]$X)
names(X) <- paste0("dos_basis_", seq_len(ncol(X)))
cat <- cbind(cat, X)


m5 = sdmTMB(WEIGHT_KG~ s(bcT),
            offset = 'leffort',
	   spatial_varying=~dos_basis_1+dos_basis_2+dos_basis_3+dos_basis_4,
            data=cat,
            family = nbinom2(link='log') ,
            mesh = bspde,
            spatial='on',
            time='SYEAR',
            spatiotemporal='IID'
)

#cAIC(m4)
#cAIC(m5)

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
X <- as.data.frame(smoothCon(   s(DOS, k = 5),  data = t1)[[1]]$X)
names(X) <- paste0("dos_basis_", seq_len(ncol(X)))
t1 <- cbind(t1, X)

 pre = merge(gtos,t1)


 require(purrr)
base_subsets <- map(1:3000, function(i) {
			pre %>%
			group_by(SYEAR) %>%
			slice_sample(n = 1, replace = TRUE) %>%
			ungroup()
																							                      })
sampled_ids <- bind_rows(base_subsets) %>% distinct()
remaining_df <- anti_join(pre, sampled_ids)

# Step 3: Randomly distribute remaining rows across the 200 subsets
remaining_split <- split(remaining_df, rep(1:3000, length.out = nrow(remaining_df)))

# Step 4: Combine base samples with remaining rows
final_subsets <- map2(base_subsets, remaining_split, bind_rows)

years = unique(cat$SYEAR)
for(i in 1:length(final_subsets)) {
	        fs = final_subsets[[i]]
          fs = subset(fs, SYEAR %in% years)
      	  g = predict(m5,newdata=fs,se_fit=F,offset=fs$leffort)
	        fs$pred = m5$family$linkinv(g$est)
		        final_subsets[[i]] = fs
		        saveRDS(fs, file=paste0('cpue_predictions',i,'.rds'))
			        rm(fs,g)
gc(reset=T)
}

fin = bind_rows(final_subsets)

saveRDS(fin,'compiled_cpue_predictions.rds')

fin = readRDS('compiled_cpue_predictions.rds')

#marginal temps
fin$temp = round(fin$bcT*2)/2
fi = aggregate(pred~temp,data=fin,FUN=function(x) quantile(x,c(0.25,0.5,0.75)))
ggplot(fi,aes(x=temp,y=pred[,2]/100,ymin=pred[,1]/100,ymax=pred[,3]/100))+
  geom_point(color='steelblue4')+
  geom_line(color='steelblue4')+
  geom_ribbon(fill='steelblue',alpha=.3)+
  labs(x='Temperature',y='Marginal CPUE')

#marginal dos
fi = aggregate(pred~DOS,data=fin,FUN=function(x) quantile(x,c(0.25,0.5,0.75)))
ggplot(fi,aes(x=DOS,y=pred[,2]/100,ymin=pred[,1]/100,ymax=pred[,3]/100))+
  geom_point(color='steelblue4')+
  geom_line(color='steelblue4')+
  geom_ribbon(fill='steelblue',alpha=.3)+
  labs(x='Temperature',y='Marginal CPUE')


