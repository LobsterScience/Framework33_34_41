require(bio.lobster)
require(devtools)
require(sf)
require(bio.utilities)
require(dplyr)
require(ggplot2)

la()

p = list()
p$yrs = 2026

#survey data


#adding in shell hardness
lobster.db('survey')
il = subset(surveyMeasurements,SPECCD_ID==2550 & !is.na(SHELL) ,select=c(TRIP_ID,SET_NO,FISHSET_ID,SET_DATE,SET_LAT,SET_LONG,SHELL, FISH_LENGTH,SEX,CALC_WT_G))
il$CALC_WT_G = lobLW(CL=il$FISH_LENGTH,sex = il$SEX)
library(dplyr)
library(tidyr)
il$id = paste(il$TRIP_ID,il$SET_NO,sep="_")
ilts_shell <- il %>%                 ``
group_by(
  id,
  SHELL
) %>%
  summarise(
    shell_wt = sum(CALC_WT_G, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    id) %>%
  mutate(
    prop_wt = shell_wt / sum(shell_wt)
  ) %>%
  ungroup() %>%
  select(
    id,SHELL, prop_wt
  ) %>%
  pivot_wider(
    names_from = SHELL,
    values_from = prop_wt,
    names_prefix = "shell_",
    values_fill = 0
  )


#vx = merge(v,shell_props_wide,by.x=c('Date','Y'),by.y=c('SET_DATE','SET_LAT'),all.x=T)

x = groundfish.db('special.lobster.sampling')
xx = groundfish.db('gs_trawl_conversions')
set = xx$gsinf
cas = xx$gscat
de = xx$gsdet
de = subset(de,spec==2550)

me = subset(x,select=c(flen,mission,setno,fsex,specimen_id,molt_stage))
mei = subset(me, !is.na(molt_stage) & molt_stage %ni% c(0,50))
mei$calwt = lobLW(CL=mei$flen,sex=mei$fsex)
mei$id = paste(mei$mission,mei$setno,sep='_')

library(dplyr)
library(tidyr)

molt_props <- mei %>%
  group_by(
    id,
    molt_stage
  ) %>%
  summarise(
    molt_wt = sum(calwt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    id) %>%
  mutate(
    prop_wt = molt_wt / sum(molt_wt)
  ) %>%
  ungroup() %>%
  select(
    id,
    molt_stage, prop_wt
  ) %>%
  pivot_wider(
    names_from = molt_stage,
    values_from = prop_wt,
    names_prefix = "shell_",
    values_fill = 0
  )

cxx = rbind(molt_props,ilts_shell)

v1 = merge(v,cxx,by='id',all.x=T)
saveRDS(v1,file = file.path(project.datadirectory('Framework_LFA33_34_41'),'Outputs','SURVEYS','SurveyOnlyData_shellhard.rds')) #sent to raph sept 9 2026

#does hardness affect catch, doing it local
v1$SOFT = v1$shell_1+v1$shell_2+v1$shell_7
v1$HARD = v1$shell_3+v1$shell_4+v1$shell_5+v1$shell_6

#cluster to 40km areas
v1 = st_as_sf(v1)
d <- st_distance(v1)
hc <- hclust(as.dist(d))
v1$cluster40 <- cutree(hc, h = 40)

library(lme4)

mod <- lmer(
  Lobster ~ SOFT + (1|cluster40),
  data = v1
)
#clusters
library(lwgeom)

hulls2 <- v1 %>%
  group_by(cluster40) %>%
  summarise() %>%
  st_convex_hull()

cc <- v1 %>%
  group_by(cluster40) %>%
  summarise(
    mean_soft = mean(SOFT, na.rm = TRUE),
    mean_lob = mean(Lobster, na.rm = TRUE),
    n = n()
  )
ggplot() +
  geom_sf(data = hulls2,
          fill ='red',
          colour = "black",
          linewidth = 0.3) +
  geom_sf(data = subset(hulls2,cluster40 %in% subset(cc,mean_soft>.1)$cluster40),
          fill ='blue',
          colour = "red",
          linewidth = 0.3) +
  geom_sf(data = v1,
          aes(),
          size = .2)


require(ggeffects)
library(ggeffects)
library(ggplot2)

p <- ggpredict(mod, terms = "SOFT")

ggplot(p, aes(x, predicted)) +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1.2,
            colour = "blue") +
  theme_bw() +
  labs(
    x = "Proportion Soft-Shelled",
    y = "Predicted Lobster Density",
    title = "Mixed Model Estimated Effect of Shell Softness"
  )
