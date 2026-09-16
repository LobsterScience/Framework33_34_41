###########################
#start with an overview of all fishing year by year
require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(devtools)
require(sf)
require(ggplot2)
require(mgcv)
require(ggeffects)
require(ggforce)

require(sdmTMB)
require(purrr)
la()

fd = file.path(project.datadirectory('Framework_LFA33_34_41'))
setwd(fd)

m = readRDS(file='ExtendedTS_LFA34_logs.rds')
ca = m
 dir.create('L34CPUE')
ca$fYear = as.factor(ca$SYEAR)
ca$SOURCE = as.factor(ca$SOURCE)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca = subset(ca,!is.na(NUM_OF_TRAPS) & !is.na(WEIGHT_KG) & NUM_OF_TRAPS<=1200 & WEIGHT_KG>0& !is.na(bcT)& SYEAR<2026&DOS<186)
ca$CPUE = ca$WEIGHT_KG/ca$NUM_OF_TRAPS

l34 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca),family = tw(link='log'),method='REML')
saveRDS(l34,file='L34CPUE/l34.rds')

l34b = bam(WEIGHT_KG~s(DOS)+fYear+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE)
saveRDS(l34b,file='L34CPUE/l34b.rds')

l34c = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l34c,file='L34CPUE/l34c.rds')

l34dt = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(bcT)+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l34dt,file='L34CPUE/l34dt.rds')

l34dlt = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+(bcT)+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l34dt,file='L34CPUE/l34dlt.rds')


l34d = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l34d,file='L34CPUE/l34d.rds')

l34bt = bam(WEIGHT_KG~s(DOS)+s(bcT)+fYear+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE)
saveRDS(l34bt,file='L34CPUE/l34bt.rds')

l34ct = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(bcT)+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l34ct,file='L34CPUE/l34ct.rds')

#l34 = lb[[1]]
#l34b = lb[[2]]
#l34c = lb[[3]]
#saveRDS(list(l34d),'CPUEmodels34re_grid.rds')
#l34d=readRDS('CPUEmodels34re_grid.rds')[[1]]


library(dplyr)
library(purrr)

mods <- list(Base = l34,  BSD = l34b,BSDxY=l34c,BSDxYrG=l34d, BSDt = l34bt, BSDxYt = l34ct, BSDxYrGt = l34dt, BSDxYrGlt = l34dlt)
gam_table <- imap_dfr(mods, function(mod, name){
  
  s <- summary(mod)
  
  tibble(
    Model = name,
    Terms = paste(
      c(attr(terms(mod), "term.labels"),
        vapply(mod$smooth, `[[`, "", "label")),
      collapse = " + "
    ),
    EDF = sum(s$edf),
    Dev_Explained = round(100 * s$dev.expl, 1),
    AIC = round(AIC(mod), 1)
  )
})

library(officer)
library(flextable)

# Create a formatted table
ft <- flextable(gam_table)
ft <- autofit(ft)

# Create Word document
doc <- read_docx()
doc <- body_add_par(doc, "GAM Model Summary", style = "heading 1")
doc <- body_add_flextable(doc, ft)

# Save
print(doc, target = "GAM_Model_Summary.docx")

#w = readRDS('first4CPUEmodels34.rds')
#l34 = w[[1]]
#l34a = w[[2]]
#l34b = w[[3]]
#l34c = w[[4]]


#################trip weighted marginal means
aT=ca
ind = aggregate(SD_LOG_ID~DOS+SYEAR,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
names(ind1)[2] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fYear=as.factor(ind$SYEAR)


# Standardization weights based on trip distribution across all years
        wts <- aggregate(SD_LOG_ID ~ DOS, data = ind, sum)
        wts <- subset(wts,SD_LOG_ID>500)
        wts$wt <- wts$SD_LOG_ID / sum(wts$SD_LOG_ID)


# Standardization weights based on trip distribution within each years
        wts2 <- ind[,c('SYEAR','DOS','prop')]
        wts2$wt = ind$prop        
        
#wts by grid if included in the model        
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
        ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
        names(ind1)[2] = 'SumTrips'
        ind = merge(ind,ind1)
        ind$prop = ind$SD_LOG_ID/ind$SumTrips
        ind$fYear=as.factor(ind$SYEAR)
        
        wts3 <- aggregate(SD_LOG_ID ~ DOS+GRID_NUM, data = ind, sum)
        wts3 <- subset(wts3,SD_LOG_ID>500)
        wts3$wt <- wts3$SD_LOG_ID / sum(wts3$SD_LOG_ID)

        
#grids commonly fished
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
        grid_years <- aggregate(SYEAR ~ GRID_NUM,  data = unique(ind[, c("SYEAR", "GRID_NUM")]),  FUN = length)
        names(grid_years)[2] <- "n_years"
        grid_years <- grid_years[order(-grid_years$n_years), ]
        
        ind = subset(ind,GRID_NUM %in% c(92,103,139,140,157,156))
        ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
        names(ind1)[2] = 'SumTrips'
        ind = merge(ind,ind1)
        ind$prop = ind$SD_LOG_ID/ind$SumTrips
        ind$fYear=as.factor(ind$SYEAR)
        
        wts4 <- aggregate(SD_LOG_ID ~ DOS+GRID_NUM, data = ind, sum)
        wts4 <- subset(wts4,SD_LOG_ID>500)
        wts4$wt <- wts4$SD_LOG_ID / sum(wts4$SD_LOG_ID)
        
        
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
        grid_years <- aggregate(SYEAR ~ GRID_NUM,  data = unique(ind[, c("SYEAR", "GRID_NUM")]),  FUN = length)
        names(grid_years)[2] <- "n_years"
        grid_years <- grid_years[order(-grid_years$n_years), ]
        
        ind = subset(ind,GRID_NUM %ni% c(92,103,139,140,157,156))
        ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
        names(ind1)[2] = 'SumTrips'
        ind = merge(ind,ind1)
        ind$prop = ind$SD_LOG_ID/ind$SumTrips
        ind$fYear=as.factor(ind$SYEAR)
        
        wts5 <- aggregate(SD_LOG_ID ~ DOS+GRID_NUM, data = ind, sum)
        wts5 <- subset(wts5,SD_LOG_ID>500)
        wts5$wt <- wts5$SD_LOG_ID / sum(wts5$SD_LOG_ID)
        
 	tdg = aggregate(bcT~DOS+GRID_NUM,data=ca,FUN=mean)       
	wt3 = merge(tdg,wt3)

	td = aggregate(bcT~DOS,data=ca,FUN=mean)
	merge(
        
####annual index with observed effort by year
  get_index_effort_wt <- function(model, ind, n_sim=1000, effort=1,name){
          
          mf  <- model.frame(model)
          yrs <- levels(mf$fYear)
          
          Xp <- predict(model, type="lpmatrix")
          b  <- coef(model)
          V  <- vcov(model)
          
          bsim <- MASS::mvrnorm(n=n_sim, mu=b, Sigma=V)
          
          res <- list()
          
          for(y in yrs){
            
            nd=ind
            
            if(any(names(ind)=='SYEAR'))nd <- subset(ind, SYEAR == as.numeric(as.character(y)))
            
            nd$fYear  <- factor(y, levels=yrs)
            nd$SOURCE <- "marfis"
            nd$leffort <- log(effort)
            
            Xp <- predict(model,
                          newdata=nd,
                          type="lpmatrix")
            
            eta <- Xp %*% t(bsim)
            mu  <- exp(eta)
            
            idx <- colSums(mu * nd$wt)
            
            res[[y]] <- data.frame(
              SYEAR = as.numeric(y),
              mean  = mean(idx),
              median= median(idx),
              se    = sd(idx),
              lwr   = quantile(idx,.025),
              upr   = quantile(idx,.975),
              model = name
            )
          }
          
          bind_rows(res)
        }
        
l34p = get_index_effort_wt(l34,name='Base',ind = wts3)
l34bp = get_index_effort_wt(l34b,name='BSD',ind = wts3)
l34cp = get_index_effort_wt(l34c,name='BSDxY',ind = wts3)
l34dp = get_index_effort_wt(l34d,name='BSDxYrG',ind = wts3)

l34dp_core_grids = get_index_effort_wt(l34d,name='BSDxYrG_co',ind = wts4)

l34dp_non_core_grids = get_index_effort_wt(l34d,name='BSDxYrG_nc',ind = wts5)

com = bind_rows(list(l34p,l34bp,l34cp,l34dp,l34dp_core_grids,l34dp_non_core_grids))
####CPUE 

aa = split(ca,f=list(ca$LFA,ca$SYEAR))
cpue.lst<-list()
cpue.ann<- list()
for(i in 1:length(aa)){
  tmp<-aa[[i]]
  if(nrow(tmp)==0) next
  tmp = tmp[,c('DATE_FISHED','WEIGHT_KG','NUM_OF_TRAPS')]
  names(tmp)<-c('time','catch','effort')
  tmp$date<-as.Date(tmp$time)
  first.day<-min(tmp$date)
  tmp$time<-round(julian(tmp$date,origin=first.day-1))
  tmp = tmp[order(tmp$time),]
  g<-as.data.frame(biasCorrCPUE(tmp,by.time=T,min.sample.size = 3))
  g$lfa=unique(aa[[i]]$LFA)
  g$yr = unique(aa[[i]]$SYEAR)
  gl = aggregate(effort~time, data=tmp, FUN=sum)
  g = merge(g,gl,by.x='t',by.y='time')
  cpue.lst[[i]] <- g
  
  g<-as.data.frame(t(biasCorrCPUE(tmp,by.time=F)))
  g=data.frame(g,LFA=as.numeric(unique(aa[[i]]$LFA)),SYEAR=as.numeric(unique(aa[[i]]$SYEAR)))
  cpue.ann[[i]]=g
}
cc = bind_rows(cpue.lst)
library(ggplot2)

ggplot(
  subset(cc, yr %in% 1989:2000),
  aes(t, unBCPUE)
) +
  geom_point(
    colour = "grey50",
    alpha = 0.3,
    size = 0.5
  ) +
  geom_smooth(
    colour = "#2C7FB8",
    linewidth = 1.1,
    se = FALSE
  ) +
  facet_wrap(~yr, ncol = 4) +
  scale_x_continuous(
    breaks = c(0, 50, 100, 150),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Day of Season",
    y = "Catch Per Unit Effort"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "grey50"
    ),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA34_cpue89-00.png')

ggplot(
  subset(cc, yr %in% 2001:2012),
  aes(t, unBCPUE)
) +
  geom_point(
    colour = "grey50",
    alpha = 0.3,
    size = 0.5
  ) +
  geom_smooth(
    colour = "#2C7FB8",
    linewidth = 1.1,
    se = FALSE
  ) +
  facet_wrap(~yr, ncol = 4) +
  scale_x_continuous(
    breaks = c(0, 50, 100, 150),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Day of Season",
    y = "Catch Per Unit Effort"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "grey50"
    ),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA34_cpue00-12.png')

ggplot(
  subset(cc, yr %in% 2013:2025),
  aes(t, unBCPUE)
) +
  geom_point(
    colour = "grey50",
    alpha = 0.3,
    size = 0.5
  ) +
  geom_smooth(
    colour = "#2C7FB8",
    linewidth = 1.1,
    se = FALSE
  ) +
  facet_wrap(~yr, ncol = 4) +
  scale_x_continuous(
    breaks = c(0, 50, 100, 150),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Day of Season",
    y = "Catch Per Unit Effort"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "grey50"
    ),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA34_cpue13-25.png')


cap =as.data.frame(do.call(rbind,cpue.ann))
cap$Model='biased_corrected'
cap = subset(cap,select=c(SYEAR, Model, unBCPUE, l95,u95))
names(cap) = names(com)[c(1,7,2,5,6)]

com = bind_rows(com,cap)


write.csv(com,'LFA34MultipleModelsCPUE.csv')
#marginal mean, internally consistent

ggplot(subset(com,SYEAR<2026 & model %ni% c('BS','BSDxYrG_co','BSDxYrG_nc')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA34_marginalCPUE_commonseasonwts_main.png')


ggplot(subset(com,SYEAR<2026 & model %in% c('BSDxYrG','BSDxYrG_co','BSDxYrG_nc')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA34_marginalCPUE_commonseasonwts_regrid.png')
