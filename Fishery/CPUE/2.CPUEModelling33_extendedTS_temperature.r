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

m = readRDS(file='ExtendedTS_LFA33_logs.rds')
ca = m
 dir.create('L33CPUE')
ca$fYear = as.factor(ca$SYEAR)
ca$SOURCE = as.factor(ca$SOURCE)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca = subset(ca,!is.na(NUM_OF_TRAPS) & !is.na(WEIGHT_KG) & NUM_OF_TRAPS<=1200 & WEIGHT_KG>0& !is.na(bcT)& SYEAR<2026&DOS<186)
ca$CPUE = ca$WEIGHT_KG/ca$NUM_OF_TRAPS
redo.models=F
if(redo.models){
l33 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca),family = tw(link='log'),method='REML')
saveRDS(l33,file='L33CPUE/l33.rds')

l33b = bam(WEIGHT_KG~s(DOS)+fYear+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE)
saveRDS(l33b,file='L33CPUE/l33b.rds')

l33c = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l33c,file='L33CPUE/l33c.rds')

l33dt = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(bcT)+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l33dt,file='L33CPUE/l33dt.rds')

l33d = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l33d,file='L33CPUE/l33d.rds')

l33bt = bam(WEIGHT_KG~s(DOS)+s(bcT)+fYear+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE)
saveRDS(l33bt,file='L33CPUE/l33bt.rds')

l33ct = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+s(bcT)+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l33ct,file='L33CPUE/l33ct.rds')

l33dlt = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+(bcT)+s(GRID_NUM,bs='re')+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
saveRDS(l33dlt,file='L33CPUE/l33dlt.rds')


} else {
  
  l33 = readRDS(file='L33CPUE/l33.rds')
  l33b = readRDS(file='L33CPUE/l33b.rds')
  l33c = readRDS(file='L33CPUE/l33c.rds')
  l33dt = readRDS(file='L33CPUE/l33dt.rds')
  l33d = readRDS(file='L33CPUE/l33d.rds')
  l33bt = readRDS(file='L33CPUE/l33bt.rds')
  l33ct = readRDS(file='L33CPUE/l33ct.rds')
  l33dlt = readRDS(file='L33CPUE/l33dlt.rds')
  
}

library(dplyr)
library(purrr)

mods <- list(Base = l33,  BSD = l33b,BSDxY=l33c,BSDxYrG=l33d, BSDt = l33bt, BSDxYt = l33ct, BSDxYrGt = l33dt,BSDxYrGLt = l33dlt)
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
print(doc, target = "GAM_Model_SummaryL33.docx")

#w = readRDS('first4CPUEmodels33.rds')
#l33 = w[[1]]
#l33a = w[[2]]
#l33b = w[[3]]
#l33c = w[[4]]

library(gratia)
bcT_grid <- seq(min(ca$bcT), max(ca$bcT), length.out = 100)

mm <- sapply(bcT_grid, function(x){
            nd <- ca
            nd$bcT <- x
            nd$leffort = 0
        mean(predict( l33dt,newdata = nd, type = "response",exclude = "s(GRID_NUM)"))
        })


mml <- sapply(bcT_grid, function(x){
  nd <- ca
  nd$bcT <- x
  nd$leffort = 0
  mean(predict( l33dlt,newdata = nd, type = "response",exclude = "s(GRID_NUM)"))
})


plot(bcT_grid, mm,
     type = "l",
     xlab = "Bottom temperature",
     ylab = "Marginal predicted CPUE")

lines(bcT_grid, mml,
     type = "l",col='red')


#################trip weighted marginal means
aT=ca
ind = aggregate(SD_LOG_ID~DOS+SYEAR,data=subset(aT,LFA ==33),FUN=function(x) length(unique(x)))
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
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==33),FUN=function(x) length(unique(x)))
        ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
        names(ind1)[2] = 'SumTrips'
        ind = merge(ind,ind1)
        ind$prop = ind$SD_LOG_ID/ind$SumTrips
        ind$fYear=as.factor(ind$SYEAR)
        
        wts3 <- aggregate(SD_LOG_ID ~ DOS+GRID_NUM, data = ind, sum)
        wts3 <- subset(wts3,SD_LOG_ID>500)
        wts3$wt <- wts3$SD_LOG_ID / sum(wts3$SD_LOG_ID)

        
#grids commonly fished
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==33),FUN=function(x) length(unique(x)))
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
        
        
        ind = aggregate(SD_LOG_ID~DOS+SYEAR+GRID_NUM,data=subset(aT,LFA ==33),FUN=function(x) length(unique(x)))
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
	wt3 = merge(tdg,wts3)
	wts4 = merge(tdg,wts4)
	wts5 = merge(tdg,wts5)
	
	

####annual index with observed effort by year
get_index_effort <- function(model, ind, n_sim=1000, effort=1,name){
          
          mf  <- model.frame(model)
          yrs <- levels(mf$fYear)
          
          Xp <- predict(model, type="lpmatrix")
          b  <- coef(model)
          V  <- vcov(model)
          
          bsim <- MASS::mvrnorm(n=n_sim, mu=b, Sigma=V)
          
          res <- list()
          head(ind)
          for(y in yrs){
            
            nd=ind
            
            if(any(names(ind)=='SYEAR'))nd <- subset(ind, SYEAR == as.numeric(as.character(y)))
            
            nd$fYear  <- factor(y, levels=yrs)
            nd$SOURCE <- "marfis"
            nd$leffort <- log(effort)
            Xp <- predict(model,
                          newdata=nd,
                          type="lpmatrix",
                          exclude='s(GRID_NUM') #turn off the grid specific deviations and let the weighting scheme take over
            
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
        
l33p = get_index_effort(l33,name='Base',ind = wts3)
l33bp = get_index_effort(l33b,name='BSD',ind = wts3)
l33cp = get_index_effort(l33c,name='BSDxY',ind = wts3)
l33dp = get_index_effort(l33d,name='BSDxYrG',ind = wts3)
l33dpl = get_index_effort(l33dlt,name='BSDxYrGt',ind = wt3)

l33dp_core_grids = get_index_effort(l33dlt,name='BSDxYrG_co',ind = wts4)

l33dp_non_core_grids = get_index_effort(l33dlt,name='BSDxYrG_nc',ind = wts5)

com = bind_rows(list(l33p,l33bp,l33cp,l33dp,l33dpl,l33dp_core_grids,l33dp_non_core_grids))
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
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_cpue89-00.png')

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
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_cpue00-12.png')

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

ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_cpue13-25.png')


cap =as.data.frame(do.call(rbind,cpue.ann))
cap$Model='biased_corrected'
cap = subset(cap,select=c(SYEAR, Model, unBCPUE, l95,u95))
names(cap) = names(com)[c(1,7,2,5,6)]

com = bind_rows(com,cap)


write.csv(com,'LFA33MultipleModelsCPUE.csv')
#marginal mean, internally consistent



ggplot(subset(com,SYEAR<2026 & SYEAR>1994 & model %in% c('BSDxYrG','BSDxYrGt')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)



ggplot(subset(com,SYEAR<2026 & SYEAR>1994 & model %ni% c('BS','BSDxYrG_co','BSDxYrG_nc')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)


ggplot(subset(com,SYEAR<2026 & SYEAR>1994 & model %ni% c('BS','BSDxYrG_co','BSDxYrG_nc')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_marginalCPUE_commonseasonwts_main.png')


ggplot(subset(com,SYEAR<2026 & SYEAR>1994 & model %in% c('BSDxYrG','BSDxYrG_co','BSDxYrG_nc')), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_marginalCPUE_commonseasonwts_regrid.png')
