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

ca$fYear = as.factor(ca$SYEAR)
ca$SOURCE = as.factor(ca$SOURCE)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca = subset(ca,!is.na(NUM_OF_TRAPS) & !is.na(WEIGHT_KG) & NUM_OF_TRAPS<=500 & WEIGHT_KG>0 & NUM_OF_TRAPS>20)
ca$CPUE = ca$WEIGHT_KG/ca$NUM_OF_TRAPS
l33 = gam((WEIGHT_KG)~fYear+offset(leffort),data=subset(ca),family = tw(link='log'),method='REML')
l33a = bam((WEIGHT_KG)~fYear+SOURCE+offset(leffort),data=(ca),family = tw(link='log'),method='fREML',discrete=TRUE)
l33b = bam((WEIGHT_KG)~s(DOS)+fYear+(SOURCE)+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE)
l33c = bam((WEIGHT_KG)~s(DOS,fYear,bs='fs')+SOURCE+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))
l33d = bam((WEIGHT_KG)~s(DOS,fYear,bs='fs')+s(GRID_NUM,bs='re')+SOURCE+offset(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))


#l33e = bam(WEIGHT_KG~s(DOS,fYear,bs='fs')+SOURCE+s(leffort),data=ca,family = tw(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))

#l33d = gam(WEIGHT_KG~s(DOS)+s(DOS,by=fYear)+fYear+s(bcT)+offset(leffort),data=subset(ca,LFA==33),family = Gamma(link='log'),method='REML')

#saveRDS(list(l33,l33a,l33b,l33c),'first4CPUEmodels33.rds')
#w = readRDS('first4CPUEmodels33.rds')
#l33 = w[[1]]
#l33a = w[[2]]
#l33b = w[[3]]
#l33c = w[[4]]

li = list(l33,l33b,l33c)
saveRDS(li,'first4CPUEmodels33.rds')


library(dplyr)
library(purrr)

mods <- list(Base = l33,  BSD = l33b,BSDxY=l33c)
library(dplyr)
library(purrr)

mods <- list(
  Base  = l33,
  BSD   = l33b,
  BSDxY = l33c
)

# First collect statistics
gam_table <- imap_dfr(mods, function(mod, name) {
  
  s <- summary(mod)
  
  p_est <- tryCatch(
    mod$family$getTheta(TRUE),
    error = function(e) NA_real_
  )
  
  tibble(
    Model = name,
    EDF = round(sum(s$edf), 1),
    Dev_Explained = round(100 * s$dev.expl, 1),
    Tweedie_p = round(p_est, 3),
    AIC = round(AIC(mod), 1)
  )
})

# Calculate delta AIC
gam_table <- gam_table %>%
  mutate(
    Delta_AIC = round(AIC - min(AIC), 1)
  ) %>%
  arrange(AIC)

gam_table

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

#w = readRDS('first4CPUEmodels33.rds')
#l33 = w[[1]]
#l33a = w[[2]]
#l33b = w[[3]]
#l33c = w[[4]]


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



####annual index with observed effort by year
get_index_annual_effort <- function(model, ind=wts2, n_sim=1000, effort=1,name){
  
  mf  <- model.frame(model)
  yrs <- levels(mf$fYear)
  
  Xp <- predict(model, type="lpmatrix")
  b  <- coef(model)
  V  <- vcov(model)
  
  bsim <- MASS::mvrnorm(n=n_sim, mu=b, Sigma=V)
  
  res <- list()
  
  for(y in yrs){
    
    nd <- subset(ind, SYEAR == as.numeric(as.character(y)))
    
    nd$fYear  <- factor(y, levels=yrs)
    nd$SOURCE <- "marfis"
    nd$leffort <- log(effort)
    
    Xp <- predict(model,
                  newdata=nd,
                  type="lpmatrix")
    
    eta <- Xp %*% t(bsim)
    mu  <- exp(eta)
    
    idx <- colSums(mu * nd$prop)
    
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
l33p = get_index_annual_effort(l33,name='Base')
l33bp = get_index_annual_effort(l33b,name='BSD')
l33cp = get_index_annual_effort(l33c,name='BSDxY')
com = bind_rows(list(l33p,l33bp,l33cp))
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

cap =as.data.frame(do.call(rbind,cpue.ann))
cap$Model='biased_corrected'
cap = subset(cap,select=c(SYEAR, Model, unBCPUE, l95,u95))
names(cap) = names(com)[c(1,7,2,5,6)]

com = bind_rows(com,cap)


write.csv(com,'LFA33MultipleModelsCPUE.csv')
#marginal mean, internally consistent

ggplot(subset(com,SYEAR<2026 & model !='BS'), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_marginalCPUE_individualseasonwts.png')

#######same weighting across years

####annual index using the same weighing index by year
get_index_common_effort <- function(model, ind=wts, n_sim=1000, effort=1,name){
  
  mf  <- model.frame(model)
  yrs <- levels(mf$fYear)
  
  Xp <- predict(model, type="lpmatrix")
  b  <- coef(model)
  V  <- vcov(model)
  
  bsim <- MASS::mvrnorm(n=n_sim, mu=b, Sigma=V)
  
  res <- list()
  
  for(y in yrs){
    
    nd <- ind
    
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

l33pco = get_index_common_effort(l33,name='Base_common')
l33bpco = get_index_common_effort(l33b,name='BSD_common')
l33cpco = get_index_common_effort(l33c,name='BSDxY_common')
com_com = bind_rows(list(l33pco,l33bpco,l33cpco))

com_com = bind_rows(com_com,cap)


ggplot(subset(com_com,SYEAR<2026 & model !='BS'), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)

ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/LFA33_marginalCPUE_commonseasonwts.png')
