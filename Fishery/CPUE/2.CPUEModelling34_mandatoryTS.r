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

ca$fYear = as.factor(ca$SYEAR)
ca$SOURCE = as.factor(ca$SOURCE)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca = subset(ca,!is.na(NUM_OF_TRAPS) & !is.na(WEIGHT_KG) & NUM_OF_TRAPS<=1200 & WEIGHT_KG>0)
ca$CPUE = ca$WEIGHT_KG/ca$NUM_OF_TRAPS
ca = subset(l34,SYEAR>2000)
l34 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca),family = Gamma(link='log'),method='REML')
l34b = bam(WEIGHT_KG~s(DOS)+fYear+offset(leffort),data=ca,family = Gamma(link='log'),method='fREML',discrete=TRUE)
l34c = bam(WEIGHT_KG~fYear+s(DOS,fYear,bs='fs')+offset(leffort),data=ca,family = Gamma(link='log'),method='fREML',discrete=TRUE,nthreads = (parallel::detectCores()-1))

#l34d = gam(WEIGHT_KG~s(DOS)+s(DOS,by=fYear)+fYear+s(bcT)+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')

saveRDS(list(l34,l34b,l34c),'mandatoryReportingCPUEmodels34.rds')


library(dplyr)
library(purrr)

mods <- list(Base = l34,  BSD = l34b,BSDxY=l34c)

  gam_table <- imap_dfr(mods, function(mod, name){
    
    s <- summary(mod)
    
    tibble(
      Model = name,
      Terms = paste(attr(terms(mod), "term.labels"),
                    collapse = " + "),
      EDF = sum(s$edf),
      Dev_Explained = round(100 * s$dev.expl, 1),
      AIC = round(AIC(mod), 1)
    )
  })
  

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
        wts <- ind[,c('SYEAR','DOS','prop')]
        
        
        
####annual index       
        get_index <- function(model, ind=wts, n_sim=1000, effort=1,name){
          
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
l34p = get_index(l34,name='Base')
l34bp = get_index(l34b,name='BSD')
l34cp = get_index(l34c,name='BSDxY')
com = bind_rows(list(l34p,l34bp,l34cp))
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


write.csv(com,'LFA34MultipleModelsCPUE.csv')
#marginal mean, internally consistent

ggplot(subset(com,SYEAR<2026 & model !='BS'), aes(x = SYEAR, y = mean ,colour=model,fill=model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)





#now just the later year models as the early logs cant be used for full standardizatino