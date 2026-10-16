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

fd = file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUE')
setwd(fd)

m = readRDS(file='data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
ca = m[[3]]

#aT = lobster.db('process.logs')
ca = subset(ca,SYEAR>2002 & SYEAR<2026 & LFA %in% c(33,34))
ca$mn = lubridate::month(ca$DATE_FISHED)
ca$SYEAR = ifelse(ca$mn %in% c(11,12), ca$yr+1, ca$yr)
ca$fYear = as.factor(ca$SYEAR)
ca$leffort = log(ca$NUM_OF_TRAPS)
ca$DOS = round(ca$DOS)

require(parallel)

l34 = gam(WEIGHT_KG~fYear+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='REML')
l34a = bam(WEIGHT_KG~s(DOS)+fYear+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='fREML',discrete = T,samfrac=.1,nthreads=detectCores()-1)
l34b = bam(WEIGHT_KG~s(DOS)+fYear+s(bcT)+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='fREML',discrete = T,samfrac=.1,nthreads=detectCores()-1)

###most of the interannual signal is getting 'eaten' up by the interaction term
#l34b = bam(WEIGHT_KG~s(DOS,k=10)+s(DOS,fYear,bs='sz',k=10)+fYear+offset(leffort),data=subset(ca,LFA==34),select=T, family = Gamma(link='log'),method='fREML',discrete = T,samfrac=.1,nthreads=detectCores()-1)
#l34c = bam(WEIGHT_KG~s(DOS,k=10)+s(DOS,fYear,bs='sz',k=10)+fYear+s(bcT)+offset(leffort),data=subset(ca,LFA==34),family = Gamma(link='log'),method='fREML',discrete = T,samfrac=.1,nthreads=detectCores()-1)


# simple predictions
###1
        yrs <- data.frame(
          fYear = levels(subset(ca, LFA == 34)$fYear),
          leffort = 0   # offset = log(1)
        )
        
        # Predict annual means
        p <- predict(
          l34,
          newdata = yrs,
          type = "response"
        )
        
        annual_means1 <- data.frame(
          Year = yrs$fYear,
          Mean = p)

#predict 2
          library(dplyr)
          
          d <- subset(ca, LFA == 34)
          
          nd <- expand.grid(
            fYear = levels(d$fYear),
            DOS = sort(unique(d$DOS))
          )
          
          nd$leffort <- 0   # offset = log(1)
          
          # Predict on link scale
          p <- predict(l34a,
                       newdata = nd,
                       type = "link",
                       se.fit = TRUE)
          
          nd$fit <- p$fit
          
          # Average over DOS distribution
          annual_means2 <- nd %>%
            group_by(fYear) %>%
            summarize(
              Mean = mean(exp(fit)),
              .groups = "drop"
            )


#predict 3
          library(dplyr)
          
          d <- subset(ca, LFA == 34)
          
          # Use all observed DOS-bcT combinations
          nd <- d[, c("DOS", "bcT")]
          nd$leffort <- 0
          
          # Replicate for each year
          nd <- merge(
            nd,
            data.frame(fYear = levels(d$fYear))
          )
          
          # Predictions
          nd$pred <- predict(
            l34b,
            newdata = nd,
            type = "response"
          )
          
          annual_me <- nd %>%
            group_by(fYear) %>%
            summarise(
              Mean = mean(pred),
              .groups = "drop"
            )
          
      

          
saveRDS(list(l34,l34a,l34b),'first3CPUEmodels34.rds')
w = readRDS('first3CPUEmodels34.rds')
l34 = w[[1]]
l34a = w[[2]]
l34b = w[[3]]



#################trip weighted marginal means
aT=ca
ind = aggregate(SD_LOG_ID~DOS+SYEAR,data=subset(aT,LFA ==34),FUN=function(x) length(unique(x)))
ind1 = aggregate(SD_LOG_ID~SYEAR,data=ind,FUN=sum)
names(ind1)[2] = 'SumTrips'
ind = merge(ind,ind1)
ind$prop = ind$SD_LOG_ID/ind$SumTrips
ind$fYear=as.factor(ind$SYEAR)


##l34c
###updating the prediction grids
tmps <- aggregate(
  bcT ~ round(DOS),
  data = subset(ca, LFA == 34),
  FUN = mean
)
names(tmps)[1] = 'DOS'

nd <- merge(
  tmps,
  data.frame(
    fYear = levels(ca$fYear),
    leffort = 0
  )
)

annual_index <- function(model, nd, ind, nsim = 1000) {
  
  # Design matrix and parameter uncertainty
  Xp <- predict(model, nd, type = "lpmatrix")
  Vp <- vcov(model)
  b  <- coef(model)
  
  # Simulate coefficients
  bsim <- MASS::mvrnorm(
    nsim,
    mu = b,
    Sigma = Vp
  )
  
  # Predictions for all draws
  eta <- Xp %*% t(bsim)
  mu  <- exp(eta)
  
  # Add weighting information
  nd <- merge(nd, ind)
  
  yrs <- unique(nd$fYear)
  
  res <- do.call(
    rbind,
    lapply(yrs, function(y) {
      
      idx <- which(nd$fYear == y)
      
      w <- nd$SD_LOG_ID[idx]
      w <- w / sum(w)
      
      annual_draws <- colSums(mu[idx, ] * w)
      
      data.frame(
        Year = y,
        Mean = mean(annual_draws),
        SE   = sd(annual_draws),
        LCL  = quantile(annual_draws, 0.025),
        UCL  = quantile(annual_draws, 0.975)
      )
    })
  )
  
  rownames(res) <- NULL
  res
}

mods <- list(
  Base = l34,
  BD = l34a,
  BDT = l34b
)


results <- do.call(
  rbind,
  lapply(names(mods), function(nm) {
    out <- annual_index(mods[[nm]], nd, ind)
    out$Model <- nm
    out
  })
)


write.csv(resutls,'LFA34MultipleModelsCPUE.csv')
#marginal mean, internally consistent
results$Year = as.numeric(results$Year)
ggplot(results, aes(x = Year, y = Mean ,colour=Model,fill=Model, group=Model)) +
  geom_line(linewidth=1) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width=0.2)+
  #scale_color_discrete(begin = 0, end = 1, option = 'viridis')+
  xlab('Fishing Season')+
  ylab('Weighted Marginal Mean CPUE')+
  theme_test(base_size = 14)

