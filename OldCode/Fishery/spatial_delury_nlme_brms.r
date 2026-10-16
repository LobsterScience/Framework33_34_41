library(nlme)

require(bio.lobster)
require(bio.utilities)
require(dplyr)
require(devtools)
require(sf)
require(ggplot2)
require(mgcv)
require(ggeffects)
require(ggforce)
require(nlme)
require(sdmTMB)
require(purrr)
la()

fd = file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUE')
setwd(fd)

m = readRDS(file='data_33_34_41_sdmtmb.rds')
gto = m[[1]]
ns_coast = m[[2]]
ca = m[[3]]
ca$WOS = floor((ca$DOS-1)/7)+1
caT = aggregate(bcT~LFA+SYEAR+GRID_NO+WOS,data=ca,FUN=mean)
dat = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+WOS+GRID_NO,data=subset(ca,LFA %in% c(33)),FUN=sum)
dat$CPUE = dat$WEIGHT_KG/dat$NUM_OF_TRAPS

dl = split(dat,f=list(dat$LFA, dat$SYEAR, dat$GRID_NO))
o = list()
m=0
for(i in 1:length(dl)){
  b = dl[[i]]
  if(nrow(b)>3){
    m = m+1
    b$cumland = cumsum(b$WEIGHT_KG)
  o[[m]] = b
    }
}
dat = do.call(rbind,o)  
yr = unique(dat$SYEAR)
dat = merge(dat, caT)
#dat = subset(dat,WOS<11)
#delury

library(dplyr)
library(tidyr)
library(purrr)



################################################################################################################################################################################
#########all years  with temperature effects, allowing b0 to vary year to year and 

m_all <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(b0 +N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -7,
    b1 = 0,
    N0 = 100000
  )
)


m_all_2 <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(b0 + b1+N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -7,
    b1 = 0,
    N0 = 100000
  )
)


m_all_3 <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag( b1+N0 ~ 1),
    GRID_NO = pdDiag(N0 + b1 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -7,
    b1 = 0,
    N0 = 100000
  )
)



###all years variable q but not temperature
m_all_nt <- nlme(
  CPUE ~ q * (N0 - cumland),
  fixed = q + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(q + N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    q = 1e-05,
    N0 = 100000
  )
)

####all years same q
m_all_nt_1q <- nlme(
  CPUE ~ q * (N0 - cumland),
  fixed = q + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    q = 1e-05,
    N0 = 100000
  )
)

anova(m_all_2, m_all_3,m_all,m_all_nt,m_all_nt_1q)


################
library(tidyr)
library(dplyr)

grid_n0 <- coef(m_all_3) |>
  as.data.frame() |>
  tibble::rownames_to_column("grp") |>
  separate(grp, c("SYEAR", "GRID_NO"), sep = "/")

o = grid_n0 %>%
  group_by(SYEAR) %>%
  summarise(total_N0 = sum(N0)/1000)

o$SYEAR = as.numeric(o$SYEAR)

ln = lobster.db('seasonal.landings')
ln$SYEAR = as.numeric(substr(ln$SYEAR,6,9))
ggplot(o,aes(x=SYEAR,y=total_N0))+geom_line(colour='red',linewidth=1)+geom_segment(data=ln,aes(x=SYEAR,xend=SYEAR,y=0,yend=LFA33),linewidth=1)
#
######################################################################

###34

dat = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+WOS+GRID_NO,data=subset(ca,LFA %in% c(34)),FUN=sum)
dat$CPUE = dat$WEIGHT_KG/dat$NUM_OF_TRAPS

dl = split(dat,f=list(dat$LFA, dat$SYEAR, dat$GRID_NO))
o = list()
m=0
for(i in 1:length(dl)){
  b = dl[[i]]
  if(nrow(b)>3){
    m = m+1
    b$cumland = cumsum(b$WEIGHT_KG)
    o[[m]] = b
  }
}
dat = do.call(rbind,o)  
yr = unique(dat$SYEAR)
dat = merge(dat, caT)
#dat = subset(dat,WOS<11)
#delury

library(dplyr)
library(tidyr)
library(purrr)



################################################################################################################################################################################
#########all years  with temperature effects, allowing b0 to vary year to year and 

m_all <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(b0 +N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -7,
    b1 = 0,
    N0 = 100000
  )
)


m_all_2 <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(b0 + b1+N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -7,
    b1 = 0,
    N0 = 100000
  )
)


m_all_3v <- nlme(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  fixed = b0 + b1 + N0 ~ 1,
  random = list(
    SYEAR = pdDiag( N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    b0 = -13.5,
    b1 = 0.33,
    N0 = 276000
  ),
  control = nlmeControl(
    maxIter = 200,
    msMaxIter = 200,
    pnlsMaxIter = 30
  )
)

mm = update(
  m_all_3v,
  weights = varPower(form = ~ fitted(.))
)

###all years variable q but not temperature
m_all_nt <- nlme(
  CPUE ~ q * (N0 - cumland),
  fixed = q + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(q + N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    q = 1e-05,
    N0 = 100000
  )
)

####all years same q
m_all_nt_1q <- nlme(
  CPUE ~ q * (N0 - cumland),
  fixed = q + N0 ~ 1,
  random = list(
    SYEAR = pdDiag(N0 ~ 1),
    GRID_NO = pdDiag(N0 ~ 1)
  ),
  data = dat,
  start = c(
    q = 1e-05,
    N0 = 100000
  )
)

anova(m_all_2, m_all_3,m_all,m_all_nt,m_all_nt_1q)


################
library(tidyr)
library(dplyr)

grid_n0 <- coef(m_all_3) |>
  as.data.frame() |>
  tibble::rownames_to_column("grp") |>
  separate(grp, c("SYEAR", "GRID_NO"), sep = "/")

o = grid_n0 %>%
  group_by(SYEAR) %>%
  summarise(total_N0 = sum(N0)/1000)

o$SYEAR = as.numeric(o$SYEAR)

ln = lobster.db('seasonal.landings')
ln$SYEAR = as.numeric(substr(ln$SYEAR,6,9))
ggplot(o,aes(x=SYEAR,y=total_N0))+geom_line(colour='red',linewidth=1)+geom_segment(data=ln,aes(x=SYEAR,xend=SYEAR,y=0,yend=LFA34),linewidth=1)

plot(m_all_2)
qqnorm(resid(m_all_3))
qqline(resid(m_all_3))
plot(fitted(m_all_3), resid(m_all_3))
abline(h = 0)


m_all_3v <- update(
  m_all_3,
  weights = varExp()
)



###trying m_all_3v in brms as a gamma model 
bf0 <- bf(
  CPUE ~ mu,
  
  mu ~ b0 +
    b1 * bcT +
    log(exp(logN0) - cumland),
  
  b0 + b1 + logN0 ~ 1,
  
  nl = TRUE
)

m0 <- brm(
  bf0,
  data = dat,
  family = Gamma(link = "log"),
  prior = c(
    prior(normal(-13.5, 1), nlpar = "b0"),
    prior(normal(0.33, 0.1), nlpar = "b1"),
    prior(normal(log(500000), 0.2), nlpar = "logN0")
  ),
  chains = 1,
  cores = 1
)

##tester 
bf0 <- bf(
  CPUE ~ exp(b0 + b1 * bcT) *
  exp(logN0) *(1 - cumland / exp(logN0)),
  
  b0 + b1 + logN0 ~ 1,
  nl = TRUE
)

init_fun <- function(){
    list(
      b_b0_Intercept = -13.56,
      b_b1_Intercept = 0.33,
      b_logN0_Intercept = log(max(dat$cumland)*1.5)
    )
  }
  
m0 <- brm(
  bf0,
  data = dat,
  init = init_fun,
  family = Gamma(link = "log"),
  chains = 1
)

require(brms)
form <- bf(
  CPUE ~ exp(b0 + b1 * bcT) * (N0 - cumland),
  
  b0 ~ 1,
  
  b1 ~ 1, #+
    #(1 | SYEAR) +
    #(1 | SYEAR:GRID_NO),
  
  N0 ~ 1 +
    (1 | SYEAR) +
    (1 | SYEAR:GRID_NO),
  
  nl = TRUE
)

priors <- c(
  prior(normal(-13.5, 2), nlpar = "b0"),
  prior(normal(0.33, 0.2), nlpar = "b1"),
  prior(normal((2000000), 1),
        nlpar = "N0"),
  
 # prior(exponential(1), class = "sd",
#        nlpar = "b1"),
  prior(exponential(1), class = "sd",
        nlpar = "N0")
)

m_brms <- brm(
  formula = form,
  data = dat,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  cores = 4,
  iter = 4000,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

library(brms)
library(dplyr)

dat2 <- dat %>%
  group_by(SYEAR, GRID_NO) %>%
  mutate(max_cumland = max(cumland)) %>%
  ungroup()

bf_gamma <- bf(
  CPUE ~ b0 +
    b1 * bcT +
    log(max_cumland + exp(theta) - cumland),
  
  b0 ~ 1,
  b1 ~ 1,
  
  theta ~ 1 +
    (1 | SYEAR) +
    (1 | SYEAR:GRID_NO),
  
  nl = TRUE
)
init_fun <- function() {
  list(
    b_b0_Intercept = -14.36,
    b_b1_Intercept = 0.21,
    b_theta_Intercept = 12.86
  )
}
m_gamma <- brm(
  bf_gamma,
  data = dat2,
  #family = Gamma(link = "log"),
  family=lognormal(),
  prior = c(
    prior(normal(-14.7, 0.5), nlpar = "b0"),
    prior(normal(0.34, 0.1), nlpar = "b1"),
    prior(normal(log(50000), 1), nlpar = "theta"),
    
   prior(exponential(1),
        class = "sd",
        nlpar = "theta",
          group = "SYEAR"),
    
    prior(exponential(1),
          class = "sd",
          nlpar = "theta",
          group = "SYEAR:GRID_NO")
#        prior(gamma(2, 0.1), class = "shape")
  ),
  init = init_fun,
  chains = 3,
  cores = 3
)





v = as.data.frame(exp(coef(m_gamma)$`SYEAR:GRID_NO`))
v$yr = unlist(substr(rownames(v),1,4))
aggregate(Estimate.theta_Intercept~yr,data=v,FUN=sum)
####extract
library(dplyr)
library(tibble)

library(brms)
library(dplyr)
library(tidyr)
library(posterior)

# Extract posterior draws for fixed and random effects
draws <- as_draws_df(m_brms)

# Fixed effect
beta0 <- draws$b_logN0_Intercept

# Year random effects
year_cols <- grep("^r_SYEAR\\[.*,Intercept\\]$", names(draws), value = TRUE)

# Grid-within-year random effects
grid_cols <- grep("^r_SYEAR:GRID_NO\\[.*,Intercept\\]$", names(draws),
                  value = TRUE)

# Build posterior N0 estimates for each year-grid combination
N0_post <- lapply(grid_cols, function(gcol) {
  
  grp <- sub("^r_SYEAR:GRID_NO\\[(.*),Intercept\\]$", "\\1", gcol)
  
  yr <- strsplit(grp, ":")[[1]][1]
  
  ycol <- paste0("r_SYEAR[", yr, ",Intercept]")
  
  if(!ycol %in% names(draws))
    return(NULL)
  
  tibble(
    group = grp,
    draw = seq_len(nrow(draws)),
    N0 = exp(
      beta0 +
        draws[[ycol]] +
        draws[[gcol]]
    )
  )
}) |>
  bind_rows()

# Posterior summaries
N0_summary <- N0_post |>
  group_by(group) |>
  summarise(
    mean_N0  = mean(N0),
    median_N0 = median(N0),
    lwr95 = quantile(N0, 0.025),
    upr95 = quantile(N0, 0.975),
    .groups = "drop"
  ) |>
  separate(group,
           into = c("SYEAR", "GRID_NO"),
           sep = ":")

N0_summary

year_totals <- N0_post |>
  separate(group, c("SYEAR", "GRID_NO"), sep = ":") |>
  group_by(draw, SYEAR) |>
  summarise(total_N0 = sum(N0), .groups = "drop")

year_totals |>
  group_by(SYEAR) |>
  summarise(
    mean = mean(total_N0),
    lwr95 = quantile(total_N0, 0.025),
    upr95 = quantile(total_N0, 0.975)
  )


