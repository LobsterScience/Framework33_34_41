#functional data analysis

require(bio.lobster)
require(ggplot2)
require(bio.utilities)
require(devtools)
require(ggpubr)
require(MARSS)
library(tidyr)
library(dplyr)

load_all('~/git/bio.utilities')

theme_set(theme_test(base_size = 14))
fig_dir = file.path('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/')
#cpue index
co = readRDS(file.path(project.datadirectory('Framework_LFA33_34_41'),'CPUE','unBIASED_CPUE.rds'))
cpu = co[[2]]

cu = subset(cpu,lfa==34)

wdf <- cu %>%
  dplyr::select(t,unBCPUE, yr) %>%
  pivot_wider(
    names_from = t,
    values_from = unBCPUE,
    names_prefix = "t_"
  )
wdf[,1] <- NULL
wdf = as.matrix(wdf)

cR = (log(wdf))

cRR = zscore(as.matrix(cR))


N.ts = nrow(cRR)
nT = ncol(cRR)



cntl.list = list(minit=200, maxit=5000, allow.degen=FALSE)
find.best=F
if(find.best){
  # set up forms of R matrices
  levels.R = c("diagonal and equal","diagonal and unequal","equalvarcov","unconstrained")
  levels.Q = c("diagonal and equal","diagonal and unequal")
  model.data = data.frame()
  for(R in levels.R) {
    for(Q in levels.Q){
      for(m in 1:(N.ts-1)) {
        dfa.model = list(A="zero", R=R,Q=Q, m=m)
        ff = MARSS(cRR, model=dfa.model, control=cntl.list,
                   form="dfa", z.score=TRUE)
        
        model.data = rbind(model.data,
                           data.frame(R=R,Q=Q,m=m,logLik=ff$logLik,K=ff$num.params,AICc=ff$AICc,stringsAsFactors=FALSE))
        assign(paste("ff", m, R, Q, sep="."), ff)
      } # end m loop
    } # end Q loop
  } # end R loop
  #best model is diag and unequal R and Q with three trends
  write.csv(model.data,'~/tmp/model.data.csv')
  
}
