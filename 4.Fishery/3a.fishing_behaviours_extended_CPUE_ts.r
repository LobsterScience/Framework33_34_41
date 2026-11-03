#vessel and fishing characteristics

require(bio.lobster)
require(ggplot2)
require(bio.utilities)
require(devtools)
require(ggpubr)
load_all('~/git/bio.utilities')
la()
theme_set(theme_test(base_size = 14))
fig_dir = file.path('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/')


fd = file.path(project.datadirectory('Framework_LFA33_34_41'))
setwd(fd)

m = readRDS(file='ExtendedTS_LFA34_logs.rds')
ca = m



x1 = aggregate(WEIGHT_KG~SYEAR+LFA+WOS,data=subset(ca,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=sum)
x1$P=x1$WEIGHT_KG


xxx = split(x1,f=list(x1$LFA,x1$SYEAR))
junk = list()

for(i in 1:length(xxx)){
  o = xxx[[i]]
  o$prp = cumsum(o$P) / sum(o$P)
  junk[[i]] = o  
}
x1a = do.call(rbind,junk)
x1a$Fishing_Season =x1a$SYEAR
max_week <- max(x1a$WOS)
ggplot(subset(x1a),aes(x=WOS,y=prp,group=Fishing_Season,colour=Fishing_Season))+
  scale_colour_viridis_c(option='inferno')+geom_line()+facet_wrap(~LFA)+xlab('Week of Season')+ylab('Proportion of Total Landings')+
  geom_abline(
    slope = 1/max_week,
    intercept = 0,
    colour = "white",
    linewidth = 2.5
  ) +
  geom_abline(
    slope = 1/max_week,
    intercept = 0,
    colour = "black",
    linewidth = 1.2,
    linetype = "longdash"
  )+
  theme_test(base_size = 14)

ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/FiguresFishingCharact/Fig.season.2.34_timing_of_landings.png')

ca$ID = paste(ca$UNIQUE_HARVESTER,ca$DATE_FISHED,sep="-")
x1 = aggregate(ID~SYEAR+UNIQUE_HARVESTER+LFA+WOS,data=subset(ca,WEIGHT_KG>0 & NUM_OF_TRAPS>0),FUN=function(x) length(unique(x)))
x1$P=x1$ID

x1a = aggregate(P~SYEAR+LFA+WOS,data=x1,FUN=sum)

xxx = split(x1a,f=list(x1a$LFA,x1a$SYEAR))
junk = list()

for(i in 1:length(xxx)){
  o = xxx[[i]]
  o$prp = cumsum(o$P) / sum(o$P)
  junk[[i]] = o  
}
x1a = do.call(rbind,junk)
x1a$Fishing_Season =x1a$SYEAR

max_week <- max(x1a$WOS)

ggplot(x1a,aes(x=WOS,y=prp,group=Fishing_Season,colour=Fishing_Season))+scale_colour_viridis_c(option='inferno')+geom_line()+facet_wrap(~LFA)+xlab('Week of Season')+ylab('Proportion of Total Trips')+  
  geom_abline(
    slope = 1/max_week,
    intercept = 0,
    colour = "white",
    linewidth = 2.5
  ) +
  geom_abline(
    slope = 1/max_week,
    intercept = 0,
    colour = "black",
    linewidth = 1.2,
    linetype = "longdash"
  )+
  theme_test(base_size = 14)
ggsave('C:/Users/cooka/OneDrive - DFO-MPO/LFA33_34_41_Framework/Documents/Figures/FiguresFishingCharact/Fig.season.2.34_timing_of_trips.png')
