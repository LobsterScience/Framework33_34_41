require(bio.lobster)
require(devtools)
la()

lS<-lobster.db('process.logs')
lS = subset(lS,SYEAR<2025)
H = lobster.db('historic.cpue')
H$CPUE = H$LBSPTRAP/2.2046

V =  lobster.db('process.vlog')
V$SYEAR = as.numeric(year(V$FDATE))
V$SYEAR = year(V$FDATE)
V$MONTH = month(V$FDATE)
ii = which(V$MONTH>8)
V$SYEAR[ii] = V$SYEAR[ii]+1 




##LFA 34

H34 = aggregate(CPUE~LFA+SYEAR+SDATE,data=subset(H, LFA==34 & SYEAR<1960),FUN=mean)
H34 = H34[order(H34$SDATE),]
names(H34) = c('LFA','SYEAR','SDATE','CPUE')
aH34 = aggregate(cbind(CPUE,SDATE)~LFA+SYEAR,data=subset(H, LFA==34 & SYEAR<1960),FUN=mean)

V34 = aggregate(cbind(W_KG,N_TRP)~LFA+SYEAR+FDATE, data=subset(V,LFA==34), FUN=sum)
V34$CPUE = V34$W_KG / V34$N_TRP
V34$W_KG = V34$N_TRP = NULL
names(V34) = c('LFA','SYEAR','SDATE','CPUE')

#kludge

V34$SYEAR[which(abs(V34$SYEAR-year(V34$SDATE))< -1)] = 1991



aV34 = aggregate(cbind(W_KG,N_TRP)~LFA+SYEAR, data=subset(V,LFA==34), FUN=sum)
aV34$CPUE = aV34$W_KG / aV34$N_TRP
aV34$W_KG = aV34$N_TRP = NULL
names(aV34) = c('LFA','SYEAR','CPUE')
aaV34 = aggregate(FDATE~LFA+SYEAR, data=subset(V,LFA==34), FUN=mean)
names(aaV34)[3] = 'SDATE'
aV34 = merge(aV34,aaV34)

L34 = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+DATE_FISHED, data=subset(lS,LFA==34), FUN=sum)
L34$CPUE = L34$WEIGHT_KG / L34$NUM_OF_TRAPS
L34$WEIGHT_KG = L34$NUM_OF_TRAPS = NULL
names(L34) = c('LFA','SYEAR','SDATE','CPUE')

aL34 = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR, data=subset(lS,LFA==34), FUN=sum)
aL34$CPUE = aL34$WEIGHT_KG / aL34$NUM_OF_TRAPS
aL34$WEIGHT_KG = aL34$NUM_OF_TRAPS = NULL
names(aL34) = c('LFA','SYEAR','CPUE')
aaL34 = aggregate(DATE_FISHED~LFA+SYEAR, data=subset(lS,LFA==34), FUN=mean)
names(aaL34)[3] = 'SDATE'

aL34 = merge(aL34,aaL34)

b34 = as.data.frame(rbind(rbind(H34,V34),L34))
bb34 = as.data.frame(rbind(rbind(aH34,aV34),aL34))

yy = unique(b34$SYEAR)

bb34$SDATE = as.POSIXct(paste(bb34$SYEAR,'01','01',sep="-"))

with(b34,plot(CPUE~SDATE,type='n',ylim=c(0,max(CPUE,na.rm=T)),col=rgb(0,0,0,0.5)))
for(i in yy){
  with(subset(b34,SYEAR==i),lines(CPUE~SDATE, col=rgb(0,0,0,0.5)))
}
with(subset(bb34,SYEAR<1970),points(CPUE~SDATE, col='red',pch=16,type='b'))
with(subset(bb34,SYEAR>1970),points(CPUE~SDATE, col='red',pch=16,type='b'))


#################
##LFA 35 not enough data to make it worth it

##LFA 38

H34 = aggregate(CPUE~LFA+SYEAR+SDATE,data=subset(H, LFA==38 & SYEAR<1960),FUN=mean)
H34 = H34[order(H34$SDATE),]
names(H34) = c('LFA','SYEAR','SDATE','CPUE')
aH34 = aggregate(cbind(CPUE,SDATE)~LFA+SYEAR,data=subset(H, LFA==38 & SYEAR<1960),FUN=mean)

#kludge



L34 = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR+DATE_FISHED, data=subset(lS,LFA==38), FUN=sum)
L34$CPUE = L34$WEIGHT_KG / L34$NUM_OF_TRAPS
L34$WEIGHT_KG = L34$NUM_OF_TRAPS = NULL
names(L34) = c('LFA','SYEAR','SDATE','CPUE')

aL34 = aggregate(cbind(WEIGHT_KG,NUM_OF_TRAPS)~LFA+SYEAR, data=subset(lS,LFA==38), FUN=sum)
aL34$CPUE = aL34$WEIGHT_KG / aL34$NUM_OF_TRAPS
aL34$WEIGHT_KG = aL34$NUM_OF_TRAPS = NULL
names(aL34) = c('LFA','SYEAR','CPUE')
aaL34 = aggregate(DATE_FISHED~LFA+SYEAR, data=subset(lS,LFA==38), FUN=mean)
names(aaL34)[3] = 'SDATE'

aL34 = merge(aL34,aaL34)

a34 = as.data.frame(rbind(H34,L34))
aa34 = as.data.frame(rbind(aH34,aL34))

yy = unique(a34$SYEAR)

aa34$SDATE = as.POSIXct(paste(aa34$SYEAR,'01','01',sep="-"))

with(a34,plot(CPUE~SDATE,type='n',ylim=c(0,max(CPUE,na.rm=T)),col=rgb(0,0,0,0.5)))
for(i in yy){
  with(subset(a34,SYEAR==i),lines(CPUE~SDATE, col=rgb(0,0,0,0.5)))
}
with(subset(aa34,SYEAR<1970),points(CPUE~SDATE, col='red',pch=16,type='b'))
with(subset(aa34,SYEAR>1970),points(CPUE~SDATE, col='red',pch=16,type='b'))


