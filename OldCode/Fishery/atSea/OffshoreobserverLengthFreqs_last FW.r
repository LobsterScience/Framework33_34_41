#at sea observer Length frequs
options(stringsAsFactors=F)
#combining mature females abundance at length and fecundity at length
require(bio.lobster)
require(bio.utilities)

la()


		LFA41 = readRDS(file.path( git.repo,'bio.lobster.data','mapping_data',"LFA41Offareas.rds"))
		
#ISDB tables 
			lobster.db(DS = 'atSea') #object called obs.samp
			obs.samp = subset(atSea,LFA==41& DESCRIPTION %ni% c('Data not collected by trap') & !is.na(LONGITUDE) & LONGITUDE< -64)
			obs.samp$EID = 1:nrow(obs.samp)
			obs.samp$yr = year(obs.samp$STARTDATE)
			obs.samp$qt = quarter(obs.samp$STARTDATE)
			obs.samp$mn = month(obs.samp$STARTDATE)
			obs.samp$X = obs.samp$LONGITUDE
			obs.samp$Y = obs.samp$LATITUDE
			
			o_sf = st_as_sf(obs.samp,coords=c('LONGITUDE','LATITUDE'),crs=4326)
			o_sf = st_join(o_sf,LFA41)
			o = subset(o_sf,!is.na(OFFAREA) & SPECIESCODE %in% c(2550) | is.na(SPECIESCODE))
      o$UID = paste(o$TRIPNO,o$TRAPNO,o$STARTDATE,o$STRINGNO,o$X,sep="-")
			
      ox = aggregate(STARTDATE~UID,data=o,FUN=min)
      
			
						
			ggplot()+geom_sf(data=o)+geom_sf(data=LFA41,fill=NA)
				

			a = aggregate(CARLENGTH~OFFAREA+yr+mn,data=o,FUN=length)
######################################################################################
#matching observed trips to logs
#how good is teh sampling, join o with sfdat from CPUE/41/tripSummaries.r
			sfdat$EID = 1:nrow(sfdat)
			log_dates <- sort(unique(sfdat$DATE_FISHED))
			obs_dates <- sort(unique(as.Date(o$STARTDATE)))
			
			nearest_obs_date <- sapply(obs_dates, function(d) {
			  log_dates[which.min(abs(as.numeric(log_dates - d)))]
			})

			
			date_match <- data.frame(
			  ob_date = obs_dates,
			  log_date = 			as.Date(nearest_obs_date, origin = "1970-01-01"),
			  day_diff = as.numeric(abs(obs_dates - 			as.Date(nearest_obs_date, origin = "1970-01-01")))
			)
			
			date_match = subset(date_match,day_diff<10)
			
			
			matches <- list()
			
			for(i in seq_len(nrow(date_match))) {
			  
			  log_dat <- date_match$log_date[i]
			  obs_dat <- date_match$ob_date[i]
			  
			  log <- sfdat[sfdat$DATE_FISHED == log_dat, ]
			  obs <- o[as.Date(o$STARTDATE) == obs_dat, ]
			  
			  if(nrow(obs) == 0 || nrow(log) == 0)
			    next
			  
			  obs_utm <- st_transform(obs, 32620)
			  log_utm <- st_transform(log, 32620)
			  
			  nn <- st_nearest_feature(obs_utm, log_utm)
			  
			  dist_m <- as.numeric(
			    st_distance(
			      obs_utm,
			      log_utm[nn, ],
			      by_element = TRUE
			    )
			  )
			  
			  matches[[i]] <- data.frame(
			    UID          = obs$UID,
			    obs_date     = obs_dat,
			    log_date     = log_dat,
			    log_EID      = log$EID[nn],
			    distance_m   = dist_m
			  )
			}
			
			matches <- do.call(rbind, matches)
			matches <- matches[!duplicated(matches),]
			ml = split(matches,f=matches$log_EID)
### match cpues from obs and trips
			df = data.frame(DateFished=NA,ObsDate=NA,Area=NA,SumFishTraps=NA,SumCatch=NA, SumObsTraps=NA,SumObsLegalCatch=NA,logsEID = NA)
			for(i in 1:length(ml)){
			    x = ml[[i]]
			    v = subset(sfdat, EID == unique(x$log_EID))
			    w = subset(o,UID %in% x$UID[which(x$distance_m<30000)])
			    if(nrow(w)==0 | nrow(v)==0) next
			    w$Legal = ifelse(w$CARLENGTH>82 & w$SEX <3,1,0)
			    w$Legal = ifelse(w$CARLENGTH==82 & w$SEX <3,sample(size=1,c(0,1),prob=c(.5,.5)),w$Legal) #randomly assign legal or sublegal
			    j = which(is.na(w$CALWT))
			    w$CALWT[j] = 0
			    w$Legal[j] = 1
			    wl = aggregate(CALWT~UID,data=subset(w,Legal==1),FUN=sum)[,2]
			    df[i,'SumObsLegalCatch'] = sum(wl)/1000
			    df[i,'SumObsTraps'] = length(wl)
			    df[i,'Area'] = unique(w$OFFAREA)[1]
			    df[i,'DateFished'] = as.Date(unique(v$DATE_FISHED))
			    df[i,'ObsDate'] = as.Date(unique(w$STARTDATE))[1]
			    df[i,'SumCatch'] = sum(v$ADJCATCH)
			    df[i,'SumFishTraps'] = sum(v$NUM_OF_TRAPS)
			    df[i,'logsEID'] = unique(x$log_EID)
			     }
			
			df = subset(df, !is.na(Area))
			df$obscpue = df$SumObsLegalCatch/df$SumObsTraps
			df$logcpue = df$SumCatch/df$SumFishTraps
			df$DateFished = as.Date(df$DateFished)
			df$ObsDate = as.Date(df$ObsDate)
			
			df$YR = lubridate::year(df$DateFished)
			df = subset(df,obscpue<30)
			df$YR10 <- paste0(floor(df$YR / 10) * 10, "s")
			ggplot(df,aes(logcpue,obscpue))+geom_point()+geom_abline(slope=1,col='red',linewidth=1.1)+
			  geom_smooth(method = MASS::rlm)+facet_wrap(~YR10)
			
			library(dplyr)
			library(tidyr)
			library(ggplot2)
			
			cdfdat <- df %>%
			  select(YR10, obscpue, logcpue) %>%
			  pivot_longer(
			    cols = c(obscpue, logcpue),
			    names_to = "metric",
			    values_to = "value"
			  )
			
			ggplot(cdfdat, aes(x = value, colour = metric)) +
			  stat_ecdf(linewidth = 1) +
			  facet_wrap(~YR10, scales = "free_x") +
			  labs(
			    x = "CPUE",
			    y = "Cumulative probability",
			    colour = NULL
			  ) +
			  theme_bw()

	##compare to all fishing
			sfdat$sampled = ifelse(sfdat$EID %in% df$logsEID,1,0)
			sfdat = st_join(sfdat,LFA41)
			sdat = as_tibble(sfdat)
			ggplot(sdat, aes(fishingYear, fill = sampled)) +
			  geom_bar(position = "fill") +
			  ylab("Proportion")	
			
			sdat
			
			
################################################################################################
			
#seasonal and annual changes by area
	l = unique(O$Area)
	pdf('~/tmp/Olenfeq.pdf')

	for(j in l) {
			P = subset(O,Area==j)
			P = P[order(P$yr),]
			y = unique(P$yr) 
	
				for(i in y) {
					J = subset(P,yr==i)
					s = unique(J$qt)

					plot(1,1,xlim=c(min(O$carlength),max(O$carlength)),type='n',main=paste(j,i,sep="-"),ylim=c(0,1))
					legend('topleft',c("Winter",'Spring','Summer','Fall'),lty=c(1,1,1,1),col=c(1,2,3,4),bty='n')
						for(m in s) {
							K = subset(J,qt==m)
							B = hist(K$carlength,plot=F,breaks=seq(20,200,1))
							lines(B$mids,cumsum(B$density)/sum(B$density),col=m)
						}


				}
		}
dev.off()

O  = subset(O,yr < 2016)

#not all seasons and areas have adequate data the combinations to use are as determined by ACook Sept 2016

              SWBs = subset(O,Area=='SW.Browns' & mn %in% c(6,7))
              SWBf = subset(O,Area=='SW.Browns' & mn %in% c(10,11,12,1))

              SEB = subset(O,Area=='SE.Browns' & mn %in% c(4,5,6,7))
              GBs = subset(O,Area=='Georges.Basin' & mn %in% c(1,2,3,4))
              GBsu = subset(O,Area=='Georges.Basin' & mn %in% c(5,6,7))

              GBa = subset(O,Area=='Georges.Bank' & mn %in% c(5,6,7))

              #Crowell Basin dropped as no fishing observed trips since 2005
	            ll = list(SWBs,SWBf,SEB,GBs,GBsu,GBa)
	      		ln=c('SW.Browns.Summer','SW.Browns.Autumn','SE.Browns.Spring','Georges.Basin.Spring','Georges.Basin.Summer','Georges.Bank.Summer')

for( i in 1:length(ll)) {
	r = ll[[i]]

			lm = aggregate(carlength~yr,data=r, FUN=function(x) quantile(x,c(0.5,0.25,0.75,0.95)))
				lm = data.frame(lm[,1],lm$carlength[,1],lm$carlength[,2],lm$carlength[,3],lm$carlength[,4])

			
			names(lm) = c('yr','medL','medLlower','medLupper','upper95')
			
			r$ff = round(r$carlength/3)*3
			
			aa = aggregate(qt~ff+yr,data = r,FUN=length)
			sH=c()
			dk = unique(aa$yr)
			for (j in dk) {
				gh = subset(aa,yr==j)
				gh$p = gh$qt / sum(gh$qt)
				sH = c(sH,-1*(sum(gh$p * log(gh$p)) / log(nrow(gh))))
			}
			sH = data.frame(yr=  dk, ShannonEquitability = sH)
			lm = merge(lm,sH)
			af = aggregate(carlength~yr,data=r, FUN=length)
			lm = merge(lm,af,by='yr')
			names(lm)[7] = 'ObsLobs'
			names(af) = c('x','y')
			write.csv(lm,file=file.path(project.datadirectory('bio.lobster'),'analysis','indicators',paste(ln[i],'obslength.csv',sep="_")))


			#grouped year length freq plots

					YG = 5 # year grouping
					y = unique(aa$yr)
					yL = y[length(y)] #last year
					yLL = length(y)-1
					yLm = yLL %% YG
					yLr = yLL %/% YG
					yw = y[which(y %in% y[1:yLm])] #add the early years to the first histogram and keep the rest at 5 years
					
					yLw = c(rep(1,yLm),rep(1:yLr,each = YG),yLr+1)
					grps = data.frame(yr = y,ry = yLw)
					aa = merge(aa,grps,by='yr',all.x=T)
					yll = max(aggregate(qt~ff+ry,data=aa,FUN=mean)$qt)
					yll = 1
					h = split(aa,f=aa$ry)
					
					for(j in 1:length(h)) {
							g = h[[j]]
							
							y = unique(g$yr)
							u = aggregate(qt~ff,data=g,FUN=mean)              
							u$qt = u$qt / max(u$qt)
							fn = paste(ln[i],min(y),max(y),'pdf',sep=".")
				            nn = sum(g$qt)
				            pdf(file.path(project.figuredirectory('bio.lobster'),fn))
							plot(u$ff,u$qt,lwd=3,xlab='Carapace Length',ylab = 'Scaled Number Measured',type='h',ylim=c(0,yll),xlim=c(50,200))
							abline(v=82.5,lty=2,col='red',lwd=3)
							legend('topleft',bty='n',pch="", legend=c(paste(min(y),max(y),sep="-"),paste('N=',nn,sep=" ")),cex=1.5)
							dev.off()
							print(fn)
							}
						


			fn = paste('max95',ln[i],'pdf',sep=".")
            nn = sum(r$ObsLobs)
            pdf(file.path(project.figuredirectory('bio.lobster'),fn))
			plot(lm$yr,lm$upper95,lwd=1,xlab='Year',ylab = 'Maximum Length (mm)',type='b',ylim=c(115,195),pch=16)
			lines(lm$yr,runmed(lm$upper95,k=3,endrule='median'),col='salmon',lwd=2)
			dev.off()
			
			fn = paste('shannon',ln[i],'pdf',sep=".")
            nn = sum(g$ObsLobs)
            oo = lm[,c('yr','ShannonEquitability')]
            ii = which(is.na(oo$ShannonEquitability))
            oo$ShannonEquitability[ii] <- oo$ShannonEquitability[ii-1]

            pdf(file.path(project.figuredirectory('bio.lobster'),fn))
			plot(oo$yr,oo$ShannonEquitability,lwd=1,xlab='Year',ylab = 'Shannon Equitability',type='b',pch=16,ylim=c(0.60,1))
			lines(oo$yr,runmed(oo$ShannonEquitability,k=3,endrule='median'),col='salmon',lwd=2)
			dev.off()




				p=list()
			                  p$add.reference.lines = F
                              p$time.series.start.year = min(r$yr)
                              p$time.series.end.year = max(r$yr)
                              p$metric = 'medianL' #weights
                              p$measure = 'stratified.mean' #'stratified.total'
                              p$figure.title = ""
                              p$reference.measure = 'median' # mean, geomean
                              p$file.name = paste('medianL',ln[i],'png',sep='.')
    		                  print(p$file.name)

                        p$y.maximum = NULL # NULL # if ymax is too high for one year
                        p$show.truncated.numbers = F #if using ymax and want to show the numbers that are cut off as values on figure
                        		p$ylim = c(60,155)
                                p$legend = FALSE
                                p$running.median = T
                                p$running.length = 3
                                p$running.mean = F #can only have rmedian or rmean
                               p$error.polygon=T
                              p$error.bars=F
                              

                              p$ylim2 = c(0,5000)
                             
                       figure.stratified.analysis(x=lm,out.dir = 'bio.lobster', x2 = af, p=p,sampleSizes=T,save=T)
		}



