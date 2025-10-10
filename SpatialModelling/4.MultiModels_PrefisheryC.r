#Create Mesh Survey Data only

require(tidyr)
require(sdmTMB)
require(sdmTMBextra)
require(splines)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(PBSmapping)
require(SpatialHub)
require(sf)
la()
fd=file.path(project.datadirectory('Framework_LFA33_34_41'),'outputs','SURVEYS')
setwd(fd)
crs_utm20 <- 32620
sf_use_s2(FALSE) #needed for cropping


a = readRDS('SurveyOnlyPreFishery_Data.rds')
aT = a[[1]]
bB = a[[2]]
cC = a[[3]]


sf_use_s2(FALSE) #needed for cropping
aT$X1000 <- st_coordinates(aT)[,1]
aT$Y1000 <- st_coordinates(aT)[,2]

survey = as_tibble(subset(aT,z>0))
survey$of = log(survey$OFFSET)

spde <- make_mesh(survey, xy_cols = c("X1000", "Y1000"),
                  cutoff = 21) #12 is for final
#n vertices
spde$mesh$n

#Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, cC, range_fraction = .1,
  proj_scaling = 1, plot = TRUE
)


survey$Legal = survey$PreF_C
k <- 3 # 4 basis splines
k3 <- Hmisc::wtd.quantile(survey$z, weights=survey$Legal, probs = seq(0, 1, len = k)[-c(1,k)]) ##quantiles of depth weighted by catch, to have splines at depth where catch located

k <- 4 # 5 basis splines
k4 <- Hmisc::wtd.quantile(survey$z, weights=survey$Legal, probs = seq(0, 1, len = k)[-c(1,k)]) ##quantiles of depth weighted by catch, to have splines at depth where catch located
survey$IDS = "I"
survey = cv_SpaceTimeFolds(survey,idCol = 'IDS', nfolds=5)
survey$n = survey$Legal
path=file.path('../Model_outputs/CommercialN')

##model selection table
##set up a blank table

columns <- c("Model", "Formula", "Time-varying", "Family", "AIC", "cAIC", "EDF_ep","EDF_om","rho", "Matern range", "Spatial SD", "Hessian_positive", "Sum loglik", "MAE_train","MAE_test", "RMSE_train", "RMSE_test" )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns

##model selection table function

mod.select.fn <- function (){
  
  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
          colnames(c) <- columns
          c$Model <- mod.label
          c$Formula <-m$formula [1]
          c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
          c$"Family" <- ifelse(m$family[1]=="tweedie", paste0(m$family[1], "(link = ", m$family[2], ")"), m$family["clean_name"])
          c$"cAIC" = cAIC(m,what='cAIC')
          ee = cAIC(m,what='EDF')
          c$"EDF_ep" = ee[1]
          c$"EDF_om" = ee[2]
          
          ##spatial model 
          c$AIC <- AIC (m)
          c$rho <- m$sd_report[[1]]["rho"]
          c$`Matern range` <- m$sd_report[[1]]["range"]
          c$`Spatial SD` <- m$sd_report[[1]]["sigma_O"]
          c$"Hessian_positive" <- m$pos_def_hessian
          
          ##model validation 
          c$"Sum loglik" <- m_cv$sum_loglik
          m_cvTT = sdmTMBcv_tntpreds(m_cv)
          fitTT = dplyr::bind_rows(m_cvTT)
          fitTT$n = fitTT$Legal
          fitTT$sqR = fitTT$n - fitTT$pred
          c$MAE_test<-  with(fitTT[fitTT$tt=='test',],mae(as.numeric(n),as.numeric(pred)))
          c$MAE_train<-  with(fitTT[fitTT$tt=='train',],mae(as.numeric(n),as.numeric(pred)))
          c$RMSE_test <- with(fitTT[fitTT$tt=='test',],rmse(as.numeric(n),as.numeric(pred)))
          c$RMSE_train <- with(fitTT[fitTT$tt=='train',],rmse(as.numeric(n),as.numeric(pred)))
          return(c)
          }

# ###MODELS

models <-c( "m1", "m2",'m3','m4') 
survey$fYear = as.factor(survey$YEAR)
survey$lz = log(survey$z)
survey$depth_scaled = (survey$z - mean(survey$z))/sd(survey$z)

if ("m1" %in% models) {
              mod.label <- "m1" 
              m <- sdmTMB(
                    data = survey,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = survey$of,
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1"
                    )
              m_cv <- sdmTMB_cv(
                    data = survey,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = 'of',
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1",
                    fold_ids='fold_id',
                    k_folds = 5
                    )

              ca <-mod.select.fn()
              mod.select <- rbind(mod.select, ca)
              
              # Save model results:
               dd = format(Sys.Date(), "%Y%m%d")
              saveRDS(m, file = paste0(path, "MAR_Lobster_CommercialN", "_", mod.label,"_",dd,  ".rds"))
              
              m1 <- m
              m1_cv <- m_cv
        }


if ("m2" %in% models) {
              mod.label <- "m2" 
              m <- sdmTMB(
                    data = survey,
                    formula = Legal ~ 0+SOURCE+fYear+s(lz), 
                    offset = survey$of,
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1"
                    )
              m_cv <- sdmTMB_cv(
                    data = survey,
                    formula = Legal ~ 0+SOURCE+fYear+s(lz), 
                    offset = 'of',
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1",
                    fold_ids='fold_id',
                    k_folds = 5
                    )

              ca <-mod.select.fn()
              mod.select <- rbind(mod.select, ca)
              
              # Save model results:
                     dd = format(Sys.Date(), "%Y%m%d")
              saveRDS(m, file = paste0(path, "MAR_Lobster_CommercialN", "_", mod.label,"_",dd,  ".rds"))
        
              m2 <- m
              m2_cv <- m_cv
        }


        if ("m3" %in% models) {
              mod.label <- "m3" 
                # # Build B-spline basis functions 
                    k <- 4
                    knots <- quantile(survey$depth_scaled, p = seq(0, 1, len = k)[-c(1,k)])
                    bs <- bs(survey$depth_scaled, knots = knots, intercept = FALSE)
                    bs <- as.data.frame(bs)
                    names(bs) <- paste0("bs", names(bs))
                    survey <- survey[, setdiff(names(survey), names(survey)[grep("^bs[0-9]+", names(survey))])]
                    survey_bs <- cbind(survey, bs)

              m <- sdmTMB(
                    data = survey_bs,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = survey_bs$of,
                    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
                    time_varying_type = 'rw0',
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1"
                    )
              m_cv <- sdmTMB_cv(
                    data = survey_bs,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = 'of',
                    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 + bs5,
                    time_varying_type = 'rw0',
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1",
                    fold_ids='fold_id',
                    k_folds = 5
                    )

              ca <-mod.select.fn()
              mod.select <- rbind(mod.select, ca)
              
              # Save model results: 
                    dd = format(Sys.Date(), "%Y%m%d")
              saveRDS(m, file = paste0(path, "MAR_Lobster_CommercialN", "_", mod.label,"_",dd,  ".rds"))
        
              m3 <- m
              m3_cv <- m_cv
        }


        if ("m4" %in% models) {
                mod.label <- "m4" 
                # # Build B-spline basis functions 
                    k <- 3
                    knots <- quantile(survey$depth_scaled, p = seq(0, 1, len = k)[-c(1,k)])
                    bs <- bs(survey$depth_scaled, knots = knots, intercept = FALSE)
                    bs <- as.data.frame(bs)
                    names(bs) <- paste0("bs", names(bs))
                    survey <- survey[, setdiff(names(survey), names(survey)[grep("^bs[0-9]+", names(survey))])]
                    survey_bs <- cbind(survey, bs)

              m <- sdmTMB(
                    data = survey_bs,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = survey_bs$of,
                    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4,
                    time_varying_type = 'rw0',
                    mesh = bspde,
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1"
                    )
              m_cv <- sdmTMB_cv(
                    data = survey_bs,
                    formula = Legal ~ 0+SOURCE+fYear, 
                    offset = 'of',
		mesh = bspde,
		    time_varying = ~ 1 + bs1 + bs2 + bs3 + bs4 ,
                    time_varying_type = 'rw0',
                    spatial = "on",
                    family =  tweedie(link = "log"),
                    time = "YEAR",
                    spatiotemporal = "ar1",
                    fold_ids='fold_id',
                    k_folds = 5
                    )

              ca <-mod.select.fn()
              mod.select <- rbind(mod.select, ca)
              
              # Save model results:
               dd = format(Sys.Date(), "%Y%m%d")
              saveRDS(m, file = paste0(path, "MAR_Lobster_CommercialN", "_", mod.label,"_",dd,  ".rds"))
              
              m4 <- m
              m4_cv <- m_cv
        }
