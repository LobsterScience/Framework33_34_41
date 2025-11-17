#


columns <- c("Model", "Formula", "Time-varying", "Family", "AIC", "cAIC", "EDF_ep","EDF_om","rho",
             "Matern range", "Spatial SD", "Hessian_positive", "Sum loglik", "MAE_train","MAE_test", "RMSE_train", "RMSE_test" )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns

##model selection table function

mod.select.fn <- function (){
  
  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
  colnames(c) <- columns
  c$Model <- mod.label
  c$Formula <-m$formula [1]
  c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
  c$"Family" <- paste0(m$family[1]$family, "(link = ", m$family$link[1], ")")
  c$"cAIC" = cAIC(m,what='cAIC')
 #ee = cAIC(m,what='EDF')
  #c$"EDF_ep" = ee[1]
  #c$"EDF_om" = ee[2]
  
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
  fitTT$n = fitTT$diff
  fitTT$sqR = fitTT$n - fitTT$pred
  c$MAE_test<-  with(fitTT[fitTT$tt=='test',],mae(as.numeric(n),as.numeric(pred)))
  c$MAE_train<-  with(fitTT[fitTT$tt=='train',],mae(as.numeric(n),as.numeric(pred)))
  c$RMSE_test <- with(fitTT[fitTT$tt=='test',],rmse(as.numeric(n),as.numeric(pred)))
  c$RMSE_train <- with(fitTT[fitTT$tt=='train',],rmse(as.numeric(n),as.numeric(pred)))
  return(c)
}
