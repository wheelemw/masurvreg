#########################################
#File MA.R
#Creator: Matthew Wheeler
#Modified Date: 9/4/2017
##########################################


###############################################################
#Fit Individual Exponentially Modified Normal Model
#Inputs:
# data   - survival data
# frange - range of the survival function interested in
#Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
expmNorm_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit <- logExpModNorm_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                                 seed = seed)
  ##################################################
  #extract the parameters
  ##################################################
  est.parms <- extract(fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  nu <- est.parms$nu
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind,vexpmn_surv(t,nu, l+reff, b,length(reff),length(t)))
  
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ##############################################################################
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
  
}

###############################################################
#Fit Individual double exponential  Model
#Inputs:
# data   - survival data
# frange - range of the survival function interested in
#Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
dexp_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- logDoubleExp_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                                 seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind,vdexp_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
  
}


###############################################################
#Fit Individual Weibull  Model
#Inputs:
# data   - survival data
# frange - range of the survival function interested in
#Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
weib_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- weib_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                         seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind, vweib_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
  
}


###############################################################
#Fit Individual Log-Gaussian
#Inputs:
# data   - survival data
# frange - range of the survival function interested in
#Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
logGaussian_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- lognorm_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                            seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind, vlnorm_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
  
}

###############################################################
#Fit Individual Inv-Gaussian
# Inputs:
# data   - survival data
# frange - range of the survival function interested in
# Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
invGaussian_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- invGaussian_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                                seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind,vigauss_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
  
}

###############################################################
#Fit Individual log-Gumbel
# Inputs:
# data   - survival data
# frange - range of the survival function interested in
# Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
logGumbel_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- logGumbel_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                              seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind,vlgum_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
}

###############################################################
#Fit Individual log-Logistic
# Inputs:
# data   - survival data
# frange - range of the survival function interested in
# Output: 
# fit  - Fitted survival model 
# mean - mean estimated survival curve - integrated over the random effects
# ub   - 100*(1-alpha)% upper bound    - integrated over the random effects
# lb   - 100*(alpha)%   lower bound    - integrated over the random effects
###############################################################
logLogistic_model_Analyze <- function(data,frange,alpha = 0.05, sample = 30000,mcmc_warmup=1000,adapt_delta=0.9,seed=8675309){
  
  t <- seq(frange[1],frange[2],(frange[2]-frange[1])/200)
  ####################################################
  #fit the actual model
  ###################################################
  fit  <- logLogistic_model_fit(data_surv,mcmc_iter = sample,mcmc_warmup = mcmc_warmup,adapt_delta = adapt_delta,
                              seed = seed)
  est.parms <- extract(fit$modelFit)
  
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  
  surv <- do.call(rbind,vlogistic_surv(t, l+reff, b,length(reff),length(t)))
  ###############################################################################
  #Create a survival function (mean), the upper bound of that survival function (ub)
  #and the lower bound of that survival functio (lb), this is done over the range of t
  ###############################################################################
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  mean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))
  ############################################################################## 
  return(list(fit=fit ,ub = ub,mean=mean,lb=lb))
}


##################################################################################
#Function: MA_fit
#
#
#
#
#
##################################################################################
MA_fit <- function(data,method="stacking",alpha = 0.05, sample = 30000,mcmc_warmup=2000,adapt_delta=0.9,seed=8675309){
  
    frange = c(0,max(data$t))
    
    fit.llogit   <- logLogistic_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
    fit.lGauss   <- logGaussian_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
    fit.weibul   <- weib_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
    fit.dexp     <- dexp_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
    fit.lGumbel  <- logGumbel_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
    fit.iGauss   <- invGaussian_model_Analyze(data,frange,alpha=alpha,sample=sample,mcmc_warmup = mcmc_warmup,adapt_delta=adapt_delta,
                                            seed = seed)
  
    
    log_lik_list <- list()
    log_lik_list[[1]] = extract(fit.llogit$fit$modelFit)[["log_lik"]]
    log_lik_list[[2]] = extract(fit.lGauss$fit$modelFit)[["log_lik"]]
    log_lik_list[[3]] = extract(fit.weibul$fit$modelFit)[["log_lik"]]
    log_lik_list[[4]] = extract(fit.dexp$fit$modelFit)[["log_lik"]]
    log_lik_list[[5]] = extract(fit.lGumbel$fit$modelFit)[["log_lik"]]
    log_lik_list[[6]] = extract(fit.iGauss$fit$modelFit)[["log_lik"]]
    
    model_weights <- model_weights(log_lik_list,method=method,optim_method = "BFGS",BB=TRUE)
 
    ##############################################################################
    #
    ##############################################################################
    # find the max and adjust for numerical reasons
    ##############################################################################
    names(model_weights) <- c("Log-Logistic","Log-Gaussian","Weibull","Log-DoubleExponential",
                                "Log-Gumbel","Log-invGauss")
    
    ##############################################################################
    #compute the functions for the model average
    ##############################################################################
    return( MA = list(fit.llogit=fit.llogit,fit.lGauss = fit.lGauss,fit.Weibul=fit.weibul,
                      fit.dexp = fit.dexp,  fit.lGumbel = fit.lGumbel,
                      fit.iGauss = fit.iGauss, posterior.probs = model_weights) )
}

##############################################################
# This function returns a function with the estimated mean model 
# average survival curve
#
###########################################################
est.meanMA_surv <- function(MAfits){
names(MAfits$posterior.probs) <- NULL
  return(meanMA <- function(x){
                   return( as.numeric(MAfits$posterior.probs[1]*MAfits$fit.llogit$mean(x) + MAfits$posterior.probs[2]*MAfits$fit.lGauss$mean(x) +
                           MAfits$posterior.probs[3]*MAfits$fit.Weib$mean(x)   + MAfits$posterior.probs[4]*MAfits$fit.dexp$mean(x) + 
                           MAfits$posterior.probs[5]*MAfits$fit.expmNorm$mean(x) + MAfits$posterior.probs[6]*MAfits$fit.lGumbel$mean(x) +  
                           MAfits$posterior.probs[7]*MAfits$fit.iGauss$mean(x)))
              })
}


##############################################################
# This function returns a function with the estimated mean model 
# average survival curve along with the (1-alpha)% and alpha 
# upper and lower bounds respectively
###########################################################
est.MAsurv_functions <- function(MAfits,frange,alpha = 0.05){
  
  t <- seq(frange[1],frange[2],0.1)
  ###################################################################################
  #double exponential  

  est.parms <- extract(MAfits$fit.dexp$fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  dexp.surv <- MAfits$posterior.probs[4]*do.call(rbind,vdexp_surv(t, l+reff, b,length(reff),length(t)))
  #Weibull
  est.parms <- extract(MAfits$fit.Weibul$fit$modelFit)
    b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  weib.surv <- MAfits$posterior.probs[3]*do.call(rbind, vweib_surv(t, l+reff, b,length(reff),length(t)))
  #log-normal
  est.parms <- extract(MAfits$fit.lGauss$fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  lnorm.surv <- MAfits$posterior.probs[2]*do.call(rbind, vlnorm_surv(t, l+reff, b,length(reff),length(t)))
  #log-Gumbel
  est.parms <- extract(MAfits$fit.lGumbel$fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  lgum.surv <-  MAfits$posterior.probs[5]*do.call(rbind,vlgum_surv(t, l+reff, b,length(reff),length(t)))
  ## log-logistic
  est.parms <- extract(MAfits$fit.llogit$fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  llogistic.surv <- MAfits$posterior.probs[1]*do.call(rbind,vlogistic_surv(t, l+reff, b,length(reff),length(t)))
  ## inverse-Gaussian
  est.parms <- extract(MAfits$fit.iGauss$fit$modelFit)
  b <- est.parms$b
  l <- est.parms$l
  reff <- rnorm(length(est.parms$lsig_sq),0,est.parms$lsig_sq)
  igauss.surv <-  MAfits$posterior.probs[6]*do.call(rbind,vigauss_surv(t, l+reff, b,length(reff),length(t)))
  ################################################################
  surv <- igauss.surv + llogistic.surv + lgum.surv + lnorm.surv + weib.surv + dexp.surv 
  
  
  ub <- splinefun(t,apply(surv,2,quantile,alpha,na.rm = T))
  tempmean <- splinefun(t,apply(surv,2,mean,na.rm = T))
  lb <- splinefun(t,apply(surv,2,quantile,1-alpha,na.rm = T))

  return(functions = list(mean=tempmean,ub=ub,lb=lb))
}
