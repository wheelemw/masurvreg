####################################################
# Fit the Weibul Model
#
####################################################
weib_model_fit <- function(data_surv,mcmc_iter = 30000, mcmc_warmup = 1000,
                           seed=8675309, adapt_delta = 0.9)
{

    fitr =    sampling(stanmodels$wei,data=data_surv,
                         iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,chains=1,
                         control=list(adapt_delta = adapt_delta))
    expose_stan_functions(stanmodels$wei)
    ################################################################
    #### compute the log prior probability WITH integrating factors
    ################################################################
    est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
    b <- est.parms$b
    l <- est.parms$l
    lsig_sq <- est.parms$lsig_sq
    
    l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
    target = 0; 
    for(i in 1:data_surv$N) { 
      if (data_surv$CENC[i] == 0) 
        target = target +   log_weib_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
      if (data_surv$CENC[i] == 1) 
        target = target +  log_weib_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
      if (data_surv$CENC[i] == 2) 
        target =  target + log_weib_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
      if (data_surv$CENC[i] == 3) 
        target = target +  log_weib_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
    }
    
    target = target + dnorm(l, log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE) - pgamma(4,1,1,log.p = T)
    
    for (i in 1:length(l_reff)){
      target = target + dnorm(l_reff[i],0,lsig_sq)
    }
    #####################################################
    #####################################################
    
    nparms <- ncol(as.matrix(fitr)) -1
    log.wei.p = target
    cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
    
    bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
    names(bayesFactor) <- "LogBF"
    return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the LogNormal
#####################################################
lognorm_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                              seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$lognormal,data=data_surv,chains=1,
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,
                     control=list(adapt_delta = adapt_delta))
  expose_stan_functions(stanmodels$lognormal)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
  b <- est.parms$b
  l <- est.parms$l
  lsig_sq <- est.parms$lsig_sq
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target = target +   log_lnorm_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 1) 
      target = target +  log_lnorm_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 2) 
      target =  target + log_lnorm_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 3) 
      target = target +  log_lnorm_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
  }
  
  target = target + dnorm(l, log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE)  -pgamma(4,1,1,log.p = T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p =target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the inverse Gaussian
#####################################################
invGaussian_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                                  seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$inv_gaussian,data=data_surv,
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,chains=1,
                     control=list(adapt_delta = adapt_delta,max_treedepth=12))
                     
  expose_stan_functions(stanmodels$inv_gaussian)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
  b <- est.parms$b
  l <- est.parms$l
  lsig_sq <- est.parms$lsig_sq
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target = target +   log_igauss_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 1) 
      target = target +  log_igauss_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 2) 
      target =  target + log_igauss_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 3) 
      target = target +  log_igauss_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
  }
  
  target = target + dnorm(l,0,10, log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE)  -pgamma(4,1,1,log.p = T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p = target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the LogGumbel
#####################################################
logGumbel_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                                seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$loggumbel,data=data_surv, chains = 1, 
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,
                     control=list(adapt_delta = adapt_delta))
  expose_stan_functions(stanmodels$loggumbel)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
  b <- est.parms$b
  l <- est.parms$l
  lsig_sq <- est.parms$lsig_sq
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
  target = 0; 
  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target = target +   log_lgum_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 1) 
      target = target +  log_lgum_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 2) 
      target =  target + log_lgum_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 3) 
      target = target +  log_lgum_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
  }
  
  target = target + dnorm(l, log=TRUE) + dnorm(b,0,10,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE)  -pgamma(4,1,1,log.p = T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p = target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the LogExpModNormal
#####################################################
logExpModNorm_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                                    seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$expmodn,data=data_surv,
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,chains=1,
                     control=list(adapt_delta = adapt_delta))
  expose_stan_functions(stanmodels$expmodn)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:4])
  nu <- est.parms$nu
  b <- est.parms$b
  l <- est.parms$l
  lsig_sq <- est.parms$lsig_sq
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:4,ncol(as.matrix(fitr)))]

  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target = target +   log_expmn_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(nu),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 1) 
      target = target +  log_expmn_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(nu),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 2) 
      target =  target + log_expmn_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(nu),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 3) 
      target = target +  log_expmn_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(nu),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
  }
  
  
  target = target + dnorm(l,0,10,log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE) -pgamma(4,1,1,log.p=T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p = target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the LogDoubleExponential
#####################################################
logDoubleExp_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                                   seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$logdexp,data=data_surv,chains=1,
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,
                     control=list(adapt_delta = adapt_delta))
  expose_stan_functions(stanmodels$logdexp)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
 
  b <- as.numeric(est.parms$b)
  l <- as.numeric(est.parms$l)
  lsig_sq <- as.numeric(est.parms$lsig_sq)
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
  
  
  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target =  target+ log_ldexp_right_censor(as.numeric(data_surv$t[i,1]), l+l_reff[ID[i]], b ); 
    if (data_surv$CENC[i] == 1) 
      target =  target+ log_ldexp_exact_lifetime(as.numeric(data_surv$t[i,1]),l+l_reff[ID[i]], b); 
    if (data_surv$CENC[i] == 2) 
      target =  target+ log_ldexp_left_censor(as.numeric(data_surv$t[i,1]),l+l_reff[ID[i]], b); 
    if (data_surv$CENC[i] == 3) 
      target =  target+ log_ldexp_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),l+l_reff[ID[i]], b);
  }
  
  
  target = target + dnorm(l,0,10,log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE)-pgamma(4,1,1,log.p=T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p = target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}

####################################################
# Fit the LogLogistic
#####################################################
logLogistic_model_fit <- function(data_surv,mcmc_iter = 10000, mcmc_warmup = 1000,
                                  seed=8675309, adapt_delta = 0.9)
{
  
  fitr =    sampling(stanmodels$logistic,data=data_surv,
                     iter=mcmc_iter,warmup=mcmc_warmup,seed=seed,chains=1,
                     control=list(adapt_delta = adapt_delta))
  expose_stan_functions(stanmodels$logistic)
  ################################################################
  #### compute the log prior probability WITH integrating factors
  ################################################################
  est.parms <- as.list(colMeans(as.matrix(fitr))[1:3])
  
  b <- as.numeric(est.parms$b)
  l <- as.numeric(est.parms$l)
  lsig_sq <- as.numeric(est.parms$lsig_sq)
  
  l_reff   <-  colMeans(as.matrix(fitr))[-c(1:3,ncol(as.matrix(fitr)))]
  
  target = 0; 
  for(i in 1:data_surv$N) { 
    if (data_surv$CENC[i] == 0) 
      target = target +   log_logistic_right_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 1) 
      target = target +  log_logistic_exact_lifetime(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 2) 
      target =  target + log_logistic_left_censor(as.numeric(data_surv$t[i,1]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) ); 
    if (data_surv$CENC[i] == 3) 
      target = target +  log_logistic_interval_censor(as.numeric(data_surv$t[i,1]),as.numeric(data_surv$t[i,2]),as.numeric(l+l_reff[ID[i]]), as.numeric(b) );
  }
  
  
  target = target + dnorm(l,log=TRUE) + dnorm(b,log=TRUE) + dgamma(lsig_sq,1,1,log=TRUE) -pgamma(4,1,1,log.p=T)
  
  for (i in 1:length(l_reff)){
    target = target + dnorm(l_reff[i],0,lsig_sq)
  }
  #####################################################
  #####################################################
  
  nparms <- ncol(as.matrix(fitr)) -1
  log.wei.p = target
  cov <- cov(as.matrix(fitr)[,-ncol(as.matrix(fitr))]) 
  
  bayesFactor = log((2*pi)^(nparms)) + 0.5*log(det(cov)) + log.wei.p
  names(bayesFactor) <- "LogBF"
  return( list(modelFit =fitr,bayesFactor = bayesFactor))
}
