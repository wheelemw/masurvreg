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
  
    return( list(modelFit =fitr))
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
  
  return( list(modelFit =fitr))
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
 
  return( list(modelFit =fitr))
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

  return( list(modelFit =fitr))
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
 
  return( list(modelFit =fitr))
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
 
  return( list(modelFit =fitr))
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
 
  return( list(modelFit =fitr))
}
