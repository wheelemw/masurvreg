

functions { 
 
    real[,] vigauss_surv(real[] t, real[] mu, real[] alpha, int ROW,int COL) {
      real output[ROW,COL]; 
     
     // this is a matrix -> best traversed in row major order
     for (j in 1:ROW){ 
        for (i in 1:COL){
              output[j,i] =  (normal_cdf( sqrt(exp(alpha[j])/t[i])*(1-t[i]/exp(mu[j])) , 0, 1)  - exp(2*exp(alpha[j]-mu[j]))*normal_cdf( -sqrt(exp(alpha[j])/t[i])*(1+t[i]/exp(mu[j])), 0, 1)); 
        }
     }  
     
     return output; 
    }
 
   real igauss_surv(real t, real mu, real alpha) { 
        return  (normal_cdf( sqrt(exp(alpha)/t)*(1-t/exp(mu)) , 0, 1)  - exp(2*exp(alpha-mu))*normal_cdf( -sqrt(exp(alpha)/t)*(1+t/exp(mu)), 0, 1)); 
   } 

   real log_igauss_right_censor(real t, real mu, real alpha) { 
        return log(igauss_surv(t,mu,alpha)); 
   } 
  
   real log_igauss_left_censor(real t, real mu, real alpha) { 
        return log(1 - igauss_surv(t,mu,alpha));
   }

   real log_igauss_interval_censor(real tl,real tr, real mu, real alpha) { 
        return log(igauss_surv(tl,mu,alpha)-igauss_surv(tr,mu,alpha)); 
   }

   real log_igauss_exact_lifetime(real t, real mu, real alpha) { 
        return (0.5)*(log(exp(alpha))-log(2*3.14159265*t^3)) - alpha*(t-exp(mu))^2/(2*exp(mu)^2*t);
   }  
} 
 

data { 
   int<lower=0> N; 
   int<lower=1> N_GROUPS; 
   real t[N,2];
   real CENC[N];
   int<lower=1> ID[N];
} 
 

parameters { 
   //real l  
   real l; 
   real b; 
   real <lower=0, upper=2> lsig_sq;
  // real<lower = 0>  lsig_sq;	
   real l_reff[N_GROUPS]; 
} 

model { 
     l ~ normal(0,10); 
     b ~ normal(0,1);
     l_reff ~ normal(0,sqrt(lsig_sq));
     lsig_sq ~ gamma(1,1);


      for(i in 1:N) { 
            if (CENC[i] == 0) 
            target +=   log_igauss_right_censor(t[i,1], l+l_reff[ID[i]], b); 
            if (CENC[i] == 1) 
            target +=   log_igauss_exact_lifetime(t[i,1],l+l_reff[ID[i]], b); 
            if (CENC[i] == 2) 
            target +=   log_igauss_left_censor(t[i,1],l+l_reff[ID[i]], b); 
            if (CENC[i] == 3) 
            target +=  log_igauss_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b);
      } 
}  

generated quantities{
  vector[N] log_lik; 
  
 for(i in 1:N) { 
        if (CENC[i] == 0) 
            log_lik[i] =   log_igauss_right_censor(t[i,1], l+l_reff[ID[i]], b); 
        if (CENC[i] == 1) 
            log_lik[i] =   log_igauss_exact_lifetime(t[i,1],l+l_reff[ID[i]], b); 
        if (CENC[i] == 2) 
            log_lik[i] =   log_igauss_left_censor(t[i,1],l+l_reff[ID[i]], b); 
        if (CENC[i] == 3) 
            log_lik[i] =  log_igauss_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b);
      }  

}
