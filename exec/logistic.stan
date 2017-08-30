

functions { 
 
    real[,] vlogistic_surv(real[] t, real[] mu, real[] alpha, int ROW,int COL) {
       real output[ROW,COL]; 
     
       // this is a matrix -> best traversed in row major order
       for (j in 1:ROW){ 
          for (i in 1:COL){
                output[j,i] =  1/(1+exp(mu[j])*t[i]^(exp(alpha[j]))); 
          }
       }  
     
       return output; 
    }
    
   real logistic_surv(real t, real mu, real alpha) { 
        return  1/(1+exp(mu)*t^(exp(alpha))); 
   } 

   real log_logistic_right_censor(real t, real mu, real alpha) { 
        return log(logistic_surv(t,mu,alpha)); 
   } 
  
   real log_logistic_left_censor(real t, real mu, real alpha) { 
        return log(1 - logistic_surv(t,mu,alpha));
   }

   real log_logistic_interval_censor(real tl,real tr, real mu, real alpha) { 
        return log(logistic_surv(tl,mu,alpha)-logistic_surv(tr,mu,alpha)); 
   }

   real log_logistic_exact_lifetime(real t, real mu, real alpha) { 
        return  (log(exp(alpha+mu))+(exp(alpha)-1)*log(t)) - 2*log(1+exp(mu)*t^(exp(alpha))); 
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
   real l; 
   real b; 
//   real<lower = 0>,upper =1.44>  lsig_sq;
   real<lower = 0>  lsig_sq;	
   real l_reff[N_GROUPS]; 
} 

model { 
     l ~ normal(0,1); 
     b ~ normal(0,1);
     l_reff ~ normal(0,sqrt(lsig_sq));
     lsig_sq ~ gamma(1,1);


      for(i in 1:N) { 
            if (CENC[i] == 0) 
            target +=   log_logistic_right_censor(t[i,1],l+l_reff[ID[i]], b ); 
            if (CENC[i] == 1) 
            target +=   log_logistic_exact_lifetime(t[i,1],l+l_reff[ID[i]], b ); 
            if (CENC[i] == 2) 
            target +=   log_logistic_left_censor(t[i,1],l+l_reff[ID[i]], b ); 
            if (CENC[i] == 3) 
            target +=  log_logistic_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b );
      } 
}  
