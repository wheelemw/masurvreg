

functions { 
   
     real[,] vdexp_surv(real[] t, real[] mu, real[] sigma, int ROW, int COL) {
          real output[ROW,COL]; 
         
         // this is a matrix -> best traversed in row major order
         for (j in 1:ROW){ 
            for (i in 1:COL){
                  output[j,i] = exp(double_exponential_lccdf(log(t[i]) | mu[j], exp(sigma[j])));
            }
         }  
         
         return output; 
    }
   
     real ldexp_surv(real t, real mu, real sigma) { 
            return exp(double_exponential_lccdf(log(t) | mu, exp(sigma))); 
    } 

    real log_ldexp_right_censor(real t, real mu, real sigma) { 
        return double_exponential_lccdf(log(t) | mu, exp(sigma)); 
    } 

    real log_ldexp_left_censor(real t, real mu, real sigma) { 
        return double_exponential_lcdf(log(t) | mu, exp(sigma)); 
    }

    real log_ldexp_interval_censor(real lt, real rt,  real mu, real sigma) { 
        return log(exp(double_exponential_lccdf(log(lt) | mu, exp(sigma))) - exp(double_exponential_lccdf(log(rt) | mu, exp(sigma)))); 
    }

    real log_ldexp_exact_lifetime(real t, real mu, real sigma) { 
      return double_exponential_lpdf(log(t) | mu, exp(sigma)); 
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
   real<lower = 0,upper =4>  lsig_sq;
   real l_reff[N_GROUPS]; 
} 

model { 
     l_reff ~ normal(0,lsig_sq);
     b ~ normal(0,10);
     lsig_sq ~ gamma(1,1);

     l ~ normal(0,1); 
   

     for(i in 1:N) { 
         if (CENC[i] == 0) 
               target +=   log_ldexp_right_censor(t[i,1], l+l_reff[ID[i]], b ); 
         if (CENC[i] == 1) 
               target +=   log_ldexp_exact_lifetime(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 2) 
               target +=   log_ldexp_left_censor(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 3) 
               target +=  log_ldexp_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b);
     } 
} 
