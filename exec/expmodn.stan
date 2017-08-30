

functions { 
  
     real[,] vexpmn_surv(real[] t,real[] nu, real[] mu, real[] sigma, int ROW, int COL) {
      real output[ROW,COL]; 
     
     // this is a matrix -> best traversed in row major order
     for (j in 1:ROW){ 
        for (i in 1:COL){
              output[j,i] =  exp(exp_mod_normal_lccdf(log(t[i]) | mu[j], exp(sigma[j]),nu[j]));
        }
     }  
     
     return output; 
    }
  
     real expmn_surv(real t,real nu, real mu, real sigma) { 
            return exp(exp_mod_normal_lccdf(log(t) | mu, exp(sigma),nu)); 
    } 

    real log_expmn_right_censor(real t,real nu, real mu, real sigma) { 
        return exp_mod_normal_lccdf(log(t) | mu, exp(sigma),nu); 
    } 

    real log_expmn_left_censor(real t,real nu, real mu, real sigma) { 
        return exp_mod_normal_lcdf(log(t) | mu, exp(sigma),nu); 
    }

    real log_expmn_interval_censor(real lt, real rt, real nu, real mu, real sigma) { 
        return log(exp(exp_mod_normal_lccdf(log(lt) | mu, exp(sigma),nu)) - exp(exp_mod_normal_lccdf(log(rt) | mu, exp(sigma),nu))); 
    }

    real log_expmn_exact_lifetime(real t,real nu, real mu, real sigma) { 
      return exp_mod_normal_lpdf(log(t) |mu, exp(sigma),nu); 
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
   real nu; 
   real l; 
   real b; 
   real<lower = 0,upper = 4>  lsig_sq;
  // real<lower = 0>  lsig_sq;	
   real l_reff[N_GROUPS]; 
} 

model { 
     l ~ normal(0,10); 
     b ~ normal(0,1);
     l_reff ~ normal(0,sqrt(lsig_sq));
     lsig_sq ~ gamma(1,1);
     nu ~ gamma(1,1);
   
     
     for(i in 1:N) { 
         if (CENC[i] == 0) 
               target +=   log_expmn_right_censor(t[i,1],nu, l+l_reff[ID[i]], b ); 
         if (CENC[i] == 1) 
               target +=   log_expmn_exact_lifetime(t[i,1],nu,l+l_reff[ID[i]], b); 
         if (CENC[i] == 2) 
               target +=   log_expmn_left_censor(t[i,1],nu,l+l_reff[ID[i]], b); 
         if (CENC[i] == 3) 
               target +=  log_expmn_interval_censor(t[i,1],t[i,2],nu,l+l_reff[ID[i]], b);
     } 
} 
