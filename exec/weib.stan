

functions { 
     
     
     real[,] vweib_surv(real[] t, real[] l, real[] b, int ROW,int COL) {
       real output[ROW,COL]; 
     
       // this is a matrix -> best traversed in row major order
       for (j in 1:ROW){ 
          for (i in 1:COL){
                output[j,i] =  exp(-exp(-l[j])*t[i]^exp(b[j])); 
          }
       }  
     
       return output; 
    }
   
   
   real weib_surv(real t, real l, real b) { 
        return exp(-exp(-l)*t^exp(b)); 
   } 

   real log_weib_right_censor(real t, real l, real b) { 
        return -exp(-l)*t^exp(b); 
   } 
  
   real log_weib_left_censor(real t, real l, real b) { 
        return log(1-exp(-exp(-l)*t^exp(b))); 
   }

   real log_weib_interval_censor(real tl,real tr, real l, real b) { 
        return log(exp(-exp(-l)*tl^exp(b)) - exp(-exp(-l)*tr^exp(b))); 
   }

   real log_weib_exact_lifetime(real t, real l, real b) { 
        return log(exp(b)*exp(-l))+(exp(b)-1)*log(t)-exp(-l)*t^exp(b); 
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
   real <lower=0, upper=2> lsig_sq;
  // real<lower = 0>  lsig_sq;	
   real l_reff[N_GROUPS]; 
} 

model { 
     l ~ normal(0,1); 
     b ~ normal(0,0.5);
     l_reff ~ normal(0,sqrt(lsig_sq));
     lsig_sq ~ gamma(1,1);
      

     for(i in 1:N) { 
         if (CENC[i] == 0) 
               target +=   log_weib_right_censor(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 1) 
               target +=   log_weib_exact_lifetime(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 2) 
               target +=   log_weib_left_censor(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 3) 
               target +=  log_weib_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b);

     } 
} 

generated quantities{
  vector[N] log_lik;
  for (i in 1:N){
         if (CENC[i] == 0) 
               log_lik[i] =   log_weib_right_censor(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 1) 
               log_lik[i] =   log_weib_exact_lifetime(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 2) 
               log_lik[i] =   log_weib_left_censor(t[i,1],l+l_reff[ID[i]], b); 
         if (CENC[i] == 3) 
               log_lik[i] =  log_weib_interval_censor(t[i,1],t[i,2],l+l_reff[ID[i]], b);
    
  }
}
