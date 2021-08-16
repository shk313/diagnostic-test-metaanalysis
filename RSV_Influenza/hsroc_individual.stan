//
// HSROC model for DTA meta-analysis


data {
  int<lower=0> J; // Number of unique studies
  int<lower=0> l; //Number of studies, > J 
  int<lower=0> N; // Total number of people
  int<lower=0> ind[N]; //Indicator for which study a person belong to and their id within that study
  int<lower=0> study[l]; //Indicator for which unique study a person belong to and their id within that study
  int T1[N]; // index test
  int T2[N]; 
  int<lower=1> ref_test[N]; // indicator for which reference test a person belongs to
  int<lower=0> RT; // Number of reference tests
}




parameters {
  real<lower=0.5,upper=1> Se_Ref[RT];
  real<lower=0.5,upper=1> Sp_Ref[RT];
  real THETA; // global mean of theta (A+B)
  real LAMBDA; // global mean of alpha (A+B)
  real<lower=-0.75, upper=0.75> beta; // common scale param
  real<lower=0> sd_theta; // changed to vector from real
  real<lower=0> sd_alpha;
  vector[J] theta_; // for reparameterisation of theta
  vector[J] alpha_; // for reparameterisation of alpha
  real<lower=0,upper=1> prevalence[l];
  vector[N] N_RE; // e.g. infection intensity
  real<lower=0> sd_re[RT+1]; // sd for infection intensity between positives in each test
  
}

transformed parameters {
  
  vector<lower=0,upper=1>[N] prob[2,2]; // I want it to be N rows, two columns(for Se and Sp), with two shelves(one for each test)
  simplex[2] prev[l];
  vector<lower=0, upper=1>[l] Se; // index test
  vector<lower=0, upper=1>[l] Sp; // index test
  vector[J] theta;
  vector[J] alpha;

  theta = THETA + theta_*sd_theta; // implies theta~N(THETA,prec[1])
  alpha = LAMBDA + alpha_*sd_alpha; //implies alpha~N(LAMBDA,prec[2])


for(j in 1:l){ 
  
  prev[j,1] = 1-prevalence[j];
  prev[j,2] = prevalence[j];
  
   // Equation from Dendukuri et al
  Se[j] = inv_logit(-(theta[study[j]] - alpha[study[j]]/2)/exp(beta/2));
  Sp[j] = inv_logit((theta[study[j]] + alpha[study[j]]/2)/exp(beta/2));

}
  

  // Conditional independence
 // for(n in 1:N){
 //   prob[2,1,n] = Se[ind[n]];   //person, se or sp, test
 //   prob[1,1,n] = 1-Sp[ind[n]];
 //   prob[2,2,n] = Se_Ref[ref_test[n]];
 //   prob[1,2,n] = 1-Sp_Ref[ref_test[n]];
 // }

  // Conditional dependence structure
  for(n in 1:N){
  prob[2,1,n] = inv_logit(logit(Se[ind[n]]) + (N_RE[n]*sd_re[7]));
  prob[1,1,n] = inv_logit(logit(1-Sp[ind[n]]));
  prob[2,2,n] = inv_logit(logit(Se_Ref[ref_test[n]]) + N_RE[n]*sd_re[ref_test[n]]);
  prob[1,2,n] = inv_logit(logit(1-Sp_Ref[ref_test[n]]));
  }
  
}

model {
  real lp[2]; // 2d vector with each entry corresponding
  // to the un-normalised posterior prob for diseased and not diseased
  

  THETA~normal(0,1.5);
  LAMBDA~normal(0,2);

  sd_re~normal(0,1);
  N_RE~std_normal();

  theta_ ~ std_normal(); // implies theta~normal(THETA, prec[1])
  alpha_ ~ std_normal(); // implies alpha~normal(LAMBDA, prec[1])

  sd_theta ~ normal(0,1);
  sd_alpha ~ normal(0,1);
  
  prevalence~beta(1,1);
  
  Se_Ref~beta(5,1);
  Sp_Ref~beta(10,1);


  for(n in 1:N){
    for(k in 1:2){  
      // log posterior probability
      
      lp[k] = log(prev[ind[n],k]);
      lp[k] = lp[k] + binomial_lpmf(T1[n] | 1, prob[k,1,n]) + binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }
    // marginalise out dependence on k and increment the overall log prob. by this amount
    target += log_sum_exp(lp);
  }

}  

generated quantities {
 real<lower=0,upper=1> Se_pooled;
 real<lower=0,upper=1> Sp_pooled;
 real<lower=0,upper=1> Se_pred;
 real<lower=0,upper=1> Sp_pred;
  

 vector[N] log_lik;
 real ll[2];
 real theta_new;
 real alpha_new;

 theta_new = normal_rng(THETA,sd_theta);
 alpha_new = normal_rng(LAMBDA,sd_alpha);
 
 Se_pooled = inv_logit(-(THETA - LAMBDA/2) / exp(beta/2));
 Sp_pooled = inv_logit((THETA + LAMBDA/2) / exp(beta/2));

 Se_pred= inv_logit(-(theta_new - alpha_new/2) / exp(beta/2));
 Sp_pred= inv_logit((theta_new + alpha_new/2) / exp(beta/2));


  for(n in 1:N){
    for(k in 1:2){
      ll[k] = log(prev[ind[n],k]) +  binomial_lpmf(T1[n] | 1, prob[k,1,n]) +  binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }

  log_lik[n] = log_sum_exp(ll);
  }
 

}
