## HSROC model with covariate for study sub-type



data {
  int<lower=0> J; // Number of unique studies
  int<lower=0> l; //Number of studies, > J 
  int<lower=0> N; // Total number of people
  int<lower=0> ind[N]; //Indicator for which study a person belong to and their id within that study
  int<lower=0> study[l]; //Indicator for which unique study a person belong to and their id within that study
  int T1[N]; // index test
  int T2[N]; 
  vector<lower=0, upper=1>[l] pathogenA; // 1 = yes to pathogen A, 0 = no
  vector<lower=0, upper=1>[l] pathogenB;
  //vector<lower=0, upper=1>[l] pathogenAB; // For RSV
  int<lower=1, upper=2> path[l];
  int<lower=1> ref_test[N]; // indicator for which reference test a person belongs to
  int<lower=0> RT; // Number of reference tests
}


parameters {
  real<lower=0.5,upper=1> Se_Ref[RT];
  real<lower=0.5,upper=1> Sp_Ref[RT];

  real THETA_pathogen[2]; // global mean of theta, 2 is the number of sub-types  
  real LAMBDA_pathogen[2]; // global mean of alpha  

  real<lower=-0.75, upper=0.75>beta_pathogen[2];
  
  real<lower=0> sd_theta[2]; // changed to vector from real
  real<lower=0> sd_alpha[2];
  vector[l] theta_; // for reparameterisation of theta
  vector[l] alpha_; // for reparameterisation of alpha
 
  real<lower=0,upper=1> prevalence[l];
  vector[N] N_RE; // e.g. infection intensity
  real<lower=0> sd_re[RT+1]; // sd for infection intensity between positive [1] and negatives [2]
  
}

transformed parameters {
  
  vector<lower=0,upper=1>[N] prob[2,2]; 
  simplex[2] prev[l];
  real<lower=0, upper=1> Se[l]; // index test
  real<lower=0, upper=1> Sp[l];
  vector[l] theta;
  vector[l] alpha;
  vector[l] beta;
  vector[l] T;
  vector[l] L;


  // theta = THETA + theta_*sd_theta; // implies theta~N(THETA,prec[1])
  // alpha = LAMBDA + alpha_*sd_alpha; //implies alpha~N(LAMBDA,prec[2])


for(j in 1:l){ 
  
  prev[j,1] = 1-prevalence[j];
  prev[j,2] = prevalence[j];
  
  // beta, theta and alpha mean are now a linear combination of the covariate for sub-type:
  beta[j] = exp((beta_pathogen[1]*pathogenA[j] + beta_pathogen[2]*pathogenB[j])/2);
  T[j] = THETA_pathogen[1]*pathogenA[j] + THETA_pathogen[2]*pathogenB[j];
  L[j] = LAMBDA_pathogen[1]*pathogenA[j] + LAMBDA_pathogen[2]*pathogenB[j];
 
  theta[j] = T[j]+ theta_[j]* sd_theta[path[j]]; // implies theta~N(THETA,prec[1])
  alpha[j] = L[j] + alpha_[j]* sd_alpha[path[j]]; //implies alpha~N(LAMBDA,prec[2])

  
// Equation from Dendukuri et al

  Se[j] = inv_logit(-(theta[j] - alpha[j]/2)/beta[j]);
  Sp[j] = inv_logit((theta[j] + alpha[j]/2)/beta[j]);

  
}
  

// Conditional Independence  
 // for(n in 1:N){
 //   prob[2,1,n] = Se[ind[n]];   //person, se or sp, test
 //   prob[1,1,n] = 1-Sp[ind[n]];
 //   prob[2,2,n] = Se_Ref[ref_test[n]];
 //   prob[1,2,n] = 1-Sp_Ref[ref_test[n]];
 // }
  
  

  //Conditional dependence
  for(n in 1:N){
    prob[2,1,n] = inv_logit(logit(Se[ind[n]]) + (N_RE[n]*sd_re[6]));
    prob[1,1,n] = inv_logit(logit(1-Sp[ind[n]]));
    prob[2,2,n] = inv_logit(logit(Se_Ref[ref_test[n]]) + N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = inv_logit(logit(1-Sp_Ref[ref_test[n]]));
  }
   
}

model {
  real lp[2];
  
  //THETA~normal(0,1.5);
  //LAMBDA~normal(0,2);
  THETA_pathogen~normal(0,1.5);
  LAMBDA_pathogen~normal(0,2);
  // beta with or without covariate
  //beta_pathogen~normal(0,0.75);
  //beta~normal(0,0.75);
  
  // non-centred parameterisation
  // helps sampler explore more efficiently
  // to standardise typically subtract mean and divide by sd
  // do this in reverse
  // multiply standardised variable N_RE by sd
  // both means are 0 in this case so don't need to add the mean back on
  sd_re~normal(0,1);
  N_RE~std_normal(); // implies N_RE ~ normal(0, sd_re)

  theta_ ~ std_normal(); // implies theta~normal(THETA, prec[1])
  alpha_ ~ std_normal(); // implies alpha~normal(LAMBDA, prec[1])


 sd_theta ~ normal(0,1);
 sd_alpha ~ normal(0,1);
  // alternative priors
  // sd_theta ~ gamma(2.1,1);
  // sd_alpha ~ gamma(2.1,1);
  
  prevalence~beta(1,1);
  
  Se_Ref~beta(5,1);
  Sp_Ref~beta(10,1);


  for(n in 1:N){
    for(k in 1:2){                                   
      lp[k] = log(prev[ind[n],k]);
      lp[k] = lp[k] + binomial_lpmf(T1[n] | 1, prob[k,1,n]) + binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }
    target += log_sum_exp(lp);
  }


}  

generated quantities {

 // real<lower=0,upper=1> Se_pooled_AB;
 // real<lower=0,upper=1> Sp_pooled_AB;
 // real<lower=0,upper=1> Se_pred_AB;
 // real<lower=0,upper=1> Sp_pred_AB;
 real<lower=0,upper=1> Se_pooled_A;
 real<lower=0,upper=1> Sp_pooled_A;
 real<lower=0,upper=1> Se_pred_A;
 real<lower=0,upper=1> Sp_pred_A;
 real<lower=0,upper=1> Se_pooled_B;
 real<lower=0,upper=1> Sp_pooled_B;
 real<lower=0,upper=1> Se_pred_B;
 real<lower=0,upper=1> Sp_pred_B;

 vector[N] log_lik;
 real ll[2];
 real alpha_new[2];
 real theta_new[2];
 
 theta_new = normal_rng(THETA_pathogen, sd_theta);
 alpha_new= normal_rng(LAMBDA_pathogen, sd_alpha);

 // Se_pooled_AB = inv_logit(-(THETA_pathogen[1] - LAMBDA_pathogen[1]/2) / exp(beta_pathogen[1]/2));
 // Sp_pooled_AB = inv_logit((THETA_pathogen[1] + LAMBDA_pathogen[1]/2) / exp(beta_pathogen[1]/2));
 Se_pooled_A = inv_logit(-(THETA_pathogen[1] - LAMBDA_pathogen[1]/2) / exp(beta_pathogen[1]/2));
 Sp_pooled_A = inv_logit((THETA_pathogen[1] + LAMBDA_pathogen[1]/2) / exp(beta_pathogen[1]/2));
 Se_pooled_B = inv_logit(-(THETA_pathogen[2] - LAMBDA_pathogen[2]/2) / exp(beta_pathogen[2]/2));
 Sp_pooled_B = inv_logit((THETA_pathogen[2] + LAMBDA_pathogen[2]/2) / exp(beta_pathogen[2]/2));

 // Se_pred_AB= inv_logit(-(theta_new[1] - alpha_new[1]/2) / exp(beta_pathogen[1]/2));
 // Sp_pred_AB= inv_logit((theta_new[1] + alpha_new[1]/2) / exp(beta_pathogen[1]/2));
 Se_pred_A= inv_logit(-(theta_new[1] - alpha_new[1]/2) / exp(beta_pathogen[1]/2));
 Sp_pred_A= inv_logit((theta_new[1] + alpha_new[1]/2) / exp(beta_pathogen[1]/2));
 Se_pred_B= inv_logit(-(theta_new[2] - alpha_new[2]/2) / exp(beta_pathogen[2]/2));
 Sp_pred_B= inv_logit((theta_new[2] + alpha_new[2]/2) / exp(beta_pathogen[2]/2));

 
 // Likelihood for use in LOO-CV
  for(n in 1:N){
    for(k in 1:2){
      ll[k] = log(prev[ind[n],k]) +  binomial_lpmf(T1[n] | 1, prob[k,1,n]) +  binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }

  log_lik[n] = log_sum_exp(ll);
  }
 

}
