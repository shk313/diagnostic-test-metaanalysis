//
// HSROC not individual level model for DTA meta-analysis

data {
  int<lower=0> J; // Number of unique studies
  int<lower=0> l; //Number of studies, > J 
  int<lower=0> N; // Total number of people
  int<lower=0> ind_unique[N]; //Indicator for which study a person belong to and their id within that study
  //int<lower=0> ind_study[N]; //Indicator for which study a person belong to and their id within that study
  int<lower=0> study[l]; //Indicator for which unique study a person belong to and their id within that study
  int T1[N]; // index test
  int T2[N]; 
  int<lower=1> ref_test[N]; // indicator for which reference test a person belongs to
  int<lower=1> RT; // Number of reference tests

  
}




parameters {
  real THETA; // global mean of theta (A+B)
  real LAMBDA; // global mean of alpha (A+B)
  real<lower=-0.75, upper=0.75> beta; // common scale param
  real<lower=0> sd_theta; // changed to vector from real
  real<lower=0> sd_alpha;
  vector[1] theta_; // for reparameterisation of theta
  vector[1] alpha_; // for reparameterisation of alpha
  real<lower=0,upper=1> prevalence[J];


  real<lower=0,upper=1> Sp_Ref1;
  real<lower=inv_logit(logit(1-Sp_Ref1)),upper=1> Se_Ref1;
  real<lower=0,upper=1> Sp_Ref2;
  real<lower=inv_logit(logit(1-Sp_Ref2)),upper=1> Se_Ref2;
  real<lower=0,upper=1> Sp_Ref3;
  real<lower=inv_logit(logit(1-Sp_Ref3)),upper=1> Se_Ref3;
  real<lower=0,upper=1> Sp_Ref4;
  real<lower=inv_logit(logit(1-Sp_Ref4)),upper=1> Se_Ref4;
  real<lower=0,upper=1> Sp_Ref5;
  real<lower=inv_logit(logit(1-Sp_Ref5)),upper=1> Se_Ref5;
  // real<lower=0,upper=1> Sp_Ref6;
  // real<lower=inv_logit(logit(1-Sp_Ref6)),upper=1> Se_Ref6;
  // real<lower=0,upper=1> Sp_Ref7;
  // real<lower=inv_logit(logit(1-Sp_Ref7)),upper=1> Se_Ref7;
  //real<lower=0,upper=1> Se_Ref[RT];
  real<lower=0> fe_pos[RT];
}

transformed parameters {
  
  vector<lower=0,upper=1>[N] prob[2,2]; // I want it to be N rows, two columns(for Se and Sp), with two shelves(one for each test)
  simplex[2] prev[J];
  vector<lower=0, upper=1>[l] Sp; // index test
  vector<lower=0, upper=1>[l] Se; // index test
  vector[1] theta;
  vector[1] alpha;
  real y;
  

  theta = THETA + theta_*sd_theta; // implies theta~N(THETA,prec[1])
  alpha = LAMBDA + alpha_*sd_alpha; // Just for a single J.
  
for(j in 1:J){

  
  prev[j,1] = 1-prevalence[j];
  prev[j,2] = prevalence[j];
}


for(j in 1:l){ 
  
   // Equation from Dendukuri et al
  Se[j] = inv_logit(-(theta[1] - alpha[1]/2)/exp(beta/2));
  Sp[j] = inv_logit((theta[1] + alpha[1]/2)/exp(beta/2));

}
  

  // Conditional independence
 // for(n in 1:N){
 //   prob[2,1,n] = Se[ind_unique[n]];   //person, se or sp, test
 //   prob[1,1,n] = 1-Sp[ind_unique[n]];
 //   prob[2,2,n] = Se_Ref[ref_test[n]];
 //   prob[1,2,n] = 1-Sp_Ref[ref_test[n]];
 // }
 
 
for(n in 1:N){
  //x = ind_unique[n]; // grab study id of person
  y = ref_test[n]; // grab reference test id of person
  
  //prob[2,1,n] = Se[ind_unique[n]]; //person, se or sp, test
  prob[2,1,n] = inv_logit(logit(Se[ind_unique[n]])+fe_pos[ref_test[n]]);
  prob[1,1,n] = 1-Sp[ind_unique[n]];
  
 if(y == 1){
    //prob[2,2,n] = Se_Ref1;
    prob[2,2,n] = inv_logit(logit(Se_Ref1)+fe_pos[1]);
    prob[1,2,n] = 1-Sp_Ref1;
    }
     if(y == 2){
    //prob[2,2,n] = Se_Ref2;
    prob[2,2,n] = inv_logit(logit(Se_Ref2)+fe_pos[2]);
    prob[1,2,n] = 1-Sp_Ref2;
    }
      if(y == 3){
    //prob[2,2,n] = Se_Ref3;
    prob[2,2,n] = inv_logit(logit(Se_Ref3)+fe_pos[3]);
    prob[1,2,n] = 1-Sp_Ref3;
    }
      if(y == 4){
    //prob[2,2,n] = Se_Ref4;
    prob[2,2,n] = inv_logit(logit(Se_Ref4)+fe_pos[4]);
    prob[1,2,n] = 1-Sp_Ref4;
    }
    if(y == 5){
    //prob[2,2,n] = Se_Ref5;
    prob[2,2,n] = inv_logit(logit(Se_Ref5)+fe_pos[5]);
    prob[1,2,n] = 1-Sp_Ref5;
    }
    // if(y == 6){
    // //prob[2,2,n] = Se_Ref6;
    // prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // prob[1,2,n] = 1-Sp_Ref6;
    // }



 }

 

  
}

model {
  real lp[2]; // 2d vector with each entry corresponding
  // to the un-normalised posterior prob for diseased and not diseased
  

  THETA~normal(0,1);
  LAMBDA~normal(0,1.5); //0,2
  
  fe_pos~gamma(1,1);

  theta_ ~ std_normal(); // implies theta~normal(THETA, prec[1])
  alpha_ ~ std_normal(); // implies alpha~normal(LAMBDA, prec[1])

  sd_theta ~ normal(0,1);
  sd_alpha ~ normal(0,1);
  
  prevalence~beta(1,1);
  
  // Serum grouped LA and CALAS
  Se_Ref1~beta(1,1); // Culture
  Se_Ref2~beta(1,1); // LA/CALAS
  Se_Ref3~beta(5,1); // Clinical
  Se_Ref4~beta(1,1); // Lung Biopsy
  Se_Ref5~beta(5,1); // Clinical
  // Se_Ref6~beta(5,1); 

  Sp_Ref1~beta(50,1); // Culture
  Sp_Ref2~beta(1,1); // LA/CALAS
  Sp_Ref3~beta(1,1); // Clinical
  Sp_Ref4~beta(5,1); // Lung Biopsy
  Sp_Ref5~beta(1,1); // Clinical
  // Sp_Ref6~beta(1,1); // Clinical

// CSF
  // Se_Ref1~beta(1,1); // Culture
  // Se_Ref2~beta(1,1); // LA
  // Se_Ref3~beta(1,1); // Clinical
  // 
  // Sp_Ref1~beta(50,1); // Culture
  // Sp_Ref2~beta(1,1); // LA
  // Sp_Ref3~beta(1,1); // Clinical
  


  for(n in 1:N){
    for(k in 1:2){  
      // log posterior probability
      
      lp[k] = log(prev[ind_unique[n],k]);
      lp[k] = lp[k] + binomial_lpmf(T1[n] | 1, prob[k,1,n]) + binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }
    // marginalise out dependence on k and increment the overall log prob. by this amount
    target += log_sum_exp(lp);
  }

}  

generated quantities {
 real<lower=0,upper=1> Se_pooled;
 real<lower=0,upper=1> Sp_pooled;
 // real<lower=0,upper=1> Se_pred;
 // real<lower=0,upper=1> Sp_pred;
  

 // vector[N] log_lik;
 // real ll[2];
 // real theta_new;
 // real alpha_new;
 //real<lower=-0.75, upper=0.75> beta_new;

 // theta_new = normal_rng(THETA,sd_theta);
 // alpha_new = normal_rng(LAMBDA,sd_alpha);
 
 // Se_pooled = inv_logit(-(THETA - LAMBDA/2) / exp(beta/2));
 // Sp_pooled = inv_logit((THETA + LAMBDA/2) / exp(beta/2));
 
 Se_pooled = inv_logit(-(theta[1] - alpha[1]/2) / exp(beta/2));
 Sp_pooled = inv_logit((theta[1] + alpha[1]/2) / exp(beta/2));
 
 

}
