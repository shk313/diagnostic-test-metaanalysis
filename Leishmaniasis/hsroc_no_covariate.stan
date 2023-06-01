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
  // real THETA_region[R];
  // real LAMBDA_region[R];
  real<lower=-0.75, upper=0.75> beta; // common scale param
  real<lower=0> sd_theta; // changed to vector from real
  real<lower=0> sd_alpha;
  vector[J] theta_; // for reparameterisation of theta
  vector[J] alpha_; // for reparameterisation of alpha
  //vector[J] alpha_;
  real<lower=0,upper=1> prevalence[J];
  vector[N] N_RE; // e.g. infection intensity
  real<lower=0> sd_re[RT+1]; // sd for infection intensity between positives in each test


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
  real<lower=0,upper=1> Sp_Ref6;
  real<lower=inv_logit(logit(1-Sp_Ref6)),upper=1> Se_Ref6;
  real<lower=0,upper=1> Sp_Ref7;
  real<lower=inv_logit(logit(1-Sp_Ref7)),upper=1> Se_Ref7;
  real<lower=0,upper=1> Sp_Ref8;
  real<lower=inv_logit(logit(1-Sp_Ref8)),upper=1> Se_Ref8;


}

transformed parameters {
  
  vector<lower=0,upper=1>[N] prob[2,2]; // I want it to be N rows, two columns(for Se and Sp), with two shelves(one for each test)
  simplex[2] prev[J];
  vector<lower=0, upper=1>[l] Sp; // index test
  vector<lower=0, upper=1>[l] Se; // index test

  vector[J] theta;
  //vector[1] theta;
  vector[J] alpha;
  //vector[J] alpha;
  //real x;
  real y;
  // vector[l] T;
  // vector[l] L;
  

  theta = THETA + theta_*sd_theta; // implies theta~N(THETA,prec[1])
  alpha = LAMBDA + alpha_*sd_alpha; //implies alpha~N(LAMBDA,prec[2])

  
for(j in 1:J){

  
  prev[j,1] = 1-prevalence[j];
  prev[j,2] = prevalence[j];
  
  
  // T[j] = THETA_region[1]*regionA[j] + THETA_region[2]*regionB[j] + THETA_region[3]*regionC[j] + THETA_region[4]*regionD[j];
  // L[j] = LAMBDA_region[1]*regionA[j] + LAMBDA_region[2]*regionB[j] + LAMBDA_region[3]*regionC[j] + LAMBDA_region[4]*regionD[j];
  // 
  // theta[j] = T[j]+ theta_[j]* sd_theta[region[j]]; // implies theta~N(THETA,prec[1])
  // alpha[j] = L[j] + alpha_[j]* sd_alpha[region[j]]; 
}


for(j in 1:l){ 
  
   // Equation from Dendukuri et al
  Se[j] = inv_logit(-(theta[study[j]] - alpha[study[j]]/2)/exp(beta/2));
  Sp[j] = inv_logit((theta[study[j]] + alpha[study[j]]/2)/exp(beta/2));
}
  


 
// Conditional in/dependence, re-write this in a more efficient way
for(n in 1:N){

  y = ref_test[n]; // grab reference test id of person
  
  prob[2,1,n] = Se[ind_unique[n]];
  //prob[2,1,n] = inv_logit(logit(Se[ind_unique[n]])+N_RE[n]*sd_re[RT+1]);
  prob[1,1,n] = 1-Sp[ind_unique[n]];

 if(y == 1){
    prob[2,2,n] = Se_Ref1
    //prob[2,2,n] = inv_logit(logit(Se_Ref1)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref1;
    }
     if(y == 2){
    prob[2,2,n] = Se_Ref2;
    //prob[2,2,n] = inv_logit(logit(Se_Ref2)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref2;
    }
      if(y == 3){
    prob[2,2,n] = Se_Ref3;   
    //prob[2,2,n] = inv_logit(logit(Se_Ref3)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref3;
    }
      if(y == 4){
    prob[2,2,n] = Se_Ref4;
    //prob[2,2,n] = inv_logit(logit(Se_Ref4)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref4;
    }
    if(y == 5){
    prob[2,2,n] = Se_Ref5;
    //prob[2,2,n] = inv_logit(logit(Se_Ref5)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref5;
    }
    if(y == 6){
    prob[2,2,n] = Se_Ref6;  
    //prob[2,2,n] = inv_logit(logit(Se_Ref6)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref6;
    }
    if(y == 7){
    prob[2,2,n] = Se_Ref7; 
    //prob[2,2,n] = inv_logit(logit(Se_Ref7)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref7;
    }
        if(y == 8){
    prob[2,2,n] = Se_Ref8;     
    //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref8;
    }


 }


  
}

model {
  real lp[2]; // 2d vector with each entry corresponding
  // to the un-normalised posterior prob for diseased and not diseased
  
  
  THETA~normal(0,1);
  LAMBDA~normal(0,1.5); //0,2

  sd_re~normal(0,1);
  N_RE~std_normal();
  
  theta_ ~ std_normal(); // implies theta~normal(THETA, prec[1])
  alpha_ ~ std_normal(); // implies alpha~normal(LAMBDA, prec[1])

  sd_theta ~ normal(0,1);
  sd_alpha ~ normal(0,1);
  
  prevalence~beta(1,1);
  
  // LQ
  // Se_Ref1~beta(1,1); // Micro
  // Se_Ref2~beta(1,1); // ELISA
  // Se_Ref3~beta(5,1); // PCR
  // Se_Ref4~beta(1,1); // Leish Skin test
  // Se_Ref5~beta(5,1.5); // DAT
  // Se_Ref6~beta(1,1); // rk39
  // Se_Ref7~beta(1,1); // Latex
  // Se_Ref8~beta(1,1); // IFAT
  // // Se_Ref9~beta(1,1); // Katex
  // // 
  // Sp_Ref1~beta(50,1);
  // Sp_Ref2~beta(1,1);
  // Sp_Ref3~beta(10,1);
  // Sp_Ref4~beta(1,1);
  // Sp_Ref5~beta(5,1);
  // Sp_Ref6~beta(1,1);
  // Sp_Ref7~beta(1,1);
  // Sp_Ref8~beta(1,1);
  // // Sp_Ref9~beta(1,1);

  // FD
  // Se_Ref1~beta(1,1); // Parasitology
  // Se_Ref2~beta(5,1); // LQ DAT
  // //Se_Ref3~beta(5,1); // FAST DAT
  // Se_Ref3~beta(1,1); // ELISA
  // Se_Ref4~beta(5,1.5); // rK39
  // Se_Ref5~beta(1,1); // Leish Skin Test
  // Se_Ref6~beta(1,1); // IFAT
  // // 
  // Sp_Ref1~beta(50,1);
  // Sp_Ref2~beta(5,1);
  // //Sp_Ref3~beta(5,1);
  // Sp_Ref3~beta(1,1);
  // Sp_Ref4~beta(5,1);
  // Sp_Ref5~beta(1,1);
  // Sp_Ref6~beta(1,1);
  
  // Both
  Se_Ref1~beta(1,1); // Micro
  Se_Ref2~beta(1,1); // ELISA
  Se_Ref3~beta(5,1); // PCR
  Se_Ref4~beta(5,1.5); // DAT
  Se_Ref5~beta(1,1); // Leish Skin test
  Se_Ref6~beta(1,1); // rk39
  Se_Ref7~beta(1,1); // IFAT
  Se_Ref8~beta(1,1); // Latex /Katex
  //
  Sp_Ref1~beta(50,1);
  Sp_Ref2~beta(1,1);
  Sp_Ref3~beta(10,1);
  Sp_Ref4~beta(5,1);
  Sp_Ref5~beta(1,1);
  Sp_Ref6~beta(1,1);
  Sp_Ref7~beta(1,1);
  Sp_Ref8~beta(1,1);

  
  
  // Symptomatic only
  // Se_Ref1~beta(1,1); // rk39
  // Sp_Ref1~beta(1,1);
  // 
  // Se_Ref2~beta(5,1.5); //DAT
  // Sp_Ref2~beta(5,1);
  // 
  // Se_Ref3~beta(1,1); // ELISA
  // Sp_Ref3~beta(1,1);
  // 
  // Se_Ref4~beta(1,1); // Micro
  // Sp_Ref4~beta(50,1);
  // 
  // Se_Ref5~beta(1,1); // IFAT
  // Sp_Ref5~beta(1,1);
  // 
  // Se_Ref6~beta(5,1); // PCR
  // Sp_Ref6~beta(10,1);
  
    // FD DAT
  // Se_Ref1~beta(1,1); // Micro
  // Sp_Ref1~beta(50,1);
  // 
  // Se_Ref2~beta(5,1.5); //DAT
  // Sp_Ref2~beta(5,1);
  // 
  // Se_Ref3~beta(1,1); // ELISA
  // Sp_Ref3~beta(1,1);
  // 
  // Se_Ref4~beta(5,1); // rK39
  // Sp_Ref4~beta(5,1);
  // 
  // Se_Ref5~beta(1,1); // Leish Skin test
  // Sp_Ref5~beta(1,1);
  // 
  // Se_Ref6~beta(1,1); // IFAT
  // Sp_Ref6~beta(1,1);
  // 
     // LQ DAT
  // Se_Ref1~beta(1,1); // Micro
  // Sp_Ref1~beta(50,1);
  // 
  // Se_Ref2~beta(5,1); //ELISA
  // Sp_Ref2~beta(5,1);
  // 
  // Se_Ref3~beta(1,1); // PCR
  // Sp_Ref3~beta(5,1);
  // 
  // Se_Ref4~beta(1,1); // Leish Skin Test
  // Sp_Ref4~beta(1,1);
  // 
  // Se_Ref5~beta(5,1.5); // DAT
  // Sp_Ref5~beta(5,1);
  // 
  // Se_Ref6~beta(5,1); // rk39
  // Sp_Ref6~beta(5,1);
  // 
  // Se_Ref7~beta(1,1); // Latex/Katex
  // Sp_Ref7~beta(1,1);
  // 
  // Se_Ref8~beta(1,1); // IFAT
  // Sp_Ref8~beta(1,1);
  
  
  ### HIV only
  // Se_Ref1~beta(1,1); // MICRO
  // Sp_Ref1~beta(50,1);
  // 
  // Se_Ref2~beta(5,1); // rk39
  // Sp_Ref2~beta(5,1);
  // 
  // Se_Ref3~beta(1,1); // IFAT
  // Sp_Ref3~beta(1,1);
  // 
  // Se_Ref4~beta(1,1); // KAtex
  // Sp_Ref4~beta(1,1);
  
  
  
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
 real<lower=0,upper=1> Se_pred;
 real<lower=0,upper=1> Sp_pred;

 real alpha_new;
 real theta_new;
 
 // theta_new = normal_rng(THETA_region, sd_theta);
 // alpha_new = normal_rng(LAMBDA_region, sd_alpha);

 theta_new = normal_rng(THETA,sd_theta);
 alpha_new = normal_rng(LAMBDA,sd_alpha);

 Se_pooled = inv_logit(-(THETA - LAMBDA/2) / exp(beta/2));
 Sp_pooled = inv_logit((THETA + LAMBDA/2) / exp(beta/2));
 
 Se_pred= inv_logit(-(theta_new - alpha_new/2) / exp(beta/2));
 Sp_pred= inv_logit((theta_new + alpha_new/2) / exp(beta/2));

 
}
