//
//' Final Model , Model3

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=0,upper=1> t1[N]; // results of T1
  int<lower=0,upper=1> t2[N]; 
  int<lower=0,upper=1> t3[N];
  int<lower=0,upper=1> t4[N];
  int<lower=0,upper=1> t5[N];
  // real alpha_1;
  // real beta_1;
  
}


parameters {
  real<lower=0,upper=1> a1; 
  //real<lower=0.99,upper=1> Sp1; 
  real<lower=0,upper=1> a2;
  real<lower=0,upper=1> Sp2;
  real<lower=0,upper=1> a3;
  real<lower=0,upper=1> Sp3;
  real<lower=0,upper=1> a4;
  real<lower=0,upper=1> Sp4;
  real<lower=0,upper=1> a5;
  real<lower=0,upper=1> Sp5;
  
  real<lower=0,upper=1> prev; 
  
  vector[N] RE;
  real<lower=0> tau;
  real<lower=0> sd_re[4]; 

}

transformed parameters {
  
  simplex[2] theta; // prob infected or not infected
  vector[N] prob[M,2];   // probabilities of being TN or TP
  real Sp1;
  // real logtau;
  // 
  // logtau = log(tau);
  
  Sp1 = 1;

  theta[1] = 1-prev;
  theta[2] = prev;
 
// correlation among infected people  
  // prob[1,1] = rep_vector(1-Sp1,N); // Test 1, individual n, not infected
  // prob[1,2] = rep_vector(a1,N); // not correalted with the other tests
  // prob[2,1] = rep_vector(1-Sp2,N);
  // prob[2,2] = inv_logit(logit(a2)+sd_re[1]*RE);
  // prob[3,1] = rep_vector(1-Sp3,N);
  // prob[3,2] = inv_logit(logit(a3)+sd_re[2]*RE);
  // prob[4,1] = rep_vector(1-Sp4,N);
  // prob[4,2] = inv_logit(logit(a4)+sd_re[3]*RE);
  // prob[5,1] = rep_vector(1-Sp5,N);
  // prob[5,2] = inv_logit(logit(a5)+sd_re[4]*RE);
  //  
  
// correlation among non infected people  
  prob[1,1] = rep_vector(1-Sp1,N); // Test 1, individual n, not infected
  prob[1,2] = rep_vector(a1,N); // not correalted with the other tests
  prob[2,1] = inv_logit(logit(1-Sp2)+sd_re[1]*RE);
  prob[2,2] = rep_vector(a2, N);
  prob[3,1] = inv_logit(logit(1-Sp3)+sd_re[2]*RE);
  prob[3,2] = rep_vector(a3, N);
  prob[4,1] = inv_logit(logit(1-Sp4)+sd_re[3]*RE);
  prob[4,2] = rep_vector(a4,N);
  prob[5,1] = inv_logit(logit(1-Sp5)+sd_re[4]*RE);
  prob[5,2] = rep_vector(a5,N);  
  
}


model {
  real ps[2];
  // priors
  a1~beta(1,1);
  a2~beta(1,1);
  a3~beta(1,1);
  a4~beta(1,1);
  a5~beta(1,1);
  //Sp1~beta(100,1);
  Sp2~beta(1,1);
  Sp3~beta(1,1);
  Sp4~beta(1,1);
  Sp5~beta(1,1);
  prev~beta(1,1); 

  RE~std_normal(); // implies RE ~normal(0,sd_re)
  tau ~ gamma(1,1); // 'base' precision
  sd_re~normal(tau,1); //'precision' parameters for each test


  
  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]) + binomial_lpmf(t4[n]| 1, prob[4,k,n]) + binomial_lpmf(t5[n]| 1, prob[5,k,n]);
    }

  target += log_sum_exp(ps);
  }


}

generated quantities {
  real Se_mean[M];
  real Sp_mean[M];
  
  // vector[N] log_lik;
  // real ll[2];
  
  // real ppv1;
  // real npv1;
  // real ppv2;
  // real npv2;
  // real ppv3;
  // real npv3;
  // real ppv4;
  // real npv4;
  // real ppv5;
  // real npv5;
  
  int<lower=0> y_pred[N,M];
  int<lower=0,upper=1> inf[N];

  real<lower=0,upper=1> p[N,M];

// mean se and sp
for(m in 1:M){
  Se_mean[m] = mean(prob[m,2,]);
  Sp_mean[m] = mean(1-prob[m,1,]);
}

  // ppv1 = Se_mean[1]*prev / (Se_mean[1]*prev+(1-Sp1)*(1-prev));
  // npv1 = Sp1*(1-prev) / (Sp1*(1-prev)+(1-Se_mean[1])*prev);
  // ppv2 = Se_mean[2]*prev / (Se_mean[2]*prev+(1-Sp2)*(1-prev));
  // npv2 = Sp2*(1-prev) / (Sp2*(1-prev)+(1-Se_mean[2])*prev);
  // ppv3 = Se_mean[3]*prev / (Se_mean[3]*prev+(1-Sp3)*(1-prev));
  // npv3 = Sp3*(1-prev) / (Sp3*(1-prev)+(1-Se_mean[3])*prev);
  // ppv4 = Se_mean[4]*prev / (Se_mean[4]*prev+(1-Sp4)*(1-prev));
  // npv4 = Sp4*(1-prev) / (Sp4*(1-prev)+(1-Se_mean[4])*prev);
  // ppv5 = Se_mean[5]*prev / (Se_mean[5]*prev+(1-Sp5)*(1-prev));
  // npv5 = Sp5*(1-prev) / (Sp5*(1-prev)+(1-Se_mean[5])*prev);
  // 
// prediction  
for(n in 1:N){
   inf[n] = binomial_rng(1,theta[2]);
 }

for(n in 1:N){
  for(m in 1:M){
      p[n,m] = (inf[n]*prob[m,2,n])+((1-inf[n])*(prob[m,1,n])); // probaility of person N being positive for test m
   }
}

for(n in 1:N){
  for(m in 1:M){
    y_pred[n,m] = binomial_rng(1, p[n,m]); // test result for person N on test M
  }
}

// Likelihood for use in LOO-CV
  // for(n in 1:N){
  //   for(k in 1:2){
  //     ll[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]) + binomial_lpmf(t4[n]| 1, prob[4,k,n]) + binomial_lpmf(t5[n]| 1, prob[5,k,n]);
  //   }
  // 
  // log_lik[n] = log_sum_exp(ll);
  // }
  
}

