//
//Three tets random effect model

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=0,upper=1> t1[N]; // results of T1
  int<lower=0,upper=1> t2[N]; 
  int<lower=0,upper=1> t3[N];
  int<lower=0,upper=1> t4[N];
  real alpha_1;
  real beta_1;
  real alpha_2;
  real beta_2;
  
}


parameters {
  real<lower=0,upper=1> a1; 
  real<lower=0.8,upper=1> Sp1; 
  real<lower=0,upper=1> a2;
  real<lower=0,upper=1> Sp2;
  real<lower=0,upper=1> a3;
  real<lower=0,upper=1> Sp3;
  real<lower=0,upper=1> a4;
  real<lower=0,upper=1> Sp4;
  
  real<lower=0,upper=1> prev; 
  
  vector[N] RE;
  real<lower=0> tau;
  real<lower=0> sd_re[3];
}

transformed parameters {
  
  simplex[2] theta; // prob infected or not infected
  vector[N] prob[M,2];   // probabilities of being TN or TP
  
  theta[1] = 1-prev;
  theta[2] = prev;
 
// for(n in 1:N){
//   prob[1,1,n] = inv_logit(logit(1-Sp1)); // Test 1, individual n, not infected, FP
//   prob[1,2,n] = inv_logit(logit(a1));  // TP
//   prob[2,1,n] = inv_logit(logit(1-Sp2));
//   prob[2,2,n] = inv_logit(logit(a2)+b1*RE[n]);
//   prob[3,1,n] = inv_logit(logit(1-Sp3));
//   prob[3,2,n] = inv_logit(logit(a3)+b2*RE[n]);
//   prob[4,1,n] = inv_logit(logit(1-Sp4));
//   prob[4,2,n] = inv_logit(logit(a4)+b3*RE[n]); 
//   }
  
// vectorised version of above, faster
  prob[1,1] = rep_vector(1-Sp1,N); // Test 1, individual n, not infected
  prob[1,2] = rep_vector(a1,N); // not correalted with the other tests
  prob[2,1] = rep_vector(1-Sp2,N);
  prob[2,2] = inv_logit(logit(a2)+sd_re[1]*RE);
  prob[3,1] = rep_vector(1-Sp3,N);
  prob[3,2] = inv_logit(logit(a3)+sd_re[2]*RE);
  prob[4,1] = rep_vector(1-Sp4,N);
  prob[4,2] = inv_logit(logit(a4)+sd_re[3]*RE);
  //  
}


model {
  real ps[2];
  // priors
  a1~beta(1,1);
  a2~beta(1,1);
  a3~beta(1,1);
  a4~beta(1,1);
  Sp1~beta(alpha_2,beta_2);
  Sp2~beta(alpha_1,beta_1);
  Sp3~beta(alpha_1,beta_1);
  Sp4~beta(alpha_1,beta_1);
  prev~beta(1,1); 

  RE~normal(0,1); 
  tau~gamma(1,1); // 'base' precision
  sd_re~normal(tau,1); //'precision' parameters for each test
  
  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]) + binomial_lpmf(t4[n]| 1, prob[4,k,n]);
    }

  target += log_sum_exp(ps);
  }


}

generated quantities {
  real Se_mean[M];
  real Sp_mean[M];
  // real Se_all_mean;
  // real Sp_all_mean;

for(m in 1:M){
  Se_mean[m] = mean(prob[m,2,]);
  Sp_mean[m] = mean(1-prob[m,1,]);
}

}


