//
// Model 0 individual

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
  //vector<lower=1,upper=1>[32] profiles[M];
  
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
  
  // vector[N] RE;
  // real<lower=0,upper=5> b1; 
  // real<lower=0,upper=5> b2;
  // real<lower=0,upper=5> b3;
  // real<lower=0,upper=5> b4;
}

transformed parameters {
  
  simplex[2] theta; // prob infected or not infected
  vector[N] prob[M,2];   // probabilities of being TN or TP
  real Sp1;

  Sp1 = 1;
  
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
  prob[2,2] = rep_vector(a2, N);
  prob[3,1] = rep_vector(1-Sp3,N);
  prob[3,2] = rep_vector(a3, N);
  prob[4,1] = rep_vector(1-Sp4,N);
  prob[4,2] = rep_vector(a4, N);
  prob[5,1] = rep_vector(1-Sp5,N);
  prob[5,2] = rep_vector(a5, N);
  //  
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

  // RE~normal(0,1); 
  // b1~gamma(1,1);
  // b2~gamma(1,1);
  // b3~gamma(1,1);
  // b4~gamma(1,1);
  
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
  
  vector[N] log_lik;
  real ll[2];
  
  int<lower=0> y_pred[N,M];
  int<lower=0,upper=1> inf[N];
  
  real<lower=0,upper=1> p[N,M];

for(m in 1:M){
  Se_mean[m] = mean(prob[m,2,]);
  Sp_mean[m] = mean(1-prob[m,1,]);
}
 
  
for(n in 1:N){
  inf[n] = binomial_rng(1,theta[2]);
}  

for(n in 1:N){
  for(m in 1:M){
    p[n,m] = inf[n]*prob[m,2,n]+(1-inf[n])*(prob[m,1,n]);
  }
}

for(n in 1:N){
  for(m in 1:M){
    y_pred[n,m] = binomial_rng(1,p[n,m]);
  }
}  

// Likelihood for use in LOO-CV
  for(n in 1:N){
    for(k in 1:2){
      ll[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]) + binomial_lpmf(t4[n]| 1, prob[4,k,n]) + binomial_lpmf(t5[n]| 1, prob[5,k,n]);
    }

  log_lik[n] = log_sum_exp(ll);
  }

}

