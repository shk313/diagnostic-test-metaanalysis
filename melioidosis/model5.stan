//
// Model 5 (model 3 + model 4)

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=0,upper=1> t1[N]; // results of T1
  int<lower=0,upper=1> t2[N]; 
  int<lower=0,upper=1> t3[N];
  int<lower=0,upper=1> t4[N];
  int<lower=0,upper=1> t5[N];

}


parameters {
  real<lower=0,upper=1> a1; 
  //real<lower=0.99,upper=1> b1; 
  real<lower=0,upper=1> a2;
  real<lower=0,upper=1> b2;
  real<lower=0,upper=1> a3;
  real<lower=0,upper=1> b3;
  real<lower=0,upper=1> a4;
  real<lower=0,upper=1> b4;
  real<lower=0,upper=1> a5;
  real<lower=0,upper=1> b5;
  
  real<lower=0,upper=1> prev; 
  
  vector[N] RE;
  real<lower=0,upper=5> sd_re1; 
  real<lower=0,upper=5> sd_re2;

}

transformed parameters {
  
  simplex[2] theta; // prob infected or not infected
  vector[N] prob[M,2];   // probabilities of being TN or TP
  real b1; 
  
  b1 = 1; // culture specificty assumed perfect
  
  theta[1] = 1-prev;
  theta[2] = prev;
  
// vectorised version of above, faster
  prob[1,1] = rep_vector(1-b1,N); // Test 1, individual n, not infected
  prob[1,2] = rep_vector(a1,N); // not correalted with the other tests
  prob[2,1] = inv_logit(logit(1-b2)+sd_re1*RE);
  prob[2,2] = inv_logit(logit(a2)+sd_re2*RE);
  prob[3,1] = inv_logit(logit(1-b3)+sd_re1*RE);
  prob[3,2] = inv_logit(logit(a3)+sd_re2*RE);
  prob[4,1] = inv_logit(logit(1-b4)+sd_re1*RE);
  prob[4,2] = inv_logit(logit(a4)+sd_re2*RE);
  prob[5,1] = inv_logit(logit(1-b5)+sd_re1*RE);
  prob[5,2] = inv_logit(logit(a5)+sd_re2*RE);
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
  //b1~beta(10,1);
  b2~beta(1,1);
  b3~beta(1,1);
  b4~beta(1,1);
  b5~beta(1,1);
  prev~beta(1,1); 

  RE~normal(0,1); 
  sd_re1~gamma(1,1);
  sd_re2~gamma(1,1);


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
  
    // for loo-cv
  vector[N] log_lik;
  real ll[2];
  
  int<lower=0> y_pred[N,M];
  int<lower=0,upper=1> inf[N];

  real<lower=0,upper=1> p[N,M];

// mean se and sp
for(m in 1:M){
  Se_mean[m] = mean(prob[m,2,]);
  Sp_mean[m] = mean(1-prob[m,1,]);
}
  
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


//Likelihood for use in LOO-CV
for(n in 1:N){
  for(k in 1:2){
    ll[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]) + binomial_lpmf(t4[n]| 1, prob[4,k,n]) + binomial_lpmf(t5[n]| 1, prob[5,k,n]);
  }

log_lik[n] = log_sum_exp(ll);
}
  
}

