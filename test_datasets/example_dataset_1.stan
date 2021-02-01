
data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=0,upper=1> t1[N]; // results of T1
  int<lower=0,upper=1> t2[N]; 
  int<lower=0,upper=1> t3[N];
  real alpha_1;
  real beta_1;
  
}


parameters {
  real<lower=0,upper=1> a1; 
  real<lower=0.6,upper=1> Sp1; 
  real<lower=0,upper=1> a2; 
  real<lower=0.6,upper=1> Sp2; 
  real<lower=0,upper=1> a3; 
  real<lower=0.6,upper=1> Sp3; 
  real<lower=0,upper=1> prev; 
  
  vector[N] RE;
  real<lower=0> b1; 
  real<lower=0> b2; 
  real<lower=0> b3; 
}

transformed parameters {
  
  real theta[2];
  vector[N] prob[M,2];   // probabilities of being TN or TP
  
  theta[1] = 1-prev;
  theta[2] = prev;


// probit
  prob[1,1] = rep_vector(1-Sp1,N);
  prob[1,2] = Phi_approx(a1+b1*RE);
  prob[2,1] = rep_vector(1-Sp2,N);
  prob[2,2] = Phi_approx(a2+b2*RE);
  prob[3,1] = rep_vector(1-Sp3,N);
  prob[3,2] = Phi_approx(a3+b3*RE); 
//

// // logit
//   prob[1,1] = rep_vector(1-Sp1,N);
//   prob[1,2] = inv_logit(logit(a1)+b1*RE);
//   prob[2,1] = rep_vector(1-Sp2,N);
//   prob[2,2] = inv_logit(logit(a2)+b2*RE);
//   prob[3,1] = rep_vector(1-Sp3,N);
//   prob[3,2] = inv_logit(logit(a3)+b3*RE); 
// //

}


model {
  real ps[2];
  // priors
  a1~beta(1,1);
  a2~beta(1,1);
  a3~beta(1,1);
  Sp1~beta(alpha_1,beta_1);
  Sp2~beta(alpha_1,beta_1);
  Sp3~beta(alpha_1,beta_1);
  prev~beta(1,1); 

  RE~normal(0,1); 
  b1~gamma(1,1);
  b2~gamma(1,1);
  b3~gamma(1,1);
  
  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]) + binomial_lpmf(t3[n]| 1, prob[3,k,n]);
    }

  target += log_sum_exp(ps);
  }


}

// generated quantities {
//   real Se_mean[M];
//   real Sp_mean[M];
// 
// for(m in 1:M){
//   Se_mean[m] = mean(prob[m,2,]);
//   Sp_mean[m] = mean(1-prob[m,1,]);
// }
// 
// }

