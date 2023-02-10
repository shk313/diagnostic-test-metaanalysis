//
// PCR and ELISA

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=0,upper=1> t1[N]; // results of T1
  int<lower=0,upper=1> t2[N]; 
  
}


parameters {
  //real<lower=0.5,upper=1> a11; // index test specificity 
  real<lower=0,upper=1> a11;
  //real<lower=0,upper=1> a12; // index test sensitivity
  //real<lower=0.5,upper=1> a21; // comparator test specificity
  real<lower=0,upper=1> a21; // comparator test specificity
  //real<lower=0,upper=1> a22; // comparator test sensitivity 

  real<lower=0,upper=1> prev; 
  
  vector[N] RE;
  real<lower=0> b1; 
  
  real<lower=1-inv_logit(logit(a11)-b1*2), upper=1> a12;
  real<lower=1-inv_logit(logit(a21)-b1*2), upper=1> a22;
  //real<lower=0> b2; 
  //real<lower=0> b3; 
}

transformed parameters {
  
  simplex[2] theta; // prob infected or not infected
  vector[N] prob[M,2];   // probabilities of being TN or TP
  
  theta[1] = 1-prev;
  theta[2] = prev;
 
for(n in 1:N){
  prob[1,1,n] = inv_logit(logit(1-a11)); // Test 1, individual n, not infected, FP
  prob[1,2,n] = inv_logit(logit(a12)+b1*RE[n]);  // TP
  prob[2,1,n] = inv_logit(logit(1-a21));
  prob[2,2,n] = inv_logit(logit(a22)+b1*RE[n]);
  }
  
  // vectorised version of above
  // prob[1,1] = rep_vector(inv_logit(logit(1-Sp1)),N); // Test 1, individual n, not infected
  // prob[1,2] = inv_logit(logit(a1)+b1*RE);
  // prob[2,1] = rep_vector(inv_logit(logit(1-Sp2)),N);
  // prob[2,2] = inv_logit(logit(a2)+b2*RE);
  // prob[3,1] = rep_vector(inv_logit(logit(1-Sp3)),N);
  // prob[3,2] = inv_logit(logit(a3)+b3*RE);
  // 
}


model {
  real ps[2];
  // priors
  a11~beta(1,1);
  a12~beta(5,1);
  a21~beta(1,1);
  a22~beta(5,1);
  prev~beta(1,1); 

  RE~normal(0,1); 
  b1~gamma(1,1);

  
  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) +  binomial_lpmf(t1[n]| 1, prob[1,k,n]) +  binomial_lpmf(t2[n]| 1, prob[2,k,n]);
    }

  target += log_sum_exp(ps);
  }


}


generated quantities {
 real Se_pooled;
 real Sp_pooled;
 real Se_Ref;
 real Sp_Ref;
 
 Se_pooled = mean(prob[1,2,]) ;
 Sp_pooled = mean(1-prob[1,1,]) ;
 Se_Ref = mean(prob[2,2,]) ;
 Sp_Ref = mean(1-prob[2,1,]) ;
 
}

