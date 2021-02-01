//
// ' Model 0'

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=1> C; // Number of test combinations
  int<lower=0> y[C]; // frequency of each combination of test results
  
}


parameters {
  real<lower=0,upper=1> Se1; // sensitivity
  // real<lower=0,upper=1> Sp1; // specifcity
  real<lower=0,upper=1> Se2; // sensitivity
  real<lower=0,upper=1> Sp2; // specifcity
  real<lower=0,upper=1> Se3; // sensitivity
  real<lower=0,upper=1> Sp3; // specifcity
  real<lower=0,upper=1> Se4; // sensitivity
  real<lower=0,upper=1> Sp4; // specifcity
  real<lower=0,upper=1> Se5; // sensitivity
  real<lower=0,upper=1> Sp5; // specifcity
  real<lower=0,upper=1> prev; // prevalence
  

}

transformed parameters {
  
  vector<lower=0,upper=1>[C] theta2;   // 
  real Sp1;
  Sp1 = 1;
  
  theta2[1] = prev*(Se1*Se2*Se3*Se4*Se5) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 1,1,1,1,1
  theta2[2] = prev*(Se1*Se2*Se3*Se4*(1-Se5)) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*(1-Sp4)*Sp5); // 1,1,1,1,0
  theta2[3] = prev*(Se1*Se2*Se3*(1-Se4)*Se5) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*Sp4*(1-Sp5)); // 1,1,1,0,1
  theta2[4] = prev*(Se1*Se2*Se3*(1-Se4)*(1-Se5)) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*Sp4*Sp5); // 1,1,1,0,0
  theta2[5] = prev*(Se1*Se2*(1-Se3)*Se4*Se5) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*(1-Sp4)*(1-Sp5)); // 1,1,0,1,1
  theta2[6] = prev*(Se1*Se2*(1-Se3)*Se4*(1-Se5)) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*(1-Sp4)*Sp5); // 1,1,0,1,0
  theta2[7] = prev*(Se1*Se2*(1-Se3)*(1-Se4)*Se5) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*Sp4*(1-Sp5)); // 1,1,0,0,1
  theta2[8] = prev*(Se1*Se2*(1-Se3)*(1-Se4)*(1-Se5)) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*Sp4*Sp5); // 1,1,0,0,0
  theta2[9] = prev*(Se1*(1-Se2)*Se3*Se4*Se5) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 1,0,1,1,1
  theta2[10] = prev*(Se1*(1-Se2)*Se3*Se4*(1-Se5)) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*(1-Sp4)*Sp5); // 1,0,1,1,0
  theta2[11] = prev*(Se1*(1-Se2)*Se3*(1-Se4)*Se5) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*Sp4*(1-Sp5)); // 1,0,1,0,1
  theta2[12] = prev*(Se1*(1-Se2)*Se3*(1-Se4)*(1-Se5)) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*Sp4*Sp5); // 1,0,1,0,0
  theta2[13] = prev*(Se1*(1-Se2)*(1-Se3)*Se4*Se5) + (1-prev)*((1-Sp1)*Sp2*Sp3*(1-Sp4)*(1-Sp5)); // 1,0,0,1,1
  theta2[14] = prev*(Se1*(1-Se2)*(1-Se3)*Se4*(1-Se5)) + (1-prev)*((1-Sp1)*Sp2*Sp3*(1-Sp4)*Sp5); // 1,0,0,1,0
  theta2[15] = prev*(Se1*(1-Se2)*(1-Se3)*(1-Se4)*Se5) + (1-prev)*((1-Sp1)*Sp2*Sp3*Sp4*(1-Sp5)); // 1,0,0,0,1
  theta2[16] = prev*(Se1*(1-Se2)*(1-Se3)*(1-Se4)*(1-Se5)) + (1-prev)*((1-Sp1)*Sp2*Sp3*Sp4*Sp5); // 1,0,0,0,0
  theta2[17] = prev*((1-Se1)*Se2*Se3*Se4*Se5) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 0,1,1,1,1
  theta2[18] = prev*((1-Se1)*Se2*Se3*Se4*(1-Se5)) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*(1-Sp4)*Sp5); // 0,1,1,1,0
  theta2[19] = prev*((1-Se1)*Se2*Se3*(1-Se4)*Se5) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*Sp4*(1-Sp5)); // 0,1,1,0,1
  theta2[20] = prev*((1-Se1)*Se2*Se3*(1-Se4)*(1-Se5)) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*Sp4*Sp5); // 0,1,1,0,0
  theta2[21] = prev*((1-Se1)*Se2*(1-Se3)*Se4*Se5) + (1-prev)*(Sp1*(1-Sp2)*Sp3*(1-Sp4)*(1-Sp5)); // 0,1,0,1,1
  theta2[22] = prev*((1-Se1)*Se2*(1-Se3)*Se4*(1-Se5)) + (1-prev)*(Sp1*(1-Sp2)*Sp3*(1-Sp4)*Sp5); // 0,1,0,1,0
  theta2[23] = prev*((1-Se1)*Se2*(1-Se3)*(1-Se4)*Se5) + (1-prev)*(Sp1*(1-Sp2)*Sp3*Sp4*(1-Sp5)); // 0,1,0,0,1
  theta2[24] = prev*((1-Se1)*Se2*(1-Se3)*(1-Se4)*(1-Se5)) + (1-prev)*(Sp1*(1-Sp2)*Sp3*Sp4*Sp5); // 0,1,0,0,0
  theta2[25] = prev*((1-Se1)*(1-Se2)*Se3*Se4*Se5) + (1-prev)*(Sp1*Sp2*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 0,0,1,1,1
  theta2[26] = prev*((1-Se1)*(1-Se2)*Se3*Se4*(1-Se5)) + (1-prev)*(Sp1*Sp2*(1-Sp3)*(1-Sp4)*Sp5); // 0,0,1,1,0
  theta2[27] = prev*((1-Se1)*(1-Se2)*Se3*(1-Se4)*Se5) + (1-prev)*(Sp1*Sp2*(1-Sp3)*Sp4*(1-Sp5)); // 0,0,1,0,1
  theta2[28] = prev*((1-Se1)*(1-Se2)*Se3*(1-Se4)*(1-Se5)) + (1-prev)*(Sp1*Sp2*(1-Sp3)*Sp4*Sp5); // 0,0,1,0,0
  theta2[29] = prev*((1-Se1)*(1-Se2)*(1-Se3)*Se4*Se5) + (1-prev)*(Sp1*Sp2*Sp3*(1-Sp4)*(1-Sp5)); // 0,0,0,1,1
  theta2[30] = prev*((1-Se1)*(1-Se2)*(1-Se3)*Se4*(1-Se5)) + (1-prev)*(Sp1*Sp2*Sp3*(1-Sp4)*Sp5); // 0,0,0,1,0
  theta2[31] = prev*((1-Se1)*(1-Se2)*(1-Se3)*(1-Se4)*Se5) + (1-prev)*(Sp1*Sp2*Sp3*Sp4*(1-Sp5)); // 0,0,0,0,1
  theta2[32] = prev*((1-Se1)*(1-Se2)*(1-Se3)*(1-Se4)*(1-Se5)) + (1-prev)*(Sp1*Sp2*Sp3*Sp4*Sp5); // 0,0,0,0,0
}


model {
  real ps[2];
  // priors
  Se1~beta(1,1);
  Se2~beta(1,1);
  Se3~beta(1,1);
  Se4~beta(1,1);
  Se5~beta(1,1);
  Sp2~beta(1,1);
  Sp3~beta(1,1);
  Sp4~beta(1,1);
  Sp5~beta(1,1);
  prev~beta(1,1); 

  y~multinomial(theta2);

}

generated quantities {
  
  int<lower=0> y_pred[C];
  real ppv1;
  real npv1;
  real ppv2;
  real npv2;
  real ppv3;
  real npv3;
  real ppv4;
  real npv4;
  real ppv5;
  real npv5;
  
  vector[N] log_lik;
  
  y_pred = multinomial_rng(theta2, N);
  
  for(n in 1:N){
    log_lik[n] = multinomial_lpmf(y | theta2);
  }
  
  ppv1 = Se1*prev / (Se1*prev+(1-Sp1)*(1-prev));
  npv1 = Sp1*(1-prev) / (Sp1*(1-prev)+(1-Se1)*prev);
  ppv2 = Se2*prev / (Se2*prev+(1-Sp2)*(1-prev));
  npv2 = Sp2*(1-prev) / (Sp2*(1-prev)+(1-Se2)*prev);
  ppv3 = Se3*prev / (Se3*prev+(1-Sp3)*(1-prev));
  npv3 = Sp3*(1-prev) / (Sp3*(1-prev)+(1-Se3)*prev);
  ppv4 = Se4*prev / (Se4*prev+(1-Sp4)*(1-prev));
  npv4 = Sp4*(1-prev) / (Sp4*(1-prev)+(1-Se4)*prev);
  ppv5 = Se5*prev / (Se5*prev+(1-Sp5)*(1-prev));
  npv5 = Sp5*(1-prev) / (Sp5*(1-prev)+(1-Se5)*prev);
  
}
