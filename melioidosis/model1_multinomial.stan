//
// ' Model 1'

data {
  
  int<lower=1> M; // Number of different tests
  int<lower=1> N; // Number of people in study
  int<lower=1> C; // Number of test combinations
  int<lower=0> y[C]; // frequency of each combination of test results
  
}


parameters {
  real<lower=0,upper=1> Se1; // sensitivity
  //real<lower=0,upper=1> Sp1; // specifcity
  real<lower=0,upper=1> Se2; // sensitivity
  real<lower=0,upper=1> Sp2; // specifcity
  real<lower=0,upper=1> Se3; // sensitivity
  real<lower=0,upper=1> Sp3; // specifcity
  real<lower=0,upper=1> Se4; // sensitivity
  real<lower=0,upper=1> Sp4; // specifcity
  real<lower=0,upper=1> Se5; // sensitivity
  real<lower=0,upper=1> Sp5; // specifcity
  real<lower=0,upper=1> prev; // prevalence
  
  real<lower=-(Se2*Se3),upper=fmin(Se2,Se3)-Se2*Se3> cov23;// covariance between tests two and three
  //real<lower=0,upper=1> cov23;// covariance between tests two and three

}

transformed parameters {
  //vector<lower=0,upper=1>[C] theta2;   // 
  simplex[C] theta2;   // 
  //real ub23;
  real Sp1;
  Sp1 = 1;
  
  //ub23 = fmin(Se2,Se3)-Se2*Se3;
  
  theta2[1] = prev*(Se1*Se2*Se3*Se4*Se5+cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 1,1,1,1,1
  theta2[2] = prev*(Se1*Se2*Se3*Se4*(1-Se5)+cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*(1-Sp4)*Sp5); // 1,1,1,1,0
  theta2[3] = prev*(Se1*Se2*Se3*(1-Se4)*Se5+cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*Sp4*(1-Sp5)); // 1,1,1,0,1
  theta2[4] = prev*(Se1*Se2*Se3*(1-Se4)*(1-Se5)+cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*(1-Sp3)*Sp4*Sp5); // 1,1,1,0,0
  theta2[5] = prev*(Se1*Se2*(1-Se3)*Se4*Se5-cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*(1-Sp4)*(1-Sp5)); // 1,1,0,1,1
  theta2[6] = prev*(Se1*Se2*(1-Se3)*Se4*(1-Se5)-cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*(1-Sp4)*Sp5); // 1,1,0,1,0
  theta2[7] = prev*(Se1*Se2*(1-Se3)*(1-Se4)*Se5-cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*Sp4*(1-Sp5)); // 1,1,0,0,1
  theta2[8] = prev*(Se1*Se2*(1-Se3)*(1-Se4)*(1-Se5)-cov23) + (1-prev)*((1-Sp1)*(1-Sp2)*Sp3*Sp4*Sp5); // 1,1,0,0,0
  theta2[9] = prev*(Se1*(1-Se2)*Se3*Se4*Se5-cov23) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 1,0,1,1,1
  theta2[10] = prev*(Se1*(1-Se2)*Se3*Se4*(1-Se5)-cov23) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*(1-Sp4)*Sp5); // 1,0,1,1,0
  theta2[11] = prev*(Se1*(1-Se2)*Se3*(1-Se4)*Se5-cov23) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*Sp4*(1-Sp5)); // 1,0,1,0,1
  theta2[12] = prev*(Se1*(1-Se2)*Se3*(1-Se4)*(1-Se5)-cov23) + (1-prev)*((1-Sp1)*Sp2*(1-Sp3)*Sp4*Sp5); // 1,0,1,0,0
  theta2[13] = prev*(Se1*(1-Se2)*(1-Se3)*Se4*Se5+cov23) + (1-prev)*((1-Sp1)*Sp2*Sp3*(1-Sp4)*(1-Sp5)); // 1,0,0,1,1
  theta2[14] = prev*(Se1*(1-Se2)*(1-Se3)*Se4*(1-Se5)+cov23) + (1-prev)*((1-Sp1)*Sp2*Sp3*(1-Sp4)*Sp5); // 1,0,0,1,0
  theta2[15] = prev*(Se1*(1-Se2)*(1-Se3)*(1-Se4)*Se5+cov23) + (1-prev)*((1-Sp1)*Sp2*Sp3*Sp4*(1-Sp5)); // 1,0,0,0,1
  theta2[16] = prev*(Se1*(1-Se2)*(1-Se3)*(1-Se4)*(1-Se5)+cov23) + (1-prev)*((1-Sp1)*Sp2*Sp3*Sp4*Sp5); // 1,0,0,0,0
  theta2[17] = prev*((1-Se1)*Se2*Se3*Se4*Se5+cov23) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 0,1,1,1,1
  theta2[18] = prev*((1-Se1)*Se2*Se3*Se4*(1-Se5)+cov23) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*(1-Sp4)*Sp5); // 0,1,1,1,0
  theta2[19] = prev*((1-Se1)*Se2*Se3*(1-Se4)*Se5+cov23) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*Sp4*(1-Sp5)); // 0,1,1,0,1
  theta2[20] = prev*((1-Se1)*Se2*Se3*(1-Se4)*(1-Se5)+cov23) + (1-prev)*(Sp1*(1-Sp2)*(1-Sp3)*Sp4*Sp5); // 0,1,1,0,0
  theta2[21] = prev*((1-Se1)*Se2*(1-Se3)*Se4*Se5-cov23) + (1-prev)*(Sp1*(1-Sp2)*Sp3*(1-Sp4)*(1-Sp5)); // 0,1,0,1,1
  theta2[22] = prev*((1-Se1)*Se2*(1-Se3)*Se4*(1-Se5)-cov23) + (1-prev)*(Sp1*(1-Sp2)*Sp3*(1-Sp4)*Sp5); // 0,1,0,1,0
  theta2[23] = prev*((1-Se1)*Se2*(1-Se3)*(1-Se4)*Se5-cov23) + (1-prev)*(Sp1*(1-Sp2)*Sp3*Sp4*(1-Sp5)); // 0,1,0,0,1
  theta2[24] = prev*((1-Se1)*Se2*(1-Se3)*(1-Se4)*(1-Se5)-cov23) + (1-prev)*(Sp1*(1-Sp2)*Sp3*Sp4*Sp5); // 0,1,0,0,0
  theta2[25] = prev*((1-Se1)*(1-Se2)*Se3*Se4*Se5-cov23) + (1-prev)*(Sp1*Sp2*(1-Sp3)*(1-Sp4)*(1-Sp5)); // 0,0,1,1,1
  theta2[26] = prev*((1-Se1)*(1-Se2)*Se3*Se4*(1-Se5)-cov23) + (1-prev)*(Sp1*Sp2*(1-Sp3)*(1-Sp4)*Sp5); // 0,0,1,1,0
  theta2[27] = prev*((1-Se1)*(1-Se2)*Se3*(1-Se4)*Se5-cov23) + (1-prev)*(Sp1*Sp2*(1-Sp3)*Sp4*(1-Sp5)); // 0,0,1,0,1
  theta2[28] = prev*((1-Se1)*(1-Se2)*Se3*(1-Se4)*(1-Se5)-cov23) + (1-prev)*(Sp1*Sp2*(1-Sp3)*Sp4*Sp5); // 0,0,1,0,0
  theta2[29] = prev*((1-Se1)*(1-Se2)*(1-Se3)*Se4*Se5+cov23) + (1-prev)*(Sp1*Sp2*Sp3*(1-Sp4)*(1-Sp5)); // 0,0,0,1,1
  theta2[30] = prev*((1-Se1)*(1-Se2)*(1-Se3)*Se4*(1-Se5)+cov23) + (1-prev)*(Sp1*Sp2*Sp3*(1-Sp4)*Sp5); // 0,0,0,1,0
  theta2[31] = prev*((1-Se1)*(1-Se2)*(1-Se3)*(1-Se4)*Se5+cov23) + (1-prev)*(Sp1*Sp2*Sp3*Sp4*(1-Sp5)); // 0,0,0,0,1
  theta2[32] = prev*((1-Se1)*(1-Se2)*(1-Se3)*(1-Se4)*(1-Se5)+cov23) + (1-prev)*(Sp1*Sp2*Sp3*Sp4*Sp5); // 0,0,0,0,0
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

  //cov23~uniform(0,ub23);
  
  y~multinomial(theta2);

}

generated quantities {
  
  int<lower=0> y_pred[C];
  y_pred = multinomial_rng(theta2, N);
  
}

