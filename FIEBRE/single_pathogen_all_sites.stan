//
// This Stan program defines a simple model, with a


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> Nstrata;  //number of strata to consider for regression
  int<lower=0> N_Cases;
  int<lower=0> J; // Number of possible causes/pathogens
  int<lower=0> M; // Number of tests

  vector[N] I1; //indicator
  vector[N] I2; //indicator

  int t1[N];
  int t2[N];

  // int Region[N]; // Laos or African site
  int Strata[N]; // which strata does each individual belong to
  
}


parameters {

  
  real<lower=0,upper=1.0> Sp1[Nstrata];
  real<lower=0,upper=1.0> Sp2[Nstrata];

  real<lower=inv_logit(logit(1-min(Sp1))),upper=1.0> Se1;
  real<lower=inv_logit(logit(1-Sp2)),upper=1.0> Se2;

  simplex[J] aeti[Nstrata];// simple of length number of causes (here a single cause or other) plus i want it for 4 sites
  // 
  // vector[N] REp; // e.g. infection intensity
  // real<lower=0> b;
  
}

transformed parameters {
  /// Probability of a positive test for each individual that is a case across both tests
  vector<lower=0,upper=1>[N] prob[2,M]; 

  for(n in 1:N){
    
    prob[1,1,n] = 1-Sp1[Strata[n]];
    //prob[1,1,n] = 1-Sp1;
    //prob[2,1,n] = inv_logit(logit(Se1) + REp[n]*b);
    //prob[2,1,n] = Se1[Region[n]];
    prob[2,1,n] = Se1;
    prob[1,2,n] = 1-Sp2[Strata[n]];
    //prob[2,2,n] = inv_logit(logit(Se2) + REp[n]*b);
    prob[2,2,n] = Se2;
    
  }

  

}


model {
  real lp[J]; // 2d vector with each entry corresponding
  //real lp2;// controls
  // to the un-normalised posterior prob for each cause
  // priors
  aeti ~ dirichlet(rep_vector(1.0,J)); //Aetiology priors

  Se1~beta(34.68921,9.609493); // malaria RDT
  //Se1~beta(71.08486,3.242059);// 90-99% resp pathogens
  //Se1~beta(13.48064, 1.630844);//70-99% resp pathogens
  // Se1~beta(13.48064, 1.630844); //brucella capt 90-99%
  // Se2~beta(7.124981, 1.329325); //brucella mag acute 55-99%
  // Se3~beta(71.08486,3.242059); // brucella mag paired 90-99%
  // Se1~beta(12.21711, 1.577343); //dengue pcr 68-99%
  // //Se2~beta(1.261544, 5.967161); //dengue acute IgM 0-50%
  // Se2~beta(2.551923, 3.419105); //dengue acute IgM 10-80%
  // //Se3~beta(71.08486,3.242059); // dengue IgG paired 90-99%
  // Se3~beta(13.48064, 1.630844); // dengue IgG paired 70-99%
  // Se1~beta(7.124981, 1.329325); // Lepto ELISA 55-99
  // Se2~beta(1.261544, 5.967161); // Lepto MAT acute 0-50%
  // Se3~beta(13.48064, 1.630844); // Lepto MAT paired 70-99%
  // Se4~beta(13.48064, 1.630844); // Lepto PCR 55-99%
  // Se1~beta(7.124981, 1.329325); // Lepto ELISA 55-99
  // Se2~beta(1.261544, 5.967161); // Lepto MAT acute 0-50%
  // Se3~beta(24.75475, 2.037352); // Lepto MAT paired 80-99
  // Se4~beta(24.75475, 2.037352); // Lepto PCR 80-99
  // 
  //Se1~beta(36.03087,54.53019); // burkholderia pseudo, Typhoidal Sal
  //Se1~beta(31.4709,28.38045); // NTS
  //Se1~beta(16.86126,12.19816); // 40-75%
  //Se1~beta(36.03087, 54.53019); //30-50%
  //Se1~beta(12.60666, 24.26818); // 20-50%
  // Se2~beta(70.13171,95.22233); //35-50%
  //Se1~beta(15.52342,23.78298); //25-55%
  //Se1~beta(15.69303,14.10203); //35-70%
  //Se1~beta(9.875229,7.020875); //35-80%
  Se2~beta(1,1);
  // // //malaria microscopy
  // Se3~beta(0.45,22.9);
  //Se3[2]~beta(1,1);
  //Se4~beta(5,1);// Luminex Influenza
  //Se1~beta(62.86,0.52); // Luminex RSV
  //Se1~beta(3.54,0.39); // BrucellaCApt
  //Se1~beta(214.34,70.8); // PCR
  // Se2~beta(2.5,7.4); // Culture
  // Se3~beta(162.15,27.82);
  // Se4~beta(162.15,27.82);
  // Se5~beta(194.04,47.8);
  
  Sp1~beta(1,1);//malaria RDT
  //Sp1~beta(87.41117,1.404222); // blood cultures 95-99.9%
  // Sp2~beta(15.56221,0.8936577); //80-99.9%
  // Sp3~beta(71.08486,3.242059); // 90-99% 
  // Sp4~beta(1,1);
  // Sp4[1]~beta(20,1);
  // Sp4[2]~beta(1,1);
  // Sp4[3]~beta(1,1);
  // Sp4[4]~beta(20,1);
  Sp2[1]~beta(1,1); // Microscopy
  Sp2[2]~beta(1,1); // MLW update when chaging between 8 and 12 strata
  Sp2[3]~beta(1,1);
  Sp2[4]~beta(50,0.5);
  Sp2[5]~beta(50,0.5);
  Sp2[6]~beta(50,0.5);
  Sp2[7]~beta(1,1);
  Sp2[8]~beta(1,1);
  // Sp2[9]~beta(1,1);
  // Sp2[10]~beta(1,1);
  // Sp2[11]~beta(1,1);
  // Sp2[12]~beta(1,1);
  
  // malaria microscopy
  // Sp3~beta(1,1);
  // Sp3[2]~beta(1,1);
  // Sp3[3]~beta(1,1);
  // Sp3[4]~beta(50,1);// Luminex Influenza
  //Sp1~beta(134.42,0.58); // Luminex RSV
  //Sp1~beta(1.93,0.34); // BrucellaCapt
  // Sp1~beta(5.2,0.5);
  // Sp2~beta(60.8,0.74);
  // Sp3~beta(224.65,95.74);
  // Sp4~beta(214.34,70.8);
  // Sp5~beta(194.04,47.8);
  // 
  // REp~normal(0,1);
  // b~gamma(1,1);
  // 
  

  for(n in 1:N_Cases){
     for(j in 1:2){

       lp[j] = log(aeti[Strata[n], j]);
       //lp[j] = log(aeti[j]);
     } 
       
       lp[1] = lp[1] + binomial_lpmf(t1[n] | 1, prob[2,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[2,2,n])*I2[n];// + binomial_lpmf(t3[n] | 1, prob[2,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[2,4,n])*I4[n];
       lp[2] = lp[2] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n];// + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n];

    target += log_sum_exp(lp);
    //target += (lp);
  }

    // Controls
    for(n in (N_Cases+1):N){

        lp2 = binomial_lpmf(t1[n]|1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n]|1, prob[1,2,n])*I2[n];// + binomial_lpmf(t3[n]|1,prob[1,3,n])*I3[n] + binomial_lpmf(t4[n]|1,prob[1,4,n])*I4[n];
   //+ binomial_lpmf(y[n,2] | 1, prob[1,2,n])+ binomial_lpmf(y[n,3] | 1, prob[1,3,n])+ binomial_lpmf(y[n,4] | 1, prob[1,4,n])+ binomial_lpmf(y[n,5]|1,prob[1,5,n]);
        target += (lp2);
}
  
}
