//
// Code for the full model.


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // Total number of individuals (ordered all cases first followed by all controls
  int<lower=0> Nstrata;  //number of strata to consider for regression
  int<lower=0> N_Cases; //Number of cases
  int<lower=0> J; // Number of possible causes/pathogens
  int<lower=0> M; // Number of tests
  vector[N] I1; // Indicator for each test
  vector[N] I2; 
  vector[N] I3; 
  vector[N] I4;
  vector[N] I5;
  vector[N] I6;
  vector[N] I7; 
  vector[N] I8;
  vector[N] I9;
  vector[N] I10;
  vector[N] I11;
  vector[N] I12; 
  vector[N] I13; 
  vector[N] I14;
  vector[N] I15;
  vector[N] I16;
  vector[N] I17; 
  vector[N] I18;
  vector[N] I19;
  vector[N] I20;
  vector[N] I21; 
  vector[N] I22; 
  vector[N] I23; 
  vector[N] I24;
  vector[N] I25;
  vector[N] I26;
  vector[N] I27; 
  vector[N] I28;
  vector[N] I29;
  vector[N] I30;
  vector[N] I31; 
  vector[N] I32; 
  vector[N] I33; 
  vector[N] I34;
  vector[N] I35;
  vector[N] I36;
  vector[N] I37; 

  array[N] int t1; //results for each test
  array[N] int t2;
  array[N] int t3;
  array[N] int t4;
  array[N] int t5;
  array[N] int t6;
  array[N] int t7;
  array[N] int t8;
  array[N] int t9;
  array[N] int t10;
  array[N] int t11;
  array[N] int t12;
  array[N] int t13;
  array[N] int t14;
  array[N] int t15;
  array[N] int t16;
  array[N] int t17;
  array[N] int t18;
  array[N] int t19;
  array[N] int t20;
  array[N] int t21;
  array[N] int t22;
  array[N] int t23;
  array[N] int t24;
  array[N] int t25;
  array[N] int t26;
  array[N] int t27;
  array[N] int t28;
  array[N] int t29;
  array[N] int t30;
  array[N] int t31;
  array[N] int t32;
  array[N] int t33;
  array[N] int t34;
  array[N] int t35;
  array[N] int t36;
  array[N] int t37;

  array[N] int Strata; // which strata (age/site) does each individual belong to
  
}


parameters {
  
  array[Nstrata] real<lower=0,upper=1.0> Sp1;
  array[Nstrata] real<lower=0,upper=1.0> Sp2;
  array[Nstrata] real<lower=0,upper=1.0> Sp3;
  array[Nstrata] real<lower=0,upper=1.0> Sp4;
  array[Nstrata] real<lower=0,upper=1.0> Sp5;
  array[Nstrata] real<lower=0,upper=1.0> Sp6;
  array[Nstrata] real<lower=0,upper=1.0> Sp7;
  array[Nstrata] real<lower=0,upper=1.0> Sp8;
  array[Nstrata] real<lower=0,upper=1.0> Sp9;
  array[Nstrata] real<lower=0,upper=1.0> Sp10;
  array[Nstrata] real<lower=0,upper=1.0> Sp11;
  array[Nstrata] real<lower=0,upper=1.0> Sp12;
  array[Nstrata] real<lower=0,upper=1.0> Sp13;
  array[Nstrata] real<lower=0,upper=1.0> Sp14;
  array[Nstrata] real<lower=0,upper=1.0> Sp15;
  array[Nstrata] real<lower=0,upper=1.0> Sp16;
  array[Nstrata] real<lower=0,upper=1.0> Sp17;
  array[Nstrata] real<lower=0,upper=1.0> Sp18;
  array[Nstrata] real<lower=0,upper=1.0> Sp19;
  array[Nstrata] real<lower=0,upper=1.0> Sp20;
  array[Nstrata] real<lower=0,upper=1.0> Sp21;
  array[Nstrata] real<lower=0,upper=1.0> Sp22;
  array[Nstrata] real<lower=0,upper=1.0> Sp23;
  array[Nstrata] real<lower=0,upper=1.0> Sp24;
  array[Nstrata] real<lower=0,upper=1.0> Sp25;
  array[Nstrata] real<lower=0,upper=1.0> Sp26;
  array[Nstrata] real<lower=0,upper=1.0> Sp27;
  array[Nstrata] real<lower=0,upper=1.0> Sp28;
  array[Nstrata] real<lower=0,upper=1.0> Sp29;
  array[Nstrata] real<lower=0,upper=1.0> Sp30;
  array[Nstrata] real<lower=0,upper=1.0> Sp31;
  array[Nstrata] real<lower=0,upper=1.0> Sp32;
  array[Nstrata] real<lower=0,upper=1.0> Sp33;
  array[Nstrata] real<lower=0,upper=1.0> Sp34;
  array[Nstrata] real<lower=0,upper=1.0> Sp35;
  array[Nstrata] real<lower=0,upper=1.0> Sp36;
  array[Nstrata] real<lower=0,upper=1.0> Sp37;

  real<lower=inv_logit(logit(1-min(Sp1))),upper=1.0> Se1;
  real<lower=inv_logit(logit(1-min(Sp2))),upper=1.0> Se2;
  real<lower=inv_logit(logit(1-min(Sp3))),upper=1.0> Se3;
  real<lower=inv_logit(logit(1-min(Sp4))),upper=1.0> Se4;
  real<lower=inv_logit(logit(1-min(Sp5))),upper=1.0> Se5;
  real<lower=inv_logit(logit(1-min(Sp6))),upper=1.0> Se6;
  real<lower=inv_logit(logit(1-min(Sp7))),upper=1.0> Se7;
  real<lower=inv_logit(logit(1-min(Sp8))),upper=1.0> Se8;
  real<lower=inv_logit(logit(1-min(Sp9))),upper=1.0> Se9;
  real<lower=inv_logit(logit(1-min(Sp10))),upper=1.0> Se10;
  real<lower=inv_logit(logit(1-min(Sp11))),upper=1.0> Se11;
  real<lower=inv_logit(logit(1-min(Sp12))),upper=1.0> Se12;
  real<lower=inv_logit(logit(1-min(Sp13))),upper=1.0> Se13;
  real<lower=inv_logit(logit(1-min(Sp14))),upper=1.0> Se14;
  real<lower=inv_logit(logit(1-min(Sp15))),upper=1.0> Se15;
  real<lower=inv_logit(logit(1-min(Sp16))),upper=1.0> Se16;
  real<lower=inv_logit(logit(1-min(Sp17))),upper=1.0> Se17;
  real<lower=0,upper=1.0> Se18;
  real<lower=inv_logit(logit(1-min(Sp19))),upper=1.0> Se19;
  real<lower=inv_logit(logit(1-min(Sp20))),upper=1.0> Se20;
  real<lower=0,upper=1> Se21;
  real<lower=inv_logit(logit(1-min(Sp22))),upper=1.0> Se22;
  real<lower=inv_logit(logit(1-min(Sp23))),upper=1.0> Se23;
  real<lower=0,upper=1> Se24;
  real<lower=inv_logit(logit(1-min(Sp25))),upper=1.0> Se25;
  real<lower=inv_logit(logit(1-min(Sp26))),upper=1.0> Se26;
  real<lower=0,upper=1> Se27;
  real<lower=inv_logit(logit(1-min(Sp28))),upper=1.0> Se28;
  real<lower=inv_logit(logit(1-min(Sp29))),upper=1.0> Se29;
  real<lower=0,upper=1.0> Se30;
  real<lower=inv_logit(logit(1-min(Sp31))),upper=1.0> Se31;
  real<lower=inv_logit(logit(1-min(Sp32))),upper=1.0> Se32;
  real<lower=inv_logit(logit(1-min(Sp33))),upper=1.0> Se33;
  real<lower=inv_logit(logit(1-min(Sp34))),upper=1.0> Se34;
  real<lower=inv_logit(logit(1-min(Sp35))),upper=1.0> Se35;
  real<lower=inv_logit(logit(1-min(Sp36))),upper=1.0> Se36;
  real<lower=inv_logit(logit(1-min(Sp37))),upper=1.0> Se37;
 
  array[Nstrata] simplex[J] aeti;
  // vector[N] REp; // e.g. infection intensity
  // real<lower=0> b;
  
}

transformed parameters {
  /// Probability of a positive test for each individual that is a case across both tests
  array[2,M] vector<lower=0,upper=1>[N] prob;
  
  for(n in 1:N){
        
    prob[1,1,n] = 1-Sp1[Strata[n]];
    //prob[1,1,n] = 1-Sp1;
    //prob[2,1,n] = inv_logit(logit(Se1) + REp[n]*b);
    prob[2,1,n] = Se1;
    //prob[1,2,n] = 1-Sp2;
    prob[1,2,n] = 1-Sp2[Strata[n]];
    // //prob[2,2,n] = inv_logit(logit(Se2) + REp[n]*b);
    prob[2,2,n] = Se2;
    
    prob[1,3,n] = 1-Sp3[Strata[n]];
    prob[2,3,n] = Se3;
    prob[1,4,n] = 1-Sp4[Strata[n]];
    prob[2,4,n] = Se4;
    prob[1,5,n] = 1-Sp5[Strata[n]];
    prob[2,5,n] = Se5;
    prob[1,6,n] = 1-Sp6[Strata[n]];
    prob[2,6,n] = Se6;
    prob[1,7,n] = 1-Sp7[Strata[n]];
    prob[2,7,n] = Se7;
    prob[1,8,n] = 1-Sp8[Strata[n]];
    prob[2,8,n] = Se8;
    prob[1,9,n] = 1-Sp9[Strata[n]];
    prob[2,9,n] = Se9;
    prob[1,10,n] = 1-Sp10[Strata[n]];
    prob[2,10,n] = Se10;
    prob[1,11,n] = 1-Sp11[Strata[n]];
    prob[2,11,n] = Se11;
    prob[1,12,n] = 1-Sp12[Strata[n]];
    prob[2,12,n] = Se12;
    prob[1,13,n] = 1-Sp13[Strata[n]];
    prob[2,13,n] = Se13;
    prob[1,14,n] = 1-Sp14[Strata[n]];
    prob[2,14,n] = Se14;
    prob[1,15,n] = 1-Sp15[Strata[n]];
    prob[2,15,n] = Se15;
    prob[1,16,n] = 1-Sp16[Strata[n]];
    prob[2,16,n] = Se16;
    prob[1,17,n] = 1-Sp17[Strata[n]];
    prob[2,17,n] = Se17;
    prob[1,18,n] = 1-Sp18[Strata[n]];
    prob[2,18,n] = Se18;
    prob[1,19,n] = 1-Sp19[Strata[n]];
    prob[2,19,n] = Se19;
    prob[1,20,n] = 1-Sp20[Strata[n]];
    prob[2,20,n] = Se20;
    prob[1,21,n] = 1-Sp21[Strata[n]];
    prob[2,21,n] = Se21;
    prob[1,22,n] = 1-Sp22[Strata[n]];
    prob[2,22,n] = Se22;
    prob[1,23,n] = 1-Sp23[Strata[n]];
    prob[2,23,n] = Se23;
    prob[1,24,n] = 1-Sp24[Strata[n]];
    prob[2,24,n] = Se24;
    prob[1,25,n] = 1-Sp25[Strata[n]];
    prob[2,25,n] = Se25;
    prob[1,26,n] = 1-Sp26[Strata[n]];
    prob[2,26,n] = Se26;
    prob[1,27,n] = 1-Sp27[Strata[n]];
    prob[2,27,n] = Se27;
    prob[1,28,n] = 1-Sp28[Strata[n]];
    prob[2,28,n] = Se28;
    prob[1,29,n] = 1-Sp29[Strata[n]];
    prob[2,29,n] = Se29;
    prob[1,30,n] = 1-Sp30[Strata[n]];
    prob[2,30,n] = Se30;
    prob[1,31,n] = 1-Sp31[Strata[n]];
    prob[2,31,n] = Se31;
    prob[1,32,n] = 1-Sp32[Strata[n]];
    prob[2,32,n] = Se32;
    prob[1,33,n] = 1-Sp33[Strata[n]];
    prob[2,33,n] = Se33;
    prob[1,34,n] = 1-Sp34[Strata[n]];
    prob[2,34,n] = Se34;
    prob[1,35,n] = 1-Sp35[Strata[n]];
    prob[2,35,n] = Se35;
    prob[1,36,n] = 1-Sp36[Strata[n]];
    prob[2,36,n] = Se36;
    prob[1,37,n] = 1-Sp37[Strata[n]];
    prob[2,37,n] = Se37;
  }

}


model {
  array[J] real lp;

  real lp2;// controls
  // to the un-normalised posterior prob for each cause
  // priors

  aeti ~ dirichlet(rep_vector(1.0,J)); //Aetiology priors

// Malaria test priors 
  Se1~beta(34.68921,9.609493); // RDT 65-2.89.0%
  Se2~beta(1,1); // Microscopy

  Sp1~beta(1,1);//malaria RDT
  Sp2[1]~beta(1,1);// microscopy
  Sp2[2]~beta(1,1);
  Sp2[3]~beta(50,0.5); // microscopy malawi
  Sp2[4]~beta(50,0.5); // microscopy malawi
  Sp2[5]~beta(1,1);
  Sp2[6]~beta(1,1);
  Sp2[7]~beta(1,1);
  Sp2[8]~beta(1,1);
  
// 11 Respiratory pathogens
  Se3~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se4~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se5~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se6~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se7~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se8~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se9~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se10~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se11~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se12~beta(71.08486,3.242059); //Luminex for all resp 90-99%
  Se13~beta(71.08486,3.242059); //Luminex for all resp 90-99%

  Sp3~beta(1,1); //Luminex for all resp 90-99%
  Sp4~beta(1,1); //Luminex for all resp 90-99%
  Sp5~beta(1,1); //Luminex for all resp 90-99%
  Sp6~beta(1,1); //Luminex for all resp 90-99%
  Sp7~beta(1,1); //Luminex for all resp 90-99%
  Sp8~beta(1,1); //Luminex for all resp 90-99%
  Sp9~beta(1,1); //Luminex for all resp 90-99%
  Sp10~beta(1,1); //Luminex for all resp 90-99%
  Sp11~beta(1,1); //Luminex for all resp 90-99%
  Sp12~beta(1,1); //Luminex for all resp 90-99%
  Sp13~beta(1,1);
  
// Brucella
  Se14~beta(13.48064, 1.630844); //brucella capt 90-99%
  Se15~beta(7.124981, 1.329325); //brucella mag acute 55-99%
  Se16~beta(71.08486,3.242059); // brucella mag paired 90-99%
  
  Sp14~beta(1,1);
  Sp15~beta(1,1);
  Sp16~beta(71.08486,3.242059); // 90-99% 
  
  // Dengue
  Se17~beta(12.21711, 1.577343); //dengue pcr 68-99%
  Se18~beta(1.261544, 5.967161); //dengue acute IgM 0-50%
  Se19~beta(71.08486,3.242059); // dengue IgG paired 90-99%
  
  Sp17~beta(1,1);
  Sp18~beta(1,1);
  Sp19~beta(71.08486,3.242059);
  
  // Zika
  Se20~beta(12.21711, 1.577343); //dengue pcr 68-99%
  Se21~beta(1.261544, 5.967161); //dengue acute IgM 0-50%
  Se22~beta(71.08486,3.242059); // dengue IgG paired 90-99%
  
  Sp20~beta(1,1);
  Sp21~beta(1,1);
  Sp22~beta(71.08486,3.242059);
  
  // CHIK
  Se23~beta(12.21711, 1.577343); //dengue pcr 68-99%
  Se24~beta(1.261544, 5.967161); //dengue acute IgM 0-50%
  Se25~beta(71.08486,3.242059); // dengue IgG paired 90-99%
  
  Sp23~beta(1,1);
  Sp24~beta(1,1);
  Sp25~beta(71.08486,3.242059);
    
  // Lepto
  Se26~beta(7.124981, 1.329325); // Lepto ELISA 55-99
  Se27~beta(1.261544, 5.967161); // Lepto MAT acute 0-50%
  Se28~beta(24.75475, 2.037352); // Lepto MAT paired 80-99
  Se29~beta(24.75475, 2.037352); // Lepto PCR 80-99
  
  Sp26~beta(1,1);
  Sp27~beta(1,1);
  Sp28~beta(71.08486,3.242059);
  Sp29~beta(1,1);
  
  // 8 Blood cultures
  Se30~beta(36.03087,54.53019); // Burkholderia 30-50%
  Se31~beta(36.03087,54.53019); // TS 30-50%
  Se32~beta(31.4709,28.38045); // NTS 40-65%
  Se33~beta(16.86126,12.19816); // Staph aure 40-75%
  Se34~beta(16.86126,12.19816); // E.coli40-75%
  Se35~beta(16.86126,12.19816); // Other Entero 40-75$
  Se36~beta(16.86126,12.19816); // Kelb pneu 40-75%
  Se37~beta(36.03087,54.53019); // Strep penu 30-50

 Sp30~beta(87.41117, 1.404222);
 Sp31~beta(87.41117, 1.404222);
 Sp32~beta(87.41117, 1.404222);
 Sp33~beta(87.41117, 1.404222);
 Sp34~beta(87.41117, 1.404222);
 Sp35~beta(87.41117, 1.404222);
 Sp36~beta(87.41117, 1.404222);
 Sp37~beta(87.41117, 1.404222);
  
  // REp~normal(0,1);
  // b~gamma(1,1);
  

  for(n in 1:N_Cases){
     for(j in 1:J){

       lp[j] = log(aeti[Strata[n], j]);
       //lp[j] = log(aeti[j]);
     }
       
       // malaria
       lp[1] = lp[1] + binomial_lpmf(t1[n] | 1, prob[2,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[2,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
       // Influenza A
        lp[2] = lp[2] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[2,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
         // Influenza B
        lp[3] = lp[3] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[2,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // RSV
        lp[4] = lp[4] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[2,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Adeno
        lp[5] = lp[5] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[2,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // RhinoEntero
        lp[6] = lp[6] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[2,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Human Boca
        lp[7] = lp[7] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[2,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // Para
        lp[8] = lp[8] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[2,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // Corona
        lp[9] = lp[9] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[2,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Human meta
        lp[10] = lp[10] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[2,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Chlamy
        lp[11] = lp[11] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[2,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Myco
        lp[12] = lp[12] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[2,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Brucella
        lp[13] = lp[13] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[2,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[2,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[2,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // Dengue
        lp[14] = lp[14] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[2,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[2,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[2,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
           // Zika
        lp[15] = lp[15] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[2,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[2,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[2,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Chik
        lp[16] = lp[16] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[2,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[2,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[2,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
   
          // Lepto
        lp[18] = lp[18] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[2,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[2,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[2,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[2,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // Burk
        lp[19] = lp[19] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[2,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // TS
        lp[20] = lp[20] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[2,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
           // NTS
        lp[21] = lp[21] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[2,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Staph
        lp[22] = lp[22] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[2,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // E.coli
        lp[23] = lp[23] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[2,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
            // OtherEntero
        lp[24] = lp[24] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[2,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
           // Kleb
        lp[25] = lp[25] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[2,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
          // Strep pneu
        lp[26] = lp[26] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[2,37,n])*I37[n];
         // Other not tested for
        lp[27] = lp[27] + binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
       
    target += log_sum_exp(lp);
    //target += (lp);
  }

//Controls
    for(n in (N_Cases+1):N){

      //lp2 = binomial_lpmf(t1[n]|1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n]|1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n]|1,prob[1,3,n])*I3[n] + binomial_lpmf(t4[n]|1,prob[1,4,n])*I4[n] + binomial_lpmf(t5[n]|1,prob[1,5,n])*I5[n]+ binomial_lpmf(t6[n]|1,prob[1,6,n])*I6[n];
       lp2 = binomial_lpmf(t1[n] | 1, prob[1,1,n])*I1[n] + binomial_lpmf(t2[n] | 1, prob[1,2,n])*I2[n] + binomial_lpmf(t3[n] | 1, prob[1,3,n])*I3[n] + binomial_lpmf(t4[n] | 1, prob[1,4,n])*I4[n] + binomial_lpmf(t5[n] | 1, prob[1,5,n])*I5[n] + binomial_lpmf(t6[n] | 1, prob[1,6,n])*I6[n]+ binomial_lpmf(t7[n] | 1, prob[1,7,n])*I7[n]+ binomial_lpmf(t8[n] | 1, prob[1,8,n])*I8[n]+ binomial_lpmf(t9[n] | 1, prob[1,9,n])*I9[n]+ binomial_lpmf(t10[n] | 1, prob[1,10,n])*I10[n]+
       binomial_lpmf(t11[n] | 1, prob[1,11,n])*I11[n] + binomial_lpmf(t12[n] | 1, prob[1,12,n])*I12[n] + binomial_lpmf(t13[n] | 1, prob[1,13,n])*I13[n] + binomial_lpmf(t14[n] | 1, prob[1,14,n])*I14[n] + binomial_lpmf(t15[n] | 1, prob[1,15,n])*I15[n] + binomial_lpmf(t16[n] | 1, prob[1,16,n])*I16[n]+ binomial_lpmf(t17[n] | 1, prob[1,17,n])*I17[n]+ binomial_lpmf(t18[n] | 1, prob[1,18,n])*I18[n]+ binomial_lpmf(t19[n] | 1, prob[1,19,n])*I19[n]+ binomial_lpmf(t20[n] | 1, prob[1,20,n])*I20[n]+
       binomial_lpmf(t21[n] | 1, prob[1,21,n])*I21[n] + binomial_lpmf(t22[n] | 1, prob[1,22,n])*I22[n] + binomial_lpmf(t23[n] | 1, prob[1,23,n])*I23[n] + binomial_lpmf(t24[n] | 1, prob[1,24,n])*I24[n] + binomial_lpmf(t25[n] | 1, prob[1,25,n])*I25[n] + binomial_lpmf(t26[n] | 1, prob[1,26,n])*I26[n]+ binomial_lpmf(t27[n] | 1, prob[1,27,n])*I27[n]+ binomial_lpmf(t28[n] | 1, prob[1,28,n])*I28[n]+ binomial_lpmf(t29[n] | 1, prob[1,29,n])*I29[n]+ binomial_lpmf(t30[n] | 1, prob[1,30,n])*I30[n]+
       binomial_lpmf(t31[n] | 1, prob[1,31,n])*I31[n] + binomial_lpmf(t32[n] | 1, prob[1,32,n])*I32[n] + binomial_lpmf(t33[n] | 1, prob[1,33,n])*I33[n] + binomial_lpmf(t34[n] | 1, prob[1,34,n])*I34[n] + binomial_lpmf(t35[n] | 1, prob[1,35,n])*I35[n] + binomial_lpmf(t36[n] | 1, prob[1,36,n])*I36[n]+ binomial_lpmf(t37[n] | 1, prob[1,37,n])*I37[n];
       
   //+ binomial_lpmf(y[n,2] | 1, prob[1,2,n])+ binomial_lpmf(y[n,3] | 1, prob[1,3,n])+ binomial_lpmf(y[n,4] | 1, prob[1,4,n])+ binomial_lpmf(y[n,5]|1,prob[1,5,n]);
        target += (lp2);
}
  
}
