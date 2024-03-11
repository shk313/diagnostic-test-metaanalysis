// HSROC not individual level model for DTA meta-analysis


data {
  int<lower=0> J; // Number of unique studies
  int<lower=0> l; //Number of studies, > J 
  int<lower=0> N; // Total number of people
  int<lower=0> ind_unique[N]; //Indicator for which study a person belong to and their id within that study
  //int<lower=0> ind_study[N]; //Indicator for which study a person belong to and their id within that study
  int<lower=0> study[l]; //Indicator for which unique study a person belong to and their id within that study
  int T1[N]; // index test
  int T2[N]; 
  int<lower=1> ref_test[N]; // indicator for which reference test a person belongs to
  int<lower=1> RT; // Number of reference tests
  // int<lower=0> R; //Number of regions
  // vector<lower=0, upper=1>[l] regionA; // 1 = yes to pathogen A, 0 = no
  // vector<lower=0, upper=1>[l] regionB;
  // vector<lower=0, upper=1>[l] regionC;
  // vector<lower=0, upper=1>[l] regionD;
  // int<lower=1, upper=R> region[l];
  //vector[RT] L; // lower bounds for sensitviity of reference test
  //array[RT] real L;
  //int<lower=1> ref_test_study[J]; 
  
}




parameters {
  // real<lower=0.5,upper=1> Se_Ref[RT];
  // real<lower=0.5,upper=1> Sp_Ref[RT];
  real THETA; // global mean of theta (A+B)
  real LAMBDA; // global mean of alpha (A+B)
  // real THETA_region[R];
  // real LAMBDA_region[R];
  real<lower=-0.75, upper=0.75> beta; // common scale param
  real<lower=0> sd_theta; // changed to vector from real
  real<lower=0> sd_alpha;
  vector[l] theta_; // for reparameterisation of theta
  vector[l] alpha_; // for reparameterisation of alpha
  //vector[J] alpha_;
  real<lower=0,upper=1> prevalence[J];
  // vector[N] N_RE; // e.g. infection intensity
  // real<lower=0> sd_re[RT+1]; // sd for infection intensity between positives in each test
  //real sd_re[RT+1];
  // real<lower=0.2,upper=1> Sp_Ref[RT];
  // real<lower=0.2,upper=1> Se_Ref[RT];

  real<lower=0,upper=1> Sp_Ref1;
  real<lower=inv_logit(logit(1-Sp_Ref1)),upper=1> Se_Ref1;
  //real<lower=0,upper=1> Se_Ref1;
  real<lower=0,upper=1> Sp_Ref2;
  real<lower=inv_logit(logit(1-Sp_Ref2)),upper=1> Se_Ref2;
  //real<lower=0,upper=1> Se_Ref2;
  real<lower=0,upper=1> Sp_Ref3;
  //real<lower=inv_logit(logit(1-Sp_Ref3)),upper=1> Se_Ref3;
  real<lower=0,upper=1> Se_Ref3;
  real<lower=0,upper=1> Sp_Ref4;
  real<lower=inv_logit(logit(1-Sp_Ref4)),upper=1> Se_Ref4;
  real<lower=0,upper=1> Sp_Ref5;
  real<lower=inv_logit(logit(1-Sp_Ref5)),upper=1> Se_Ref5;
  //real<lower=0,upper=1> Se_Ref5;
  real<lower=0,upper=1> Sp_Ref6;
  real<lower=inv_logit(logit(1-Sp_Ref6)),upper=1> Se_Ref6;
  real<lower=0,upper=1> Sp_Ref7;
  real<lower=inv_logit(logit(1-Sp_Ref7)),upper=1> Se_Ref7;
  real<lower=0,upper=1> Sp_Ref8;
  real<lower=inv_logit(logit(1-Sp_Ref8)),upper=1> Se_Ref8;
  // real<lower=0,upper=1> Se_Ref8;
  real<lower=0,upper=1> Sp_Ref9;
  real<lower=inv_logit(logit(1-Sp_Ref9)),upper=1> Se_Ref9;
  // real<lower=0,upper=1> Sp_Ref10;
  // real<lower=inv_logit(logit(1-Sp_Ref10)),upper=1> Se_Ref10;
  // real<lower=0,upper=1> Sp_Ref11;
  // real<lower=inv_logit(logit(1-Sp_Ref11)),upper=1> Se_Ref11;
  // real<lower=0,upper=1> Sp_Ref12;
  // real<lower=inv_logit(logit(1-Sp_Ref12)),upper=1> Se_Ref12;
  // real<lower=0,upper=1> Sp_Ref13;
  // real<lower=inv_logit(logit(1-Sp_Ref13)),upper=1> Se_Ref13;
  // real<lower=0,upper=1> Sp_Ref14;
  // real<lower=inv_logit(logit(1-Sp_Ref14)),upper=1> Se_Ref14;
  // real<lower=0,upper=1> Sp_Ref15;
  // real<lower=inv_logit(logit(1-Sp_Ref15)),upper=1> Se_Ref15;
  // real<lower=0,upper=1> Sp_Ref16;
  // real<lower=inv_logit(logit(1-Sp_Ref16)),upper=1> Se_Ref16;
  // real<lower=0,upper=1> Sp_Ref17;
  // real<lower=inv_logit(logit(1-Sp_Ref17)),upper=1> Se_Ref17;
  // real<lower=0,upper=1> Sp_Ref18;
  // real<lower=inv_logit(logit(1-Sp_Ref18)),upper=1> Se_Ref18;
  // real<lower=0,upper=1> Sp_Ref19;
  // real<lower=inv_logit(logit(1-Sp_Ref19)),upper=1> Se_Ref19;
  // real<lower=0,upper=1> Sp_Ref20;
  // real<lower=inv_logit(logit(1-Sp_Ref20)),upper=1> Se_Ref20;
  //real<lower=0,upper=1> Se_Ref[RT];
  //real<lower=0> fe_pos[RT];

}

transformed parameters {
  
  vector<lower=0,upper=1>[N] prob[2,2]; // I want it to be N rows, two columns(for Se and Sp), with two shelves(one for each test)
  simplex[2] prev[J];
  vector<lower=0, upper=1>[l] Sp; // index test
  vector<lower=0, upper=1>[l] Se; // index test
  //vector<lower=0, upper=1>[l] Sp; // index test
  vector[l] theta;
  //vector[1] theta;
  vector[l] alpha;
  //vector[J] alpha;
  //real x;
  real y;
  // real<lower=1> sum_accuracy[J];
  // vector[l] T;
  // vector[l] L;
  
  
  

  theta = THETA + theta_*sd_theta; // implies theta~N(THETA,prec[1])
  alpha = LAMBDA + alpha_*sd_alpha; //implies alpha~N(LAMBDA,prec[2])

  
for(j in 1:J){

  
  prev[j,1] = 1-prevalence[j];
  prev[j,2] = prevalence[j];
  
  
  
  // T[j] = THETA_region[1]*regionA[j] + THETA_region[2]*regionB[j] + THETA_region[3]*regionC[j] + THETA_region[4]*regionD[j];
  // L[j] = LAMBDA_region[1]*regionA[j] + LAMBDA_region[2]*regionB[j] + LAMBDA_region[3]*regionC[j] + LAMBDA_region[4]*regionD[j];
  // 
  // theta[j] = T[j]+ theta_[j]* sd_theta[region[j]]; // implies theta~N(THETA,prec[1])
  // alpha[j] = L[j] + alpha_[j]* sd_alpha[region[j]]; 
}


for(j in 1:l){ 
  
   // Equation from Dendukuri et al
  Se[j] = inv_logit(-(theta[study[j]] - alpha[study[j]]/2)/exp(beta/2));
  Sp[j] = inv_logit((theta[study[j]] + alpha[study[j]]/2)/exp(beta/2));
  // 
  // Se[j] = inv_logit(-(theta[study[j]] - alpha[1]/2)/exp(beta/2));
  // Sp[j] = inv_logit((theta[study[j]] + alpha[1]/2)/exp(beta/2));
  // Se[j] = inv_logit(-(theta[1] - alpha[1]/2)/exp(beta/2));
  // Sp[j] = inv_logit((theta[1] + alpha[1]/2)/exp(beta/2));
  // sum_accuracy[j] = Se[j]+Sp[j];

}
  

  // Conditional independence
 // for(n in 1:N){
 //   prob[2,1,n] = Se[ind_unique[n]];   //person, se or sp, test
 //   prob[1,1,n] = 1-Sp[ind_unique[n]];
 //   prob[2,2,n] = Se_Ref[ref_test[n]];
 //   prob[1,2,n] = 1-Sp_Ref[ref_test[n]];
 // }
 
 
for(n in 1:N){
  //x = ind_unique[n]; // grab study id of person
  y = ref_test[n]; // grab reference test id of person
  
  prob[2,1,n] = Se[ind_unique[n]]; //person, se or sp, test
  //prob[2,1,n] = inv_logit(logit(Se[ind_unique[n]])+fe_pos[ref_test[n]]);
  //prob[2,1,n] = inv_logit(logit(Se[ind_unique[n]])+N_RE[n]*sd_re[RT+1]);
  prob[1,1,n] = 1-Sp[ind_unique[n]];
  
 if(y == 1){
    prob[2,2,n] = Se_Ref1;
    //prob[2,2,n] = inv_logit(logit(Se_Ref1)+fe_pos[1]);
    // prob[2,2,n] = inv_logit(logit(Se_Ref1)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref1;
    }
     if(y == 2){
    prob[2,2,n] = Se_Ref2;
    //prob[2,2,n] = inv_logit(logit(Se_Ref2)+fe_pos[2]);
    // prob[2,2,n] = inv_logit(logit(Se_Ref2)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref2;
    }
      if(y == 3){
    prob[2,2,n] = Se_Ref3;
    //prob[2,2,n] = inv_logit(logit(Se_Ref3)+fe_pos[3]);
    //prob[2,2,n] = inv_logit(logit(Se_Ref3)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref3;
    }
      if(y == 4){
    prob[2,2,n] = Se_Ref4;
    //prob[2,2,n] = inv_logit(logit(Se_Ref4)+fe_pos[4]);
    //prob[2,2,n] = inv_logit(logit(Se_Ref4)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref4;
    }
    if(y == 5){
    prob[2,2,n] = Se_Ref5;
    //prob[2,2,n] = inv_logit(logit(Se_Ref5)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref5;
    }
    if(y == 6){
    prob[2,2,n] = Se_Ref6;
    // prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // prob[2,2,n] = inv_logit(logit(Se_Ref6)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref6;
    }
    if(y == 7){
    prob[2,2,n] = Se_Ref7;
    //prob[2,2,n] = inv_logit(logit(Se_Ref7)+fe_pos[7]);
    //prob[2,2,n] = inv_logit(logit(Se_Ref7)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref7;
    }
        if(y == 8){
    prob[2,2,n] = Se_Ref8;
    //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref8;
    }
        if(y == 9){
    prob[2,2,n] = Se_Ref9;
    //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    prob[1,2,n] = 1-Sp_Ref9;
    }
    //     if(y == 10){
    // prob[2,2,n] = Se_Ref10;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref10;
    // }
    //     if(y == 11){
    // prob[2,2,n] = Se_Ref11;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref11;
    // }
    //    if(y == 12){
    // prob[2,2,n] = Se_Ref12;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref12;
    // }
    //   if(y == 13){
    // prob[2,2,n] = Se_Ref13;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref13;
    // }
    //   if(y == 14){
    // prob[2,2,n] = Se_Ref14;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref14;
    // }
    //   if(y == 15){
    // prob[2,2,n] = Se_Ref15;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref15;
    // }
    //     if(y == 16){
    // prob[2,2,n] = Se_Ref16;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref16;
    //   }
    //     if(y == 17){
    // prob[2,2,n] = Se_Ref17;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref17;
    // }
    //       if(y == 18){
    // prob[2,2,n] = Se_Ref18;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref18;
    // }
    //     if(y == 19){
    // prob[2,2,n] = Se_Ref19;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref19;
    //   }
    //   if(y == 20){
    // prob[2,2,n] = Se_Ref20;
    // //prob[2,2,n] = inv_logit(logit(Se_Ref6)+fe_pos[6]);
    // //prob[2,2,n] = inv_logit(logit(Se_Ref8)+N_RE[n]*sd_re[ref_test[n]]);
    // prob[1,2,n] = 1-Sp_Ref20;
    // }


 }

 
 

  // Conditional dependence structure
  // for(n in 1:N){
  // prob[2,1,n] = inv_logit(logit(Se[ind_unique[n]]) + (N_RE[n]*sd_re[RT+1]));
  // prob[1,1,n] = inv_logit(logit(1-Sp[ind_unique[n]]));
  // prob[2,2,n] = inv_logit(logit(Se_Ref[ref_test[n]]) + N_RE[n]*sd_re[ref_test[n]]);
  // prob[1,2,n] = inv_logit(logit(1-Sp_Ref[ref_test[n]]));
  // }
  
}

model {
  real lp[2]; // 2d vector with each entry corresponding
  // to the un-normalised posterior prob for diseased and not diseased
  
  // THETA_region~normal(0,1);
  // LAMBDA_region~normal(0,1.5); 
  
  THETA~normal(0,1);
  LAMBDA~normal(0,1.5); //0,2

  // sd_re~normal(0,1);
  // N_RE~std_normal();
  // 
  // fe_pos~gamma(1,1);

  theta_ ~ std_normal(); // implies theta~normal(THETA, prec[1])
  alpha_ ~ std_normal(); // implies alpha~normal(LAMBDA, prec[1])

  sd_theta ~ normal(0,1);
  sd_alpha ~ normal(0,1);
  
  prevalence~beta(1,1);


  //RT-PCR all the data
  // Se_Ref1~beta(49.4,19.1); // NS1 ELISA acute
  // Se_Ref2~beta(22.9,6.3); // IgM ELISA acute
  // Se_Ref3~beta(3.3,0.4); // RPA
  // Se_Ref4~beta(98.9,10.8); // NS1 RDT
  // Se_Ref5~beta(1421.4,731.8); // NS1 Ag Strip
  // Se_Ref6~beta(3.3,0.4); // Pocket dengue virus
  // Se_Ref7~beta(1.9,0.3); // Virus Isolation
  // Se_Ref8~beta(3.9,0.4); // RT-LAMP
  // Se_Ref9~beta(6.5,0.4); // RT-pCR
  // Se_Ref10~beta(4.5,0.4); // RDT+IgM
  // Se_Ref11~beta(4.5,0.4); // Combined serology
  // Se_Ref12~beta(6.5,0.4); // Neutralization
  // Se_Ref13~beta(3.3,0.4); // RT-RAA-LFT
  // Se_Ref14~beta(3.3,0.4); // STH_PAS
  // Se_Ref15~beta(26.6,19.9); // Torniquet test
  // Se_Ref16~beta(0.9,1.5); // ZIKV RT-PCR
  // Se_Ref17~beta(3.3,0.4); // Urine IgM
  // //
  // Sp_Ref1~beta(78.5,0.8);
  // Sp_Ref2~beta(115.6,15.8);
  // Sp_Ref3~beta(3.3,0.4);
  // Sp_Ref4~beta(15.6,5.4);
  // Sp_Ref5~beta(86.2,0.5);
  // Sp_Ref6~beta(3.3,0.4);
  // Sp_Ref7~beta(78.5,0.8);
  // Sp_Ref8~beta(3.9,0.4);
  // Sp_Ref9~beta(6.5,0.4);
  // Sp_Ref10~beta(4.5,0.4);
  // Sp_Ref11~beta(4.5,0.4);
  // Sp_Ref12~beta(6.5,0.4);
  // Sp_Ref13~beta(3.3,0.4);
  // Sp_Ref14~beta(3.3,0.4);
  // Sp_Ref15~beta(49.8,1);
  // Sp_Ref16~beta(0.9,1.5);
  // Sp_Ref17~beta(3.3,0.4);
  
    // RT-PCR 0-4days
  // Se_Ref1~beta(5,3); // NS1 RDT
  // Sp_Ref1~beta(5,3);
  // Se_Ref2~beta(5,2); // NS1 ELISA
  // Sp_Ref2~beta(5,2);
  // Se_Ref3~beta(8,5); // NS1 ELISA in chik patients
  // Sp_Ref3~beta(8,5);
  
     // RT-PCR 1-7days
  // Se_Ref1~beta(5,5); // NS1 RDT
  // Sp_Ref1~beta(5,5);
  // Se_Ref2~beta(5,3); // IgM ELISA
  // Sp_Ref2~beta(5,3);
  // Se_Ref3~beta(5,4); // NS1 ELISA
  // Sp_Ref3~beta(5,4);
  // Se_Ref4~beta(1,1); // Viral isolaiton
  // Sp_Ref4~beta(78.5,0.8);
  // Se_Ref5~beta(5,1); // neutralization
  // Sp_Ref5~beta(5,1);
  // Se_Ref6~beta(5,1); // IgM ELISA 5-7 days
  // Sp_Ref6~beta(5,1);
  // Se_Ref7~beta(5,1); // sequencing - PCR
  // Sp_Ref7~beta(5,1);
  // Se_Ref8~beta(5,2); // Urine DENV IgM
  // Sp_Ref8~beta(5,2);
  
       // RT-PCR ALl
  // Se_Ref1~beta(5,2); // IgM ELISA 1- 14days
  // Sp_Ref1~beta(5,2);
  // Se_Ref2~beta(5,3); // RPA?
  // Sp_Ref2~beta(5,3);
  // Se_Ref3~beta(4,5); // NS1 RDT 1-14days
  // Sp_Ref3~beta(4,5);
  // Se_Ref4~beta(5,3); // IgM ELISA 1-7days
  // Sp_Ref4~beta(5,4);
  // Se_Ref5~beta(5,1); // Pocket dengue virus??
  // Sp_Ref5~beta(5,1);
  // Se_Ref6~beta(5,4); // NS1 ELISA 1-7days
  // Sp_Ref6~beta(5,4);
  // Se_Ref7~beta(1,1); // Virus isolation
  // Sp_Ref7~beta(78.5,0.8);
  // Se_Ref8~beta(5,1); // RT_LAMP
  // Sp_Ref8~beta(5,1);
  // Se_Ref9~beta(25,5); // NS1 RDT 0-4days
  // Sp_Ref9~beta(25,5);
  // Se_Ref10~beta(12.3,0.6); // RT_PCR 1-14days
  // Sp_Ref10~beta(12.3,0.6);
  // Se_Ref11~beta(4,5); // NS1 ELISA 1-14days
  // Sp_Ref11~beta(4,5);
  // Se_Ref12~beta(25,5); // Combined serology 1-14days
  // Sp_Ref12~beta(25,5);
  // Se_Ref13~beta(25,5); // NS1 ELISA 0-4days
  // Sp_Ref13~beta(25,5);
  // Se_Ref14~beta(25,5); // neutralisation
  // Sp_Ref14~beta(25,5);
  // Se_Ref15~beta(20,5); // NS1 RDT 1-5days
  // Sp_Ref15~beta(25,5);
  // Se_Ref16~beta(1,1); // TOrniquet test? 1-14 days
  // Sp_Ref16~beta(1,1);
  // Se_Ref17~beta(5,1); //IgM ELISA 5-7 days
  // Sp_Ref17~beta(4,3);
  // Se_Ref18~beta(4,3); //IgM ELISA 0-5 days
  // Sp_Ref18~beta(4,3);
  // Se_Ref19~beta(25,5); //Sequencing
  // Sp_Ref19~beta(25,5);
  // Se_Ref20~beta(1,1); //Urine Denv IgM 1-7days
  // Sp_Ref20~beta(4,3);

// Viral Neutralization  
  // Se_Ref1~beta(1,3); // IgG ELISA - acute or combined?
  // Sp_Ref1~beta(1,1);
  // Se_Ref2~beta(5,1); // IgM ELISA acute
  // Sp_Ref2~beta(1,1);
  // Se_Ref3~beta(12.3,0.6);// RT-PCR
  // Sp_Ref3~beta(12.3,0.6);
  // 
  // 
  // Viral Neutralization without asympt 
  // Se_Ref1~beta(5,3); // RT-PCR or IgM and IgM
  // Sp_Ref1~beta(5,2);
  // Se_Ref2~beta(5,2); // IgG RDT or ELISA
  // Sp_Ref2~beta(5,3);
  // Se_Ref3~beta(1,1);// RT-PCR
  // Sp_Ref3~beta(5,1);
  // Se_Ref1~beta(5,1);// IgG - looks like paired sample results
  // Sp_Ref1~beta(5,1);
  
  // Viral neut with asympt
  // Se_Ref1~beta(5,1); // IgG RDT
  // Sp_Ref1~beta(5,1);
  // Se_Ref2~beta(2.3,6); // NS1 ELISA
  // Sp_Ref2~beta(2.3,6);
  // Se_Ref4~beta(5,1);// IgM ELISA
  // Sp_Ref4~beta(5,1);
  // Se_Ref3~beta(1,1);// IgG - looks like paired sample results
  // Sp_Ref3~beta(1,1);
  // Se_Ref5~beta(5,1);// IgM - looks like paired sample results
  // Sp_Ref5~beta(5,1);
  
  
    // NS1 ELISA All
  // Se_Ref1~beta(3,2); // Igm elisa 1-7days
  // Sp_Ref1~beta(3,2);
  // Se_Ref2~beta(3,2); // IgM ELISA 1-14days
  // Sp_Ref2~beta(3,2);
  // Se_Ref3~beta(1,1);// ns1 rdt 1-14days
  // Sp_Ref3~beta(5,1);
  // Se_Ref4~beta(12.3,0.6);// rt-pcr 1-7 or 1-14
  // Sp_Ref4~beta(12.3,0.6);
  // Se_Ref5~beta(2.3,6);// igm ELISA 0-4 days
  // Sp_Ref5~beta(2.3,6);
  // Se_Ref6~beta(1,1);// virus isolation 1-14days
  // Sp_Ref6~beta(3,1);
  // Se_Ref7~beta(3,2);// virus isolation 0-4days
  // Sp_Ref7~beta(78.5,0.8);
  // Se_Ref8~beta(3,2);// igm rdt 1-14days
  // Sp_Ref8~beta(3,2);
  // Se_Ref9~beta(12.3,0.6);// ns1 rdt 1-7 days
  // Sp_Ref9~beta(12.3,0.6);
  // Se_Ref10~beta(12.3,0.6);// pcr 0-4 days
  // Sp_Ref10~beta(12.3,0.6);
  // Se_Ref11~beta(10,4);// ns1 elisa on filter paper 0-4days
  // Sp_Ref11~beta(10,4);
  
      // NS1 ELISA 0-4days
  // Se_Ref1~beta(2.3,6); // IgM ELISA
  // Sp_Ref1~beta(2.3,6);
  // Se_Ref2~beta(12.3,0.6); // NS1 ELISA/ NS1 RDT
  // Sp_Ref2~beta(5,1);
  // Se_Ref3~beta(12.3,0.6); // PCR >80%
  // Sp_Ref3~beta(12.3,0.6);
  // Se_Ref4~beta(1,1); // Virus isolation not sensitive but specific
  // Sp_Ref4~beta(78.5,0.8);
  // Se_Ref5~beta(2.3,6); // IgM ELISA for days 1-2 only
  // Sp_Ref5~beta(2.3,6);


  
        // NS1 ELISA 1-7days
  // Se_Ref1~beta(12.3,0.6); // PCR
  // Sp_Ref1~beta(12.3,0.6);
  // Se_Ref2~beta(3,2); // IgM ELISA
  // Sp_Ref2~beta(3,2);
  // Se_Ref3~beta(5,1); // NS1 RDT
  // Sp_Ref3~beta(5,1);
  
      // NS1 ELISA 1-14days
  // Se_Ref2~beta(3,3); // IgM ELISA
  // Sp_Ref2~beta(3,3);
  // Se_Ref1~beta(3,3); // NS1 RDT
  // Sp_Ref1~beta(3,3);
  // Se_Ref3~beta(5.2,0.5); // PCR >80%
  // Sp_Ref3~beta(5.2,0.5);
  // Se_Ref4~beta(1,1); // Virus isolation not sensitive but specific
  // Sp_Ref4~beta(30,2);
  // Se_Ref5~beta(3,2); // IgM ELISA 4+ days
  // Sp_Ref5~beta(3,2);
  // Se_Ref6~beta(5,1); // IgM RDT
  // Sp_Ref6~beta(5,1);
  
  // NS1 QUADAS
  // Se_Ref1~beta(5,2); // IgM ELISA 5+ days
  // Sp_Ref1~beta(5,2);
  // Se_Ref2~beta(1,1); // Virus isolation
  // Sp_Ref2~beta(78.5,0.8);
  // Se_Ref3~beta(12.3,0.6); // PCR >80%
  // Sp_Ref3~beta(12.3,0.6);
  // Se_Ref4~beta(2.3,6); // Igm ELISA 0-5 days
  // Sp_Ref4~beta(5,2);
  // Se_Ref5~beta(5,2); // NS1 ELISA on filter paper
  // Sp_Ref5~beta(5,2);
  // Se_Ref6~beta(5,2); // NS1 RDT
  // Sp_Ref6~beta(5,2);
  // Se_Ref7~beta(1,1); // IgM RDT all dpo
  // Sp_Ref7~beta(5,2);
  
  
  // IgG QUADAS
  // Se_Ref1~beta(5,2); // IgM ELISA 5+ days
  // Sp_Ref1~beta(5,2);
  // Se_Ref2~beta(1,1); // Virus isolation
  // Sp_Ref2~beta(78.5,0.8);
  // Se_Ref3~beta(12.3,0.6); // PCR >80%
  // Sp_Ref3~beta(12.3,0.6);
  // Se_Ref4~beta(10,1.2); // PCR 0-4 days
  // Sp_Ref4~beta(10,1.2);
  // Se_Ref5~beta(10,1.2); // DENV infection
  // Sp_Ref5~beta(10,1.2);
  // Se_Ref6~beta(3,2); // ICT rapid IgG
  // Sp_Ref6~beta(5,2);
  // Se_Ref7~beta(1,1); // IgM 0-5 days
  // Sp_Ref7~beta(5,2);
  // Se_Ref8~beta(3,2); // dired serum spot IgG
  // Sp_Ref8~beta(5,2);
  
    // RT-PCR QUADAS
  Se_Ref1~beta(5,2); // IgM ELISA 1-14days
  Sp_Ref1~beta(12.3,0.6);
  Se_Ref2~beta(5,2); // NS1 RDT 1-7
  Sp_Ref2~beta(12.3,0.6);
  Se_Ref3~beta(4,2); // IgM ELISA early 1-7
  Sp_Ref3~beta(12.3,0.6);
  Se_Ref4~beta(1,1); // Virus isolation
  Sp_Ref4~beta(78.5,0.8);
  Se_Ref5~beta(10,1.2); // NS1 ELISA
  Sp_Ref5~beta(10,1.2);
  Se_Ref6~beta(12.3,0.6); // NS1 RDT 0-3 - should be good
  Sp_Ref6~beta(12.3,0.6);
  Se_Ref7~beta(12.3,0.6); // Combined
  Sp_Ref7~beta(12.3,0.6);
  Se_Ref8~beta(12.3,0.6); // Neutralization
  Sp_Ref8~beta(12.3,0.6);
  Se_Ref9~beta(12.3,0.6); // PCR
  Sp_Ref9~beta(12.3,0.6);
  
          // IgM ELISA 0-4days
  // Se_Ref1~beta(10,2); // NS1 ELISA
  // Sp_Ref1~beta(10,2);
  // Se_Ref2~beta(3,8); // IgM ELISA
  // Sp_Ref2~beta(10,2);
  // Se_Ref3~beta(12.3,0.6); // Composite (very good)
  // Sp_Ref3~beta(12.3,0.6);


      // IgM ELISA 1-7days
  // Se_Ref1~beta(5.2,0.5); // PCR
  // Sp_Ref1~beta(5.2,0.5);
  // Se_Ref2~beta(3,3); // NS1 ELISA
  // Sp_Ref2~beta(3,3);
  // Se_Ref3~beta(5,1); // IgG ELISA
  // Sp_Ref3~beta(5,1);
  // Se_Ref4~beta(1,1); // Virus isolation not sensitive but specific
  // Sp_Ref4~beta(30,2);
  // Se_Ref5~beta(3,2); // IgM RDT
  // Sp_Ref5~beta(3,2);
  // Se_Ref6~beta(3,2); // IgM ELISA
  // Sp_Ref6~beta(3,2);
  // Se_Ref7~beta(5,1); // Clincal score
  // Sp_Ref7~beta(5,1);
  // Se_Ref8~beta(3,2); // Dried serum spot IgM
  // Sp_Ref8~beta(3,2);
  // Se_Ref9~beta(3,2); // Urine IgM
  // Sp_Ref9~beta(3,2);
  
        // IgM ELISA 1-14days
  // Se_Ref3~beta(5.2,0.5); // PCR
  // Sp_Ref3~beta(5.2,0.5);
  // Se_Ref1~beta(3,3); // NS1 ELISA
  // Sp_Ref1~beta(3,3);
  // Se_Ref2~beta(5,1); // IgG ELISA
  // Sp_Ref2~beta(5,1);
  // Se_Ref4~beta(5,2); // RPA
  // Sp_Ref4~beta(5,2);
  // Se_Ref5~beta(3,2); // IgM RDT
  // Sp_Ref5~beta(3,2);
  // Se_Ref6~beta(2,10); // Zika pcr
  // Sp_Ref6~beta(2,10);
  // Se_Ref7~beta(7,2); // Composite
  // Sp_Ref7~beta(7,2);
  // Se_Ref8~beta(3,2); // Torniquet
  // Sp_Ref8~beta(3,2);
  // Se_Ref9~beta(2,10); // CHIKv IgM
  // Sp_Ref9~beta(2,10);
  // Se_Ref10~beta(5,2); // Raman Spec
  // Sp_Ref10~beta(5,2);
  // Se_Ref11~beta(5,2); // BIOCHIP
  // Sp_Ref11~beta(5,2);


      // IgG ELISA 0-7days
  // Se_Ref1~beta(12.3,0.6); // PCR
  // Sp_Ref1~beta(12.3,0.6);
  // Se_Ref2~beta(5,1); // IgM ELISA
  // Sp_Ref2~beta(5,1);
  // Se_Ref3~beta(5,1); // Rapid IgG
  // Sp_Ref3~beta(5,1);
  // Se_Ref4~beta(1,1); // Virus isolation
  // Sp_Ref4~beta(30,2);
  
  
    // IgG ELISA 0-14days
  // Se_Ref1~beta(5,3); // IgM ELISA
  // Sp_Ref1~beta(5,3);
  // Se_Ref2~beta(3,4); // NS1 ELISA
  // Sp_Ref2~beta(5,1);
  // Se_Ref3~beta(2,5); // ZIka something
  // Sp_Ref3~beta(2,5);
  // Se_Ref4~beta(1,1); // Virus isolation
  // Sp_Ref4~beta(5,1);
  // Se_Ref5~beta(12.3,0.6); // IgG ELISA
  // Sp_Ref5~beta(12.3,0.6);
  // Se_Ref6~beta(12.3,0.6); // Combination
  // Sp_Ref6~beta(12.3,0.6);
  // Se_Ref7~beta(1,1); // RT-PCR
  // Sp_Ref7~beta(12.3,0.6);
  // Se_Ref8~beta(4,2); // Special ELISA
  // Sp_Ref8~beta(4,2);
  // Se_Ref9~beta(4,2); // Raman spec
  // Sp_Ref9~beta(4,2);
  // Se_Ref10~beta(4,2); // BIOCHIP
  // Sp_Ref10~beta(4,2);
  
  // IgG ALL
  // Se_Ref1~beta(3,2); // Igm elisa 1-14days
  // Sp_Ref1~beta(3,2);
  // Se_Ref2~beta(3,2); // NS1 ELISA 1-14days
  // Sp_Ref2~beta(3,2);
  // Se_Ref3~beta(5,15);// ZIKA PRNT
  // Sp_Ref3~beta(5,15);
  // Se_Ref4~beta(1,1);// virus isolation 1-14days
  // Sp_Ref4~beta(12.3,0.6);
  // Se_Ref5~beta(12.3,0.6);// RT-PCR 1-7days
  // Sp_Ref5~beta(12.3,0.6);
  // Se_Ref6~beta(1,1);// IGg ELISA acute 1-14days 
  // Sp_Ref6~beta(3,1);
  // Se_Ref7~beta(12.3,0.6);// RT_PCR 0-4days
  // Sp_Ref7~beta(12.3,0.6);
  // Se_Ref8~beta(2,3);// igm ELISA 1-5days
  // Sp_Ref8~beta(2,3);
  // Se_Ref9~beta(12.3,0.6);// Dengue infection combined results
  // Sp_Ref9~beta(12.3,0.6);
  // Se_Ref10~beta(3,2);// rapid iGG acute
  // Sp_Ref10~beta(3,2);
  // Se_Ref11~beta(10,4);// special novel ELISA
  // Sp_Ref11~beta(10,4);
  // Se_Ref12~beta(8,4);// Raman spec, dried serum spot, biochip - what are these.
  // Sp_Ref12~beta(8,4);
  
  for(n in 1:N){
    for(k in 1:2){  
      // log posterior probability
      
      lp[k] = log(prev[ind_unique[n],k]);
      lp[k] = lp[k] + binomial_lpmf(T1[n] | 1, prob[k,1,n]) + binomial_lpmf(T2[n]| 1, prob[k,2,n]);
    }
    // marginalise out dependence on k and increment the overall log prob. by this amount
    target += log_sum_exp(lp);
  }

}  

generated quantities {
 real<lower=0,upper=1> Se_pooled;
 real<lower=0,upper=1> Sp_pooled;
 real<lower=0,upper=1> Se_pred;
 real<lower=0,upper=1> Sp_pred;
 

 real alpha_new;
 real theta_new;
 
 
 // theta_new = normal_rng(THETA_region, sd_theta);
 // alpha_new = normal_rng(LAMBDA_region, sd_alpha);

 theta_new = normal_rng(THETA,sd_theta);
 alpha_new = normal_rng(LAMBDA,sd_alpha);

 Se_pooled = inv_logit(-(THETA - LAMBDA/2) / exp(beta/2));
 Sp_pooled = inv_logit((THETA + LAMBDA/2) / exp(beta/2));
 
 Se_pred= inv_logit(-(theta_new - alpha_new/2) / exp(beta/2));
 Sp_pred= inv_logit((theta_new + alpha_new/2) / exp(beta/2));

}
