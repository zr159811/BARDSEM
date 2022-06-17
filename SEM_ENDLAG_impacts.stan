//BARDSEM from Roman and Brandt (2020)
// Zachary Roman Ph.D.
// <ZacharyJoseph.Roman@UZH.CH>
data {
 int<lower =0 > N;
 int<lower =0 > Kx;
 int<lower =0 > Ky;
 matrix [N, Kx] x;
 matrix [N, Ky] y;
 matrix<lower = 0, upper = 1>[N, N] W;
 matrix<lower = 0, upper = 1>[N, N] I;
 row_vector[N] IV;
}

parameters {

 vector [4] b1;
 vector<lower=0>[Kx] sigmax;
 vector<lower=0>[Ky] sigmay;
 //real<lower=0> sigmaeta;
 real<lower=0> sigmaeta;
 vector<lower=0>[2] sigmaxi;
 cholesky_factor_corr[2] L1;
 //test
 vector[N] zeta;
 matrix [N ,2] zi;

 vector[Kx] tx; //X Item intercepts
 vector[Ky-1] ty; //Y Item intercepts
 
 vector<lower=0>[Kx-2] lx; //X item Loadings
 vector<lower=0>[Ky-1] ly; //Y item Loadings
 
 real<lower = 0, upper = 1> rho;
}

transformed parameters {


}


model {
 matrix [N, Kx] mux ; // E [ x | xi ]
 matrix [N, Ky] muy ; // E [ y | eta ]
 vector [N] eta ;   // E [ eta | xi ]
 matrix [N ,2] xi ;     // Xi
 
 xi = zi * diag_pre_multiply(sigmaxi,L1)'; 

 eta = b1[1] + 
       b1[2]*xi[,1] +
       b1[3]*xi[,2] +
       b1[4]*(xi[,2].*xi[,1]) +
       zeta +
       rho*(W*zeta);
 
 mux [,1] = tx[1]   + xi[,1];
 mux [,2] = tx[2]   + lx[1]*xi[,1];
 mux [,3] = tx[3]   + lx[2]*xi[,1];

 mux [,4] = tx[4]   + xi[,2];
 mux [,5] = tx[5]   + lx[3]*xi[,2];
 mux [,6] = tx[6]   + lx[4]*xi[,2];
 
 muy [,1] = eta;
 muy [,2] = ty[1]  + ly[1]*eta;
 muy [,3] = ty[2]  + ly[2]*eta;


  for(z in 1:Kx){x[,z] ~ normal(mux[,z], sigmax [z]);}
  for(z in 1:Ky){y[,z] ~ normal(muy[,z], sigmay [z]);}
  
  // HB: Is the whole equation eta=alpha*(Wx*X*b1)+zeta? 
  // you generate data based on eta=X*b1+zeta + alpha*Wx*zeta
  
  //eta ~ normal(alpha*(W*mueta), sigmaeta);    // latent
  // HB: both done work
  zeta ~ normal(0, sigmaeta);    // latent
  //zeta ~ normal(0, sigmaeta);    // latent
  
  to_vector(zi) ~ normal(0,1);     // Cholesky
  b1 ~ normal(0,1);
  sigmax ~ cauchy(0,2.5);          // Prior SDs
  sigmay ~ cauchy(0,2.5);
  //sigmaeta ~ cauchy(0,2.5);
  sigmaxi ~ cauchy(0,2.5);
  sigmaeta ~ cauchy(0,2.5);
  L1 ~ lkj_corr_cholesky(2);
  tx ~ normal(0,1);
  ty ~ normal(0,1);
  lx ~ normal(0,1);
  ly ~ normal(0,1);
  rho ~ uniform(0,1); 
  //alpha ~ beta(1,1); 
}


generated quantities{
  //CPDM formulation with Gen. Quant. 
  //Provides posterior estimates for spatial spillover
  matrix[2,2] phi;
  matrix[N,N] part_deriv1;
  matrix[N,N] part_deriv2;
  matrix[N,N] part_deriv3;
  
  real direct1;
  real direct2;
  real direct3;
  
  real indirect1;
  real indirect2;
  real indirect3;
  
  real total1;
  real total2;
  real total3;
    
  phi = diag_pre_multiply(sigmaxi,L1)*diag_pre_multiply(sigmaxi,L1)';
  
  part_deriv1 = (I - rho*W)' * I * b1[2];
  part_deriv2 = (I - rho*W)' * I * b1[3];
  part_deriv3 = (I - rho*W)' * I * b1[4];
  
  direct1 = mean(diagonal(part_deriv1));
  total1 = mean(part_deriv1);
  indirect1 = (total1 - direct1)/2;
  
  direct2 = mean(diagonal(part_deriv2));
  total2 = mean(part_deriv2);
  indirect2 = (total2 - direct2)/2;
  
  direct3 = mean(diagonal(part_deriv3));
  total3 = mean(part_deriv3);
  indirect3 = (total3 - direct3)/2;
}
