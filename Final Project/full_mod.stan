
data {
   int<lower=0> n;
   int<lower=0> y[n];
   real<lower=0,upper=1> atrac[n];
   int<lower=0> vac[n];
   int<lower=0> sexo[n];
   int<lower=0> edad[n];
   int<lower=0> dpto[n];
   int<lower=0> zbs[n];
   int<lower=0> tiempo[n];
   int<lower=0> N[n];
   int<lower=0> nperiods;
   int<lower=0> nzbs;
   
   int<lower=0> N_car;
   int<lower=0> N_edges;
   int node1[N_edges];
   int node2[N_edges];
   
   int<lower=0> N_carT;
   int<lower=0> N_edgesT;
   int node1T[N_edgesT];
   int node2T[N_edgesT];
}

parameters {
   real beta_0;
   real beta_1;
   real beta_2;
   vector[3] beta_3;
   matrix[6,241] theta_S;
   real<lower=0> sd_d;
   real<lower=0> sd_inter;
   vector[6] inter;
   real<lower=0> sd_het;
   real<lower=0> sd_spat;
   real<lower=0,upper=1> ro;
   vector[24] alpha_d;
   matrix[nperiods, nzbs] het;
}

transformed parameters {
   //Make parameters a function of other parameters
   vector<lower=0, upper = 1>[n] p;
   matrix[nperiods, nzbs] BYM;
   matrix[nperiods, nzbs] theta_ST;
   
   for(z in 1:nzbs){
      BYM[1,z] = sd_spat*theta_S[1,z] + sd_het*het[1,z];
      theta_ST[1,z] = pow(1-ro*ro, -0.5)*BYM[1,z];
   }
   
   for(t in 2:nperiods){
      for(z in 1:nzbs){
         BYM[t,z] = sd_spat*theta_S[t,z] + sd_het*het[t,z];
         theta_ST[t,z] = ro*theta_ST[t-1,z]+BYM[1,z];
      }
   }
   
   for(i in 1:n){
      p[i] = inv_logit(logit(atrac[i])+beta_0+ beta_1*vac[i]+ beta_2*sexo[i]+beta_3[edad[i]]+alpha_d[dpto[i]] +inter[tiempo[i]] + theta_ST[tiempo[i],zbs[i]]);
   }
   
}

model {
   real tau_d;
   real tau_inter;

   for(i in 1:n){
      y[i] ~ binomial(N[i], p[i]);
   }
   
   for(z in 1:nzbs){
      het[1,z] ~ normal(0,1);
   }
   
   target += -0.5 * dot_self(theta_S[1, node1] - theta_S[1, node2]);
   sum(theta_S[1,:]) ~ normal(0, 0.01 * N_car);
   //theta_S[1,1:nzbs]~car.normal(map[],w[],nvec[],1) //FIXME
   
   for(t in 2:nperiods){
      for(z in 1:nzbs){
         het[t,z] ~ normal(0,1);
      }
      target += -0.5 * dot_self(theta_S[t,node1] - theta_S[t,node2]);
      sum(theta_S[t,:]) ~ normal(0, 0.01 * N_car);
      //theta_S[t,1:nzbs]~car.normal(map[],w[],nvec[],1); //FIXME
   }
   
   beta_0 ~ uniform(-10000, 10000);
   beta_1 ~ uniform(-10000, 10000);
   beta_2 ~ uniform(-10000, 10000);
   
   //beta_3[1] ~ normal(0, 0.000001);
   for(k in 2:3){beta_3[k] ~ uniform(-10000, 10000);}
   
   sd_het ~ uniform(0,5);
   sd_spat ~ uniform(0,5);
   
   tau_d = pow(sd_d, -2);
   sd_d ~ uniform(0,5);
   for(d in 1:24){
   alpha_d[d] ~ normal(0, tau_d);
   }
   
   tau_inter = pow(sd_inter, -2);
   sd_inter ~ uniform(0,5);
   
   target += -0.5 * dot_self(inter[node1T] - inter[node2T]);
   sum(inter) ~ normal(0, 0.01 * N_carT);
   //inter[1:6]~car.normal(mapT[],wT[],nvecT[],tau.inter); //FIXME
   
   ro ~ uniform(-1,1);
}
