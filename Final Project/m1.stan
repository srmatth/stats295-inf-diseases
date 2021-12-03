
data {
   int<lower=0> n;
   int<lower=0> y[n];
   int<lower=0> N[n];
   
   int<lower=0> N_car;
   int<lower=0> N_edges;
   int node1[N_edges];
   int node2[N_edges];
}

parameters {
   real m;
   real het[n];
   vector[n] sp;
   real<lower=0> sdhet;
   real<lower=0> sdsp;
   
   //vector[N_car] phi;
   real sigma;
}

transformed parameters {
   vector<lower=0, upper = 1>[n] p;
   for(i in 1:n){
      p[i] = inv_logit(m+het[i]+sp[i]);
   }
}

model{

 real prechet;
 real precsp;
 
 prechet=pow(sdhet,-2);
 precsp=pow(sdsp,-2);

 for(i in 1:n)
	{
 		het[i]~normal(0,prechet);
 		
 		y[i]~binomial(N[i],p[i]);
	}
	
 target += -0.5 * dot_self(sp[node1] - sp[node2]);
 sum(sp) ~ normal(0, 0.01 * N_car);
 //sp[1:n]~car.normal(adj[],w[],num[],precsp);

 m~uniform(-10000, 10000);
 sdhet~uniform(0,10);
 sdsp~uniform(0,10);
}
