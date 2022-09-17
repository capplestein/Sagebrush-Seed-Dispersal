functions{
real wald_lpdf(real x, real mu, real lambda){ //defining wald pdf
	real prob;
	real lprob;
	prob = ((lambda/(2.00*pi()*(x^3.00)))^0.5)*exp(-(lambda*((x-mu)^2.00))/(2.00*(mu^2.00)*x));
	lprob  = log(prob);
	return lprob;
}
}

data {
int <lower = 0> N;      //number of observations
int numseeds[N]; //count seed trapped within 2 height bands
vector[N] nrem; //number of remnants in patch
real <lower = 0> dist[N]; //trap distance
vector[N] fecundity; //average fecundity of individuals
real <lower=0> area[N]; //area offset
real <lower=0> stddays[N]; //number of days offset
real <lower=0> httrap[N]; //height caught
int <lower=0> TR; // number of transects
int transect[N];

//for WALD kernel
vector[N] MaxWind; //wind speed
vector[N] h; //vegetation canopy height
real tvelo; //terminal velocity
real<lower=0.3,upper=0.4> kappa; //canopy density
real <lower = 0> mu_wald[N]; //pre-calculated mu_wald
real <lower = 0> lambda[N]; //pre-calculated lambda
real <lower = 0> sigma2[N]; //pre-calculated sigmna
}

parameters {				
real<lower=0> p; //parameter for empirical dispersal kernel         
real<lower=0> u; //parameter for empirical dispersal kernel
real<lower=0> f1;
real<lower=0.00001> phi;
real<lower=0> sigma_transect;// transect level variance
real<lower=0> sigma2_transect; //transect level variance 2
real omega[TR]; //deviation between p and each transect p
real gamma[TR]; //deviation between u and each transect u
real<lower=0,upper=50> distest[N];   //estimated ground distance from WALD kernal
}

transformed parameters{
vector[N] disp; //dispersal kernel
vector[N] mu; //mean seed count
vector[N] fecund; //total seed availability
real p_group[TR]; //transect-level p parameter
real u_group[TR]; //transect-level u parameter
real <lower=0> dist2[N]; //latent estimated ground distance 

//calculate individual transect slopes
for (i in 1:TR){
p_group[i] = p + omega[i]*sigma_transect;
u_group[i] = u + gamma[i]*sigma2_transect;
}

for(s in 1:N){
if(httrap[s] > 0.2 && mu_wald[s] > 0 && lambda[s] >0){
dist2[s] = dist[s] + distest[s]; //latent ground distance from wald kernel
}
else{
dist2[s] = dist[s]; //if height <= 0.2m, est distance is measured distance
}
}

for (g in 1:N) {
disp[g] = p_group[transect[g]]*(1/(pi()*u_group[transect[g]] * 
(1 + ((dist2[g]^2)/u_group[transect[g]]))^(p_group[transect[g]]+1)));
fecund[g] = f1*(fecundity[g]/1000)*nrem[g];
mu[g] = area[g] * disp[g] * fecund[g] * stddays[g];
}
}

model{
p ~ normal(0,1);
u ~ normal(0,1);
f1 ~ normal(0,1);
sigma_transect ~ exponential(4);
sigma2_transect ~ exponential(4);
phi ~ exponential(1);
omega ~ normal(0,0.5);
gamma ~ normal(0,0.5);

for(n in 1:N){
if (mu_wald[n] > 0 && lambda[n] >0) {
target += wald_lpdf(distest[n]|mu_wald[n],lambda[n]);
}
numseeds[n] ~ neg_binomial_2(mu[n], phi);
}
}

generated quantities{
vector[N] log_lik;

for (n in 1:N){
log_lik[n] = neg_binomial_2_lpmf(numseeds[n]|mu[n], phi);
}
}
