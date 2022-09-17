##fit STAN model
library(rstan)
rstan_options(auto_write = TRUE)
library(pbmcapply)
library(loo)
library(matrixStats)
library(Metrics)

conts <- list(adapt_delta = 0.95, max_treedepth = 14)
load("SeedTrapData_STAN.Rdata")

##This data file is organized to run in STAN models
##Full data with metadata is available at https://www.fs.usda.gov/rds/archive/Catalog/RDS-2021-0073
##in the SeedTrapData_STAN.Rdata file, seed counts are aggregated
##across 2 height bands on a trap and the associated height of capture is the mid-point

##calculations for mu and lambda prior to model run:
##set parameters
kappa <-0.38
seedtrapdat$sigma2 <- kappa*seedtrapdat$h*(2*(1/2))

##terminal velocities
sagetv <-c(2.1189591,0.1944729,0.2548055,0.4086022,0.2024867,1.6332378,0.5084746)
tvelo <- median(sagetv)

##calculating lambda and mu parameters for WALD kernel
seedtrapdat$lambda <-(seedtrapdat$h/sqrt(seedtrapdat$sigma2))^2
seedtrapdat$mu <-(seedtrapdat$h*seedtrapdat$MaxWind)/tvelo

##Run STAN model
modelsim_fit <-stan(file ="Transect_LatentWALD_SimModel.stan", 
        control = conts,
        data = seedtrapdat, 
        iter = 5000, warmup = 500, 
        cores = 4, chains = 4, seed = 786, 
      pars = c("p", "u", "disp", "phi", "fecund",
      "f1", "p_group", "u_group", 
      "mu", "log_lik", "preds", "omega",
      "gamma", "sigma_transect", "sigma2_transect",
      "distest", "dist2", "sim_mu", "sim_disp"))

save(modelsim_fit, file = "ModelSim_TransLatentWALD.Rdata")

##look at model parameters
parms <-data.frame(extract(modelsim_fit,
                           par = c("p", "u", "disp", "phi", "fecund",
                                   "f1", "p_group", "u_group", 
                                   "mu", "log_lik", "preds", "omega",
                                   "gamma", "sigma_transect", "sigma2_transect",
                                   "distest", "dist2", "sim_mu", "sim_disp"))) 

# Function for simulating seed rain based on model posteriors
gen_quant_r <- function(x) {
  n <-nrow(x)
  p <-sample(parms$p, size = 1)
  u <-sample(parms$u, size = 1)
  disp <- p*(1/(pi*u *(1 +(x$dist^2/u))^(p+1)))
  fecund <-sample(parms$f1, size = 1)*(x$fecundity)*x$nrem
  mu <-x$area * disp * fecund * x$stddays
  phi <-sample(parms$phi, size = 1)
  predictions <-rnbinom(n, mu = mu, size = phi)
  names(predictions) <-x$dist
  return(predictions)
}

dist = seq(from=0.01,to=100,by=0.5)
newdat <-data.frame(dist = dist,
                    fecundity = rep(30000, 200),
                    stddays = rep(2,200), nrem = rep(25, 200), 
                    area = 1)
set.seed(35)
pred <- replicate(10000, {
  preds <-gen_quant_r(newdat)
})
newdat2 <-data.frame(rowQuantiles(pred, probs = c(0.05,0.5,0.95)))
newdat2 <-cbind(newdat2, dist)
