library(brms)
library(tidyr)
library(plyr)
library(dplyr)
library(loo)
library(bayesplot)
library(ggplot2)
library(gridExtra)

load("SeedTrapData_STAN.Rdata")

##STAN data

##create a dataframe from seedtrapdat
datsum <-data.frame(trapID = as.factor(seedtrapdat$trapID),
                patchID = as.factor(seedtrapdat$patchID),
                site = as.factor(seedtrapdat$site),
                transect = as.factor(seedtrapdat$transect),
                httrap = seedtrapdat$httrap,
                dist = seedtrapdat$dist,
                fecundity = seedtrapdat$fecundity,
                numseeds = seedtrapdat$numseeds,
                area = seedtrapdat$area,
                nrem = seedtrapdat$nrem,
                stddays = seedtrapdat$stddays,
                year = seedtrapdat$year,
                winddir = seedtrapdat$winddir,
                windbin = seedtrapdat$windbin,
                siteyr = seedtrapdat$siteyr)

##total seed availability as a combo of ind fecundity * nrem
datsum$fecun <- datsum$fecundity*datsum$nrem

##scale data
scaledat <-datsum
scaledat$httrap <-scale(scaledat$httrap)
scaledat$dist <-scale(scaledat$dist)
scaledat$fecun <-scale(scaledat$fecun)
scaledat$winddir <- scale(scaledat$winddir)

##No random effects
conts = list(max_treedepth = 13)
model1 <-brm(numseeds ~ httrap*dist + fecun +  
          offset(log(stddays)) + 
          offset(log(area)), 
          data = scaledat, family = negbinomial(),
          iter = 2000, cores = 4, chains = 4,
          control = conts)
model1 <-add_criterion(model1, c("loo", "waic"), reloo = TRUE,
                       seed = TRUE)

save(model1, file = "Model1brms_NoRE.Rdata")

##site random effect
conts = list(max_treedepth = 13, adapt_delta = 0.96)
model2 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|siteyr) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 2000, cores = 4, chains = 4,
              control = conts)
model2 <-add_criterion(model2, c("loo", "waic"), reloo = TRUE)
save(model2, file = "Model2brms_Site.Rdata")

##site x patch random effect
model3 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|siteyr/patchID) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4, 
              control = conts)
model3 <-add_criterion(model3, c("loo", "waic"), reloo = TRUE)
save(model3, file = "Model3brms_SitePatch.Rdata")

##patch only
model4 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|patchID) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4,
              control = conts)
model4 <-add_criterion(model4, c("loo", "waic"), reloo = TRUE)
save(model4, file = "Model4brms_Patch.Rdata")

##patch and transect
model5 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|patchID/transect) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4,
              control = conts)
model5 <-add_criterion(model5, c("loo", "waic"), reloo = TRUE)
save(model5, file = "Model5brms_PatchTrans.Rdata")

##transect only
model6 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|transect) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4, 
              control = conts)
model6 <-add_criterion(model6, c("loo", "waic"), reloo = TRUE)
save(model6, file = "Model6brms_Trans.Rdata")

##all nested effects only
model7 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|siteyr/patchID/transect) + 
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4, 
              control = conts)
model7 <-add_criterion(model7, c("loo", "waic"), reloo = TRUE)
save(model7, file = "Model7brms_SitePatchTrans.Rdata")

##model6 with wind direction
##all nested effects only
model8 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|transect) + winddir +
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4, 
              control = conts)
model8 <-add_criterion(model8, c("loo", "waic"), reloo = TRUE)
save(model8, file = "Model8brms_TransWindDir.Rdata")

##wind bin
model9 <- brm(numseeds ~ httrap*dist + fecun +  
                (httrap*dist|transect) + windbin +
                offset(log(stddays)) + 
                offset(log(area)), 
              data = scaledat, family = negbinomial(),
              iter = 3000, cores = 4, 
              control = conts)
model9 <-add_criterion(model9, c("loo", "waic"), reloo = TRUE)
save(model9, file = "Model9brms_TransWindbin.Rdata")
