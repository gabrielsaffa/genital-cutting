###############################################################################################################################
### R script  for: Global phylogenetic analysis reveals multiple origins and correlates of genital cutting ###
###############################################################################################################################

### fitting phylogenetic logistic regression ###

library (phylolm)
library (MuMIn)

setwd ("")


##########
### EA ###

### load the EA data and the supertree

GM_EA_data <- read.csv ("GM_EA_imputed.csv", header=TRUE, row.names=1)
GM_EA_tree <- read.tree ("SPT.EA.tre")


### FGC ###
### NOTE: some models gave a warning about increasing 'btol' (bound on the parameter space), which, when increased, in some cases produced unrealistically high estimates and standard errors. Therefore the warning was ignored, which, after the discussion with one of the package authors, should be OK.

{
  
  # run a global model (i.e. all predictors of FGC considered)
  
  FGC_mg <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  # run an intercept-only model
  
  FGC_m0 <- phyloglm (FGC ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  # rum a model based on the mate guarding/paternity uncertainty hypothesis
  
  FGC_m1 <- phyloglm (FGC ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  # run a model based on the fraternal interest groups hypothesis
  
  FGC_m2 <- phyloglm (FGC ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  # run a model based on the male mating variance hypothesis
  
  FGC_m3 <- phyloglm (FGC ~ class,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  # fit model with pastoralism based on the mate guarding/paternity uncertainty hypothesis
  
  FGC_m4 <- phyloglm (FGC ~ past,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  # fit model with bride-price based on the market value hypothesis
  
  FGC_m5 <- phyloglm (FGC ~ bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  # we fit our own model
  
  FGC_m6 <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6")
aic <- round(c(FGC_mg$aic,FGC_m0$aic,FGC_m1$aic,FGC_m2$aic,FGC_m3$aic,FGC_m4$aic,FGC_m5$aic,FGC_m6$aic),2) # get AIC
ll <- round(c(FGC_mg$logLik,FGC_m0$logLik,FGC_m1$logLik,FGC_m2$logLik,FGC_m3$logLik,FGC_m4$logLik,FGC_m5$logLik,FGC_m6$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(FGC_mg$aic,FGC_m0$aic,FGC_m1$aic,FGC_m2$aic,FGC_m3$aic,FGC_m4$aic,FGC_m5$aic,FGC_m6$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

FGC_m4 <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)

FGC_mg <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- FGC_m4$coefficients # coefficients
sd <- FGC_m4$sd # SD
b_est <- as.data.frame (FGC_m4$bootmean) # bootstrapped coefficients
b_est <- b_est[-9, ] # remove alpha
b_ci <- as.data.frame (t(FGC_m4$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-9, ] # remove alpha
fgc_m4 <- cbind (est, sd, b_est, b_ci)
fgc_m4 <- fgc_m4[-1, ] # remove intercept

# repeat for the other model

est <- FGC_mg$coefficients 
sd <- FGC_mg$sd # SD
b_est <- as.data.frame (FGC_mg$bootmean) 
b_est <- b_est[-9, ] 
b_ci <- as.data.frame (t(FGC_mg$bootconfint95)) 
b_ci <- b_ci[-9, ] 
fgc_mg <- cbind (est, sd, b_est, b_ci) 
fgc_mg <- fgc_mg[-1, ] 
fgc_all <- round(rbind (fgc_m4, fgc_mg),2) # combine

write.csv (fgc_all, file="fgc_all_EA.csv") # store the values


### clitoridectomy ###

{
  
  clit_mg <- phyloglm (clit ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0,
                       log.alpha.bound=5
  )
  
  
  clit_m0 <- phyloglm (clit ~ 1, phy=GM_EA_tree, data=GM_EA_data, log.alpha.bound=5)
  
  
  clit_m1 <- phyloglm (clit ~ dist,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0,
                       log.alpha.bound=5
  )
  
  
  clit_m2 <- phyloglm (clit ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0,
                       log.alpha.bound=5
  )
  
  
  clit_m3 <- phyloglm (clit ~ class,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0,
                       log.alpha.bound=5
  )
  
  
  clit_m4 <- phyloglm (clit ~ past,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  clit_m5 <- phyloglm (clit ~ bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0,
                      log.alpha.bound=5
  )
  
  
  clit_m6 <- phyloglm (clit ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0,
                       log.alpha.bound=5
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6") 
aic <- round(c(clit_mg$aic,clit_m0$aic,clit_m1$aic,clit_m2$aic,clit_m3$aic,clit_m4$aic,clit_m5$aic,clit_m6$aic),2) # get AIC
ll <- round(c(clit_mg$logLik,clit_m0$logLik,clit_m1$logLik,clit_m2$logLik,clit_m3$logLik,clit_m4$logLik,clit_m5$logLik,clit_m6$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(clit_mg$aic,clit_m0$aic,clit_m1$aic,clit_m2$aic,clit_m3$aic,clit_m4$aic,clit_m5$aic,clit_m6$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

clit_m6 <- phyloglm (clit ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                     phy=GM_EA_tree,
                     data=GM_EA_data,
                     boot=2000,
                     log.alpha.bound=5
)


### get bootstrapped values

est <- clit_m6$coefficients # coefficients
sd <- clit_m6$sd
b_est <- as.data.frame (clit_m6$bootmean) # bootstrapped coefficients
b_est <- b_est[-9, ] # remove alpha
b_ci <- as.data.frame (t(clit_m6$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-9, ] # remove alpha
cl_m6 <- cbind (est, sd, b_est, b_ci)
cl_m6 <- cl_m6[-1, ] # remove intercept
cl_all <- round(rbind (cl_m6),2) # combine

write.csv (cl_all, file="cl_all_EA.csv") # store the values


### excision ###

{
  
  exc_mg <- phyloglm (exc ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  exc_m0 <- phyloglm (exc ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  exc_m1 <- phyloglm (exc ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  exc_m2 <- phyloglm (exc ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  exc_m3 <- phyloglm (exc ~ class,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  exc_m4 <- phyloglm (exc ~ past,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0
  )
  
  
  exc_m5 <- phyloglm (exc ~ bride_pr,
                       phy=GM_EA_tree,
                       data=GM_EA_data,
                       boot=0
  )
  
  
  exc_m6 <- phyloglm (exc ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6") 
aic <- round(c(exc_mg$aic,exc_m0$aic,exc_m1$aic,exc_m2$aic,exc_m3$aic,exc_m4$aic,exc_m5$aic,exc_m6$aic),2) # get AIC
ll <- round(c(exc_mg$logLik,exc_m0$logLik,exc_m1$logLik,exc_m2$logLik,exc_m3$logLik,exc_m4$logLik,exc_m5$logLik,exc_m6$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(exc_mg$aic,exc_m0$aic,exc_m1$aic,exc_m2$aic,exc_m3$aic,exc_m4$aic,exc_m5$aic,exc_m6$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

exc_mg <- phyloglm (exc ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- exc_mg$coefficients # coefficients
sd <- exc_mg$sd # SD
b_est <- as.data.frame (exc_mg$bootmean) # bootstrapped coefficients
b_est <- b_est[-12, ] # remove alpha
b_ci <- as.data.frame (t(exc_mg$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-12, ] # remove alpha
ex_mg <- cbind (est, sd, b_est, b_ci)
ex_mg <- ex_mg[-1, ] # remove intercept
ex_all <- round(rbind (ex_mg),2) # combine

write.csv (ex_all, file="ex_all_EA.csv") # store the values


### infibulation ###

{
  
  inf_mg <- phyloglm (inf ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      
  )
  
  
  inf_m0 <- phyloglm (inf ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  inf_m1 <- phyloglm (inf ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      
  )
  
  
  inf_m2 <- phyloglm (inf ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      
  )
  
  
  inf_m3 <- phyloglm (inf ~ class,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      
  )
  
  
  inf_m4 <- phyloglm (inf ~ past,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  inf_m5 <- phyloglm (inf ~ bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  inf_m6 <- phyloglm (inf ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6")
aic <- round(c(inf_mg$aic,inf_m0$aic,inf_m1$aic,inf_m2$aic,inf_m3$aic,inf_m4$aic,inf_m5$aic,inf_m6$aic),2) # get AIC
ll <- round(c(inf_mg$logLik,inf_m0$logLik,inf_m1$logLik,inf_m2$logLik,inf_m3$logLik,inf_m4$logLik,inf_m5$logLik,inf_m6$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(inf_mg$aic,inf_m0$aic,inf_m1$aic,inf_m2$aic,inf_m3$aic,inf_m4$aic,inf_m5$aic,inf_m6$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

inf_m3 <- phyloglm (inf ~ class,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)

inf_m4 <- phyloglm (inf ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- inf_m3$coefficients # coefficients
sd <- inf_m3$sd
b_est <- as.data.frame (inf_m3$bootmean) # bootstrapped coefficients
b_est <- b_est[-3, ] # remove alpha
b_ci <- as.data.frame (t(inf_m3$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-3, ] # remove alpha
infib_m3 <- cbind (est, sd, b_est, b_ci)
infib_m3 <- infib_m3[-1, ] # remove intercept

# repeat for the other model

est <- inf_m4$coefficients 
sd <- inf_m4$sd
b_est <- as.data.frame (inf_m4$bootmean) 
b_est <- b_est[-9, ] 
b_ci <- as.data.frame (t(inf_m4$bootconfint95)) 
b_ci <- b_ci[-9, ] 
infib_m4 <- cbind (est, sd, b_est, b_ci) 
infib_m4 <- infib_m4[-1, ] 
infib_all <- round(rbind (infib_m3, infib_m4),2) # combine

write.csv (infib_all, file="infib_all_EA.csv") # store the values


### MGC ###

{
  
  MGC_mg <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )

  
  MGC_m0 <- phyloglm (MGC ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  MGC_m1 <- phyloglm (MGC ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  MGC_m2 <- phyloglm (MGC ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  MGC_m3 <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  MGC_m4 <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4")
aic <- round(c(MGC_mg$aic,MGC_m0$aic,MGC_m1$aic,MGC_m2$aic,MGC_m3$aic,MGC_m4$aic),2) # get AIC
ll <- round(c(MGC_mg$logLik,MGC_m0$logLik,MGC_m1$logLik,MGC_m2$logLik,MGC_m3$logLik,MGC_m4$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(MGC_mg$aic,MGC_m0$aic,MGC_m1$aic,MGC_m2$aic,MGC_m3$aic,MGC_m4$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

MGC_mg <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- MGC_mg$coefficients 
sd <- MGC_mg$sd
b_est <- as.data.frame (MGC_mg$bootmean) 
b_est <- b_est[-14, ] 
b_ci <- as.data.frame (t(MGC_mg$bootconfint95)) 
b_ci <- b_ci[-14, ] 
mgc_mg <- cbind (est, sd, b_est, b_ci) 
mgc_mg <- mgc_mg[-1, ] 
mgc_mg <- round(rbind (mgc_mg),2)

write.csv (mgc_mg, file="mgc_mg_EA.csv") # store the values


### circumcision ###

{
  
  cir_mg <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )

  
  cir_m0 <- phyloglm (cir ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  cir_m1 <- phyloglm (cir ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  cir_m2 <- phyloglm (cir ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  cir_m3 <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  cir_m4 <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4")
aic <- round(c(cir_mg$aic,cir_m0$aic,cir_m1$aic,cir_m2$aic,cir_m3$aic,cir_m4$aic),2) # get AIC
ll <- round(c(cir_mg$logLik,cir_m0$logLik,cir_m1$logLik,cir_m2$logLik,cir_m3$logLik,cir_m4$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(cir_mg$aic,cir_m0$aic,cir_m1$aic,cir_m2$aic,cir_m3$aic,cir_m4$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

cir_mg <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- cir_mg$coefficients 
sd <- cir_mg$sd
b_est <- as.data.frame (cir_mg$bootmean) 
b_est <- b_est[-14, ] 
b_ci <- as.data.frame (t(cir_mg$bootconfint95)) 
b_ci <- b_ci[-14, ] 
circum_mg <- cbind (est, sd, b_est, b_ci) 
circum_mg <- circum_mg[-1, ] 
circum_mg <- round(rbind (circum_mg),2)

write.csv (circum_mg, file="circum_mg_EA.csv") # store the values 


### superincision ###

{
  
  sup_mg <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  sup_m0 <- phyloglm (sup ~ 1, phy=GM_EA_tree, data=GM_EA_data)
  
  
  sup_m1 <- phyloglm (sup ~ dist,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  sup_m2 <- phyloglm (sup ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  sup_m3 <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
  
  sup_m4 <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                      phy=GM_EA_tree,
                      data=GM_EA_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4")
aic <- round(c(sup_mg$aic,sup_m0$aic,sup_m1$aic,sup_m2$aic,sup_m3$aic,sup_m4$aic),2) # get AIC
ll <- round(c(sup_mg$logLik,sup_m0$logLik,sup_m1$logLik,sup_m2$logLik,sup_m3$logLik,sup_m4$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(sup_mg$aic,sup_m0$aic,sup_m1$aic,sup_m2$aic,sup_m3$aic,sup_m4$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

sup_m4 <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)

sup_mg <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=2000
)


### get bootstrapped values

est <- sup_m4$coefficients # coefficients
sd <- sup_m4$sd
b_est <- as.data.frame (sup_m4$bootmean) # bootstrapped coefficients
b_est <- b_est[-9, ] # remove alpha
b_ci <- as.data.frame (t(sup_m4$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-9, ] # remove alpha
super_m4 <- cbind (est, sd, b_est, b_ci)
super_m4 <- super_m4[-1, ] # remove intercept

# repeat for the other model

est <- sup_mg$coefficients 
sd <- sup_mg$sd
b_est <- as.data.frame (sup_mg$bootmean) 
b_est <- b_est[-14, ] 
b_ci <- as.data.frame (t(sup_mg$bootconfint95)) 
b_ci <- b_ci[-14, ] 
super_mg <- cbind (est, sd, b_est, b_ci) 
super_mg <- super_mg[-1, ] 
super_all <- round(rbind (super_m4, super_mg),2) # combine

write.csv (super_all, file="super_all_EA.csv") # store the values


############
### SCCS ###

GM_SCCS_data <- read.csv ("GM_SCCS_imputed.csv", header=TRUE, row.names=1)
GM_SCCS_tree <- read.tree ("SPT.SCCS.tre")


### FGC ###

{
  
  FGC_mg <- phyloglm (FGC ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )

  
  FGC_m0 <- phyloglm (FGC ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  FGC_m1 <- phyloglm (FGC ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  FGC_m2 <- phyloglm (FGC ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  FGC_m3 <- phyloglm (FGC ~ class,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  FGC_m4 <- phyloglm (FGC ~ past,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  FGC_m5 <- phyloglm (FGC ~ bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  # here we fit also a model with female scarifications
  
  FGC_m6 <- phyloglm (FGC ~ scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )

  
  FGC_m7 <- phyloglm (FGC ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6","m7")
aic <- round(c(FGC_mg$aic,FGC_m0$aic,FGC_m1$aic,FGC_m2$aic,FGC_m3$aic,FGC_m4$aic,FGC_m5$aic,FGC_m6$aic,FGC_m7$aic),2) # get AIC
ll <- round(c(FGC_mg$logLik,FGC_m0$logLik,FGC_m1$logLik,FGC_m2$logLik,FGC_m3$logLik,FGC_m4$logLik,FGC_m5$logLik,FGC_m6$logLik,FGC_m7$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(FGC_mg$aic,FGC_m0$aic,FGC_m1$aic,FGC_m2$aic,FGC_m3$aic,FGC_m4$aic,FGC_m5$aic,FGC_m6$aic,FGC_m7$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

FGC_m5 <- phyloglm (FGC ~ bride_pr,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)

FGC_m7 <- phyloglm (FGC ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)


### get bootstrapped values

est <- FGC_m5$coefficients # coefficients
sd <- FGC_m5$sd # SD
b_est <- as.data.frame (FGC_m5$bootmean) # bootstrapped coefficients
b_est <- b_est[-3, ] # remove alpha
b_ci <- as.data.frame (t(FGC_m5$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-3, ] # remove alpha
fgc_m5 <- cbind (est, sd, b_est, b_ci)
fgc_m5 <- fgc_m5[-1, ] # remove intercept

#

est <- FGC_m7$coefficients 
sd <- FGC_m7$sd # SD
b_est <- as.data.frame (FGC_m7$bootmean) 
b_est <- b_est[-10, ] 
b_ci <- as.data.frame (t(FGC_m7$bootconfint95)) 
b_ci <- b_ci[-10, ] 
fgc_m7 <- cbind (est, sd, b_est, b_ci) 
fgc_m7 <- fgc_m7[-1, ] 
fgc_all <- round(rbind (fgc_m5,fgc_m7),2) # combine

write.csv (fgm_all, file="fgm_all_SCCS.csv") # store the values


### clitoridectomy ###

{
  
  clit_mg <- phyloglm (clit ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f + islam,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0
  )

  
  clit_m0 <- phyloglm (clit ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  clit_m1 <- phyloglm (clit ~ dist + ext_aff,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0, 
                       log.alpha.bound=10,
                       btol=20
  )
  
  
  clit_m2 <- phyloglm (clit ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0
  )
  
  
  clit_m3 <- phyloglm (clit ~ class,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0
  )
  
  
  clit_m4 <- phyloglm (clit ~ past,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0
  )
  
  
  clit_m5 <- phyloglm (clit ~ bride_pr,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0,
                       log.alpha.bound=10
  )
  
  
  clit_m6 <- phyloglm (clit ~ scars_f,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0,
                       log.alpha.bound=10
  )
  
  
  # we fit also a model with Islam
  
  clit_m7 <- phyloglm (clit ~ islam,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0,
                       log.alpha.bound=10
  )

  
  clit_m8 <- phyloglm (clit ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste + islam,
                       phy=GM_SCCS_tree,
                       data=GM_SCCS_data,
                       boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6","m7","m8")
aic <- round(c(clit_mg$aic,clit_m0$aic,clit_m1$aic,clit_m2$aic,clit_m3$aic,clit_m4$aic,clit_m5$aic,clit_m6$aic,clit_m7$aic,clit_m8$aic),2) # get AIC
ll <- round(c(clit_mg$logLik,clit_m0$logLik,clit_m1$logLik,clit_m2$logLik,clit_m3$logLik,clit_m4$logLik,clit_m5$logLik,clit_m6$logLik,clit_m7$logLik,clit_m8$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(clit_mg$aic,clit_m0$aic,clit_m1$aic,clit_m2$aic,clit_m3$aic,clit_m4$aic,clit_m5$aic,clit_m6$aic,clit_m7$aic,clit_m8$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

clit_m1 <- phyloglm (clit ~ dist + ext_aff,
                     phy=GM_SCCS_tree,
                     data=GM_SCCS_data,
                     boot=2000,
                     btol=20,
                     log.alpha.bound=10
)

clit_m5 <- phyloglm (clit ~ bride_pr,
                     phy=GM_SCCS_tree,
                     data=GM_SCCS_data,
                     boot=2000,
                     log.alpha.bound=10
)

clit_m7 <- phyloglm (clit ~ islam,
                     phy=GM_SCCS_tree,
                     data=GM_SCCS_data,
                     boot=2000,
                     log.alpha.bound=10
)


### get bootstrapped values

est <- clit_m1$coefficients # coefficients
sd <- clit_m1$sd
b_est <- as.data.frame (clit_m1$bootmean) # bootstrapped coefficients
b_est <- b_est[-4, ] # remove alpha
b_ci <- as.data.frame (t(clit_m1$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-4, ] # remove alpha
cl_m1 <- cbind (est, sd, b_est, b_ci)
cl_m1 <- cl_m1[-1, ] # remove intercept

#

est <- clit_m5$coefficients 
sd <- clit_m5$sd
b_est <- as.data.frame (clit_m5$bootmean) 
b_est <- b_est[-3, ] 
b_ci <- as.data.frame (t(clit_m5$bootconfint95)) 
b_ci <- b_ci[-3, ] 
cl_m5 <- cbind (est, sd, b_est, b_ci) 
cl_m5 <- cl_m5[-1, ] 

#

est <- clit_m7$coefficients 
sd <- clit_m7$sd
b_est <- as.data.frame (clit_m7$bootmean) 
b_est <- b_est[-3, ] 
b_ci <- as.data.frame (t(clit_m7$bootconfint95)) 
b_ci <- b_ci[-3, ] 
cl_m7 <- cbind (est, sd, b_est, b_ci) 
cl_m7 <- cl_m7[-1, ] 
cl_all <- round(rbind (cl_m1,cl_m5,cl_m7),2) # combine

write.csv (cl_all, file="cl_all_SCCS.csv") # store the values


### excision ###

{
  
  exc_mg <- phyloglm (exc ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      btol=20
  )
  
  
  exc_m0 <- phyloglm (exc ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  exc_m1 <- phyloglm (exc ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0, 
                      log.alpha.bound=10
  )
  
  
  exc_m2 <- phyloglm (exc ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  exc_m3 <- phyloglm (exc ~ class,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )

  
  exc_m4 <- phyloglm (exc ~ past,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  exc_m5 <- phyloglm (exc ~ bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      log.alpha.bound=10
  )
  
  
  exc_m6 <- phyloglm (exc ~ scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      log.alpha.bound=10
  )
  
  
  exc_m7 <- phyloglm (exc ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      btol=20
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6","m7")
aic <- round(c(exc_mg$aic,exc_m0$aic,exc_m1$aic,exc_m2$aic,exc_m3$aic,exc_m4$aic,exc_m5$aic,exc_m6$aic,exc_m7$aic),2) # get AIC
ll <- round(c(exc_mg$logLik,exc_m0$logLik,exc_m1$logLik,exc_m2$logLik,exc_m3$logLik,exc_m4$logLik,exc_m5$logLik,exc_m6$logLik,exc_m7$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(exc_mg$aic,exc_m0$aic,exc_m1$aic,exc_m2$aic,exc_m3$aic,exc_m4$aic,exc_m5$aic,exc_m6$aic,exc_m7$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

exc_mg <- phyloglm (exc ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000,
                    btol=20
)


### get bootstrapped values

est <- exc_mg$coefficients # coefficients
sd <- exc_mg$sd
b_est <- as.data.frame (exc_mg$bootmean) # bootstrapped coefficients
b_est <- b_est[-14, ] # remove alpha
b_ci <- as.data.frame (t(exc_mg$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-14, ] # remove alpha
ex_mg <- cbind (est, sd, b_est, b_ci)
ex_mg <- ex_mg[-1, ] # remove intercept
ex_all <- round(rbind (ex_mg),2) # combine

write.csv (ex_all, file="ex_all_SCCS.csv") # store the values


### infibulation ###

{
  
  inf_mg <- phyloglm (inf ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      btol=20
  )
  
  
  inf_m0 <- phyloglm (inf ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  inf_m1 <- phyloglm (inf ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  inf_m2 <- phyloglm (inf ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  inf_m3 <- phyloglm (inf ~ class,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  inf_m4 <- phyloglm (inf ~ past,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  inf_m5 <- phyloglm (inf ~ bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      log.alpha.bound=10
  )
  
  
  inf_m6 <- phyloglm (inf ~ scars_f,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      log.alpha.bound=10
  )
  
  
  inf_m7 <- phyloglm (inf ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6","m7")
aic <- round(c(inf_mg$aic,inf_m0$aic,inf_m1$aic,inf_m2$aic,inf_m3$aic,inf_m4$aic,inf_m5$aic,inf_m6$aic,inf_m7$aic),2) # get AIC
ll <- round(c(inf_mg$logLik,inf_m0$logLik,inf_m1$logLik,inf_m2$logLik,inf_m3$logLik,inf_m4$logLik,inf_m5$logLik,inf_m6$logLik,inf_m7$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(inf_mg$aic,inf_m0$aic,inf_m1$aic,inf_m2$aic,inf_m3$aic,inf_m4$aic,inf_m5$aic,inf_m6$aic,inf_m7$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

inf_m5 <- phyloglm (inf ~ bride_pr,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)

inf_m1 <- phyloglm (inf ~ dist + ext_aff,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)

inf_m4 <- phyloglm (inf ~ past,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)

inf_m7 <- phyloglm (inf ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)


### get bootstrapped values

est <- inf_m5$coefficients # coefficients
sd <- inf_m5$sd
b_est <- as.data.frame (inf_m5$bootmean) # bootstrapped coefficients
b_est <- b_est[-3, ] # remove alpha
b_ci <- as.data.frame (t(inf_m5$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-3, ] # remove alpha
infib_m5 <- cbind (est, sd, b_est, b_ci)
infib_m5 <- infib_m5[-1, ] # remove intercept

#

est <- inf_m1$coefficients # coefficients
sd <- inf_m1$sd
b_est <- as.data.frame (inf_m1$bootmean) # bootstrapped coefficients
b_est <- b_est[-4, ] # remove alpha
b_ci <- as.data.frame (t(inf_m1$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-4, ] # remove alpha
infib_m1 <- cbind (est, sd, b_est, b_ci)
infib_m1 <- infib_m1[-1, ] # remove intercept

#

est <- inf_m4$coefficients 
sd <- inf_m4$sd
b_est <- as.data.frame (inf_m4$bootmean) 
b_est <- b_est[-3, ] 
b_ci <- as.data.frame (t(inf_m4$bootconfint95)) 
b_ci <- b_ci[-3, ] 
infib_m4 <- cbind (est, sd, b_est, b_ci) 
infib_m4 <- infib_m4[-1, ] 

#

est <- inf_m7$coefficients 
sd <- inf_m7$sd
b_est <- as.data.frame (inf_m7$bootmean) 
b_est <- b_est[-10, ] 
b_ci <- as.data.frame (t(inf_m7$bootconfint95)) 
b_ci <- b_ci[-10, ] 
infib_m7 <- cbind (est, sd, b_est, b_ci) 
infib_m7 <- infib_m7[-1, ] 
infib_all <- round(rbind (infib_m5,infib_m1,infib_m4,infib_m7),2) # combine

write.csv (infib_all, file="infib_all_SCCS.csv") # store the values


### MGC ###

{
  
  MGC_mg <- phyloglm (MGC ~ dist + ext_aff + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + ext_war + forag + class + chiefdoms + states + caste + scars_m,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  MGC_m0 <- phyloglm (MGC ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  MGC_m1 <- phyloglm (MGC ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  MGC_m2 <- phyloglm (MGC ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  # we fit a model based on Sosis et al. (2007) and add chiefdoms and states as they suggest
  
  MGC_m3 <- phyloglm (MGC ~ ext_war + class + forag + chiefdoms + states,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  # we also fit a model with male scarifications
  
  MGC_m4 <- phyloglm (MGC ~ scars_m,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  MGC_m5 <- phyloglm (MGC ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5")
aic <- round(c(MGC_mg$aic,MGC_m0$aic,MGC_m1$aic,MGC_m2$aic,MGC_m3$aic,MGC_m4$aic,MGC_m5$aic),2) # get AIC
ll <- round(c(MGC_mg$logLik,MGC_m0$logLik,MGC_m1$logLik,MGC_m2$logLik,MGC_m3$logLik,MGC_m4$logLik,MGC_m5$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(MGC_mg$aic,MGC_m0$aic,MGC_m1$aic,MGC_m2$aic,MGC_m3$aic,MGC_m4$aic,MGC_m5$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

MGC_m5 <- phyloglm (MGC ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)


### get bootstrapped values

est <- MGC_m5$coefficients # coefficients
sd <- MGC_m5$sd
b_est <- as.data.frame (MGC_m5$bootmean) # bootstrapped coefficients
b_est <- b_est[-9, ] # remove alpha
b_ci <- as.data.frame (t(MGC_m5$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-9, ] # remove alpha
mgc_m5 <- cbind (est, sd, b_est, b_ci)
mgc_m5 <- mgc_m5[-1, ] # remove intercept
mgc_m5 <- round(rbind (mgc_m5),2)

write.csv (mgc_m5, file="mgc_m5_SCCS.csv") # store the values


### circumcision ###

{
  
  cir_mg <- phyloglm (cir ~ dist + ext_aff + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + ext_war + class + forag + chiefdoms + states + caste + scars_m + islam,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0,
                      btol=20
  )
  
  
  cir_m0 <- phyloglm (cir ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  cir_m1 <- phyloglm (cir ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  cir_m2 <- phyloglm (cir ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  cir_m3 <- phyloglm (cir ~ ext_war + class + forag + chiefdoms + states,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  # we fit a model with Islam
  
  cir_m4 <- phyloglm (cir ~ islam,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  cir_m5 <- phyloglm (cir ~ scars_m,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  cir_m6 <- phyloglm (cir ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste + islam,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5","m6")
aic <- round(c(cir_mg$aic,cir_m0$aic,cir_m1$aic,cir_m2$aic,cir_m3$aic,cir_m4$aic,cir_m5$aic,cir_m6$aic),2) # get AIC
ll <- round(c(cir_mg$logLik,cir_m0$logLik,cir_m1$logLik,cir_m2$logLik,cir_m3$logLik,cir_m4$logLik,cir_m5$logLik,cir_m6$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(cir_mg$aic,cir_m0$aic,cir_m1$aic,cir_m2$aic,cir_m3$aic,cir_m4$aic,cir_m5$aic,cir_m6$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

cir_m6 <- phyloglm (cir ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste + islam,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)


### get bootstrapped values

est <- cir_m6$coefficients # coefficients
sd <- cir_m6$sd
b_est <- as.data.frame (cir_m6$bootmean) # bootstrapped coefficients
b_est <- b_est[-10, ] # remove alpha
b_ci <- as.data.frame (t(cir_m6$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-10, ] # remove alpha
circum_m6 <- cbind (est, sd, b_est, b_ci)
circum_m6 <- circum_m6[-1, ] # remove intercept
circum_all <- round(rbind (circum_m6),2)

write.csv (circum_all, file="circum_all_SCCS.csv") # store the values


### superincision ###

{
  
  sup_mg <- phyloglm (sup ~ dist + ext_aff + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + ext_war + class + forag + chiefdoms + states + caste + scars_m,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  sup_m0 <- phyloglm (sup ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)
  
  
  sup_m1 <- phyloglm (sup ~ dist + ext_aff,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  sup_m2 <- phyloglm (sup ~ patriloc + patrilin + past + ext_agric + int_agric + bride_pr,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  sup_m3 <- phyloglm (sup ~ ext_war + class + forag + chiefdoms + states,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
  
  sup_m4 <- phyloglm (sup ~ scars_m,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )

  
  sup_m5 <- phyloglm (sup ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste,
                      phy=GM_SCCS_tree,
                      data=GM_SCCS_data,
                      boot=0
  )
  
}


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("mg","m0","m1","m2","m3","m4","m5") 
aic <- round(c(sup_mg$aic,sup_m0$aic,sup_m1$aic,sup_m2$aic,sup_m3$aic,sup_m4$aic,sup_m5$aic),2) # get AIC
ll <- round(c(sup_mg$logLik,sup_m0$logLik,sup_m1$logLik,sup_m2$logLik,sup_m3$logLik,sup_m4$logLik,sup_m5$logLik),2) # get log-likelihood values
w_aic <- round (Weights(c(sup_mg$aic,sup_m0$aic,sup_m1$aic,sup_m2$aic,sup_m3$aic,sup_m4$aic,sup_m5$aic)),2) # get AIC weights

# get delta AIC 

d_aic <- 0
for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order the data by AIC

aic_all <- data.frame (model,ll,aic,d_aic,w_aic)
aic_all <- aic_all[order(aic),] 
aic_all


### compute bootstrapped parameter estimates and confidence intervals for the model/s selected by delta AIC and AIC weight

sup_m1 <- phyloglm (sup ~ dist + ext_aff,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=2000
)


### get bootstrapped values

est <- sup_m1$coefficients # coefficients
b_est <- as.data.frame (sup_m1$bootmean) # bootstrapped coefficients
b_est <- b_est[-4, ] # remove alpha
b_ci <- as.data.frame (t(sup_m1$bootconfint95)) # bootstrapped confidence intervals
b_ci <- b_ci[-4, ] # remove alpha
super_m1 <- cbind (est, b_est, b_ci)
super_m1 <- super_m1[-1, ] # remove intercept
super_m1 <- round(rbind (super_m1),2)

write.csv (super_m1, file="super_m1_SCCS.csv") # store the values


#############################################
#############################################
