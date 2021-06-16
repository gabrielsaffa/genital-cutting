################################################################################################################
### R script for: Global phylogenetic analysis reveals multiple origins and functions of genital mutilations ###
################################################################################################################

### phylogenetic signal, ancestral state reconstructions and correlated evolution ###

library (geiger)
library (caper)   
library (phytools) 

setwd ("")


##########
### EA ###

GM_EA_data <- read.csv ("GM_EA_imputed.csv", header=TRUE, row.names=1)
GM_EA_tree <- read.tree ("SPT.EA.tre")

name.check (GM_EA_tree,GM_EA_data) # check which populations are in the tree but not in the data and remove them; we remove those populations that have missing data for any of the outcomes

GM_EA_tree <- drop.tip (GM_EA_tree, tip=c("Abron","Amis","Aweti","Batanga","Bunun","Burarra","Bwaidoka","Creek","Dobuan","Guaymi",              
                                          "Ho","Huastec","Kachari","Kaqchikel","Karbi","Kariyarra","Kim_Mun","Koronadal_Blaan",
                                          "Kui","Kunda","Lele","Mazatec","Mbala_Ivukuna","Mbugwe","Mixtec","Miyako","Molima",
                                          "Murrinh_Patha","Naga_Southern_Rengma","Ndobo","Nung_China","Opata","Paiwan","Puyuma",
                                          "Regi","Sanga","Selaru","Shuwa","Sumbwa","Taabwa","Tepehuan","Tzetzal","Waropen","Wayana",
                                          "Yaeyama","Yami"))

name.check (GM_EA_tree,GM_EA_data) # OK
write.tree (GM_EA_tree, "SPT.EA.tre") # write the reduced tree


###########################
### phylogenetic signal ###

# FGM

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=FGM, permut=1000)
GM_EA_phylo_sig

# clitoridectomy

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=clit, permut=1000)
GM_EA_phylo_sig

# excision

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=exc, permut=1000)
GM_EA_phylo_sig

# infibulation 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=inf, permut=1000) 
GM_EA_phylo_sig

# MGM 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=MGM, permut=1000) 
GM_EA_phylo_sig

# circumcision 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=cir, permut=1000) 
GM_EA_phylo_sig

# superincision 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=sup, permut=1000) 
GM_EA_phylo_sig

# co-wives separate 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=dist, permut=1000) 
GM_EA_phylo_sig

# sex norms 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=sex_norms, permut=1000) 
GM_EA_phylo_sig

# patrilocality 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=patriloc, permut=1000) 
GM_EA_phylo_sig

# patrilineality 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=patrilin, permut=1000) 
GM_EA_phylo_sig

# pastoralism 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=past, permut=1000) 
GM_EA_phylo_sig

# bride-price 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=bride_pr, permut=1000) 
GM_EA_phylo_sig

# chiefdoms 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=chiefdoms, permut=1000) 
GM_EA_phylo_sig

# classes 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=class, permut=1000) 
GM_EA_phylo_sig

# castes 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=caste, permut=1000) 
GM_EA_phylo_sig

# high gods 

GM_EA_phylo_sig <- phylo.d (GM_EA_data, GM_EA_tree, names.col=pop, binvar=high_gods, permut=1000) 
GM_EA_phylo_sig


#########################################################################
### ancestral state reconstruction using stochastic character mapping ###

### NOTE: only the EA sample was used for this type of analysis


### FGM ###

FGM <- setNames (GM_EA_data[,3], rownames(GM_EA_data))
sim_FGM <- make.simmap (GM_EA_tree, FGM, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_FGM_pd <- summary (sim_FGM)

cols <- setNames (c("gray85","blue"), levels(as.factor(FGM)))

plot (sim_FGM_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("blue", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### clitoridectomy ###

clit <- setNames (GM_EA_data[,4], rownames(GM_EA_data))
sim_clit <- make.simmap (GM_EA_tree, clit, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_clit_pd <- summary (sim_clit)

cols <- setNames (c("gray85","turquoise2"), levels(as.factor(clit)))

plot (sim_clit_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("turquoise2", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### excision ###

exc <- setNames (GM_EA_data[,5], rownames(GM_EA_data))
sim_exc <- make.simmap (GM_EA_tree, exc, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_exc_pd <- summary (sim_exc)

cols <- setNames (c("gray85","turquoise4"), levels(as.factor(exc)))

plot (sim_exc_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("turquoise4", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### infibulation ###

inf <- setNames (GM_EA_data[,6], rownames(GM_EA_data))
sim_inf <- make.simmap (GM_EA_tree, inf, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_inf_pd <- summary (sim_inf)

cols <- setNames (c("gray85","dodgerblue4"), levels(as.factor(inf)))

plot (sim_inf_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("dodgerblue4", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### MGM ###

MGM <- setNames (GM_EA_data[,7], rownames(GM_EA_data))
sim_MGM <- make.simmap (GM_EA_tree, MGM, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_MGM_pd <- summary (sim_MGM)

cols <- setNames (c("gray85","red1"), levels(as.factor(MGM)))

plot (sim_MGM_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("red1", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### circumcision ###

cir <- setNames (GM_EA_data[,8], rownames(GM_EA_data))
sim_cir <- make.simmap (GM_EA_tree, cir, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_cir_pd <- summary (sim_cir)

cols <- setNames (c("gray85","gold"), levels(as.factor(cir)))

plot (sim_cir_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("gold", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


### superincision ###

sup <- setNames (GM_EA_data[,9], rownames(GM_EA_data))
sim_sup <- make.simmap (GM_EA_tree, sup, model="ARD", Q="empirical", nsim=1000, pi=c(1,0))
sim_sup_pd <- summary (sim_sup)

cols <- setNames (c("gray85","red4"), levels(as.factor(sup)))

plot (sim_sup_pd,
      offset=1.0,
      fsize=0.2,
      lwd=1,
      cex=c(0.1,0.1),
      colors=cols
)

add.simmap.legend (x=0,
                   y=570,
                   colors=setNames(c("red4", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)


############################
### correlated evolution ###


##########
### EA ###

GM_EA_data_phylo <- GM_EA_data[GM_EA_tree$tip.label, ] # match the data with the phylogeny

# assign row names to the selected variables

### FGM ~ co-wives separate ###

FGM <- GM_EA_data_phylo$FGM
names (FGM) <- row.names (GM_EA_data_phylo) 

dist <- GM_EA_data_phylo$dist
names (dist) <- row.names (GM_EA_data_phylo)

# ARD model was used based on ancestral state reconstructions
# perform all three types of test for each pair of selected variables

fit1 <- fitPagel (GM_EA_tree, x=FGM, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=FGM, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=FGM, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE) # plot the selected model(s)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ sex norms ###

sex_norms <- GM_EA_data_phylo$sex_norms
names (sex_norms) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=FGM, y=sex_norms, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=FGM, y=sex_norms, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=FGM, y=sex_norms, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ MGM ###

MGM <- GM_EA_data_phylo$MGM
names (MGM) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=FGM, y=MGM, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=FGM, y=MGM, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=FGM, y=MGM, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ bride-price ###

bride_pr <- GM_EA_data_phylo$bride_pr
names (bride_pr) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=FGM, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=FGM, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=FGM, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### FGM ~ castes ###

castes <- GM_EA_data_phylo$caste
names (castes) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=FGM, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=FGM, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=FGM, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### clit ~ dist ###

clit <- GM_EA_data_phylo$clit
names (clit) <- row.names (GM_EA_data_phylo)

dist <- GM_EA_data_phylo$dist
names (dist) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=clit, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=clit, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=clit, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### clit ~ sex norms ###

sex_norms <- GM_EA_data_phylo$sex_norms
names (sex_norms) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=clit, y=sex_norms, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=clit, y=sex_norms, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=clit, y=sex_norms, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### clit ~ bride-price ###

bride_pr <- GM_EA_data_phylo$bride_pr
names (bride_pr) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=clit, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=clit, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=clit, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### clit ~ classes ###

class <- GM_EA_data_phylo$class
names (class) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=clit, y=class, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=clit, y=class, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=clit, y=class, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### exc ~ dist ###

exc <- GM_EA_data_phylo$exc
names (exc) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=exc, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=exc, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=exc, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### exc ~ patriloc ###

patriloc <- GM_EA_data_phylo$patriloc
names (patriloc) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=exc, y=patriloc, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=exc, y=patriloc, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=exc, y=patriloc, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### exc ~ castes ###

fit1 <- fitPagel (GM_EA_tree, x=exc, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=exc, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=exc, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### clit ~ classes ###

fit1 <- fitPagel (GM_EA_tree, x=exc, y=class, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=exc, y=class, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=exc, y=class, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### inf ~ dist ###

inf <- GM_EA_data_phylo$inf
names (inf) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=inf, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=inf, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=inf, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### inf ~ patrilin ###

fit1 <- fitPagel (GM_EA_tree, x=inf, y=patrilin, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=inf, y=patrilin, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=inf, y=patrilin, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### inf ~ pastoralism ###

past <- GM_EA_data_phylo$past
names (past) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=inf, y=past, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=inf, y=past, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=inf, y=past, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### inf ~ bride-price ###

bride_pr <- GM_EA_data_phylo$bride_pr
names (bride_pr) <- row.names (GM_EA_data_phylo) 

fit1 <- fitPagel (GM_EA_tree, x=inf, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=inf, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=inf, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### inf ~ class ###

class <- GM_EA_data_phylo$class
names (class) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=inf, y=class, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=inf, y=class, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=inf, y=class, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### MGM ~ dist ###

MGM <- GM_EA_data_phylo$MGM
names (MGM) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=MGM, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=MGM, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=MGM, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### MGM ~ patrilin ###

fit1 <- fitPagel (GM_EA_tree, x=MGM, y=patrilin, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=MGM, y=patrilin, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=MGM, y=patrilin, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### MGM ~ castes ###

fit1 <- fitPagel (GM_EA_tree, x=MGM, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=MGM, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=MGM, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### MGM ~ high gods ###

high_gods <- GM_EA_data_phylo$high_gods
names (high_gods) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=MGM, y=high_gods, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=MGM, y=high_gods, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=MGM, y=high_gods, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### cir ~ dist ###

cir <- GM_EA_data_phylo$cir
names (cir) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=cir, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=cir, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=cir, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### cir ~ castes ###

fit1 <- fitPagel (GM_EA_tree, x=cir, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=cir, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=cir, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### cir ~ high gods ###

fit1 <- fitPagel (GM_EA_tree, x=cir, y=high_gods, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=cir, y=high_gods, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=cir, y=high_gods, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### sup ~ dist ###

sup <- GM_EA_data_phylo$sup
names (sup) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=sup, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=sup, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=sup, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### sup ~ male segregation ###

segr_adol <- GM_EA_data_phylo$segr_adol
names (segr_adol) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=sup, y=segr_adol, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=sup, y=segr_adol, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=sup, y=segr_adol, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### sup ~ chiefdoms ###

chiefdoms <- GM_EA_data_phylo$chiefdoms
names (chiefdoms) <- row.names (GM_EA_data_phylo)

fit1 <- fitPagel (GM_EA_tree, x=sup, y=chiefdoms, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=sup, y=chiefdoms, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=sup, y=chiefdoms, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### sup ~ high gods ###

fit1 <- fitPagel (GM_EA_tree, x=sup, y=high_gods, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_EA_tree, x=sup, y=high_gods, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_EA_tree, x=sup, y=high_gods, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


############
### SCCS ###


###########################
### phylogenetic signal ###

GM_SCCS_data <- read.csv ("GM_SCCS_imputed.csv", header=TRUE) 
GM_SCCS_tree <- read.tree ("SPT.SCCS.tre")


# FGM

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=FGM, permut=1000)
GM_SCCS_phylo_sig

# clitoridectomy

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=clit, permut=1000)
GM_SCCS_phylo_sig

# excision

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=exc, permut=1000)
GM_SCCS_phylo_sig

# infibulation

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=inf, permut=1000) 
GM_SCCS_phylo_sig

# MGM

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=MGM, permut=1000) 
GM_SCCS_phylo_sig

# circumcision

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=cir, permut=1000) 
GM_SCCS_phylo_sig

# superincision

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=sup, permut=1000) 
GM_SCCS_phylo_sig

# co-wives separate

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=dist, permut=1000) 
GM_SCCS_phylo_sig

# patrilocality 

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=patriloc, permut=1000) 
GM_SCCS_phylo_sig

# patrilineality 

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=patrilin, permut=1000) 
GM_SCCS_phylo_sig

# pastoralism

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=past, permut=1000) 
GM_SCCS_phylo_sig

# bride-price

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=bride_pr, permut=1000) 
GM_SCCS_phylo_sig

# chiefdoms

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=chiefdoms, permut=1000) 
GM_SCCS_phylo_sig

# classes

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=class, permut=1000) 
GM_SCCS_phylo_sig

# castes

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=caste, permut=1000) 
GM_SCCS_phylo_sig

# female scarification

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=scars_f, permut=1000)
GM_SCCS_phylo_sig

# male scarification

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=scars_m, permut=1000)
GM_SCCS_phylo_sig

# Islam

GM_SCCS_phylo_sig <- phylo.d (GM_SCCS_data, GM_SCCS_tree, names.col=pop, binvar=islam, permut=1000) 
GM_SCCS_phylo_sig


############################
### correlated evolution ###

GM_SCCS_data$dist_bin <- ifelse (GM_SCCS_data$dist==4,1,0) # first make a binary version of co-wives separate
GM_SCCS_data_phylo <- GM_SCCS_data[GM_SCCS_tree$tip.label, ] # match the data with the phylogeny


### FGM ~ co-wives separate ###

FGM <- GM_SCCS_data_phylo$FGM
names (FGM) <- row.names (GM_SCCS_data_phylo) 

dist <- GM_SCCS_data_phylo$dist_bin
names (dist) <- row.names (GM_SCCS_data_phylo)


fit1 <- fitPagel (GM_SCCS_tree, x=FGM, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=FGM, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=FGM, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ MGM ###

MGM <- GM_SCCS_data_phylo$MGM
names (MGM) <- row.names (GM_SCCS_data_phylo) 


fit1 <- fitPagel (GM_SCCS_tree, x=FGM, y=MGM, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=FGM, y=MGM, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=FGM, y=MGM, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)


### FGM ~ bride-price ###

bride_pr <- GM_SCCS_data_phylo$bride_pr
names (bride_pr) <- row.names (GM_SCCS_data_phylo) 


fit1 <- fitPagel (GM_SCCS_tree, x=FGM, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=FGM, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=FGM, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ castes ###

castes <- GM_SCCS_data_phylo$caste
names (castes) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=FGM, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=FGM, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=FGM, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### FGM ~ scars-female ###

scars_f <- GM_SCCS_data_phylo$scars_f
names (scars_f) <- row.names (GM_SCCS_data_phylo) 


fit1 <- fitPagel (GM_SCCS_tree, x=FGM, y=scars_f, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=FGM, y=scars_f, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=FGM, y=scars_f, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)


### clit ~ co-wives separate ###

clit <- GM_SCCS_data_phylo$clit
names (clit) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=clit, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=clit, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=clit, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### clit ~ Islam ###

islam <- GM_SCCS_data_phylo$islam
names (islam) <- row.names (GM_SCCS_data_phylo) 


fit1 <- fitPagel (GM_SCCS_tree, x=clit, y=islam, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=clit, y=islam, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=clit, y=islam, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### clit ~ bride_pr ###

fit1 <- fitPagel (GM_SCCS_tree, x=clit, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=clit, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=clit, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### clit ~ classes ###

fit1 <- fitPagel (GM_SCCS_tree, x=clit, y=class, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=clit, y=class, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=clit, y=class, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### exc ~ dist ###

exc <- GM_SCCS_data_phylo$exc
names (exc) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=exc, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=exc, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=exc, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### exc ~ patriloc ###

patriloc <- GM_SCCS_data_phylo$patriloc
names (patriloc) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=exc, y=patriloc, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=exc, y=patriloc, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=exc, y=patriloc, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### exc ~ castes ###

castes <- GM_SCCS_data_phylo$caste
names (castes) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=exc, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=exc, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=exc, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### exc ~ scars-f ###

fit1 <- fitPagel (GM_SCCS_tree, x=exc, y=scars_f, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=exc, y=scars_f, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=exc, y=scars_f, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### clit ~ classes ###

fit1 <- fitPagel (GM_SCCS_tree, x=exc, y=class, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=exc, y=class, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=exc, y=class, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### inf ~ dist ###

inf <- GM_SCCS_data_phylo$inf
names (inf) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=inf, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=inf, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=inf, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### inf ~ patrilin ###

fit1 <- fitPagel (GM_SCCS_tree, x=inf, y=patrilin, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=inf, y=patrilin, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=inf, y=patrilin, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### inf ~ past ###

past <- GM_SCCS_data_phylo$past
names (past) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=inf, y=past, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=inf, y=past, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=inf, y=past, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)


### inf ~ bride-price ###

bride_pr <- GM_SCCS_data_phylo$bride_pr
names (bride_pr) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=inf, y=bride_pr, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=inf, y=bride_pr, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=inf, y=bride_pr, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)


### inf ~ classes ###

classes <- GM_SCCS_data_phylo$class
names (classes) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=inf, y=classes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=inf, y=classes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=inf, y=classes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### MGM ~ dist ###

MGM <- GM_SCCS_data_phylo$MGM
names (MGM) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=MGM, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=MGM, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=MGM, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### MGM ~ pattrilin ###

fit1 <- fitPagel (GM_SCCS_tree, x=MGM, y=patrilin, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=MGM, y=patrilin, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=MGM, y=patrilin, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)


### MGM ~ castes ###

fit1 <- fitPagel (GM_SCCS_tree, x=MGM, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=MGM, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=MGM, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### MGM ~ scars-m ###

scars_m <- GM_SCCS_data_phylo$scars_m
names (scars_m) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=MGM, y=scars_m, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=MGM, y=scars_m, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=MGM, y=scars_m, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### cir ~ dist ###

cir <- GM_SCCS_data_phylo$cir
names (cir) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=cir, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=cir, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=cir, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### cir ~ castes ###

fit1 <- fitPagel (GM_SCCS_tree, x=cir, y=castes, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=cir, y=castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=cir, y=castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit2, lwd.by.rate=TRUE)
plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


### cir ~ Islam ###

fit1 <- fitPagel (GM_SCCS_tree, x=cir, y=islam, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=cir, y=islam, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=cir, y=islam, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### sup ~ dist ###

sup <- GM_SCCS_data_phylo$sup
names (sup) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=sup, y=dist, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=sup, y=dist, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=sup, y=dist, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### sup ~ male segregation ###

segr_adol <- GM_SCCS_data_phylo$segr_adol
names (segr_adol) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=sup, y=segr_adol, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=sup, y=segr_adol, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=sup, y=segr_adol, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit1, lwd.by.rate=TRUE)
plot.fitPagel (fit2, lwd.by.rate=TRUE)


### sup ~ chiefdoms ###

chiefdoms <- GM_SCCS_data_phylo$chiefdoms
names (chiefdoms) <- row.names (GM_SCCS_data_phylo) 

fit1 <- fitPagel (GM_SCCS_tree, x=sup, y=chiefdoms, model="ARD") # transition rate in 'x' depends on 'y' and vice versa
fit2 <- fitPagel (GM_SCCS_tree, x=sup, y=chiefdoms, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit3 <- fitPagel (GM_SCCS_tree, x=sup, y=chiefdoms, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa


### get likelihoods, AIC, delta AIC and AIC weights

model <- c("independent","x dependent","y dependent","both dependent")

aic <- round(c(fit1$independent.AIC, fit2$dependent.AIC, fit3$dependent.AIC, fit1$dependent.AIC),2) # get AIC
ll <- round(c(fit1$independent.logL, fit2$dependent.logL, fit3$dependent.logL, fit1$dependent.logL),2) # get log-likelihood values
w_aic <- round (aic.w (aic),2) # get AIC weights

# get delta AIC 

d_aic <- 0

for (i in 1:length (aic))
{
  aic_diff <- (aic[i] - min(aic))
  d_aic[i] <- round (aic_diff,2)
}

# create a data frame with all computed values and order it by AIC

aic_all <- data.frame (cbind (model,ll,aic,d_aic,w_aic))
aic_all <- aic_all[order(aic),] 
aic_all

plot.fitPagel (fit3, lwd.by.rate=TRUE)
plot.fitPagel (fit1, lwd.by.rate=TRUE)


######################################
######################################