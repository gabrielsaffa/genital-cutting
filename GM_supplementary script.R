################################################################################################################
### R script for: Global phylogenetic analysis reveals multiple origins and functions of genital mutilations ###
################################################################################################################

# load the packages

library (car)        # assessing multicollinearity among candidate variables
library (missMDA)    # imputing missing values necessary for correspondence analysis only (do not confuse with missForest imputation below) 
library (FactoMineR) # correspondence analysis
library (factoextra) # visualizing the outputs of correspondence analysis
library (missForest) # imputation of missing values
library (caper)      # estimation of the phylogenetic signal with the method of Fritz & Purvis (2010)
library (phytools)   # marginal ancestral state reconstructions (Yang et al. 1995) & correlated evolution analysis (Pagel 1999)
library (phylolm)    # phylogenetic logistic regression (PLR; Ives & Garland 2010)
library (rr2)        # calculation of the explained variance (R2) of the PLR models (Ives 2019)
library (MuMIn)      # Akaike weights of the PLR models 
library (ggplot2)    # creating coefficient plots
library (ggpubr)     # creating figure panels


#######################################################
### Data preparation & predictor variable selection ###
#######################################################


setwd ("")

GMSCCS <- read.csv ("GM-SCCS_original.csv", header=TRUE, row.names=1) # load the SCCS data with original coding

# because the EA does not contain as many variables as the SCCS, variables included in the EA were those that are included in both samples (except Polygyny variables and High Gods); while the SCCS is a subset of the EA in terms of societies, EA is a subset of SCCS in terms of predictor variables


### consider variables category by category (see categories in Supplementary Table 1) ###

# marriage system

vif (glm(MGM~MaritalComposition+DistanceCowives+CulturalBasisPolygyny+Polygamy+StandardPolygamyCode,data=GMSCCS,family="binomial")) # remove StandardPolygamyCode 
vif (glm(MGM~MaritalComposition+DistanceCowives+CulturalBasisPolygyny+Polygamy,data=GMSCCS,family="binomial")) # remove Polygamy 
vif (glm(MGM~MaritalComposition+DistanceCowives+CulturalBasisPolygyny,data=GMSCCS,family="binomial")) # values <3

gmsccsimp <- imputeCA (GMSCCS, ncp=2, threshold=1e-08, maxiter=1000) # need to impute missing values before proceeding to correspondence analysis (CA); imputed values have no influence on CA
gmsccsca <- CA (gmsccsimp[,8:10], graph=TRUE) # run CA on marriage variables
summary (gmsccsca) 
gmsccsca$col$contrib # asses the contribution of each variable

# marriage practices

vif (glm(MGM~MarriageArrangementsM+MarriageArrangementsF+TransactionsMarriage,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,13:15], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# sexual behaviour

vif (glm(MGM~ExtramaritalSex+PremaritalSexNorms+SegrAdolBoys+PremaritalSexAttitudesF+MaleSegregation,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,16:20], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# rituals

vif (glm(MGM~ScarificationM+ScarificationF,data=GMSCCS,family="binomial")) # values <3

# can't run CA on a 2x2 contingency table; however, both variables are important and would be included in the study

# residence

vif (glm(MGM~TransferResidence+MaritalResidence+MaritalResidencePP,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,23:25], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# social stratification

vif (glm(MGM~JurisHierarchy+ClassDiff+CasteDiff+SocialStratification,data=GMSCCS,family="binomial")) # remove SocialStratification
vif (glm(MGM~JurisHierarchy+ClassDiff+CasteDiff,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,26:28], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# subsistence

vif (glm(MGM~SubsistenceEconomy+PrincipalSubsistenceCategory+PrimarySubsistence+DomesticAnimals,data=GMSCCS,family="binomial")) # remove PrincipalSubsistenceCategory
vif (glm(MGM~SubsistenceEconomy+PrimarySubsistence+DomesticAnimals,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,30:32], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# warfare

vif (glm(MGM~ExternalWarfare+InternalWarfare+OverallWarfare,data=GMSCCS,family="binomial")) # remove OverallWarfare
vif (glm(MGM~ExternalWarfare+InternalWarfare,data=GMSCCS,family="binomial")) # values <3

# can't run CA on a 2x2 contingency table; however, both variables are important and would be included in the study

# descent

vif (glm(MGM~Descent+DescentRules+DescentCorporateGroups,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,37:39], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib

# religion

vif (glm(MGM~DeepIslam+HighGods+WorldReligions,data=GMSCCS,family="binomial")) # values <3

gmsccsca <- CA (gmsccsimp[,40:42], graph=TRUE)
summary (gmsccsca) 
gmsccsca$col$contrib


# after a set of predictor variables was selected, many predictors were re-coded and/or split into several binary characters (see Supplementary Table 1) and so VIFs were assessed again to avoid potential problems with model fitting

### VIFs of the recoded data ###


# SCCS

vif (glm(MGM~CulturalBasisPolygyny+DistanceCowives,data=GMSCCS,family="binomial")) # ok
vif (glm(MGM~VirginityInsistence+SegrAdolBoys+PremSexAttF,data=GMSCCS,family="binomial")) # ok
vif (glm(MGM~ScarificationM+ScarificationF,data=GMSCCS,family="binomial")) # ok
vif (glm(MGM~ExternalWarfare+InternalWarfare,data=GMSCCS,family="binomial")) # ok
vif (glm(MGM~Chiefdoms+States+ClassDiff+CasteDiff,data=GMSCCS,family="binomial")) # ok
vif (glm(MGM~Foraging+Fishing+Pastoralism+ExtAgric+IntAgric,data=GMSCCS,family="binomial")) # VIFs >3 but still <10, all retained (see Methods)
vif (glm(MGM~MarrArrM+MarrArrF+BridePrice,data=GMSCCS,family="binomial")) # ok

# EA

vif (glm(MGM~PolygynyMild+PolygynyTrue+ResidDistance,data=GMEA,family="binomial")) # ok
vif (glm(MGM~PremaritalSexNorms+VirginityInsistence+SegrAdolBoys,data=GMEA,family="binomial")) # ok
vif (glm(MGM~Chiefdoms+States+ClassDiff+CasteDiff,data=GMEA,family="binomial")) # ok
vif (glm(MGM~Foraging+Fishing+Pastoralism+ExtAgric+IntAgric,data=GMEA,family="binomial")) # VIFs >3 but still <10, all retained (see Methods)


####################################
### Imputation of missing values ###
####################################


# compare the 'performance' of imputing the original and re-coded data by comparing their 'imputation errors' (see Stekhoven & Bühlmann 2012)

# SCCS-original

GMSCCS_o <- read.csv ("GM-SCCS_original.csv", header=TRUE, row.names=1) # load the original data with the selected predictors; don't forget to set societies' names (1st column) as row names, otherwise imputation will fail

cols <- c(1:27)
GMSCCS_o[,cols] <- data.frame (apply(GMSCCS_o[cols], 2, as.factor)) # convert to factors

GMSCCSimp_o <- missForest (GMSCCS_o, variablewise=TRUE) # run the imputaiton function; argument 'variablewise' needs to be TRUE to return imputation error for each variable separately

GMSCCSimp_o$OOBerror # check the imputation error of each imputed variable

GMSCCSimp_o <- as.data.frame (GMSCCSimp_o$OOBerror)
GMSCCSimp_o <- GMSCCSimp_o [-c(1:8,25:27,32:33),] # remove values with zero imputation error, i.e. predictors with no missing values that were included only to increase the estimation accuracy of imputed values, otherwise means will be biased

# SCCS-recoded

GMSCCS_r <- read.csv ("GM-SCCS_recoded.csv", header=TRUE, row.names=1) # load the recoded data and repeat the same procedure as above

cols <- c(1:32)
GMSCCS_r[,cols] <- data.frame (apply(GMSCCS_r[cols], 2, as.factor))

GMSCCSimp_r <- missForest (GMSCCS_r, variablewise=TRUE) 

GMSCCSimp_r$OOBerror

GMSCCSimp_r <- as.data.frame (GMSCCSimp_r$OOBerror) 
GMSCCSimp_r <- GMSCCSimp_r [-c(1:8,18,27),]
GMSCCSimp_r

# means of imputation errors of the original and recoded data

mean (GMSCCSimp_o) 
mean (GMSCCSimp_r)

t.test (GMSCCSimp_o, GMSCCSimp_r, paired=FALSE) # the difference is significant (it was also significant for the EA). We choose the recoded datasets for imputation.


### imputation of the recoded data ###

GMSCCS <- read.csv ("GM-SCCS_recoded.csv", header=TRUE, row.names=1) 

cols <- c(1:31)
GMSCCS[,cols] <- data.frame (apply(GMSCCS[cols], 2, as.factor)) 

GMSCCSimp <- missForest (GMSCCS, variablewise=TRUE) 

write.csv (GMSCCSimp$ximp, file="GM-SCCS-dataset_missForest.csv") # save the imputed values

# EA

GMEA <- read.csv ("GM-EA_recoded.csv", header=TRUE, row.names=1)

cols <- c(1:26)
GMEA[,cols] <- data.frame (apply(GMEA[cols], 2, as.factor))

GMEAimp <- missForest (GMEA, variablewise=TRUE)

write.csv (GMEAimp$ximp, file="GM-EA-dataset_missForest.csv") # save the imputed values


#############################################
### Estimation of the phylogenetic signal ###
#############################################


GMdata <- read.csv ("GM-SCCS_imputed.csv", header=TRUE) # load the data; in this case without assigning the first column to the data ('row.names=1'), otherwise caper's function 'comparative.data' will not be able to match societies in the data and the tree
GMtree <- read.tree ("SPT.SCCS.tre") # load the phylogeny

# caper package requires to combine phylogeny and the data into one object

GM_SCCS <- comparative.data (phy=GMtree, 
                             data=GMdata, 
                             names.col="Population", 
                             vcv=TRUE, 
                             na.omit=FALSE, 
                             warn.dropped=TRUE) 

GMphylo <- phylo.d (GM_SCCS, binvar=MGM, permut=1000) # estimate the phylogenetic signal 
GMphylo

# repeat for each GM variable/sample (see Supplementary Table 2)


##########################################
### Reconstruction of ancestral states ###
##########################################

# note that only EA sample was used for this type of analysis (see Methods)

# match tips of the phylogeny to societies in the dataset and assign it to a single variable

GMdata_phylo <- GMdata[GMtree$tip.label, ] 
MGM <- GMdata_phylo$MGM
names (MGM) <- row.names (GMdata_phylo) 

ER <- rerootingMethod (GMtree, MGM, model="ER") # run simpler 'Equal Rates' model
ER

ARD <- rerootingMethod (GMtree, MGM, model="ARD") # run more complex 'All Rates Differ' model
ARD

1-pchisq (2*abs(ER$loglik - ARD$loglik), 1) # perform likelihood ratio test; if the p-value is <0.05, choose more complex (ARD) model

# ARD model was a better fit for each GM type

write.csv (ARD$marginal.anc, file="MGM-EA_anc_values.csv") # save the file with the reconstructed values at the nodes


########################################
### Phylogenetic logistic regression ###
########################################

# run the full model with all the variables included

MGM.SCCS.plogreg <- phyloglm (MGM~CulturalBasisPolygyny+DistanceCowives+VirginityInsistence+SegrAdolBoys+PremSexAttF+ScarificationM+ScarificationF+InternalWarfare+ExternalWarfare+Chiefdoms+States+ClassDiff+CasteDiff+Patriloc+Patrilin+Foraging+Fishing+Pastoralism+ExtAgric+IntAgric+MarrArrM+MarrArrF+BridePrice+DeepIslam,
                              phy=GMtree,
                              data=GMdata,
                              boot=100
)

summary (MGM.SCCS.plogreg) # assess the model fit (SE, CI), AIC and p-values for hypotheses tests of the estimated parameters

# continue with stepwise backward (and eventually forward) elimination procedure until a set of the three best fitting models is selected

# create bootstrapped estimates for coefficients, confidence intervals and parameter alpha by running 2000 simulations for each of the selected models

MGM.SCCS.plogreg1 <- phyloglm (MGM~DistanceCowives+PremSexAttF+ScarificationM+CasteDiff+Patriloc+ExtAgric+BridePrice,phy=GMtree,data=GMdata,boot=2000)
MGM.SCCS.plogreg2 <- phyloglm (MGM~DistanceCowives+PremSexAttF+ScarificationM+CasteDiff+Patriloc+BridePrice+DeepIslam,phy=GMtree,data=GMdata,boot=2000)
MGM.SCCS.plogreg3 <- phyloglm (MGM~DistanceCowives+PremSexAttF+ScarificationM+CasteDiff+Patriloc+BridePrice,phy=GMtree,data=GMdata,boot=2000)

# get R2 values of the best fitting models

R2.lik (MGM.SCCS.plogreg1)
R2.lik (MGM.SCCS.plogreg2)
R2.lik (MGM.SCCS.plogreg3)

# calculate Akaike model weights (see Supplementary Tables 5 & 6)

round (Weights(c(MGM.SCCS.plogreg1$aic,MGM.SCCS.plogreg2$aic,MGM.SCCS.plogreg3$aic)), 3)


### assess how much of the 'explained variance' (R2) of each selected model is accounted for by phylogeny ###

GMtree_star <- ape::compute.brlen (GMtree, method="Grafen", power=.0001) # convert the tree into a star-like phylogeny

# run a 'non-phylogenetic' model

MGM.SCCS.plogreg1_nophylo <- phyloglm (MGM~DistanceCowives+PremSexAttF+ScarificationM+CasteDiff+Patriloc+ExtAgric+BridePrice,
                                       phy=GMtree_star, # transformed phylogeny
                                       data=GMdata,
                                       boot=100 # bootstrapping doesn't have any influence on R2
)

# calculate the difference in R2 between phylogenetic and non-phylogenetic models (Supplementary Table 7)

R2.lik (MGM.SCCS.plogreg1, MGM.SCCS.plogreg1_nophylo) 

# compare R2 between phylogenetic and non-phylogenetic models (see Supplementary Table 8)

R2_all <- read.csv ("R2-all.csv", header=TRUE)         # load the file with all the R2 values
R2_FGM <- read.csv ("R2-FGM.csv", header=TRUE)         # load the file with values for all FGM models
R2_MGM <- read.csv ("R2-MGM.csv", header=TRUE)         # load the file with values for all MGM models
R2_EA <- read.csv ("R2-EA.csv", header=TRUE)           # load the file with values for all EA models
R2_SCCS <- read.csv ("R2-SCCS.csv", header=TRUE)       # load the file with values for all SCCS models
R2_FGM_MGM <- read.csv ("R2-FGM-MGM.csv", header=TRUE) # load the file with values for all FGM and all MGM models

with (R2_all, t.test (phylo, nophylo, paired=TRUE))    # run a paired t-test on all models
with (R2_FGM, t.test (phylo, nophylo, paired=TRUE))    # run a paired t-test on all FGM models
with (R2_MGM, t.test (phylo, nophylo, paired=TRUE))    # run a paired t-test on all MGM models
with (R2_EA, t.test (phylo, nophylo, paired=TRUE))     # run a paired t-test on all EA models
with (R2_SCCS, t.test (phylo, nophylo, paired=TRUE))   # run a paired t-test on all SCCS models
with (R2_FGM_MGM, t.test (MGM, FGM, paired=FALSE))     # run a two-sample t-test to compare the difference in R2 between all FGM and all MGM models 


### prepare the files with coefficient estimates and confidence intervals (ci) for plotting ###

estimate <- MGM.SCCS.plogreg1$coefficients # get coefficient estimates
estimate

bm <- as.data.frame (MGM.SCCS.plogreg1$bootmean) # get also bootstrapped mean coefficient estimates
bm

bootm <- bm[-10, ] # remove parameter alpha
bootm

df <- as.data.frame (t(MGM.SCCS.plogreg1$bootconfint95)) # convert bootstrapped ci into a data frame and use the transpose function to convert rows into columns to bind it with coefficients

ci <- df[-10, ] # remove parameter alpha
ci

MGM1 <- cbind (estimate, bootm, ci) # combine all columns (coefficients, bootstrapped coefficients and bootstrapped ci) into one object

model1 <- MGM1[-1, ] # now remove intercept
model1

MGM_SCCS_all <- rbind (model1, model2, model3) # do the same with each selected model for a given GM type and bind them all together

names (MGM_SCCS_all)[names(MGM_SCCS_all) == "2.5%"] <- "lci" # replace with a more convenient name
names (MGM_SCCS_all)[names(MGM_SCCS_all) == "97.5%"] <- "uci" # same

write.csv (MGM_SCCS_all, file="MGM-SCCS_coefs.csv") # save the file with coefficients and ci of all models

# repeat the same procedure for each dependent variable


###################################################
### Correspondence analysis of model predictors ###
###################################################


# see Methods for details of the procedure
# assess if the associations between rows (GMs) and columns (model predictors) are due to chance; files were used for creating Fig. 4

GM_models_ea <- read.csv ("GM-EA_predictors.csv", header=TRUE, row.names=1)
GM_models_sccs <- read.csv ("GM-SCCS_predictors.csv", header=TRUE, row.names=1)

gmeaca <- CA (GM_models_ea, graph=FALSE) # run correspondence analysis
summary (gmeaca)

chisq.test (gmeaca) 
fisher.test (gmeaca, simulate.p.value=2000) # run Fisher's exact test to obtain p-value as all cells in the contingency table have counts <5

# do the same for the SCCS sample

gmsccsca <- CA (GM_models_sccs, graph=FALSE)
summary (gmsccsca)

chisq.test (gmsccsca)
fisher.test (gmsccsca, simulate.p.value=2000)


#####################################################
### Analysis of correlated evolution (Pagel 1994) ###
#####################################################


GMdata_phylo <- GMdata[GMtree$tip.label, ] # match the data with the phylogeny

# assign row names to the selected variables

MGM <- GMdata_phylo$MGM 
names (MGM) <- row.names (GMdata_phylo)

Castes <- GMdata_phylo$CasteDiff
names (Castes) <- row.names (GMdata_phylo)

# ARD model used according to ancestral state reconstructions
# perform all types of test for each pair of selected variables

fit1 <- fitPagel (GMtree, x=MGM, y=Castes, model="ARD", dep.var="xy") # transition rate in 'x' depends on 'y' and vice versa
fit1 

fit2 <- fitPagel (GMtree, x=MGM, y=Castes, model="ARD", dep.var="x") # transition rate in 'x' depends on 'y' but NOT vice versa
fit2

fit3 <- fitPagel (GMtree, x=MGM, y=Castes, model="ARD", dep.var="y") # transition rate in 'y' depends on 'x' but NOT vice versa
fit3

# models 'fit2' and 'fit3' are not nested so can't assess their fit with a likelihood ratio test, simply choose the model with higher likelihood
# then compare the fit of either 'fit2' or 'fit3' to 'fit1' with a likelihood ratio test


1-pchisq (2*abs(fit1$dependent.logL - fit2$dependent.logL), 2) 

plot.fitPagel (fit2, lwd.by.rate=TRUE) # plot the selected model

# repeat for each pair of selected variables 


################################################
### Plotting ancestral state reconstructions ###
################################################

# plot the phylogeny

par(mfrow=c(1,1))

plotTree (GMtree, 
          setEnv=TRUE, 
          offset=0.5, 
          fsize=0.2, 
          lwd=1
          )

# add values at the tips

tiplabels (pie=to.matrix(MGM, sort(unique(MGM))), 
           piecol = c("gray65", "red1"), 
           cex = 0.1
           )

# add reconstructed values at internal nodes from the ARD model

nodelabels (node=as.numeric(rownames(ARD$marginal.anc)), 
            pie=ARD$marginal.anc,
            piecol = c("gray65", "red1"), 
            cex=0.1
            )

# add a simple legend

add.simmap.legend (x=0,
                   y=180,
                   colors=setNames(c("gray65", "red1"),
                                   c("Absent", "Present")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.6
                   )


#########################
### Coefficient plots ###
#########################

# load all the files with coefficient and ci estimates saved previously

MGM_SCCS <- read.csv ("MGM-SCCS_coefs.csv", header=TRUE)
FGM_SCCS <- read.csv ("FGM-SCCS_coefs.csv", header=TRUE)
Circumcision_SCCS <- read.csv ("Cir-SCCS_coefs.csv", header=TRUE)
Superincision_SCCS <- read.csv ("Sup-SCCS_coefs.csv", header=TRUE)
Excision_SCCS <- read.csv ("Exc-SCCS_coefs.csv", header=TRUE)
Infibulation_SCCS <- read.csv ("Inf-SCCS_coefs.csv", header=TRUE)
Clitoridectomy_SCCS <- read.csv ("Clit-SCCS_coefs.csv", header=TRUE)

# create coefficient plot with point estimates and error bars represented by boostrapped mean estimate and bootstrapped confidence intervals, respectively

MGM <- ggplot(MGM_SCCS, 
              aes(x=predictor,
                  y=bootm,
                  fill=factor(model),
                  shape=factor(model))
              ) + geom_point(position=position_dodge(-0.65), 
             size=2, 
             color="red1"
             ) + scale_fill_discrete(name="predictor")+ geom_errorbar(aes(ymin=lci, 
                    ymax=uci),
                position=position_dodge(-0.65),
                width=.5, 
                color="red1"
                ) + ggtitle("e") + geom_hline(yintercept=0, 
             linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank()) + theme(legend.position="none") + theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank()) 
+ theme(plot.title=element_text(size=12,
                                face="bold")
        ) + scale_x_discrete(limits=c("Islam", 
                            "Bride-price", 
                            "Extensive agr.", 
                            "Patrilocality",
                            "Castes", 
                            "Scarifications-M", 
                            "Sex attitudes-F",
                            "DistanceCowives")
)

# repeat the same for each dependent variable stored as separate object and plot in one panel

ggarrange(FGM, 
          Clitoridectomy, 
          Excision, 
          Infibulation, 
          MGM, 
          Circumcision, 
          Superincision, 
          ncol=4, 
          nrow=3, 
          legend="bottom",  
          common.legend=TRUE
) 


############################################################
### Plots of correspondence analysis of model predictors ###
############################################################


GM_models_ea <- read.csv ("GM-EA_predictors.csv", header=TRUE, row.names=1)
GM_models_sccs <- read.csv ("GM-SCCS_predictors.csv", header=TRUE, row.names=1)

gmeaca <- CA (GM_models_ea, graph=FALSE)
gmsccsca <- CA (GM_models_sccs, graph=FALSE)

ea <- fviz_ca_biplot(gmeaca, 
                     map="rowprincipal", 
                     arrow=c(FALSE, FALSE), 
                     repel=TRUE, 
                     col.col="gray65", 
                     col.row="#FFCC00", 
                     geom=c("text", 
                            "point")
                     ) + labs(title="a") + theme_classic(base_size = 12) + theme(legend.position="none") + theme(panel.background=element_rect(fill="white", 
                                      colour="black")
)

sccs <- fviz_ca_biplot(gmsccsca, 
                       map="rowprincipal", 
                       arrow=c(FALSE, FALSE), 
                       repel=TRUE, 
                       col.col="gray65", 
                       col.row="#FFCC00", 
                       geom=c("text", 
                             "point")
                       ) + labs(title="b") + theme_classic(base_size=12) + theme(legend.position="none") + theme(panel.background=element_rect(fill="white", 
                                      colour="black")
)

# combine in one panel

ggarrange(ea, 
          sccs, 
          ncol=2, 
          nrow=2, 
          common.legend=TRUE, 
          legend="bottom"
          )

# plots were further adjusted in a graphical software


#####################################################
#####################################################
#####################################################

