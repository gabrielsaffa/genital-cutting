#########################################################################################################################
### R script  for: Global phylogenetic analysis reveals multiple origins and correlates of genital mutilation/cutting ###
#########################################################################################################################

### computing the 'explained variance' or R2 ###

library (ape)
library (phylolm)
library (rr2)

setwd ("")


##########
### EA ###

GM_EA_data <- read.csv ("GM_EA_imputed.csv", header=TRUE, row.names=1)
GM_EA_tree <- read.tree ("SPT.EA.tre")
GM_EA_tree_star <- compute.brlen (GM_EA_tree, method="Grafen", power=0.0001) # convert the tree into a star-like phylogeny


### FGC ###
  
# fit the phylogenetic model

FGC_m4 <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

# fit the same model but with a star-like phylogeny

FGC_m4_nophy <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste, 
                          phy=GM_EA_tree_star,
                          data=GM_EA_data,
                          boot=0
)

round (R2.lik (FGC_m4),2) # compute the total 'variance explained' of the phylogenetic model (predictors + phylogeny)
round (R2.lik (FGC_m4, FGC_m4_nophy),2) # what is the amount of the 'variance explained' by the phylogeny alone?

#

FGC_mg <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)


FGC_mg_nophy <- phyloglm (FGC ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (FGC_mg),2) 
round (R2.lik (FGC_m4, FGC_mg_nophy),2)


### clitoridectomy ###

clit_m6 <- phyloglm (clit ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                     phy=GM_EA_tree,
                     data=GM_EA_data,
                     boot=0,
                     log.alpha.bound=5
)

clit_m6_nophy <- phyloglm (clit ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                     phy=GM_EA_tree_star,
                     data=GM_EA_data,
                     boot=0,
                     log.alpha.bound=5
)

round (R2.lik (clit_m6),2) 
round (R2.lik (clit_m6, clit_m6_nophy),2)


### excision ###

exc_mg <- phyloglm (exc ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

exc_mg_nophy <- phyloglm (exc ~ dist + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (exc_mg),2) 
round (R2.lik (exc_mg, exc_mg_nophy),2)


### infibulation ###

inf_m3 <- phyloglm (inf ~ class,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

inf_m3_nophy <- phyloglm (inf ~ class,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (inf_m3),2) 
round (R2.lik (inf_m3, inf_m3_nophy),2)

#

inf_m4 <- phyloglm (inf ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

inf_m4_nophy <- phyloglm (inf ~ dist + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (inf_m4),2) 
round (R2.lik (inf_m4, inf_m4_nophy),2)


### MGC ###

MGC_mg <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

MGC_mg_nophy <- phyloglm (MGC ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (MGC_mg),2) 
round (R2.lik (MGC_mg, MGC_mg_nophy),2)


### circumcision ###

cir_mg <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

cir_mg_nophy <- phyloglm (cir ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (cir_mg),2) 
round (R2.lik (cir_mg, cir_mg_nophy),2)


### superinision ###

sup_m4 <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

sup_m4_nophy <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + chiefdoms + caste + high_gods,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (sup_m4),2) 
round (R2.lik (sup_m4, sup_m4_nophy),2)

#

sup_mg <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree,
                    data=GM_EA_data,
                    boot=0
)

sup_mg_nophy <- phyloglm (sup ~ dist + segr_adol + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + chiefdoms + states + caste + high_gods,
                    phy=GM_EA_tree_star,
                    data=GM_EA_data,
                    boot=0
)

round (R2.lik (sup_mg),2) 
round (R2.lik (sup_mg, sup_mg_nophy),2)


############
### SCCS ###

GM_SCCS_data <- read.csv ("GM_SCCS_imputed.csv", header=TRUE, row.names=1)
GM_SCCS_tree <- read.tree ("SPT.SCCS.tre")

GM_SCCS_tree_star <- compute.brlen (GM_SCCS_tree, method="Grafen", power=0.0001) # convert the tree into a star-like phylogeny


### FGC ###

FGC_m5 <- phyloglm (FGC ~ bride_pr,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

FGC_m5_nophy <- phyloglm (FGC ~ bride_pr,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)


round (R2.lik (FGC_m5),2) 
round (R2.lik (FGC_m5, FGC_m5_nophy),2)

#

FGC_m7 <- phyloglm (FGC ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

FGC_m7_nophy <- phyloglm (FGC ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (FGC_m7),2) 
round (R2.lik (FGC_m7, FGC_m7_nophy),2)


### clitoridectomy ###

clit_m1 <- phyloglm (clit ~ dist + ext_aff,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0,
                    btol=20,
                    log.alpha.bound=10
)

clit_m1_nophy <- phyloglm (clit ~ dist + ext_aff,
                     phy=GM_SCCS_tree_star,
                     data=GM_SCCS_data,
                     boot=0,
                     btol=20,
                     log.alpha.bound=10
)

round (R2.lik (clit_m1),2) 
round (R2.lik (clit_m1, clit_m1_nophy),2)

#

clit_m5 <- phyloglm (clit ~ bride_pr,
                     phy=GM_SCCS_tree,
                     data=GM_SCCS_data,
                     boot=0,
                     log.alpha.bound=10
)

clit_m5_nophy <- phyloglm (clit ~ bride_pr,
                     phy=GM_SCCS_tree_star,
                     data=GM_SCCS_data,
                     boot=0,
                     log.alpha.bound=10
)

round (R2.lik (clit_m5),2) 
round (R2.lik (clit_m5, clit_m5_nophy),2)

#

clit_m4 <- phyloglm (clit ~ islam,
                     phy=GM_SCCS_tree,
                     data=GM_SCCS_data,
                     boot=0,
                     log.alpha.bound=10
)

clit_m4_nophy <- phyloglm (clit ~ islam,
                     phy=GM_SCCS_tree_star,
                     data=GM_SCCS_data,
                     boot=0,
                     log.alpha.bound=10
)

round (R2.lik (clit_m4),2) 
round (R2.lik (clit_m4, clit_m4_nophy),2)


### excision ###

exc_mg <- phyloglm (exc ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0,
                    btol=20
)

exc_mg_nophy <- phyloglm (exc ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + ext_agric + int_agric + bride_pr + class + caste + scars_f,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0,
                    btol=20
)

round (R2.lik (exc_mg),2) 
round (R2.lik (exc_mg, exc_mg_nophy),2)


### infibulation ###

inf_m5 <- phyloglm (inf ~ bride_pr,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

inf_m5_nophy <- phyloglm (inf ~ bride_pr,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (inf_m5),2) 
round (R2.lik (inf_m5, inf_m5_nophy),2)

#

inf_m1 <- phyloglm (inf ~ dist + ext_aff,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

inf_m1_nophy <- phyloglm (inf ~ dist + ext_aff,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (inf_m1),2) 
round (R2.lik (inf_m1, inf_m1_nophy),2)

#

inf_m4 <- phyloglm (inf ~ past,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

inf_m4_nophy <- phyloglm (inf ~ past,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (inf_m4),2) 
round (R2.lik (inf_m4, inf_m4_nophy),2)

#

inf_m7 <- phyloglm (inf ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

inf_m7_nophy <- phyloglm (inf ~ dist + ext_aff + sex_norms + patriloc + patrilin + past + bride_pr + caste,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (inf_m7),2) 
round (R2.lik (inf_m7, inf_m7_nophy),2)


### MGC ###

MGC_m5 <- phyloglm (MGC ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

MGC_m5_nophy <- phyloglm (MGC ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (MGC_m5),2) 
round (R2.lik (MGC_m5, MGC_m5_nophy),2)


### circumcision ###

cir_m6 <- phyloglm (cir ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste + islam,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

cir_m6_nophy <- phyloglm (cir ~ dist + ext_aff + segr_adol + patriloc + patrilin + chiefdoms + caste + islam,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (cir_m6),2) 
round (R2.lik (cir_m6, cir_m6_nophy),2)


### superincision ###

sup_m1 <- phyloglm (sup ~ dist + ext_aff,
                    phy=GM_SCCS_tree,
                    data=GM_SCCS_data,
                    boot=0
)

sup_m1_nophy <- phyloglm (sup ~ dist + ext_aff,
                    phy=GM_SCCS_tree_star,
                    data=GM_SCCS_data,
                    boot=0
)

round (R2.lik (sup_m1),2) 
round (R2.lik (sup_m1, sup_m1_nophy),2)

#

sup_m0 <- phyloglm (sup ~ 1, phy=GM_SCCS_tree, data=GM_SCCS_data, log.alpha.bound=10)

sup_m0_nophy <- phyloglm (sup ~ 1, phy=GM_SCCS_tree_star, data=GM_SCCS_data, log.alpha.bound=10)

round (R2.lik (sup_m0),2) 
round (R2.lik (sup_m0, sup_m0_nophy),2)


#############################################
#############################################
