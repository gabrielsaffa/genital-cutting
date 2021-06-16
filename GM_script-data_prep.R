################################################################################################################
### R script for: Global phylogenetic analysis reveals multiple origins and functions of genital mutilations ###
################################################################################################################

### predictor selection, data coding and missing values imputation ###

library (car)        
library (missMDA)    
library (FactoMineR) 
library (factoextra) 
library (missForest)

setwd ("")

### NOTE: here we'll select predictors for the SCCS only because the EA does not contain enough variables for each category to choose from

GM_SCCS_data <- read.csv ("GM_SCCS_original_full.csv", header=TRUE)
str (GM_SCCS_data)

# consider variables category by category


### marriage system ###

vif (glm(MGM ~ marit_comp + dist + cult_polyg + polyg + std_polyg, data=GM_SCCS_data, family="binomial")) # remove std_polyg 
vif (glm(MGM ~ marit_comp + dist + cult_polyg + polyg, data=GM_SCCS_data, family="binomial")) # remove polyg 
vif (glm(MGM ~ marit_comp + dist + cult_polyg, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_imp <- imputeCA (GM_SCCS_data, ncp=2, threshold=1e-08, maxiter=1000) # we need to impute missing values before proceeding to correspondence analysis (CA); imputed values have no influence on CA

GM_SCCS_ca <- CA (GM_SCCS_imp[,8:10], graph=TRUE) # run CA on marriage variables
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # assess the contribution of each variable; although both marital composition and distance between co-wives perform similarly, we choose the latter because it was used by Wilson (2008)

### sexual behavior ###

vif (glm(MGM ~ ext_sex + ext_aff + freq_ext_sex, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,13:15], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose extramarital affairs

vif (glm(MGM ~ sex_norms + segr_adol + sex_att + male_segr, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,16:19], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose segregation of adolescent boys and premarital sex norms

### residence ###

vif (glm(MGM ~ trans_resid + marit_resid + marit_resid_pp, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,20:22], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose transfer of residence

### descent ###

vif (glm(MGM ~ desc + desc_rules + desc_corp_gr, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,23:25], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose descent: major type

### subsistence ###

vif (glm(MGM ~ sub_econ + agric_inten + prin_subs, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,26:28], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose subsistence economy and agriculture: intensity

### marriage practices ###

vif (glm(MGM ~ marr_arr_m + marr_arr_f + marr_tr, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,29:31], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose marriage transactions

### social stratification ###

vif (glm(MGM ~ juris_hier + class + caste + soc_strat, data=GM_SCCS_data, family="binomial")) # remove soc_strat
vif (glm(MGM ~ juris_hier + class + caste, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,33:35], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose all three variables

### rituals ###

vif (glm(MGM ~ scars_m + scars_f, data=GM_SCCS_data, family="binomial")) # values <3

# both variables are important and will be included in the study

### religion ###

vif (glm(MGM ~ high_gods + islam + world_rel, data=GM_SCCS_data, family="binomial")) # values <3

GM_SCCS_ca <- CA (GM_SCCS_imp[,39:41], graph=TRUE)
summary (GM_SCCS_ca) 
GM_SCCS_ca$col$contrib # we choose Islam due to its high loading on the second dimension and because it is of particular interest


GM_SCCS_data <- GM_SCCS_data[,-c(8,10,11,12,13,15,18,19,21,22,24,25,28,30,31,36,39,41)] # keep only the selected variables
str (GM_SCCS_data) # OK

write.csv (GM_SCCS_data, "GM_SCCS_original.csv")


### let's first impute data with the original coding

############
### SCCS ###

GM_SCCS_data <- read.csv ("GM_SCCS_original.csv", header=TRUE)

pop <- GM_SCCS_data[,1] # store the population column
GM_SCCS_data <- GM_SCCS_data[,-1] # remove the population column


is.na(GM_SCCS_data$ext_war) <- GM_SCCS_data$ext_war=="88" # replace with NAs
is.na(GM_SCCS_data$ext_war) <- GM_SCCS_data$ext_war=="0" # replace with NAs
is.na(GM_SCCS_data$trans_resid) <- GM_SCCS_data$trans_resid=="4" # replace with NAs
is.na(GM_SCCS_data$sub_econ) <- GM_SCCS_data$sub_econ=="9" # replace with NAs

cols <- c (1:23)
GM_SCCS_data[cols] <- lapply (GM_SCCS_data[cols], factor) # convert to factors before imputing


GM_SCCS_imp <- missForest (GM_SCCS_data, variablewise=TRUE) # impute

GM_SCCS_imp$OOBerror # check the imputation error of each imputed variable
GM_SCCS_imp_e <- rbind (GM_SCCS_imp$OOBerror)


### now we will re-code the original SCCS data

GM_SCCS_data <- read.csv ("GM_SCCS_original.csv", header=TRUE)

# distance between co-wives - code as in Wilson (2008)

GM_SCCS_data$dist <- as.factor(GM_SCCS_data$dist)
levels(GM_SCCS_data$dist) <- list("1"=c("0","1"), "2"="2", "3"=c("3","4"), "4"=c("5","6"))

# extramarital affairs - invert the scale

GM_SCCS_data$ext_aff <- as.factor(GM_SCCS_data$ext_aff)
levels(GM_SCCS_data$ext_aff) <- list("1"="3", "2"="2", "3"="1")

# transfer of residence - we want only patrilocality

GM_SCCS_data$patriloc <- as.factor(GM_SCCS_data$trans_resid)
GM_SCCS_data$patriloc <- ifelse (GM_SCCS_data$patriloc==1,1,0)

# descent - we want only patrilineality

GM_SCCS_data$patrilin <- as.factor(GM_SCCS_data$desc)
GM_SCCS_data$patrilin <- ifelse (GM_SCCS_data$patrilin==1,1,0)

# subsistence economy - foraging = hunting + gathering

GM_SCCS_data$gath <- as.factor(GM_SCCS_data$sub_econ)
GM_SCCS_data$gath <- ifelse (GM_SCCS_data$gath==1,1,0)

GM_SCCS_data$hunt <- as.factor(GM_SCCS_data$sub_econ)
GM_SCCS_data$hunt <- ifelse (GM_SCCS_data$hunt==3,1,0)

GM_SCCS_data$forag <- paste(GM_SCCS_data$gath, GM_SCCS_data$hunt, sep="")
GM_SCCS_data$forag <- as.factor(GM_SCCS_data$forag)
levels(GM_SCCS_data$forag) <- list("1"=c("01","11","10"), "0"=c("00"))

# subsistence economy - pastoralism

GM_SCCS_data$past <- as.factor(GM_SCCS_data$sub_econ)
GM_SCCS_data$past <- ifelse (GM_SCCS_data$past==4,1,0)

# agriculture intensity - extensive agriculture

GM_SCCS_data$ext_agric <- as.factor(GM_SCCS_data$agric_inten)
GM_SCCS_data$ext_agric <- ifelse (GM_SCCS_data$ext_agric==3,1,0)

# agriculture intensity - intensive agriculture

GM_SCCS_data$int_agric <- as.factor(GM_SCCS_data$agric_inten)
levels(GM_SCCS_data$int_agric) <- list("1"=c("5","6"), "0"=c("1","2","3","4"))

# transactions at marriage - we want only bride-price

GM_SCCS_data$bride_pr <- as.factor(GM_SCCS_data$marr_tr)
GM_SCCS_data$bride_pr <- ifelse (GM_SCCS_data$bride_pr==1,1,0)

# premarital sex norms - code as in Ericksen (1989)

GM_SCCS_data$sex_norms <- as.factor(GM_SCCS_data$sex_norms)
levels(GM_SCCS_data$sex_norms) <- list("1"=c("1","2","3"), "0"=c("4","5","6"))

# external warfare 

is.na(GM_SCCS_data$ext_war) <- GM_SCCS_data$ext_war== "88" # replace with NAs
is.na(GM_SCCS_data$ext_war) <- GM_SCCS_data$ext_war== "0" # replace with NAs

GM_SCCS_data$ext_war <- as.factor(GM_SCCS_data$ext_war)
levels(GM_SCCS_data$ext_war) <- list("1"=c("1","2","3","4"), "2"=c("5","6","7","8"), "3"=c("9","10","11","12"), "4"=c("13","14","15","16"), "5"="17")

# classes

GM_SCCS_data$class <- as.factor(GM_SCCS_data$class)
levels(GM_SCCS_data$class) <- list("1"=c("1","2","3","4"),"0"="0")

# male segregation

GM_SCCS_data$segr_adol <- as.factor(GM_SCCS_data$segr_adol)
levels(GM_SCCS_data$segr_adol) <- list("1"=c("3","4","5"), "0"=c("1","2"))

# jurisdictional hierarchy - we want chiefdoms

GM_SCCS_data$chiefdoms <- as.factor(GM_SCCS_data$juris_hier)
levels(GM_SCCS_data$chiefdoms) <- list("1"=c("1","2"), "0"=c("0","3","4"))

# jurisdictional hierarchy - we want states

GM_SCCS_data$states <- as.factor(GM_SCCS_data$juris_hier)
levels(GM_SCCS_data$states) <- list("0"=c("0","1","2"), "1"=c("3","4"))

# castes

GM_SCCS_data$caste <- as.factor(GM_SCCS_data$caste)
levels(GM_SCCS_data$caste) <- list("0"="1", "1"=c("2","3","4"))

# scarifications - males

GM_SCCS_data$scars_m <- as.factor(GM_SCCS_data$scars_m)
levels(GM_SCCS_data$scars_m) <- list("0"=c("0","1"), "1"=c("2","3"))

# scarifications - females

GM_SCCS_data$scars_f <- as.factor(GM_SCCS_data$scars_f)
levels(GM_SCCS_data$scars_f) <- list("0"=c("0","1"), "1"=c("2","3"))

# Islam

GM_SCCS_data$islam <- as.factor(GM_SCCS_data$islam)
levels(GM_SCCS_data$islam) <- list("0"=c("0","2"), "1"="1")

str (GM_SCCS_data) 

GM_SCCS_data <- GM_SCCS_data[,-c(10,11,12,13,14,19,26,27)] # remove variables with the original coding
str (GM_SCCS_data) # OK

write.csv (GM_SCCS_data, "GM_SCCS_original_recoded.csv")


### we'll impute the re-coded data

GM_SCCS_data_rec <- read.csv ("GM_SCCS_original_recoded.csv", header=TRUE)

pop <- GM_SCCS_data_rec[,1] # store the population column
GM_SCCS_data_rec <- GM_SCCS_data_rec[,-1] # remove the population column

cols <- c (1:26)
GM_SCCS_data_rec[cols] <- lapply (GM_SCCS_data_rec[cols], factor) # convert to factors before imputing


GM_SCCS_imp_r <- missForest (GM_SCCS_data_rec, variablewise=TRUE) # impute

GM_SCCS_imp_r$OOBerror # check the imputation error of each imputed variable
GM_SCCS_imp_r_e <- rbind (GM_SCCS_imp_r$OOBerror)

### compare the imputation error between the original and the re-coded data

mean (GM_SCCS_imp_e)
mean (GM_SCCS_imp_r_e) # error is lower for the re-coded data

GM_SCCS_imp_r <- cbind (pop, GM_SCCS_imp_r$ximp) # bind the population column with the imputed data
write.csv (GM_SCCS_imp_r, file="GM_SCCS_imputed.csv") # store the imputed data


##########
### EA ###

### we'll do the same for the EA

GM_EA_data <- read.csv ("GM_EA_original.csv", header=TRUE)

GM_EA_data <- GM_EA_data[complete.cases(GM_EA_data[, 1:7]), ] # remove rows according to missing data for any GM
str (GM_EA_data) # N=575

pop <- GM_EA_data[,1] # store the population column
GM_EA_data <- GM_EA_data[,-1] # remove the population column

cols <- c (1:19)
GM_EA_data[cols] <- lapply (GM_EA_data[cols], factor) # convert to factors before imputing


GM_EA_imp <- missForest (GM_EA_data, variablewise=TRUE) # impute

GM_EA_imp$OOBerror # check the imputation error of each imputed variable
GM_EA_imp_e <- rbind (GM_EA_imp$OOBerror)


### now we'll re-code the original data and impute it as well

GM_EA_data <- read.csv ("GM_EA_original.csv", header=TRUE)
GM_EA_data <- GM_EA_data[complete.cases(GM_EA_data[, 1:7]), ] # remove rows according to missing data for any GM

# distance between co-wives

GM_EA_data$dist <- as.factor(GM_EA_data$marit_comp)
levels(GM_EA_data$dist) <- list("1"=c("4","5"), "0"=c("1","2","3","6","7"))

# premarital sex norms 

GM_EA_data$sex_norms <- as.factor(GM_EA_data$sex_norms)
levels(GM_EA_data$sex_norms) <- list("1"=c("1","2","3"), "0"=c("4","5","6"))

# male segregation

GM_EA_data$segr_adol <- as.factor(GM_EA_data$segr_adol)
levels(GM_EA_data$segr_adol) <- list("1"=c("3","4","5"), "0"=c("1","2"))

# transfer of residence - we want only patrilocality

GM_EA_data$patriloc <- as.factor(GM_EA_data$trans_resid)
GM_EA_data$patriloc <- ifelse (GM_EA_data$patriloc==1,1,0)

# descent - we want only patrilineality

GM_EA_data$patrilin <- as.factor(GM_EA_data$desc)
GM_EA_data$patrilin <- ifelse (GM_EA_data$patrilin==1,1,0)

# subsistence economy - pastoralism

GM_EA_data$past <- as.factor(GM_EA_data$sub_econ)
GM_EA_data$past <- ifelse (GM_EA_data$past==4,1,0)

# agriculture intensity - extensive agriculture

GM_EA_data$ext_agric <- as.factor(GM_EA_data$agric_inten)
GM_EA_data$ext_agric <- ifelse (GM_EA_data$ext_agric==3,1,0)

# agriculture intensity - intensive agriculture

GM_EA_data$int_agric <- as.factor(GM_EA_data$agric_inten)
levels(GM_EA_data$int_agric) <- list("1"=c("5","6"), "0"=c("1","2","3","4"))

# transactions at marriage - we want only bride-price

GM_EA_data$bride_pr <- as.factor(GM_EA_data$marr_tr)
GM_EA_data$bride_pr <- ifelse (GM_EA_data$bride_pr==1,1,0)

# jurisdictional hierarchy - we want chiefdoms

GM_EA_data$chiefdoms <- as.factor(GM_EA_data$juris_hier)
levels(GM_EA_data$chiefdoms) <- list("1"=c("2","3"), "0"=c("1","4","5"))

# jurisdictional hierarchy - we want states

GM_EA_data$states <- as.factor(GM_EA_data$juris_hier)
levels(GM_EA_data$states) <- list("0"=c("1","2","3"), "1"=c("4","5"))

# classes

GM_EA_data$class <- as.factor(GM_EA_data$class)
levels(GM_EA_data$class) <- list("1"=c("2","3","4","5"),"0"="1")

# castes

GM_EA_data$caste <- as.factor(GM_EA_data$caste)
levels(GM_EA_data$caste) <- list("0"="1", "1"=c("2","3","4"))

# high gods

GM_EA_data$high_gods <- as.factor(GM_EA_data$high_gods)
levels(GM_EA_data$high_gods) <- list("1"=c("3","4"), "0"=c("1","2"))


str (GM_EA_data)

GM_EA_data <- GM_EA_data[,-c(8,11,12,13,14,15,16)] # remove unnecessary variables
str (GM_EA_data) # OK

write.csv (GM_EA_data, "GM_EA_original_recoded.csv")


### impute

GM_EA_data_rec <- read.csv ("GM_EA_original_recoded.csv", header=TRUE)

pop <- GM_EA_data_rec[,1] # store the population column
GM_EA_data_rec <- GM_EA_data_rec[,-1] # remove the population column

cols <- c (1:21)
GM_EA_data_rec[cols] <- lapply (GM_EA_data_rec[cols], factor) # convert to factors before imputing


GM_EA_imp_r <- missForest (GM_EA_data_rec, variablewise=TRUE) # impute

GM_EA_imp_r$OOBerror # check the imputation error of each imputed variable
GM_EA_imp_r_e <- rbind (GM_EA_imp_r$OOBerror)

### compare the imputation error between the original and the re-coded data

mean (GM_EA_imp_e)
mean (GM_EA_imp_r_e) # error is lower for the re-coded data

GM_EA_imp_r <- cbind (pop, GM_EA_imp_r$ximp) 

write.csv (GM_EA_imp_r, file="GM_EA_imputed.csv") 


###################################
###################################