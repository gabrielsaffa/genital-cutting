###############################################################################################################################
### R script  for: Global phylogenetic analysis reveals multiple origins and correlates of genital cutting ###
###############################################################################################################################

### visualizing the data and the results ###

library(rgeos)
library(rworldmap)
library(cleangeo)
library (ggplot2)
library (ggpubr)

setwd ("")


### Fig 1 and Supplementary Fig 1 ###

GM_EA_data <- read.csv ("GM_EA_imputed.csv", header=TRUE) 
GM_SCCS_data <- read.csv ("GM_SCCS_imputed.csv", header=TRUE) 

sPDF <- getMap ("continent")
sPDF <- clgeo_Clean (sPDF)
poly_globe <- gUnionCascaded (sPDF)


### Fig 1 ###

par(mfrow=c(3,2), mar=c(2,2,5,2))

plot (poly_globe)
title ("a", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_EA_data$clit==1,"turquoise2","gray85")
points(GM_EA_data$long, GM_EA_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("d", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_EA_data$cir==1,"gold","gray85")
points(GM_EA_data$long, GM_EA_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("b", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_EA_data$exc==1,"turquoise4","gray85")
points(GM_EA_data$long, GM_EA_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("e", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_EA_data$sup==1,"red4","gray85")
points(GM_EA_data$long, GM_EA_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("c", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_EA_data$inf==1,"dodgerblue4","gray85")
points(GM_EA_data$long, GM_EA_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

legend (380,80, 
        legend=c("clitoridectomy","excision","infibulation","circumcision","superincision","absent"), 
        pt.bg=c("turquoise2","turquoise4","dodgerblue4","gold","red4","gray85"), 
        pch=21, bg="black", cex=2, bty="n", xpd=NA)

dev.off()


### Extended Data Fig 1 ###

par(mfrow=c(3,2), mar=c(2,2,5,2))

plot (poly_globe)
title ("a", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_SCCS_data$clit==1,"turquoise2","gray85")
points(GM_SCCS_data$long, GM_SCCS_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("d", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_SCCS_data$cir==1,"gold","gray85")
points(GM_SCCS_data$long, GM_SCCS_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("b", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_SCCS_data$exc==1,"turquoise4","gray85")
points(GM_SCCS_data$long, GM_SCCS_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("e", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_SCCS_data$sup==1,"red4","gray85")
points(GM_SCCS_data$long, GM_SCCS_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

plot (poly_globe)
title ("c", adj=0.05, line=0.1, cex.main=2)
bg_cols <- ifelse (GM_SCCS_data$inf==1,"dodgerblue4","gray85")
points(GM_SCCS_data$long, GM_SCCS_data$lat, col="black", bg=bg_cols, pch=21, cex=2)

legend (380,80, 
        legend=c("clitoridectomy","excision","infibulation","circumcision","superincision","absent"), 
        pt.bg=c("turquoise2","turquoise4","dodgerblue4","gold","red4","gray85"), 
        pch=21, bg="black", cex=2, bty="n", xpd=NA)

dev.off()


### NOTE: code for creating Supplementary Figs 3-17 is part of the scripts for phylogenetic analyses; Fig 2 and Supplementary Fig 2 were, after plotting individual phylogenies, created in a graphical software


### Fig 3 and Extended Data Fig 3 ###
### Fig 3 ###

{
  
  fgc_EA <- read.csv ("fgc_all_EA.csv", header=TRUE)
  cl_EA <- read.csv ("cl_all_EA.csv", header=TRUE)
  ex_EA <- read.csv ("ex_all_EA.csv", header=TRUE)
  infib_EA <- read.csv ("infib_all_EA.csv", header=TRUE)
  mgc_EA <- read.csv ("mgc_mg_EA.csv", header=TRUE)
  circum_EA <- read.csv ("circum_mg_EA.csv", header=TRUE)
  super_EA <- read.csv ("super_all_EA.csv", header=TRUE)
  
  
  ### FGC ###
  
  var_order <- c ("caste","class","bride_pr","int_agric","ext_agric","past","patrilin","patriloc","sex_norms","dist")
  names <- c ("castes","classes","bride-price","inten. agric.","exten. agric.","pastoralism","patrilineality","patrilocality","sex norms","co-wives separate")
  FGC_EA <- ggplot(fgc_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="blue") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="blue") + ggtitle("a") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### clit ###
  
  var_order <- c ("caste","bride_pr","past","patrilin","patriloc","sex_norms","dist")
  names <- c ("castes","bride-price","pastoralism","patrilineality","patrilocality","sex norms","co-wives separate")
  clit_EA <- ggplot(cl_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="turquoise2") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="turquoise2") + ggtitle("b") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### exc ###
  
  exc_EA <- ggplot(ex_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="turquoise4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="turquoise4") + ggtitle("c") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### inf ###
  
  var_order <- c ("caste","class","bride_pr","past","patrilin","patriloc","sex_norms","dist")
  names <- c ("castes","classes","bride-price","pastoralism","patrilineality","patrilocality","sex norms","co-wives separate")
  inf_EA <- ggplot(infib_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="dodgerblue4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="dodgerblue4") + ggtitle("d") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### MGC ###
  
  var_order <- c ("high_gods","caste","states","chiefdoms","bride_pr","int_agric","ext_agric","past","patrilin","patriloc","segr_adol","dist")
  names <- c ("high gods","castes","states","chiefdoms","bride-price","inten. agric.","exten. agric.","pastoralism","patrilineality","patrilocality","male segregation","co-wives separate")
  MGC_EA <- ggplot(mgc_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="red1") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="red1") + ggtitle("e") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### cir ###
  
  var_order <- c ("high_gods","caste","states","chiefdoms","bride_pr","int_agric","ext_agric","past","patrilin","patriloc","segr_adol","dist")
  names <- c ("high gods","castes","states","chiefdoms","bride-price","inten. agric.","exten. agric.","pastoralism","patrilineality","patrilocality","male segregation","co-wives separate")
  cir_EA <- ggplot(circum_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="gold") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="gold") + ggtitle("f") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### sup ###
  
  var_order <- c ("high_gods","caste","states","chiefdoms","bride_pr","int_agric","ext_agric","past","patrilin","patriloc","segr_adol","dist")
  names <- c ("high gods","castes","states","chiefdoms","bride-price","inten. agric.","exten. agric.","pastoralism","patrilineality","patrilocality","male segregation","co-wives separate")
  sup_EA <- ggplot(super_EA,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="red4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="red4") + ggtitle("g") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none") + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  
  # combine in a multiple plot
  
  ggarrange(FGC_EA,clit_EA,exc_EA,inf_EA,MGC_EA,cir_EA,sup_EA, 
            ncol=4, 
            nrow=2,
            common.legend=FALSE
  ) 
  
}


### Extended Data Fig 3 ###

{
  
  fgc_SCCS <- read.csv ("fgc_all_SCCS.csv", header=TRUE)
  cl_SCCS <- read.csv ("cl_all_SCCS.csv", header=TRUE)
  ex_SCCS <- read.csv ("ex_all_SCCS.csv", header=TRUE)
  infib_SCCS <- read.csv ("infib_all_SCCS.csv", header=TRUE)
  mgc_SCCS <- read.csv ("mgc_m5_SCCS.csv", header=TRUE)
  circum_SCCS <- read.csv ("circum_all_SCCS.csv", header=TRUE)
  super_SCCS <- read.csv ("super_m1_SCCS.csv", header=TRUE)
  
  
  ### FGC ###
  
  var_order <- c ("caste","bride_pr","past","patrilin","patriloc","sex_norms","ext_aff","dist")
  names <- c ("castes","bride-price","pastoralism","patrilineality","patrilocality","sex norms","extramarital sex","co-wives separate")
  FGM_SCCS <- ggplot(fgc_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="blue") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="blue") + ggtitle("a") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### clit ###
  
  var_order <- c ("islam","bride_pr","ext_aff","dist")
  names <- c ("Islam","bride-price","extramarital sex","co-wives separate")
  clit_SCCS <- ggplot(cl_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="turquoise2") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="turquoise2") + ggtitle("b") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### exc ###
  
  var_order <- c ("scars_f","caste","class","bride_pr","int_agric","ext_agric","past","patrilin","patriloc","sex_norms","ext_aff","dist")
  names <- c ("scars-female","castes","classes","bride-price","inten. agric.","exten. agric.","pastoralism","patrilineality","patrilocality","sex norms","extramarital sex","co-wives separate")
  exc_SCCS <- ggplot(ex_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="turquoise4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="turquoise4") + ggtitle("c") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### inf ###
  
  var_order <- c ("caste","bride_pr","past","patrilin","patriloc","sex_norms","ext_aff","dist")
  names <- c ("castes","bride-price","pastoralism","patrilineality","patrilocality","sex norms","extramarital sex","co-wives separate")
  inf_SCCS <- ggplot(infib_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="dodgerblue4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="dodgerblue4") + ggtitle("d") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### MGC ###
  
  var_order <- c ("caste","chiefdoms","patrilin","patriloc","segr_adol","ext_aff","dist")
  names <- c ("castes","chiefdoms","patrilineality","patrilocality","male segregation","extramarital sex","co-wives separate")
  MGM_SCCS <- ggplot(mgc_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="red1") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="red1") + ggtitle("e") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### cir ###
  
  var_order <- c ("islam","caste","chiefdoms","patrilin","patriloc","segr_adol","ext_aff","dist")
  names <- c ("Islam","castes","chiefdoms","patrilineality","patrilocality","male segregation","extramarital sex","co-wives separate")
  cir_SCCS <- ggplot(circum_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="gold") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="gold") + ggtitle("f") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  ### sup ###
  
  var_order <- c ("ext_aff","dist")
  names <- c ("extramarital sex","co-wives separate")
  sup_SCCS <- ggplot(super_SCCS,aes(x=var,y=b_est,fill=factor(mod),shape=factor(mod))) + geom_point(position=position_dodge(-0.65), size=2, color="red4") + scale_fill_discrete(name="var") + geom_errorbar(aes(ymin=lci, ymax=uci),position=position_dodge(-0.65),width=0,color="red4") + ggtitle("g") + geom_hline(yintercept=0,linetype="dashed") + theme_classic(base_size=12) + coord_flip() + guides(fill=FALSE) + theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), plot.title=element_text(size=14,face="bold"), plot.margin=unit(c(0.75,0.25,0.75,0.25),"cm")) + scale_x_discrete(limits=var_order, labels=names)
  
  
  # combine in a multiple plot
  
  ggarrange(FGC_SCCS,clit_SCCS,exc_SCCS,inf_SCCS,MGC_SCCS,cir_SCCS,sup_SCCS, 
            ncol=4, 
            nrow=2,
            common.legend=FALSE
  ) 
  
}


#########################################
#########################################
