###############################################################################################################################
### R script  for: Global phylogenetic analysis reveals multiple origins and correlates of genital cutting ###
###############################################################################################################################

### ancestral state reconstructions of GM for individual language families ###

library (phytools)

setwd ("")


### Austronesian ###

subtree_aus <- read.nexus ("AUS.nex")
subdata_aus <- read.csv ("sample-gray-et-al.---spt.csv", header=TRUE, row.names=1)
subtree_aus_cons <- ls.consensus (subdata_aus)


### clitoridectomy ###

{
   
clit <- setNames (subdata_aus[,2], rownames(subdata_aus))
aus_clit <- make.simmap (subtree_aus, clit, model="ARD", Q="empirical")
aus_clit_pd <- summary (aus_clit, consensus.tree=subtree_aus_cons)

   cols <- setNames (c("gray85","turquoise2"), levels(as.factor(clit)))
   plot (aus_clit_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=80,
                      colors=setNames(c("turquoise2", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )

}


### circumcision ###

{
   
   cir <- setNames (subdata_aus[,4], rownames(subdata_aus))
   aus_cir <- make.simmap (subtree_aus, cir, model="ARD", Q="empirical")
   aus_cir_pd <- summary (aus_cir, consensus.tree=subtree_aus_cons)
   
   cols <- setNames (c("gray85","gold"), levels(as.factor(cir)))
   plot (aus_cir_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=80,
                      colors=setNames(c("gold", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )
   
}


### superincision ###

{
   
   sup <- setNames (subdata_aus[,5], rownames(subdata_aus))
   aus_sup <- make.simmap (subtree_aus, sup, model="ARD", Q="empirical")
   aus_sup_pd <- summary (aus_sup, consensus.tree=subtree_aus_cons)
   
   cols <- setNames (c("gray85","red4"), levels(as.factor(sup)))
   plot (aus_sup_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=80,
                      colors=setNames(c("red4", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )
   
}


### Bantu ###

subtree_ban <- read.nexus ("BAN.nex")
subdata_ban <- read.csv ("sample-grollemund---spt.csv", header=TRUE, row.names=1)
subtree_ban_cons <- ls.consensus (subtree_ban)


### FGC ###

{

FGC <- setNames (subdata_ban[,2], rownames(subdata_ban))
ban_FGC <- make.simmap (subtree_ban, FGC, model="ARD", Q="empirical")
ban_FGC_pd <- summary (ban_FGC, consensus.tree=subtree_ban_cons)

cols <- setNames (c("gray85","blue"), levels(as.factor(FGC)))
plot (ban_FGC_pd,
      offset=0.5,
      fsize=0.6,
      lwd=1,
      cex=c(0.3,0.3),
      colors=cols
)

add.simmap.legend (x=0,
                   y=85,
                   colors=setNames(c("blue", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)

}


### clitoridectomy ###

{
   
clit <- setNames (subdata_ban[,3], rownames(subdata_ban))
ban_clit <- make.simmap (subtree_ban, clit, model="ARD", Q="empirical")
ban_clit_pd <- summary (ban_clit, consensus.tree=subtree_ban_cons)

cols <- setNames (c("gray85","turquoise2"), levels(as.factor(clit)))
plot (ban_clit_pd,
      offset=0.5,
      fsize=0.6,
      lwd=1,
      cex=c(0.3,0.3),
      colors=cols
)

add.simmap.legend (x=0,
                   y=85,
                   colors=setNames(c("turquoise2", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)

}


### excision ###

{
   
exc <- setNames (subdata_ban[,4], rownames(subdata_ban))
ban_exc <- make.simmap (subtree_ban, exc, model="ARD", Q="empirical")
ban_exc_pd <- summary (ban_exc, consensus.tree=subtree_ban_cons)

cols <- setNames (c("gray85","turquoise4"), levels(as.factor(exc)))
plot (ban_exc_pd,
      offset=0.5,
      fsize=0.6,
      lwd=1,
      cex=c(0.3,0.3),
      colors=cols
)

add.simmap.legend (x=0,
                   y=85,
                   colors=setNames(c("turquoise4", "gray85"),
                                   c("present", "absent")),
                   prompt=FALSE,
                   vertical=TRUE,
                   shape="circle",
                   cex=0.3
)

}


### circumcision ###

{
   
   cir <- setNames (subdata_ban[,5], rownames(subdata_ban))
   ban_cir <- make.simmap (subtree_ban, cir, model="ARD", Q="empirical")
   ban_cir_pd <- summary (ban_cir, consensus.tree=subtree_ban_cons)
   
   cols <- setNames (c("gray85","gold"), levels(as.factor(cir)))
   plot (ban_cir_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=85,
                      colors=setNames(c("gold", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )
   
}


### Indo-European ###

subtree_ie <- read.nexus ("IE.nex")
subdata_ie <- read.csv ("sample-bouckaert---spt.csv", header=TRUE, row.names=1)
subtree_ie_cons <- ls.consensus (subtree_ie)


### clitoridectomy ###

{
   
   clit <- setNames (subdata_ie[,2], rownames(subdata_ie))
   ie_clit <- make.simmap (subtree_ie, clit, model="ARD", Q="empirical")
   ie_clit_pd <- summary (ie_clit, consensus.tree=subtree_ie_cons)
   
   cols <- setNames (c("gray85","turquoise2"), levels(as.factor(clit)))
   plot (ie_clit_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=30,
                      colors=setNames(c("turquoise2", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )
   
}


### circumcision ### 

{
   
   cir <- setNames (subdata_ie[,3], rownames(subdata_ie))
   ie_cir <- make.simmap (subtree_ie, cir, model="ARD", Q="empirical")
   ie_cir_pd <- summary (ie_cir, consensus.tree=subtree_ie_cons)
   
   cols <- setNames (c("gray85","gold"), levels(as.factor(cir)))
   plot (ie_cir_pd,
         offset=0.5,
         fsize=0.6,
         lwd=1,
         cex=c(0.3,0.3),
         colors=cols
   )
   
   add.simmap.legend (x=0,
                      y=30,
                      colors=setNames(c("gold", "gray85"),
                                      c("present", "absent")),
                      prompt=FALSE,
                      vertical=TRUE,
                      shape="circle",
                      cex=0.3
   )
   
}


###################################
###################################
