# preparation -------------------------------------------------------------
rm(list = ls())
library(psych)
library(ggplot2)
library(plyr)

i <- 0
subnum <- 20
# subject <- c(7:16,18,19,21:28)
cond <- c("N_pre","N_sp","N_rest","N_bl")
ROIs <- c("r_OFC","r_DLPFC","r_preC","r_postC","r_IFG","r_TP","r_STG","r_MTG","r_pSTG","r_C","r_T",
          "l_OFC","l_DLPFC","l_preC","l_postC","l_IFG","l_TP","l_STG","l_MTG","l_pSTG","l_C","l_T")

for(icond in cond){
  i <- i+1
  
  # read data ---------------------------------------------------------------
  beta <- read.csv(paste("../nirs-beta/",icond,"_beta.csv",sep = ""),header = TRUE)
  
  ROI <- matrix(0,nrow = length(beta[,1]),ncol = length(ROIs)+1)
  colnames(ROI) <- c("subject",
                     "r_OFC","r_DLPFC","r_preC","r_postC","r_IFG","r_TP","r_STG","r_MTG","r_pSTG","r_C","r_T",
                     "l_OFC","l_DLPFC","l_preC","l_postC","l_IFG","l_TP","l_STG","l_MTG","l_pSTG","l_C","l_T")
  beta[is.na(beta)] <- 0

  # calculate ROI -----------------------------------------------------------
  attach(beta)
  
  ## Right Hemisphere
  ROI[,"r_OFC"] <- CH1 + CH5*0.829 + CH6*0.086
  ROI[,"r_DLPFC"] <- CH5*0.171 + CH10 + CH14 + CH15*0.352 + CH19 + CH20*0.799
  ROI[,"r_preC"] <- CH16*0.869 + CH20*0.030 + CH21*0.366
  ROI[,"r_postC"] <- CH12 + CH16*0.122 + CH17*0.007 + CH21*0.634
  ROI[,"r_IFG"] <- CH6*0.914 + CH7*0.541 + CH11 + CH15*0.648 + CH16*0.009 + CH20*0.172
  ROI[,"r_TP"] <- CH2 + CH7*0.173
  ROI[,"r_STG"] <- CH3*0.321 + CH7*0.286 + CH8 + CH13*0.756
  ROI[,"r_MTG"] <- CH3*0.679 + CH4 + CH9
  ROI[,"r_pSTG"] <- CH13*0.244 + CH17*0.993 + CH18 + CH22
  ROI[,"r_C"] <- CH12 + CH16*0.991 + CH17*0.007 + CH20*0.030 + CH21
  ROI[,"r_T"] <- CH2 + CH3 + CH4 + CH7*0.459 + CH8 + CH9 + CH13 + CH17*0.993 + CH18 + CH22
  
  ## Left Hemisphere
  ROI[,"l_OFC"] <- CH26 + CH31*0.829 + CH30*0.086
  ROI[,"l_DLPFC"] <- CH31*0.171 + CH35 + CH40 + CH39*0.352 + CH44 + CH43*0.799
  ROI[,"l_preC"] <- CH38*0.869 + CH43*0.030 + CH42*0.366
  ROI[,"l_postC"] <- CH33 + CH38*0.122 + CH37*0.007 + CH42*0.634
  ROI[,"l_IFG"] <- CH30*0.914 + CH29*0.541 + CH34 + CH39*0.648 + CH38*0.009 + CH43*0.172
  ROI[,"l_TP"] <- CH25 + CH29*0.173
  ROI[,"l_STG"] <- CH24*0.321 + CH29*0.286 + CH28 + CH32*0.756
  ROI[,"l_MTG"] <- CH24*0.679 + CH23 + CH27
  ROI[,"l_pSTG"] <- CH32*0.244 + CH37*0.993 + CH36 + CH41
  ROI[,"l_C"] <- CH33 + CH38*0.991 + CH37*0.007 + CH43*0.030 + CH42
  ROI[,"l_T"] <- CH25 + CH24 + CH23 + CH29*0.459 + CH28 + CH27 + CH32 + CH37*0.993 + CH36 + CH41
  
  detach(beta)
  
  ROI[,"subject"] <- beta$subject
  
  # output ------------------------------------------------------------------
  # write.csv(ROI,paste(icond,"_ROI.csv",sep = ""),row.names = FALSE)
  write.csv(ROI,paste("narm_",icond,"_ROI.csv",sep = ""),row.names = FALSE)
}
