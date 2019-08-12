# preparation -------------------------------------------------------------
rm(list = ls())
library(psych)
library(ggplot2)
library(plyr)

i <- 0
cond <- c("N_pre","N_sp","N_rest","N_bl")
all_r <- matrix(nrow = 18,ncol = length(cond))
colnames(all_r) <- cond
all_p <- matrix(nrow = 18,ncol = length(cond))
colnames(all_p) <- cond
all_p_fdr <- matrix(nrow = 18,ncol = length(cond))
colnames(all_p_fdr) <- cond

# DV <- "PRCS"
DV <- "Performance"
allsig_ROI <- vector("list",length=length(cond))

for(icond in cond){
  i <- i+1
  
  # read data ---------------------------------------------------------------
  beta <- read.csv(paste("./beta/narm_",icond,"_ROI.csv",sep = ""),header = TRUE)
  ROInum <- length(colnames(beta)) - 1
  sum <- read.csv("sum.csv",header = TRUE)
  sumdata <- merge(sum,beta,by="subject")
  rm(beta,sum)
  
  sumdata$r_C <- NULL
  sumdata$l_C <- NULL
  sumdata$r_T <- NULL
  sumdata$l_T <- NULL
  ROInum <- ROInum - 4
  
  corr_data <- sumdata[,c(colnames(sumdata)[6:(ROInum+5)],DV)]
  
  # correlation -------------------------------------------------------------
  corr_test <- corr.test(corr_data)
  corr_r <- matrix(corr_test$r, nrow = (ROInum+1), ncol = (ROInum+1))
  colnames(corr_r) <- colnames(corr_test$r)
  rownames(corr_r) <- colnames(corr_test$r)
  corr_p <- matrix(corr_test$p, nrow = (ROInum+1), ncol = (ROInum+1))
  colnames(corr_p) <- colnames(corr_test$p)
  rownames(corr_p) <- colnames(corr_test$p)
  corr_p_fdr <- matrix(p.adjust(corr_test$p,method = "fdr"), nrow = (ROInum+1), ncol = (ROInum+1))
  colnames(corr_p_fdr) <- colnames(corr_test$p)
  rownames(corr_p_fdr) <- colnames(corr_test$p)
  
  sigROI <- which(corr_test$p[DV,1:ROInum]<0.05)
  sigROIname <- colnames(corr_data)[sigROI]
  sigROI_r <- corr_test$r[DV,sigROI]
  sigROI_p <- corr_test$p[DV,sigROI]
  sigROI_p_fdr <- p.adjust(corr_test$p[DV,1:ROInum],method = "fdr",n=ROInum)[sigROI]
  sig_data <- data.frame(sigROIname,sigROI_r,sigROI_p,sigROI_p_fdr)
  colnames(sig_data) <- c("ROI","r","p","p.FDR")
  allsig_ROI[i] <- list(sig_data)
  
  rm(sigROI,sigROI_r,sigROI_p,sigROI_p_fdr)
  
  # output ------------------------------------------------------------------
  all_r[1:ROInum,i] <- corr_r[(ROInum+1),1:ROInum]
  all_p[1:ROInum,i] <- corr_p[(ROInum+1),1:ROInum]
  all_p_fdr[1:ROInum,i] <- corr_p_fdr[(ROInum+1),1:ROInum]
  rm(corr_data,corr_test,corr_r,corr_p,corr_p_fdr)
}
rm(i,icond)
write.csv(all_r,paste("./",DV,"/",DV,"_r.csv",sep = ""),row.names = FALSE)
write.csv(all_p,paste("./",DV,"/",DV,"_p.csv",sep = ""),row.names = FALSE)
write.csv(all_p_fdr,paste("./",DV,"/",DV,"_p_fdr.csv",sep = ""),row.names = FALSE)