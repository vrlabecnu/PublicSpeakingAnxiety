# preparation -------------------------------------------------------------
rm(list = ls())
library(MASS)
library(leaps)
library(car)
library(psych)
library(ggplot2)
library(plyr)
library(tidyr)
library(bootstrap)

cond <- c("N_pre","N_sp","N_rest","N_bl")
i <- 0


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

windowsFonts(TimesFont=windowsFont("Times New Roman"),ArialFont=windowsFont("Arial"),Dots=windowsFont("ZapfDingbats BT"))
theme_corr <-  theme(
  aspect.ratio = 1,
  panel.grid = element_blank(),
  axis.title = element_text(face="bold", size=12, colour="black", family="ArialFont"),
  axis.text  = element_text(size=10.5, colour="black", family="ArialFont"),
  axis.line = element_line(size=0.5),
  legend.title = element_text(face="bold", size=12, colour="white", family="ArialFont"),
  legend.text  = element_text(size=10.5, colour="black", family="ArialFont"),
  legend.position = c(0.75,0.75)
)
theme_reg <-  theme(
  aspect.ratio = 0.75,
  panel.grid = element_blank(),
  axis.title = element_text(face="bold", size=12, colour="black", family="ArialFont"),
  axis.text  = element_text(size=10.5, colour="black", family="ArialFont"),
  axis.line = element_line(size=0.5),
  legend.title = element_text(face="bold", size=12, colour="black", family="ArialFont"),
  legend.text  = element_text(size=10.5, colour="black", family="ArialFont")
)


# tidy data ---------------------------------------------------------------
pre <- read.csv("./beta/reg_narm_N_pre_ROI.csv",header = TRUE)
sp <- read.csv("./beta/reg_narm_N_sp_ROI.csv",header = TRUE)

# regression --------------------------------------------------------------
lm_pre_PRCS <- lm(PRCS ~ r_DLPFC+r_preC+r_postC+r_STG+l_IFG+l_DLPFC+l_pSTG, data = pre)
lm_sp_Performance <- lm(Performance ~ r_postC+r_TP+r_STG+r_MTG, data = sp)

# scatter plot ------------------------------------------------------------
pre$reg_PRCS <- lm_pre_PRCS$fitted.values
cor_pre_PRCS <- cor(pre$PRCS,pre$reg_PRCS)
figcor_pre_PRCS <- ggplot(pre, aes(x=reg_PRCS,y=PRCS)) +
  geom_point(shape=1,alpha=1) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "black",linetype = "dashed",size = 0.5) +
  scale_x_continuous(name="Anticipated Phase Prediction",limits=c(0,30),breaks=seq(0,30,6),expand = c(0,0)) +
  scale_y_continuous(name="PRCS",limits=c(0,30),breaks=seq(0,30,6),expand = c(0,0)) +
  theme_classic()
figcor_pre_PRCS <- figcor_pre_PRCS + theme_corr +
  annotate("text",x=Inf,y=-Inf,vjust=-1,hjust=1,family="ArialFont",
           label=paste("R2 = ",as.character(round(as.numeric(summary(lm_pre_PRCS)["r.squared"]),3)),sep=""))

sp$reg_Performance <- lm_sp_Performance$fitted.values
cor_sp_Performance <- cor(sp$Performance,sp$reg_Performance)
figcor_sp_Performance <- ggplot(sp, aes(x=reg_Performance,y=Performance)) +
  geom_point(shape=1,alpha=1) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, colour = "black",linetype = "dashed",size = 0.5) +
  scale_x_continuous(name="Speaking Phase Prediction",limits=c(10,25),breaks=seq(10,25,3),expand = c(0,0)) +
  scale_y_continuous(name="Performance",limits=c(10,25),breaks=seq(10,25,3),expand = c(0,0)) +
  theme_classic()
figcor_sp_Performance <- figcor_sp_Performance + theme_corr +
  annotate("text",x=Inf,y=-Inf,vjust=-1,hjust=1,family="ArialFont",
           label=paste("R2 = ",as.character(round(as.numeric(summary(lm_sp_Performance)["r.squared"]),3)),sep=""))


# shuffle -----------------------------------------------------------------
shrinkage <- function(fit, k) {
  require(bootstrap)
  
  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
  
  x <- fit$model[,2:ncol(fit$model)]
  y <- fit$model[,1]
  
  results <- crossval(x, y, theta.fit, theta.predict, ngroup=k)
  r2 <- cor(y, fit$fitted.values)^2
  r2cv <- cor(y, results$cv.fit)^2
  
  return(r2cv)
}
shufflenum <- 1000


# PRCS pre ----------------------------------------------------------------
reg_pre_PRCS_Rsq <- shrinkage(lm_pre_PRCS, k=10)
reg_pre_PRCS_randRsq <- array(data = NA, dim = shufflenum)
for (i in 1:shufflenum) {
  randnum <- sample.int(20,20,replace = FALSE)
  pre$randPRCS <- pre$PRCS[randnum]
  currfit <- lm(randPRCS ~ r_DLPFC+r_preC+r_postC+r_STG+l_IFG+l_DLPFC+l_pSTG, data = pre)
  reg_pre_PRCS_randRsq[i] <- shrinkage(currfit, k=10)
}
rm(i,randnum)

correlations <- data.frame(rep(0.494,shufflenum),reg_pre_PRCS_randRsq)
colnames(correlations) <- c("True", "Permutation")
corrgraph <- gather(correlations,condition,rsq, True:Permutation, factor_key=TRUE)
corrdat <- ddply(corrgraph, "condition", summarise, rsq.mean=mean(rsq))
corrgraph <- data.frame(correlations$Permutation)
colnames(corrgraph) <- 'rsq'

figreg_pre_PRCS <- ggplot(corrgraph, aes(x=rsq,y=..count../1000)) +
  geom_area(stat="bin",binwidth=0.03, position="identity", colour="black", fill="white") +
  geom_vline(data=corrdat, aes(xintercept=rsq.mean, color = condition),linetype="dashed", size=1) +
  scale_colour_grey() + guides(colour=FALSE) +
  geom_vline(data=correlations, aes(xintercept=mean(Permutation)+2*sd(Permutation)),linetype=3, size=0.5, colour = "black") +
  scale_x_continuous(name = "R-squared", limits = c(0,1), breaks = seq(0,1,0.25), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion",limits = c(0,0.5), breaks = seq(0,0.5,0.1), expand = c(0,0)) +
  theme_classic()
figreg_pre_PRCS <- figreg_pre_PRCS + theme_reg

# Performance sp ----------------------------------------------------------------
reg_sp_Performance_Rsq <- shrinkage(lm_sp_Performance, k=10)
reg_sp_Performance_randRsq <- array(data = NA, dim = shufflenum)
for (i in 1:shufflenum) {
  randnum <- sample.int(20,20,replace = FALSE)
  sp$randPerformance <- sp$Performance[randnum]
  currfit <- lm(randPerformance ~ r_postC+r_TP+r_STG+r_MTG, data = sp)
  reg_sp_Performance_randRsq[i] <- shrinkage(currfit, k=10)
}
rm(i,randnum)
print(mean(reg_sp_Performance_randRsq)+2*sd(reg_sp_Performance_randRsq))

correlations <- data.frame(rep(0.419,shufflenum),reg_sp_Performance_randRsq)
colnames(correlations) <- c("True", "Permutation")
corrgraph <- gather(correlations,condition,rsq, True:Permutation, factor_key=TRUE)
corrdat <- ddply(corrgraph, "condition", summarise, rsq.mean=mean(rsq))
corrgraph <- data.frame(correlations$Permutation)
colnames(corrgraph) <- 'rsq'

figreg_sp_Performance <- ggplot(corrgraph, aes(x=rsq,y=..count../1000)) +
  geom_area(stat="bin",binwidth=0.03, position="identity", colour="black", fill="white") +
  geom_vline(data=corrdat, aes(xintercept=rsq.mean, color = condition),linetype="dashed", size=1) +
  scale_colour_grey() + guides(colour=FALSE) +
  geom_vline(data=correlations, aes(xintercept=mean(Permutation)+2*sd(Permutation)),linetype=3, size=0.5, colour = "black") +
  scale_x_continuous(name = "R-squared", limits = c(0,1), breaks = seq(0,1,0.25), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion",limits = c(0,0.5), breaks = seq(0,0.5,0.1), expand = c(0,0)) +
  theme_classic()
figreg_sp_Performance <- figreg_sp_Performance + theme_reg

# plotting ----------------------------------------------------------------
multiplot(figcor_pre_PRCS, figcor_sp_Performance, figreg_pre_PRCS, figreg_sp_Performance, cols = 2)