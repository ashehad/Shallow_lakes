## Parameter - Description - (unit) - initial value (when applicable)
# BIOFL - Floating plant biomass - (g dw m-2)
# BIOSUBM - Submerged plant biomass - (g dw m-2)
# NINORG - Maximum nutrient conc. - (mg N liter-1)
## Setting for running a submerged-floating plant competition model
getwd()
setwd <- getwd()
dat <- read.csv("inputs.csv", header=TRUE,sep=",")
source("model_dynamics.R")

library(tidyverse)

biofl_list <- c(1)
biosunm_list <- c(1)

startcon_longData <- NULL

  nlevel <- 4
  NINORG <- nlevel
  days <- 365*2;
  nRecs<-nrow(dat); 
  REPS <- days;
  
    case <-  1
    BIOFL <- biofl_list[case]; BIOSUBM <- biosunm_list[case];
    flbio <- rep(0,days); submbio <- rep(0,days); flfraction <- rep(0,days);
	  nlF <- rep(0,days); nlS <- rep(0,days);
    llF <- rep(0,days); llS <- rep(0,days); Teff <- rep(0,days);
    
    results <- NULL
    # Running the model for each of the parameter values
    
    for(j in nrow(dat)){
      rf <- dat$RF[j]; rs <- dat$RS[j];
      lf <- dat$LF[j]; ls <- dat$LS[j];
      af <- dat$AF[j];
      as <- dat$AS[j]; hf <- dat$HF[j];
      hs <- 0; 
      qf <- dat$QF[j];
      qs <- dat$QS[j];
      bf<-dat$BF[j];
      FS <- biomass(BIOFL, BIOSUBM, NINORG, REPS, rf,rs,lf,ls,af,as,hf,hs,qf,qs,bf);
      flbio <- FS$FloatB; submbio <- FS$SubmB;
	    nlF <- FS$N_F; nlS <- FS$N_S;
      llF <- FS$L_F; llS <- FS$L_S; teff <- FS$T_eff
      flfraction <- FS$FloatB/(FS$FloatB+FS$SubmB);
      rm(rf,rs,lf,ls,af,as,hf,hs,qf,qs,bf,FS) #rm()remove object
      }
    output <- data.frame("day" = 1:days, "flfraction" = flfraction, "NlimitF" = nlF, "NlimitS" = nlS, 
                         "LightLF" = llF, "LightLS" = llS, "FLBIO" = flbio, "SUBMBIO" = submbio, "Teff"=teff,
                         "nlevel" = rep(nlevel, days), "case" = rep(case,days))
    startcon_longData <- rbind(startcon_longData, output)