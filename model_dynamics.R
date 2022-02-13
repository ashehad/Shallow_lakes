# Parameter - Description - (unit) - initial value (when applicable)
# BIOFL - Floating plant biomass - (g dw m-2)
# BIOSUBM - Submerged plant biomass - (g dw m-2)
# NINORG - Maximum nutrient conc. - (mg N liter-1)
# n - Total inorganic N conc. - (mg N liter-1)
# AF - density where r is 50% - ((g dw m-2)-1) - 0.01
# AS - density where r is 50% - (g dw m-2)-1 - 0.01
# B - scale for shading plants - ((g dw m-2)-1) - 0.02
# HF - half saturation conc. - (mg N liter-1) - 0.2
# HS - half saturation conc. - (mg N liter-1) - 0
# LF - biomass loss - (day-1) - 0.05
# LS - biomass loss - (day-1) - 0.05
# QF - effect of plants on N conc. in water - ((g dw m-2)-1) - 0.005
# QS - effect of plants on N conc. in water - ((g dw m-2)-1) - 0.075
# RF - max growth rate - (day-1) - 0.5
# RS - max growth rate - (day-1) - 0.5
# W - water column - (NA) - 0

biomass <- function(BIOFL, BIOSUBM, NINORG, REPS, rf,rs,lf,ls,af,as,hf,hs,qf,qs,bf){
  
  FLBIO <- rep(0,REPS + 1); SUBMBIO <- rep(0,REPS + 1); 
  FLBIONEW <- rep(0,REPS); SUBMBIONEW <- rep(0,REPS); n <- rep(0,REPS);
  n_f <- rep(0,REPS); n_s <- rep(0,REPS); l_f <- rep(0,REPS); l_s <- rep(0,REPS);
  t_l <- rep(0,REPS);
  
   #climate includes site-level data for years
   climate <- read.csv("climate.csv", header=TRUE,sep=",")
   Tair <- climate$TMAX[1:REPS]; Topt <- 26; Tmin <- 5;     
   Tlim <- (Tair-Tmin)/(Topt-Tmin);
   Tlim[Tlim < 0] <- 0;
	
  for(j in 1:REPS){
    
    FLBIO[1] <- BIOFL; SUBMBIO[1] <- BIOSUBM    
    n[j] <- NINORG/(1.0+qs*SUBMBIO[j]+qf*FLBIO[j])
	
	  yr <- ceiling(j/365)
    k <- j - (yr-1)*365  
	  rf1 = rf * Tlim[k]
	  rs1 = rs * Tlim[k]
    
	  FLBIONEW[j] <- rf1*FLBIO[j]*n[j]/(n[j]+hf)*af/(af+FLBIO[j]) - lf*FLBIO[j]+FLBIO[j]
	  SUBMBIONEW[j] <- rs1*SUBMBIO[j]*n[j]/(n[j]+hs)*as/(as+SUBMBIO[j])*1/(1+bf*FLBIO[j]) - ls*SUBMBIO[j] + SUBMBIO[j]
    
    FLBIO[j+1] <- FLBIONEW[j];
    SUBMBIO[j+1] <- SUBMBIONEW[j];
    
	  n_f[j] <- n[j]/(n[j]+hf); n_s[j] <- n[j]/(n[j]+hs);
    l_f[j] <- af/(af+FLBIO[j]); l_s[j] <- as/(as+SUBMBIO[j]); 
    t_l[j] <- rf1/rf
    
  }
  
  biomass <-list("FloatB"=FLBIONEW, "SubmB"=SUBMBIONEW, "N_F"=n_f, "N_S"=n_s, 
                 "L_F"=l_f, "L_S"=l_s, "T_eff"=t_l);
}