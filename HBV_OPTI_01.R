#///////////////////////////////////////////////////////////////////////////////////////////
# HBV-TEC-96 HYDROLOGICAL MODEL
#///////////////////////////////////////////////////////////////////////////////////////////
# VERSION: 1.0.2016
# Instituto Tecnologico de Costa Rica (www.tec.ac.cr)
# Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)
# Luis Alexander Calvo-V (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)
# This script is structured in R (www.r-project.org)
# General purpose: Implementation of the SMHI-HBV-96 Hydrological Model on R
# Input files: "hbvtecptq.txt", "hbvtecpar.txt", "hbvtecevap.txt", "hbvtecattri.txt" 
# Output files: "hbvtecptq_desc_hbv.csv", "hbvtecout_hbv.csv", "hbvtecout_desc_hbv.csv",
# "hbvqtecsim.csv", "hbvqteceff.csv"
#
#///////////////////////////////////////////////////////////////////////////////////////////
# REFERENCES:
#
# Aghakouchak, A., Habib, E. 2010. Application of a conceptual hydrologic model in 
# teaching hydrologic processes. Int J Engng Ed 26: 963-973.
# 
# Bergström, S. 1995. The HBV Model. In: V.P. Singh (editor) Computer Models of Watershed
# Hydrology. Water Resources Publications, Highlands Ranch, Colorado pp. 443- 470.
#
# Herman, J.D., Reed, P.M., Wagener, T. 2013. Time-varying sensitivity analysis 
# clarifies the effects of watershed model formulation on model behavior.
# Water Resour Res 49: 1400-1414.
#
# Lindström, G., Johansson, B., Persson, M., Gardelin, M., Bergström, S. 1997. Development
# and test of the distributed HBV-96 hydrological model. J Hydrol 201: 272-288.
# 
# Seibert, J. Vis, MJP. 2012. Teaching hydrological modeling with a user-friendly 
# catchment-runoff-model software package. Hydrol Earth Syst Sci 16: 3315-3325.
#
#///////////////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# CRAN libraries are loaded
require(alabama)
require(BB)
require(DEoptim)
require(DescTools)
require(dfoptim)
require(doParallel)
require(doSNOW)
require(FME)
require(foreach)
require(GA)
require(GenSA)
require(ggplot2)
require(graphics)
require(grid)
require(gridExtra)
require(hydromad)
require(hydroPSO)
require(lubridate)
require(MASS)
require(minpack.lm)
require(minqa)
require(nloptr)
require(optimx)
require(parallel)
require(pastecs)
require(PolynomF)
require(pracma)
require(Rcgmin)
require(RColorBrewer)
require(reshape)
require(rgenoud)
require(rootSolve)
require(Rsolnp)
require(Rvmmin)
require(snow)
require(visreg)

# Working directory is defined
setwd("C:/DATOS/R_ITC/HBV_TEC/HBV_TEC_96_VECTOR")

# /////////////////////////////////////////////////////////
# BLOCK: Parallel functions are defined
# /////////////////////////////////////////////////////////

# The cluster is defined with a number of cores
cores <- 8
cluster <- makeCluster(cores)  

# The cluster is registered
registerDoParallel(cluster)  

# /////////////////////////////////////////////////////////
# BLOCK: Custom functions are created
# /////////////////////////////////////////////////////////

# /////////////////////////////////////////////////////////////
f.HBV01 <- function(par) { # BLOCK: Main mass-balance loop
 # /////////////////////////////////////////////////////////////
 identity(1)
 # Parameters container is transformed from vector to list
 par <- as.list(par)
 
 # /////////////////////////////////////////////////////////
 # BLOCK: Defining model state-variables
 # /////////////////////////////////////////////////////////
 
 # Model state-variables containers are defined
 MASSIN <- vector() # mass entering the loop as precipitation [mm/T]
 Rstore <- vector() # mass storage [mm/T]
 SM <- vector() # soil moisture [mm/T]
 AET <- vector() # actual evapotranspiration [mm/T]
 PET <- vector() # potential evapotranspiration [mm/T]
 R <- vector() # recharge [mm/T]
 SUZ <- vector() # storage in the upper zone [mm/T]
 SLZ <- vector() # storage in the lower zone [mm/T]
 Q0 <- vector() # mass coming from SUZ [mm/T] (if SUZ > par$uzl)
 Q1 <- vector() # mass coming from SUZ [mm/T]
 Q2 <- vector() # mass coming from SLZ [mm/T]
 QSR <- vector() # mass coming from surface runoff [mm/T] if (SM >= par$fc)
 MASSOUT <- vector() # mass exiting the loop [mm/T]
 DELTAMASS <- vector() # delta of mass within the loop (MASSIN - MASSOUT) [mm/T]
 QOT <- vector() # total mass summation (Q0 + Q1 + Q2) [mm/T]
 QOUL <- vector() # mass summation of upper and lower zones (Q1 + Q2) [mm/T]
 QOSR <- vector() # mass summation of surface runoff and par$uzl (Q0 + QSR) [mm/T]
 
 # Model transformation-variables containers are defined
 mxrelweisum <- 0 # par$maxbas relative weigths summation [unitless]
 mxabsweisum <- 0 # par$maxbas absolute weigths summation [unitless]
 mxrelwei <- vector() # par$maxbas relative weigth [unitless]
 mxabswei <- vector() # par$maxbas absolute weigth [unitless]
 QSIM <- vector() # model simulated mass [mm/T]
 QSIM3 <- vector() # model simulated mass [m3/sec]
 
 # Length of relevant vectors is defined
 SUZ <- c(rep(0, countermain))
 SM <- c(rep(0, countermain))
 Q0 <- c(rep(0, countermain))
 
 # MASSINf is defined [mm/T] (entering the loop as precipitation/snow)
 MASSIN <- MASSINf
 
 # Main mass-balance loop is initialized
 for (i in 1 : countermain) {
  
  # ----------------------------------------------------------
  # Evapotranspiration sorting
  # ----------------------------------------------------------
  
  # PET values are selected [mm/month] (according to month of the year and date)
  if (f_month[i] == 1) {
   PET[i] <- evap.jan }
  else if (f_month[i] == 2) {
   PET[i] <- evap.feb }
  else if (f_month[i] == 3) {
   PET[i] <- evap.mar }
  else if (f_month[i] == 4) {
   PET[i] <- evap.apr }
  else if (f_month[i] == 5) {
   PET[i] <- evap.may }
  else if (f_month[i] == 6) {
   PET[i] <- evap.jun }
  else if (f_month[i] == 7) {
   PET[i] <- evap.jul }
  else if (f_month[i] == 8) {
   PET[i] <- evap.aug }
  else if (f_month[i] == 9) {
   PET[i] <- evap.set }
  else if (f_month[i] == 10) {
   PET[i] <- evap.oct }
  else if (f_month[i] == 11) {
   PET[i] <- evap.nov }
  else if (f_month[i] == 12) {
   PET[i] <- evap.dic }
  
  # ----------------------------------------------------------
  # Soil moisture routine
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) { 
   SM[i] <- SMtemp
  } else { # This is for i > 1; meaning day > 1
   SMtemp2 <- SM[i - 1]
   SM[i] <- SMtemp2 
  }
  
  # The process continues...
  if (SM[i] >= par$fc) {
   # R[i] <- PRECtemp + (SM[i] - par$fc) # original HBV statement
   QSR[i] <- SM[i] - par$fc # new HBV-TEC statement
   R[i] <- PRECtemp[i] + QSR[i] # new HBV-TEC statement
   SM[i] <- par$fc
  } else {
   QSR[i] <- 0 # new HBV-TEC statement
   Rstore[i] <- PRECtemp[i] * (1 - ((SM[i] / par$fc) ^ par$beta))
   SM[i] <- SM[i] + Rstore[i]
   R[i] <- PRECtemp[i] - Rstore[i]
   if (SM[i] > par$fc) {
    R[i] <- SM[i] - par$fc
    SM[i] <- par$fc
   }
  }
  
  # AET compensation within soil moisture routine
  if ((SM[i] / par$fc) > par$lp) {
   AET[i] <- PET[i]
  } else {
   AET[i] <- (SM[i] / (par$fc * par$lp)) * PET[i]
  }
  if (AET[i] < 0) {
   AET[i] <- 0
  }
  if (SM[i] > AET[i]) {
   SM[i] <- SM[i] - AET[i]
  } else {
   AET[i] <- SM[i]
   SM[i] <- 0
  } 
  
  # ----------------------------------------------------------
  # Responce function routine
  # ----------------------------------------------------------
  
  # Storage in the Upper Zone (SUZ)
  # This is for i = 1; meaning day = 1; meaning initialization
  if(i == 1) {
   SUZ[i] <- SUZtemp + R[i]
  } else { # This is for i > 1; meaning day > 1
   SUZ[i] <- SUZ[i - 1] + R[i]
  }
  
  # The process continues...
  if (SUZ[i] > par$uzl) {
   Q0[i] <- (par$k0 * (SUZ[i] - par$uzl))
   #SUZ[i] <- SUZ[i] - Q0[i]
  } else {
   Q0[i] <- Q0[i]
  }
  # MUCHO OJO!!!!!!!!!!!!!!!
  SUZ[i] <- SUZ[i] - Q0[i]
  
  # Percolation threshold
  if (SUZ[i] > par$perc) {
   SUZ[i] <- SUZ[i] - par$perc
   Q1[i] <- par$k1 * (SUZ[i])
   SUZ[i] <- SUZ[i] - Q1[i]
   
   # Storage in the Lower Zone (SLZ)
   if (i == 1) {
    SLZ[i] <- SLZtemp + par$perc
   } else { 
    SLZ[i] <- SLZ[i - 1] + par$perc
   }
  } else {
   Q1[i] <- 0 # New HBV-TEC statement
   if (i == 1) {
    SLZ[i] <- SLZtemp + SUZ[i]
   } else { 
    SLZ[i] <- SLZ[i - 1] + SUZ[i]
   }
   SUZ[i] <- 0
  }
  if (SLZ[i] > 0) {
   Q2[i] <- par$k2 * SLZ[i]
   SLZ[i] <- SLZ[i] - Q2[i]  
  } else {
   Q2[i] <- 0
   SLZ[i] <- SLZ[i - 1] # a warning should be issued here!!!!! 
  }
  
  # ----------------------------------------------------------
  # Check mass-balance loop 
  # ----------------------------------------------------------
  
  MASSOUT[i] <- ((SM[i] - SM[i - 1]) +
                  (SUZ[i] - SUZ[i - 1]) +
                  (SLZ[i] - SLZ[i - 1]) +
                  (Q0[i]) +
                  (Q1[i]) +
                  (Q2[i]) +
                  (AET[i]))  
  
  # Delta of mass within the loop (MASSIN - MASSOUT) is calculated [mm/T]  
  DELTAMASS[i] <- MASSIN[i] - MASSOUT[i]
  
  # Total mass summation (Q0 + Q1 + Q2) is calculated [mm/T]      
  QOT[i] <- Q0[i] + Q1[i] + Q2[i] # original HBV statement
  
  # Mass summation of upper and lower zones (Q1 + Q2) is calculated [mm/T]
  QOUL[i] <- Q1[i] + Q2[i] # new HBV-TEC statement
  
  # Mass summation of surface runoff and par$uzl (Q0 + QSR) is calculated [mm/T]
  QOSR[i] <- Q0[i] + QSR[i] # new HBV-TEC statement
  
 }
 
 # ----------------------------------------------------------
 # Main mass-balance loop is closed
 # ----------------------------------------------------------
 
 # /////////////////////////////////////////////////////////////
 # BLOCK: Transformation function routine (par$maxbas-t) Q0+Q1+Q2
 # /////////////////////////////////////////////////////////////
 
 # par$maxbas counters are defined
 mx2 <- (par$maxbas / 2) # half par$maxbas
 maxbascount <- ceiling(par$maxbas) # par$maxbas counter based on ceiling
 
 # ----------------------------------------------------------
 # par$maxbas relative weights are calculated
 # ----------------------------------------------------------
 
 # If par$maxbas < 2
 if (par$maxbas < 2) {
  for (i in 1 : (maxbascount)) {
   if (i <= par$maxbas) {
    mxrelwei[i] <- (par$maxbas / 2) 
   } else {
    mxrelwei[i] <- (((par$maxbas - i) + 1.0)) 
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum
  }
 }
 
 # If par$maxbas >= 2
 if (par$maxbas >= 2) {
  for (i in 1 : (maxbascount)) {    
   if (i <= mx2) {
    mxrelwei[i] <- (i) 
   } else {
    mxrelwei[i] <- (((par$maxbas - i) + 1.0))
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum             
  }
 }
 
 # ----------------------------------------------------------
 # par$maxbas relative weights normalization loop
 # ----------------------------------------------------------
 
 # par$maxbas relative weights normalization loop is initialized
 for (i in 1 : (maxbascount)) {
  mxabswei[i] <- (mxrelwei[i]) / mxrelweisum
  mxabsweisum <- (mxabswei[i]) + mxabsweisum
 }
 # par$maxbas relative weights normalization loop is closed
 
 # Total duration of hydraulic routing + par$maxbas is defined
 routingdur <- ((length(QOT)) + maxbascount - 1 + 1)
 
 # A "qresult[i]" vector inizialized to cero for the duration of "routingdur" is created
 qresult <- c(rep(0, routingdur))
 
 # Duration of hydraulic routing ONLY is defined
 QOTdur <- (length(QOT))
 
 # ----------------------------------------------------------
 # Flow integration over-time loop
 # ----------------------------------------------------------
 
 # Flow integration over-time loop is initialized
 for (n in 1 : (QOTdur)) { # day of the year external Loop
  for (i in 1 : (maxbascount)) { # internal par$maxbas Loop                 
   qrouting <- ((QOT[n]) * (mxabswei[i]))
   qresult[n + i - 1] <- qresult[n + i - 1] + qrouting
  }
 } 
 # Flow integration over-time loop is closed
 
 # Define QSIM based on "qresult[i]" for the duration of countermain ONLY 
 QSIM <- (qresult[1 : countermain])
 
 # Convert QSIM from mm/day to m3/sec
 QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 
 # If df.attri data.frame routing = 0; par$maxbas is IGNORED and NO tranformation is calculated
 if (routing == 0) {
  QSIM <- Q0 + Q1 + Q2 + QSR
  QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 }
 
 # ******************************************************
 # AS REQUESTED by modCost {FME}
 yy <- QSIM
 xx <- seq(1, (length(QSIM)), by = 1)
 return(data.frame(xx,yy))
 # ******************************************************
 
 # ******************************************************
 # AS REQUESTED by optim and optimx
 # return(sum((outQOBS - QSIM)^2))
 # END of f.HBV01 function
 # ******************************************************
 
 # ******************************************************
 # AS REQUESTED by nls.lm {minpack.lm}
 # return(((outQOBS - QSIM))) 
 # ******************************************************
 
}

# ============================================================================
# ============================================================================
# ============================================================================

# /////////////////////////////////////////////////////////////
f.HBV02 <- function(par) { # BLOCK: Main mass-balance loop
 # /////////////////////////////////////////////////////////////
 identity(1)  
 # /////////////////////////////////////////////////////////
 # BLOCK: Defining model state-variables
 # /////////////////////////////////////////////////////////
 
 # Model state-variables containers are defined
 MASSIN <- vector() # mass entering the loop as precipitation [mm/T]
 Rstore <- vector() # mass storage [mm/T]
 SM <- vector() # soil moisture [mm/T]
 AET <- vector() # actual evapotranspiration [mm/T]
 PET <- vector() # potential evapotranspiration [mm/T]
 R <- vector() # recharge [mm/T]
 SUZ <- vector() # storage in the upper zone [mm/T]
 SLZ <- vector() # storage in the lower zone [mm/T]
 Q0 <- vector() # mass coming from SUZ [mm/T] (if SUZ > par[5])
 Q1 <- vector() # mass coming from SUZ [mm/T]
 Q2 <- vector() # mass coming from SLZ [mm/T]
 QSR <- vector() # mass coming from surface runoff [mm/T] if (SM >= par[1])
 MASSOUT <- vector() # mass exiting the loop [mm/T]
 DELTAMASS <- vector() # delta of mass within the loop (MASSIN - MASSOUT) [mm/T]
 QOT <- vector() # total mass summation (Q0 + Q1 + Q2) [mm/T]
 QOUL <- vector() # mass summation of upper and lower zones (Q1 + Q2) [mm/T]
 QOSR <- vector() # mass summation of surface runoff and par[5] (Q0 + QSR) [mm/T]
 
 # Model transformation-variables containers are defined
 mxrelweisum <- 0 # par[9] relative weigths summation [unitless]
 mxabsweisum <- 0 # par[9] absolute weigths summation [unitless]
 mxrelwei <- vector() # par[9] relative weigth [unitless]
 mxabswei <- vector() # par[9] absolute weigth [unitless]
 QSIM <- vector() # model simulated mass [mm/T]
 QSIM3 <- vector() # model simulated mass [m3/sec]
 
 # Length of relevant vectors is defined
 SUZ <- c(rep(0, countermain))
 SM <- c(rep(0, countermain))
 Q0 <- c(rep(0, countermain))
 
 # MASSINf is defined [mm/T] (entering the loop as precipitation/snow)
 MASSIN <- MASSINf
 
 # Main mass-balance loop is initialized
 for (i in 1 : countermain) {
  
  # ----------------------------------------------------------
  # Evapotranspiration sorting
  # ----------------------------------------------------------
  
  # PET values are selected [mm/month] (according to month of the year and date)
  if (f_month[i] == 1) {
   PET[i] <- evap.jan }
  else if (f_month[i] == 2) {
   PET[i] <- evap.feb }
  else if (f_month[i] == 3) {
   PET[i] <- evap.mar }
  else if (f_month[i] == 4) {
   PET[i] <- evap.apr }
  else if (f_month[i] == 5) {
   PET[i] <- evap.may }
  else if (f_month[i] == 6) {
   PET[i] <- evap.jun }
  else if (f_month[i] == 7) {
   PET[i] <- evap.jul }
  else if (f_month[i] == 8) {
   PET[i] <- evap.aug }
  else if (f_month[i] == 9) {
   PET[i] <- evap.set }
  else if (f_month[i] == 10) {
   PET[i] <- evap.oct }
  else if (f_month[i] == 11) {
   PET[i] <- evap.nov }
  else if (f_month[i] == 12) {
   PET[i] <- evap.dic }
  
  # ----------------------------------------------------------
  # Soil moisture routine
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) { 
   SM[i] <- SMtemp
  } else { # This is for i > 1; meaning day > 1
   SMtemp2 <- SM[i - 1]
   SM[i] <- SMtemp2 
  }
  
  # The process continues...
  if (SM[i] >= par[1]) {
   # R[i] <- PRECtemp + (SM[i] - par[1]) # original HBV statement
   QSR[i] <- SM[i] - par[1] # new HBV-TEC statement
   R[i] <- PRECtemp[i] + QSR[i] # new HBV-TEC statement
   SM[i] <- par[1]
  } else {
   QSR[i] <- 0 # new HBV-TEC statement
   Rstore[i] <- PRECtemp[i] * (1 - ((SM[i] / par[1]) ^ par[3]))
   SM[i] <- SM[i] + Rstore[i]
   R[i] <- PRECtemp[i] - Rstore[i]
   if (SM[i] > par[1]) {
    R[i] <- SM[i] - par[1]
    SM[i] <- par[1]
   }
  }
  
  # AET compensation within soil moisture routine
  if ((SM[i] / par[1]) > par[2]) {
   AET[i] <- PET[i]
  } else {
   AET[i] <- (SM[i] / (par[1] * par[2])) * PET[i]
  }
  if (AET[i] < 0) {
   AET[i] <- 0
  }
  if (SM[i] > AET[i]) {
   SM[i] <- SM[i] - AET[i]
  } else {
   AET[i] <- SM[i]
   SM[i] <- 0
  } 
  
  # ----------------------------------------------------------
  # Responce function routine
  # ----------------------------------------------------------
  
  # Storage in the Upper Zone (SUZ)
  # This is for i = 1; meaning day = 1; meaning initialization
  if(i == 1) {
   SUZ[i] <- SUZtemp + R[i]
  } else { # This is for i > 1; meaning day > 1
   SUZ[i] <- SUZ[i - 1] + R[i]
  }
  
  # The process continues...
  if (SUZ[i] > par[5]) {
   Q0[i] <- (par[6] * (SUZ[i] - par[5]))
   #SUZ[i] <- SUZ[i] - Q0[i]
  } else {
   Q0[i] <- Q0[i]
  }
  # MUCHO OJO!!!!!!!!!!!!!!!
  SUZ[i] <- SUZ[i] - Q0[i]
  
  # Percolation threshold
  if (SUZ[i] > par[4]) {
   SUZ[i] <- SUZ[i] - par[4]
   Q1[i] <- par[7] * (SUZ[i])
   SUZ[i] <- SUZ[i] - Q1[i]
   
   # Storage in the Lower Zone (SLZ)
   if (i == 1) {
    SLZ[i] <- SLZtemp + par[4]
   } else { 
    SLZ[i] <- SLZ[i - 1] + par[4]
   }
  } else {
   Q1[i] <- 0 # New HBV-TEC statement
   if (i == 1) {
    SLZ[i] <- SLZtemp + SUZ[i]
   } else { 
    SLZ[i] <- SLZ[i - 1] + SUZ[i]
   }
   SUZ[i] <- 0
  }
  if (SLZ[i] > 0) {
   Q2[i] <- par[8] * SLZ[i]
   SLZ[i] <- SLZ[i] - Q2[i]  
  } else {
   Q2[i] <- 0
   SLZ[i] <- SLZ[i - 1] # a warning should be issued here!!!!! 
  }
  
  # ----------------------------------------------------------
  # Check mass-balance loop 
  # ----------------------------------------------------------
  
  MASSOUT[i] <- ((SM[i] - SM[i - 1]) +
                  (SUZ[i] - SUZ[i - 1]) +
                  (SLZ[i] - SLZ[i - 1]) +
                  (Q0[i]) +
                  (Q1[i]) +
                  (Q2[i]) +
                  (AET[i]))  
  
  # Delta of mass within the loop (MASSIN - MASSOUT) is calculated [mm/T]  
  DELTAMASS[i] <- MASSIN[i] - MASSOUT[i]
  
  # Total mass summation (Q0 + Q1 + Q2) is calculated [mm/T]      
  QOT[i] <- Q0[i] + Q1[i] + Q2[i] # original HBV statement
  
  # Mass summation of upper and lower zones (Q1 + Q2) is calculated [mm/T]
  QOUL[i] <- Q1[i] + Q2[i] # new HBV-TEC statement
  
  # Mass summation of surface runoff and par[5] (Q0 + QSR) is calculated [mm/T]
  QOSR[i] <- Q0[i] + QSR[i] # new HBV-TEC statement
  
 }
 
 # ----------------------------------------------------------
 # Main mass-balance loop is closed
 # ----------------------------------------------------------
 
 # /////////////////////////////////////////////////////////////
 # BLOCK: Transformation function routine (par[9]-t) Q0+Q1+Q2
 # /////////////////////////////////////////////////////////////
 
 # par[9] counters are defined
 mx2 <- (par[9] / 2) # half par[9]
 maxbascount <- ceiling(par[9]) # par[9] counter based on ceiling
 
 # ----------------------------------------------------------
 # par[9] relative weights are calculated
 # ----------------------------------------------------------
 
 # If par[9] < 2
 if (par[9] < 2) {
  for (i in 1 : (maxbascount)) {
   if (i <= par[9]) {
    mxrelwei[i] <- (par[9] / 2) 
   } else {
    mxrelwei[i] <- (((par[9] - i) + 1.0)) 
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum
  }
 }
 
 # If par[9] >= 2
 if (par[9] >= 2) {
  for (i in 1 : (maxbascount)) {    
   if (i <= mx2) {
    mxrelwei[i] <- (i) 
   } else {
    mxrelwei[i] <- (((par[9] - i) + 1.0))
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum             
  }
 }
 
 # ----------------------------------------------------------
 # par[9] relative weights normalization loop
 # ----------------------------------------------------------
 
 # par[9] relative weights normalization loop is initialized
 for (i in 1 : (maxbascount)) {
  mxabswei[i] <- (mxrelwei[i]) / mxrelweisum
  mxabsweisum <- (mxabswei[i]) + mxabsweisum
 }
 # par[9] relative weights normalization loop is closed
 
 # Total duration of hydraulic routing + par[9] is defined
 routingdur <- ((length(QOT)) + maxbascount - 1 + 1)
 
 # A "qresult[i]" vector inizialized to cero for the duration of "routingdur" is created
 qresult <- c(rep(0, routingdur))
 
 # Duration of hydraulic routing ONLY is defined
 QOTdur <- (length(QOT))
 
 # ----------------------------------------------------------
 # Flow integration over-time loop
 # ----------------------------------------------------------
 
 # Flow integration over-time loop is initialized
 for (n in 1 : (QOTdur)) { # day of the year external Loop
  for (i in 1 : (maxbascount)) { # internal par[9] Loop                 
   qrouting <- ((QOT[n]) * (mxabswei[i]))
   qresult[n + i - 1] <- qresult[n + i - 1] + qrouting
  }
 } 
 # Flow integration over-time loop is closed
 
 # Define QSIM based on "qresult[i]" for the duration of countermain ONLY 
 QSIM <- (qresult[1 : countermain])
 
 # Convert QSIM from mm/day to m3/sec
 QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 
 # If df.attri data.frame routing = 0; par[9] is IGNORED and NO tranformation is calculated
 if (routing == 0) {
  QSIM <- Q0 + Q1 + Q2 + QSR
  QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 }
 
 # ******************************************************
 # AS REQUESTED by optim and optimx
 return(sum((outQOBS - QSIM)^2))
 # ******************************************************
 # END of f.HBV01 function
 
 # ******************************************************
 # AS REQUESTED by ga{GA}
 # Mmean <- mean(outQOBS)
 # return(  1  - ((sum((outQOBS - QSIM)^2)) / (sum((outQOBS - Mmean)^2)))     )
 # ******************************************************
 # END of f.HBV01 function
 
}

# ============================================================================
# ============================================================================
# ============================================================================

# ============================================================================
# ============================================================================
# ============================================================================

# /////////////////////////////////////////////////////////////
f.HBV03 <- function(par) { # BLOCK: Main mass-balance loop
 # /////////////////////////////////////////////////////////////
 identity(1)
 # /////////////////////////////////////////////////////////
 # BLOCK: Defining model state-variables
 # /////////////////////////////////////////////////////////
 
 # Model state-variables containers are defined
 MASSIN <- vector() # mass entering the loop as precipitation [mm/T]
 Rstore <- vector() # mass storage [mm/T]
 SM <- vector() # soil moisture [mm/T]
 AET <- vector() # actual evapotranspiration [mm/T]
 PET <- vector() # potential evapotranspiration [mm/T]
 R <- vector() # recharge [mm/T]
 SUZ <- vector() # storage in the upper zone [mm/T]
 SLZ <- vector() # storage in the lower zone [mm/T]
 Q0 <- vector() # mass coming from SUZ [mm/T] (if SUZ > par[5])
 Q1 <- vector() # mass coming from SUZ [mm/T]
 Q2 <- vector() # mass coming from SLZ [mm/T]
 QSR <- vector() # mass coming from surface runoff [mm/T] if (SM >= par[1])
 MASSOUT <- vector() # mass exiting the loop [mm/T]
 DELTAMASS <- vector() # delta of mass within the loop (MASSIN - MASSOUT) [mm/T]
 QOT <- vector() # total mass summation (Q0 + Q1 + Q2) [mm/T]
 QOUL <- vector() # mass summation of upper and lower zones (Q1 + Q2) [mm/T]
 QOSR <- vector() # mass summation of surface runoff and par[5] (Q0 + QSR) [mm/T]
 
 # Model transformation-variables containers are defined
 mxrelweisum <- 0 # par[9] relative weigths summation [unitless]
 mxabsweisum <- 0 # par[9] absolute weigths summation [unitless]
 mxrelwei <- vector() # par[9] relative weigth [unitless]
 mxabswei <- vector() # par[9] absolute weigth [unitless]
 QSIM <- vector() # model simulated mass [mm/T]
 QSIM3 <- vector() # model simulated mass [m3/sec]
 
 # Length of relevant vectors is defined
 SUZ <- c(rep(0, countermain))
 SM <- c(rep(0, countermain))
 Q0 <- c(rep(0, countermain))
 
 # MASSINf is defined [mm/T] (entering the loop as precipitation/snow)
 MASSIN <- MASSINf
 
 # Main mass-balance loop is initialized
 for (i in 1 : countermain) {
  
  # ----------------------------------------------------------
  # Evapotranspiration sorting
  # ----------------------------------------------------------
  
  # PET values are selected [mm/month] (according to month of the year and date)
  if (f_month[i] == 1) {
   PET[i] <- evap.jan }
  else if (f_month[i] == 2) {
   PET[i] <- evap.feb }
  else if (f_month[i] == 3) {
   PET[i] <- evap.mar }
  else if (f_month[i] == 4) {
   PET[i] <- evap.apr }
  else if (f_month[i] == 5) {
   PET[i] <- evap.may }
  else if (f_month[i] == 6) {
   PET[i] <- evap.jun }
  else if (f_month[i] == 7) {
   PET[i] <- evap.jul }
  else if (f_month[i] == 8) {
   PET[i] <- evap.aug }
  else if (f_month[i] == 9) {
   PET[i] <- evap.set }
  else if (f_month[i] == 10) {
   PET[i] <- evap.oct }
  else if (f_month[i] == 11) {
   PET[i] <- evap.nov }
  else if (f_month[i] == 12) {
   PET[i] <- evap.dic }
  
  # ----------------------------------------------------------
  # Soil moisture routine
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) { 
   SM[i] <- SMtemp
  } else { # This is for i > 1; meaning day > 1
   SMtemp2 <- SM[i - 1]
   SM[i] <- SMtemp2 
  }
  
  # The process continues...
  if (SM[i] >= par[1]) {
   # R[i] <- PRECtemp + (SM[i] - par[1]) # original HBV statement
   QSR[i] <- SM[i] - par[1] # new HBV-TEC statement
   R[i] <- PRECtemp[i] + QSR[i] # new HBV-TEC statement
   SM[i] <- par[1]
  } else {
   QSR[i] <- 0 # new HBV-TEC statement
   Rstore[i] <- PRECtemp[i] * (1 - ((SM[i] / par[1]) ^ par[3]))
   SM[i] <- SM[i] + Rstore[i]
   R[i] <- PRECtemp[i] - Rstore[i]
   if (SM[i] > par[1]) {
    R[i] <- SM[i] - par[1]
    SM[i] <- par[1]
   }
  }
  
  # AET compensation within soil moisture routine
  if ((SM[i] / par[1]) > par[2]) {
   AET[i] <- PET[i]
  } else {
   AET[i] <- (SM[i] / (par[1] * par[2])) * PET[i]
  }
  if (AET[i] < 0) {
   AET[i] <- 0
  }
  if (SM[i] > AET[i]) {
   SM[i] <- SM[i] - AET[i]
  } else {
   AET[i] <- SM[i]
   SM[i] <- 0
  } 
  
  # ----------------------------------------------------------
  # Responce function routine
  # ----------------------------------------------------------
  
  # Storage in the Upper Zone (SUZ)
  # This is for i = 1; meaning day = 1; meaning initialization
  if(i == 1) {
   SUZ[i] <- SUZtemp + R[i]
  } else { # This is for i > 1; meaning day > 1
   SUZ[i] <- SUZ[i - 1] + R[i]
  }
  
  # The process continues...
  if (SUZ[i] > par[5]) {
   Q0[i] <- (par[6] * (SUZ[i] - par[5]))
   #SUZ[i] <- SUZ[i] - Q0[i]
  } else {
   Q0[i] <- Q0[i]
  }
  # MUCHO OJO!!!!!!!!!!!!!!!
  SUZ[i] <- SUZ[i] - Q0[i]
  
  # Percolation threshold
  if (SUZ[i] > par[4]) {
   SUZ[i] <- SUZ[i] - par[4]
   Q1[i] <- par[7] * (SUZ[i])
   SUZ[i] <- SUZ[i] - Q1[i]
   
   # Storage in the Lower Zone (SLZ)
   if (i == 1) {
    SLZ[i] <- SLZtemp + par[4]
   } else { 
    SLZ[i] <- SLZ[i - 1] + par[4]
   }
  } else {
   Q1[i] <- 0 # New HBV-TEC statement
   if (i == 1) {
    SLZ[i] <- SLZtemp + SUZ[i]
   } else { 
    SLZ[i] <- SLZ[i - 1] + SUZ[i]
   }
   SUZ[i] <- 0
  }
  if (SLZ[i] > 0) {
   Q2[i] <- par[8] * SLZ[i]
   SLZ[i] <- SLZ[i] - Q2[i]  
  } else {
   Q2[i] <- 0
   SLZ[i] <- SLZ[i - 1] # a warning should be issued here!!!!! 
  }
  
  # ----------------------------------------------------------
  # Check mass-balance loop 
  # ----------------------------------------------------------
  
  MASSOUT[i] <- ((SM[i] - SM[i - 1]) +
                  (SUZ[i] - SUZ[i - 1]) +
                  (SLZ[i] - SLZ[i - 1]) +
                  (Q0[i]) +
                  (Q1[i]) +
                  (Q2[i]) +
                  (AET[i]))  
  
  # Delta of mass within the loop (MASSIN - MASSOUT) is calculated [mm/T]  
  DELTAMASS[i] <- MASSIN[i] - MASSOUT[i]
  
  # Total mass summation (Q0 + Q1 + Q2) is calculated [mm/T]      
  QOT[i] <- Q0[i] + Q1[i] + Q2[i] # original HBV statement
  
  # Mass summation of upper and lower zones (Q1 + Q2) is calculated [mm/T]
  QOUL[i] <- Q1[i] + Q2[i] # new HBV-TEC statement
  
  # Mass summation of surface runoff and par[5] (Q0 + QSR) is calculated [mm/T]
  QOSR[i] <- Q0[i] + QSR[i] # new HBV-TEC statement
  
 }
 
 # ----------------------------------------------------------
 # Main mass-balance loop is closed
 # ----------------------------------------------------------
 
 # /////////////////////////////////////////////////////////////
 # BLOCK: Transformation function routine (par[9]-t) Q0+Q1+Q2
 # /////////////////////////////////////////////////////////////
 
 # par[9] counters are defined
 mx2 <- (par[9] / 2) # half par[9]
 maxbascount <- ceiling(par[9]) # par[9] counter based on ceiling
 
 # ----------------------------------------------------------
 # par[9] relative weights are calculated
 # ----------------------------------------------------------
 
 # If par[9] < 2
 if (par[9] < 2) {
  for (i in 1 : (maxbascount)) {
   if (i <= par[9]) {
    mxrelwei[i] <- (par[9] / 2) 
   } else {
    mxrelwei[i] <- (((par[9] - i) + 1.0)) 
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum
  }
 }
 
 # If par[9] >= 2
 if (par[9] >= 2) {
  for (i in 1 : (maxbascount)) {    
   if (i <= mx2) {
    mxrelwei[i] <- (i) 
   } else {
    mxrelwei[i] <- (((par[9] - i) + 1.0))
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum             
  }
 }
 
 # ----------------------------------------------------------
 # par[9] relative weights normalization loop
 # ----------------------------------------------------------
 
 # par[9] relative weights normalization loop is initialized
 for (i in 1 : (maxbascount)) {
  mxabswei[i] <- (mxrelwei[i]) / mxrelweisum
  mxabsweisum <- (mxabswei[i]) + mxabsweisum
 }
 # par[9] relative weights normalization loop is closed
 
 # Total duration of hydraulic routing + par[9] is defined
 routingdur <- ((length(QOT)) + maxbascount - 1 + 1)
 
 # A "qresult[i]" vector inizialized to cero for the duration of "routingdur" is created
 qresult <- c(rep(0, routingdur))
 
 # Duration of hydraulic routing ONLY is defined
 QOTdur <- (length(QOT))
 
 # ----------------------------------------------------------
 # Flow integration over-time loop
 # ----------------------------------------------------------
 
 # Flow integration over-time loop is initialized
 for (n in 1 : (QOTdur)) { # day of the year external Loop
  for (i in 1 : (maxbascount)) { # internal par[9] Loop                 
   qrouting <- ((QOT[n]) * (mxabswei[i]))
   qresult[n + i - 1] <- qresult[n + i - 1] + qrouting
  }
 } 
 # Flow integration over-time loop is closed
 
 # Define QSIM based on "qresult[i]" for the duration of countermain ONLY 
 QSIM <- (qresult[1 : countermain])
 
 # Convert QSIM from mm/day to m3/sec
 QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 
 # If df.attri data.frame routing = 0; par[9] is IGNORED and NO tranformation is calculated
 if (routing == 0) {
  QSIM <- Q0 + Q1 + Q2 + QSR
  QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 }
 
 # ******************************************************
 # AS REQUESTED by optim and optimx
 # return(sum((outQOBS - QSIM)^2))
 # ******************************************************
 # END of f.HBV01 function
 
 # ******************************************************
 # AS REQUESTED by ga{GA}
 Mmean <- mean(outQOBS)
 return(  1  - ((sum((outQOBS - QSIM)^2)) / (sum((outQOBS - Mmean)^2)))     )
 # ******************************************************
 # END of f.HBV01 function
 
}

# ============================================================================
# ============================================================================
# ============================================================================

# /////////////////////////////////////////////////////////////
f.HBV04 <- function(parS) { # BLOCK: Main mass-balance loop
 # /////////////////////////////////////////////////////////////
 identity(1)
 # Parameters container is transformed from vector to list
 # par <- as.list(par)
 
 # /////////////////////////////////////////////////////////
 # BLOCK: Defining model state-variables
 # /////////////////////////////////////////////////////////
 
 # Model state-variables containers are defined
 MASSIN <- vector() # mass entering the loop as precipitation [mm/T]
 Rstore <- vector() # mass storage [mm/T]
 SM <- vector() # soil moisture [mm/T]
 AET <- vector() # actual evapotranspiration [mm/T]
 PET <- vector() # potential evapotranspiration [mm/T]
 R <- vector() # recharge [mm/T]
 SUZ <- vector() # storage in the upper zone [mm/T]
 SLZ <- vector() # storage in the lower zone [mm/T]
 Q0 <- vector() # mass coming from SUZ [mm/T] (if SUZ > par$uzl)
 Q1 <- vector() # mass coming from SUZ [mm/T]
 Q2 <- vector() # mass coming from SLZ [mm/T]
 QSR <- vector() # mass coming from surface runoff [mm/T] if (SM >= par$fc)
 MASSOUT <- vector() # mass exiting the loop [mm/T]
 DELTAMASS <- vector() # delta of mass within the loop (MASSIN - MASSOUT) [mm/T]
 QOT <- vector() # total mass summation (Q0 + Q1 + Q2) [mm/T]
 QOUL <- vector() # mass summation of upper and lower zones (Q1 + Q2) [mm/T]
 QOSR <- vector() # mass summation of surface runoff and par$uzl (Q0 + QSR) [mm/T]
 
 # Model transformation-variables containers are defined
 mxrelweisum <- 0 # par$maxbas relative weigths summation [unitless]
 mxabsweisum <- 0 # par$maxbas absolute weigths summation [unitless]
 mxrelwei <- vector() # par$maxbas relative weigth [unitless]
 mxabswei <- vector() # par$maxbas absolute weigth [unitless]
 QSIM <- vector() # model simulated mass [mm/T]
 QSIM3 <- vector() # model simulated mass [m3/sec]
 
 # Length of relevant vectors is defined
 SUZ <- c(rep(0, countermain))
 SM <- c(rep(0, countermain))
 Q0 <- c(rep(0, countermain))
 
 # MASSINf is defined [mm/T] (entering the loop as precipitation/snow)
 MASSIN <- MASSINf
 
 # Main mass-balance loop is initialized
 for (i in 1 : countermain) {
  
  # ----------------------------------------------------------
  # Evapotranspiration sorting
  # ----------------------------------------------------------
  
  # PET values are selected [mm/month] (according to month of the year and date)
  if (f_month[i] == 1) {
   PET[i] <- evap.jan
   parS["fcG1"] -> fc
   parS["lpG1"] -> lp
   parS["betaG1"] -> beta
   parS["percG1"] -> perc
   parS["uzlG1"] -> uzl
   parS["k0G1"] -> k0
   parS["k1G1"] -> k1
   parS["k2G1"] -> k2
   parS["maxbasG1"] -> maxbas }
  else if (f_month[i] == 2) {
   PET[i] <- evap.feb 
   parS["fcG2"]  -> fc
   parS["lpG2"]  -> lp
   parS["betaG2"]  -> beta
   parS["percG2"]  -> perc
   parS["uzlG2"]  -> uzl
   parS["k0G2"]  -> k0
   parS["k1G2"]  -> k1
   parS["k2G2"]  -> k2
   parS["maxbasG2"]  -> maxbas }
  else if (f_month[i] == 3) {
   PET[i] <- evap.mar 
   parS["fcG34"]  -> fc
   parS["lpG34"]  -> lp
   parS["betaG34"]  -> beta
   parS["percG34"]  -> perc
   parS["uzlG34"]  -> uzl
   parS["k0G34"]  -> k0
   parS["k1G34"]  -> k1
   parS["k2G34"]  -> k2
   parS["maxbasG34"]  -> maxbas }
  else if (f_month[i] == 4) {
   PET[i] <- evap.apr
   parS["fcG34"]  -> fc
   parS["lpG34"]  -> lp
   parS["betaG34"]  -> beta
   parS["percG34"]  -> perc
   parS["uzlG34"]  -> uzl
   parS["k0G34"]  -> k0
   parS["k1G34"]  -> k1
   parS["k2G34"]  -> k2
   parS["maxbasG34"]  -> maxbas }
  else if (f_month[i] == 5) {
   PET[i] <- evap.may 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 6) {
   PET[i] <- evap.jun 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 7) {
   PET[i] <- evap.jul 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 8) {
   PET[i] <- evap.aug 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 9) {
   PET[i] <- evap.set 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 10) {
   PET[i] <- evap.oct 
   parS["fcGX"]  -> fc
   parS["lpGX"]  -> lp
   parS["betaGX"]  -> beta
   parS["percGX"]  -> perc
   parS["uzlGX"]  -> uzl
   parS["k0GX"]  -> k0
   parS["k1GX"]  -> k1
   parS["k2GX"]  -> k2
   parS["maxbasGX"]  -> maxbas }
  else if (f_month[i] == 11) {
   PET[i] <- evap.nov 
   parS["fcG11"]  -> fc
   parS["lpG11"]  -> lp
   parS["betaG11"]  -> beta
   parS["percG11"]  -> perc
   parS["uzlG11"]  -> uzl
   parS["k0G11"]  -> k0
   parS["k1G11"]  -> k1
   parS["k2G11"]  -> k2
   parS["maxbasG11"]  -> maxbas }
  else if (f_month[i] == 12) {
   PET[i] <- evap.dic
   parS["fcG12"]  -> fc
   parS["lpG12"]  -> lp
   parS["betaG12"]  -> beta
   parS["percG12"]  -> perc
   parS["uzlG12"]  -> uzl
   parS["k0G12"]  -> k0
   parS["k1G12"]  -> k1
   parS["k2G12"]  -> k2
   parS["maxbasG12"]  -> maxbas }
  
  # Parameters container is transformed from vector to list
  par <- as.list(c(fc, lp, beta, perc, uzl, k0, k1, k2, maxbas))
  
  names(par) <- c("fc", "lp", "beta", "perc", "uzl", "k0", "k1", "k2", "maxbas")
  
  # ----------------------------------------------------------
  # Soil moisture routine
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) { 
   SM[i] <- SMtemp
  } else { # This is for i > 1; meaning day > 1
   SMtemp2 <- SM[i - 1]
   SM[i] <- SMtemp2 
  }
  
  # The process continues...
  if (SM[i] >= par$fc) {
   # R[i] <- PRECtemp + (SM[i] - par$fc) # original HBV statement
   QSR[i] <- SM[i] - par$fc # new HBV-TEC statement
   R[i] <- PRECtemp[i] + QSR[i] # new HBV-TEC statement
   SM[i] <- par$fc
  } else {
   QSR[i] <- 0 # new HBV-TEC statement
   Rstore[i] <- PRECtemp[i] * (1 - ((SM[i] / par$fc) ^ par$beta))
   SM[i] <- SM[i] + Rstore[i]
   R[i] <- PRECtemp[i] - Rstore[i]
   if (SM[i] > par$fc) {
    R[i] <- SM[i] - par$fc
    SM[i] <- par$fc
   }
  }
  
  # AET compensation within soil moisture routine
  if ((SM[i] / par$fc) > par$lp) {
   AET[i] <- PET[i]
  } else {
   AET[i] <- (SM[i] / (par$fc * par$lp)) * PET[i]
  }
  if (AET[i] < 0) {
   AET[i] <- 0
  }
  if (SM[i] > AET[i]) {
   SM[i] <- SM[i] - AET[i]
  } else {
   AET[i] <- SM[i]
   SM[i] <- 0
  } 
  
  # ----------------------------------------------------------
  # Responce function routine
  # ----------------------------------------------------------
  
  # Storage in the Upper Zone (SUZ)
  # This is for i = 1; meaning day = 1; meaning initialization
  if(i == 1) {
   SUZ[i] <- SUZtemp + R[i]
  } else { # This is for i > 1; meaning day > 1
   SUZ[i] <- SUZ[i - 1] + R[i]
  }
  
  # The process continues...
  if (SUZ[i] > par$uzl) {
   Q0[i] <- (par$k0 * (SUZ[i] - par$uzl))
   #SUZ[i] <- SUZ[i] - Q0[i]
  } else {
   Q0[i] <- Q0[i]
  }
  # MUCHO OJO!!!!!!!!!!!!!!!
  SUZ[i] <- SUZ[i] - Q0[i]
  
  # Percolation threshold
  if (SUZ[i] > par$perc) {
   SUZ[i] <- SUZ[i] - par$perc
   Q1[i] <- par$k1 * (SUZ[i])
   SUZ[i] <- SUZ[i] - Q1[i]
   
   # Storage in the Lower Zone (SLZ)
   if (i == 1) {
    SLZ[i] <- SLZtemp + par$perc
   } else { 
    SLZ[i] <- SLZ[i - 1] + par$perc
   }
  } else {
   Q1[i] <- 0 # New HBV-TEC statement
   if (i == 1) {
    SLZ[i] <- SLZtemp + SUZ[i]
   } else { 
    SLZ[i] <- SLZ[i - 1] + SUZ[i]
   }
   SUZ[i] <- 0
  }
  if (SLZ[i] > 0) {
   Q2[i] <- par$k2 * SLZ[i]
   SLZ[i] <- SLZ[i] - Q2[i]  
  } else {
   Q2[i] <- 0
   SLZ[i] <- SLZ[i - 1] # a warning should be issued here!!!!! 
  }
  
  # ----------------------------------------------------------
  # Check mass-balance loop 
  # ----------------------------------------------------------
  
  MASSOUT[i] <- ((SM[i] - SM[i - 1]) +
                  (SUZ[i] - SUZ[i - 1]) +
                  (SLZ[i] - SLZ[i - 1]) +
                  (Q0[i]) +
                  (Q1[i]) +
                  (Q2[i]) +
                  (AET[i]))  
  
  # Delta of mass within the loop (MASSIN - MASSOUT) is calculated [mm/T]  
  DELTAMASS[i] <- MASSIN[i] - MASSOUT[i]
  
  # Total mass summation (Q0 + Q1 + Q2) is calculated [mm/T]      
  QOT[i] <- Q0[i] + Q1[i] + Q2[i] # original HBV statement
  
  # Mass summation of upper and lower zones (Q1 + Q2) is calculated [mm/T]
  QOUL[i] <- Q1[i] + Q2[i] # new HBV-TEC statement
  
  # Mass summation of surface runoff and par$uzl (Q0 + QSR) is calculated [mm/T]
  QOSR[i] <- Q0[i] + QSR[i] # new HBV-TEC statement
  
 }
 
 # ----------------------------------------------------------
 # Main mass-balance loop is closed
 # ----------------------------------------------------------
 
 # /////////////////////////////////////////////////////////////
 # BLOCK: Transformation function routine (par$maxbas-t) Q0+Q1+Q2
 # /////////////////////////////////////////////////////////////
 
 # par$maxbas counters are defined
 mx2 <- (par$maxbas / 2) # half par$maxbas
 maxbascount <- ceiling(par$maxbas) # par$maxbas counter based on ceiling
 
 # ----------------------------------------------------------
 # par$maxbas relative weights are calculated
 # ----------------------------------------------------------
 
 # If par$maxbas < 2
 if (par$maxbas < 2) {
  for (i in 1 : (maxbascount)) {
   if (i <= par$maxbas) {
    mxrelwei[i] <- (par$maxbas / 2) 
   } else {
    mxrelwei[i] <- (((par$maxbas - i) + 1.0)) 
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum
  }
 }
 
 # If par$maxbas >= 2
 if (par$maxbas >= 2) {
  for (i in 1 : (maxbascount)) {    
   if (i <= mx2) {
    mxrelwei[i] <- (i) 
   } else {
    mxrelwei[i] <- (((par$maxbas - i) + 1.0))
   }
   mxrelwei2 <- (mxrelwei[i])
   mxrelweisum <- (mxrelwei[i]) + mxrelweisum             
  }
 }
 
 # ----------------------------------------------------------
 # par$maxbas relative weights normalization loop
 # ----------------------------------------------------------
 
 # par$maxbas relative weights normalization loop is initialized
 for (i in 1 : (maxbascount)) {
  mxabswei[i] <- (mxrelwei[i]) / mxrelweisum
  mxabsweisum <- (mxabswei[i]) + mxabsweisum
 }
 # par$maxbas relative weights normalization loop is closed
 
 # Total duration of hydraulic routing + par$maxbas is defined
 routingdur <- ((length(QOT)) + maxbascount - 1 + 1)
 
 # A "qresult[i]" vector inizialized to cero for the duration of "routingdur" is created
 qresult <- c(rep(0, routingdur))
 
 # Duration of hydraulic routing ONLY is defined
 QOTdur <- (length(QOT))
 
 # ----------------------------------------------------------
 # Flow integration over-time loop
 # ----------------------------------------------------------
 
 # Flow integration over-time loop is initialized
 for (n in 1 : (QOTdur)) { # day of the year external Loop
  for (i in 1 : (maxbascount)) { # internal par$maxbas Loop                 
   qrouting <- ((QOT[n]) * (mxabswei[i]))
   qresult[n + i - 1] <- qresult[n + i - 1] + qrouting
  }
 } 
 # Flow integration over-time loop is closed
 
 # Define QSIM based on "qresult[i]" for the duration of countermain ONLY 
 QSIM <- (qresult[1 : countermain])
 
 # Convert QSIM from mm/day to m3/sec
 QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 
 # If df.attri data.frame routing = 0; par$maxbas is IGNORED and NO tranformation is calculated
 if (routing == 0) {
  QSIM <- Q0 + Q1 + Q2 + QSR
  QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
 }
 
 # ******************************************************
 # AS REQUESTED by modCost {FME}
 #yy <- QSIM
 #xx <- seq(1, (length(QSIM)), by = 1)
 #return(data.frame(xx,yy))
 # ******************************************************
 
 # ******************************************************
 # AS REQUESTED by optim and optimx
 #return(sum((outQOBS - QSIM)^2))
 # END of f.HBV01 function
 # ******************************************************
 
 # ******************************************************
 # AS REQUESTED by nls.lm {minpack.lm}
 return(((outQOBS - QSIM))) 
 # ******************************************************
 
}

# ============================================================================
# ============================================================================
# ============================================================================

# /////////////////////////////////////////////////////////
# BLOCK: Creating and organizing input data.frames
# /////////////////////////////////////////////////////////

# Precipitation, temperature and Q-observed data.frame is loaded
df.ptq <- read.delim ("hbvtecptq_Y16.txt", header = TRUE, sep = "\t")

# Descriptive statistics are requested and rounded to 4 decimals
df.ptq.desc <- round((as.data.frame(stat.desc(df.ptq[, 2 : 4]))),4)

# Model input-parameters data.frame is loaded
df.param <- read.delim ("hbvtecpar.txt", header = TRUE, sep = "\t")

# Monthly potential evapotranspiration data.frame is loaded
df.evap <- read.delim ("hbvtecevap.txt", header = TRUE, sep = "\t")

# Watershed attributes data.frame is loaded
df.attri <- read.delim ("hbvtecattri.txt", header = TRUE, sep = "\t")

# /////////////////////////////////////////////////////////
# BLOCK: Loading input variables and model parameters
# /////////////////////////////////////////////////////////

# Model input parameters are loaded
fc <- df.param [1, 2] # soil field capacity [mm]
lp <- df.param [2, 2] # threshold at which AET reaches PET [unitless]
beta <- df.param [3, 2] # soil shape calibration parameter [unitless]
perc <- df.param [4, 2] # percolation rate [mm/T]
uzl <- df.param [5, 2] # threshold parameter for quick flow [mm]
k0 <- df.param [6, 2] # recession coefficient for upper zone [1/T]
k1 <- df.param [7, 2] # recession coefficient for upper zone [1/T]
k2 <- df.param [8, 2] # recession coefficient for lower zone [1/T]
maxbas <- df.param [9, 2] # length of weighting  transformation function [T]

# Monthly potential evapotranspiration containers are loaded [mm/month]
evap.jan <- df.evap [1, 2]
evap.feb <- df.evap [2, 2]
evap.mar <- df.evap [3, 2]
evap.apr <- df.evap [4, 2]
evap.may <- df.evap [5, 2]
evap.jun <- df.evap [6, 2]
evap.jul <- df.evap [7, 2]
evap.aug <- df.evap [8, 2]
evap.set <- df.evap [9, 2]
evap.oct <- df.evap [10, 2]
evap.nov <- df.evap [11, 2]
evap.dic <- df.evap [12, 2]

# Watershed attributes are loaded
warea <- as.numeric(df.attri [1, 2]) # catchment area (km2)
fract_SMi <- as.numeric(df.attri [2, 2]) # initial SM fraction [unitless]
fract_SLZi <- as.numeric(df.attri [3, 2]) # initial SLZ fraction [unitless]
resolution <- as.numeric(df.attri [4, 2]) # model temporal resolution (1 = daily; 2 = hourly)
routing <- as.numeric(df.attri [5, 2])  # maxbas selection-control parameter

# MASSINf <- vector() # mass entering the loop as precipitation [mm/T]

# Numeric counters are defined
countermain <- length(df.ptq$QOBS)

# A sequence variable is created at df.ptq data.frame
df.ptq$SEQ <- seq(1, (length(df.ptq$QOBS)), by = 1)

# df.ptq$DATE factor class is converted to date class 
DATEtemp01 <- df.ptq$DATE
df.ptq$DATE <- as.Date(DATEtemp01, format = "%d/%m/%Y")

# lubridate Library functions are loaded
df.ptq$YEAR <- year(df.ptq$DATE) # year component of DATE
df.ptq$YEAR_CH <- as.character(year(df.ptq$DATE)) # year component of DATE as character
df.ptq$MONTH <- month(df.ptq$DATE, label = FALSE) # month component of DATE
df.ptq$MONTH_CH <- month(df.ptq$DATE, label = TRUE) # month component of DATE as character
df.ptq$WEEK <- week(df.ptq$DATE) # week component of DATE
df.ptq$DAY <- yday(df.ptq$DATE) # day component of a DATE
df.ptq$DAY_MONTH <- days_in_month(df.ptq$DATE) # number of days in the month of DATE

# A selection data.frame is created based on df.ptq data.frame
df.selection <- df.ptq[c("SEQ", "DATE", "YEAR", "YEAR_CH", "MONTH", "MONTH_CH", "DAY")]

# Model state variables are initialized
SMtemp <- fract_SMi * fc # temporal soil moisture [mm/T]
SLZtemp <- fract_SLZi * (perc / k2) # temporal storage in the lower zone [mm/T]
SUZtemp <- 0 # temporal storage in the upper zone [mm/T]

# "MONTH" data.frame variable is converted to vector
f_month <- df.ptq[, "MONTH"]

# ----------------------------------------------------------
# Precipitation routine
# ----------------------------------------------------------

# MASSINf is defined [mm/T] (entering the loop as precipitation/snow) 
PRECtemp <- df.ptq[, "PREC"]    
MASSINf <- df.ptq[, "PREC"]  
outQOBS <- df.ptq[, "QOBS"]
outPREC <- df.ptq[, "PREC"]
outTEMP <- df.ptq[, "TEMP"]

# A vector containing optimized PEST parameters is created
pars.OPTI <- c(fc=fc,
               lp=lp,
               beta=beta,
               perc=perc,
               uzl=uzl,
               k0=k0,
               k1=k1,
               k2=k2,
               maxbas=maxbas)

# A vector containing optimized PEST parameters is created
pars.OPTI.G <- c(fcG1=fc,
                 lpG1=lp,
                 betaG1=beta,
                 percG1=perc,
                 uzlG1=uzl,
                 k0G1=k0,
                 k1G1=k1,
                 k2G1=k2,
                 maxbasG1=maxbas,
                 
                 fcG2=fc,
                 lpG2=lp,
                 betaG2=beta,
                 percG2=perc,
                 uzlG2=uzl,
                 k0G2=k0,
                 k1G2=k1,
                 k2G2=k2,
                 maxbasG2=maxbas,
                 
                 fcG34=fc,
                 lpG34=lp,
                 betaG34=beta,
                 percG34=perc,
                 uzlG34=uzl,
                 k0G34=k0,
                 k1G34=k1,
                 k2G34=k2,
                 maxbasG34=maxbas,
                 
                 fcGX=fc,
                 lpGX=lp,
                 betaGX=beta,
                 percGX=perc,
                 uzlGX=uzl,
                 k0GX=k0,
                 k1GX=k1,
                 k2GX=k2,
                 maxbasGX=maxbas,
                 
                 fcG11=fc,
                 lpG11=lp,
                 betaG11=beta,
                 percG11=perc,
                 uzlG11=uzl,
                 k0G11=k0,
                 k1G11=k1,
                 k2G11=k2,
                 maxbasG11=maxbas,
                 
                 fcG12=fc,
                 lpG12=lp,
                 betaG12=beta,
                 percG12=perc,
                 uzlG12=uzl,
                 k0G12=k0,
                 k1G12=k1,
                 k2G12=k2,
                 maxbasG12=maxbas)

# A vector containing optimized PEST parameters is created
pars.LOW.G = c(fcG1 = 100,
               lpG1 = 0.10,
               betaG1 = 0.10,
               percG1 = 0.10,
               uzlG1 = 50,
               k0G1 = 0.001,
               k1G1 = 0.01,
               k2G1 = 0.01,
               maxbasG1 = 1,
               
               fcG2 = 100,
               lpG2 = 0.10,
               betaG2 = 0.10,
               percG2 = 0.10,
               uzlG2 = 50,
               k0G2 = 0.001,
               k1G2 = 0.01,
               k2G2 = 0.01,
               maxbasG2 = 1,
               
               fcG34 = 100,
               lpG34 = 0.10,
               betaG34 = 0.10,
               percG34 = 0.10,
               uzlG34 = 50,
               k0G34 = 0.001,
               k1G34 = 0.01,
               k2G34 = 0.01,
               maxbasG34 = 1,
               
               fcGX = 100,
               lpGX = 0.10,
               betaGX = 0.10,
               percGX = 0.10,
               uzlGX = 50,
               k0GX = 0.001,
               k1GX = 0.01,
               k2GX = 0.01,
               maxbasGX = 1,
               
               fcG11 = 100,
               lpG11 = 0.10,
               betaG11 = 0.10,
               percG11 = 0.10,
               uzlG11 = 50,
               k0G11 = 0.001,
               k1G11 = 0.01,
               k2G11 = 0.01,
               maxbasG11 = 1,
               
               fcG12 = 100,
               lpG12 = 0.10,
               betaG12 = 0.10,
               percG12 = 0.10,
               uzlG12 = 50,
               k0G12 = 0.001,
               k1G12 = 0.01,
               k2G12 = 0.01,
               maxbasG12 = 1)

# A vector containing optimized PEST parameters is created
pars.HIG.G = c(fcG1 = 800,
               lpG1 = 0.99,
               betaG1 = 4,
               percG1 = 20,
               uzlG1 = 120,
               k0G1 = 0.99,
               k1G1 = 0.99,
               k2G1 = 0.99,
               maxbasG1 = 5,
               
               fcG2 = 800,
               lpG2 = 0.99,
               betaG2 = 4,
               percG2 = 20,
               uzlG2 = 120,
               k0G2 = 0.99,
               k1G2 = 0.99,
               k2G2 = 0.99,
               maxbasG2 = 5,
               
               fcG34 = 800,
               lpG34 = 0.99,
               betaG34 = 4,
               percG34 = 20,
               uzlG34 = 120,
               k0G34 = 0.99,
               k1G34 = 0.99,
               k2G34 = 0.99,
               maxbasG34 = 5,
               
               fcGX = 800,
               lpGX = 0.99,
               betaGX = 4,
               percGX = 20,
               uzlGX = 120,
               k0GX = 0.99,
               k1GX = 0.99,
               k2GX = 0.99,
               maxbasGX = 5,
               
               fcG11 = 800,
               lpG11 = 0.99,
               betaG11 = 4,
               percG11 = 20,
               uzlG11 = 120,
               k0G11 = 0.99,
               k1G11 = 0.99,
               k2G11 = 0.99,
               maxbasG11 = 5,
               
               fcG12 = 800,
               lpG12 = 0.99,
               betaG12 = 4,
               percG12 = 20,
               uzlG12 = 120,
               k0G12 = 0.99,
               k1G12 = 0.99,
               k2G12 = 0.99,
               maxbasG12 = 5)

pars.INI <- NULL
pars.INI <- list()

for (i in 1:50) {
 set.seed(115 + i)
 pars.INI[i] <- list(c(fc=runif(1, 101, 799),
                       lp=runif(1, 0.11, 0.98),
                       beta=runif(1, 0.11, 3.9),
                       perc=runif(1, 0.11, 19),
                       uzl=runif(1, 51, 119),
                       k0=runif(1, 0.0012, 0.98),
                       k1=runif(1, 0.012, 0.98),
                       k2=runif(1, 0.012, 0.98),
                       maxbas=runif(1, 1.1, 4.9)))
}

# A vector containing low PEST parameters is created
pars.LOW <- c(fc=100, #800,
              lp=0.10, #0.25,
              beta=0.10,#2.576,
              perc=0.10,#13.815,
              uzl=50,#94.15,
              k0=0.001,#0.303963,
              k1=0.01,#0.1687,
              k2=0.01,#0.023536,
              maxbas=1.00)

# A vector containing high PEST parameters is created
pars.HIG <- c(fc=800, #800,
              lp=0.99, #0.25,
              beta=4.0,#2.576,
              perc=20.0,#13.815,
              uzl=120,#94.15,
              k0=0.99,#0.303963,
              k1=0.99,#0.1687,
              k2=0.99,#0.023536,
              maxbas=5.00)

# ============================================================================
# ============================================================================
# ============================================================================

# Start the clock!
ptm <- proc.time()

# f.HBV01 function is evaluated
f.HBV01(par =pars.OPTI)
#f.HBV04(parS =pars.OPTI.G)

# Stop the clock
proc.time() - ptm

# ============================================================================
# ============================================================================
# ============================================================================

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package nmkb {dfoptim}
# Nelder-Mead optimziation algorithm for derivative-free optimization
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.NMK <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.NMK <- foreach(i = 1:length(pars.INI),
                          .packages = c("dfoptim"),
                          .export=c("f.HBV02"),
                          .multicombine = TRUE) %dopar% {
                           # --------------------------------------------------
                           # BLOCK: tracer
                           .count <- 0
                           trace(identity, tracer=function() .count <<- .count +1, print = F)
                           # --------------------------------------------------
                           CNT <- nmkb(par = unlist(pars.INI[i]),
                                       fn = f.HBV02, 
                                       lower = pars.LOW,
                                       upper = pars.HIG,
                                       control = list(maxfeval = 20000))
                           # --------------------------------------------------
                           .count
                           untrace(identity)
                           # --------------------------------------------------
                           list(.count, CNT)
                          }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.NMK, file = "Y16_result.NMK")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package optim {stats}
# Function = L-BFGS-B
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.LBFGSB <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.LBFGSB <- foreach(i = 1:length(pars.INI),
                             #.packages = c("optim"),
                             .export=c("f.HBV01"),
                             .multicombine = TRUE) %dopar% {
                              # --------------------------------------------------
                              # BLOCK: tracer
                              .count <- 0
                              trace(identity, tracer=function() .count <<- .count +1, print = F)
                              # --------------------------------------------------
                              CNT <- optim(par = unlist(pars.INI[i]),
                                           fn = f.HBV01, 
                                           method ="L-BFGS-B",
                                           lower = pars.LOW,
                                           upper = pars.HIG,
                                           control = list(maxit = 500))
                              # --------------------------------------------------
                              .count
                              untrace(identity)
                              # --------------------------------------------------
                              list(.count, CNT)
                             }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.LBFGSB, file = "Y16_result.LBFGSB")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package hjkb {dfoptim}
# Hooke-Jeeves derivative-free minimization algorithm
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.HJKB <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.HJKB <- foreach(i = 1:length(pars.INI),
                           .packages = c("dfoptim"),
                           .export=c("f.HBV02"),
                           .multicombine = TRUE) %dopar% {
                            # --------------------------------------------------
                            # BLOCK: tracer
                            .count <- 0
                            trace(identity, tracer=function() .count <<- .count +1, print = F)
                            # --------------------------------------------------
                            CNT <- hjkb(par = unlist(pars.INI[i]),
                                        fn = f.HBV02, 
                                        lower = pars.LOW,
                                        upper = pars.HIG,
                                        control = list(maxfeval = 20000))
                            # --------------------------------------------------
                            .count
                            untrace(identity)
                            # --------------------------------------------------
                            list(.count, CNT)
                           }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.HJKB, file = "Y16_result.HJKB")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package Rvmmin {Rvmmin}
# Function =  Variable metric nonlinear function minimization, driver.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.RVMMIN <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.RVMMIN <- foreach(i = 1:length(pars.INI),
                             .packages = c("Rvmmin","optextras","numDeriv"),
                             .export=c("f.HBV02"),
                             .multicombine = TRUE) %dopar% {
                              # --------------------------------------------------
                              # BLOCK: tracer
                              .count <- 0
                              trace(identity, tracer=function() .count <<- .count +1, print = F)
                              # --------------------------------------------------
                              CNT <- Rvmmin(par = unlist(pars.INI[i]),
                                            fn=f.HBV02, 
                                            lower = pars.LOW,
                                            upper = pars.HIG,
                                            bdmsk = c(1,1,1,1,1,1,1,1,1),
                                            control = list(maxit = 500))
                              # --------------------------------------------------
                              .count
                              untrace(identity)
                              # --------------------------------------------------
                              list(.count, CNT)
                             }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.RVMMIN, file = "Y16_result.RVMMIN")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package bobyqa {minqa}
# Function = bobyqa implementation of Powell
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.BOBYQA <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.BOBYQA <- foreach(i = 1:length(pars.INI),
                             .packages = c("minqa"),
                             .export=c("f.HBV02"),
                             .multicombine = TRUE) %dopar% {
                              # --------------------------------------------------
                              # BLOCK: tracer
                              .count <- 0
                              trace(identity, tracer=function() .count <<- .count +1, print = F)
                              # --------------------------------------------------
                              CNT <- bobyqa(par = unlist(pars.INI[i]),
                                            fn=f.HBV02, 
                                            lower = pars.LOW,
                                            upper = pars.HIG,
                                            control = list(maxfun=20000))
                              # --------------------------------------------------
                              .count
                              untrace(identity)
                              # --------------------------------------------------
                              list(.count, CNT)
                             }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.BOBYQA, file = "Y16_result.BOBYQA")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package spg {BB}
# Function = Large-Scale Optimization
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.SPG <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.SPG <- foreach(i = 1:length(pars.INI),
                          .packages = c("BB"),
                          .export=c("f.HBV02"),
                          .multicombine = TRUE) %dopar% {
                           # --------------------------------------------------
                           # BLOCK: tracer
                           .count <- 0
                           trace(identity, tracer=function() .count <<- .count +1, print = F)
                           # --------------------------------------------------
                           CNT <- spg(par = unlist(pars.INI[i]),
                                      fn=f.HBV02, 
                                      lower = pars.LOW,
                                      upper = pars.HIG,
                                      method =3,
                                      control = list(maxfeval=20000))
                           # --------------------------------------------------
                           .count
                           untrace(identity)
                           # --------------------------------------------------
                           list(.count, CNT)
                          }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.SPG, file = "Y16_result.SPG")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package nlminb {stats}
# Function = nlm {stats}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.NLMINB <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.NLMINB <- foreach(i = 1:length(pars.INI),
                             #.packages = c("optimx"),
                             .export=c("f.HBV02"),
                             .multicombine = TRUE) %dopar% {
                              # --------------------------------------------------
                              # BLOCK: tracer
                              .count <- 0
                              trace(identity, tracer=function() .count <<- .count +1, print = F)
                              # --------------------------------------------------
                              CNT <- nlminb(start = unlist(pars.INI[i]),
                                            objective=f.HBV02, 
                                            lower = pars.LOW,
                                            upper = pars.HIG,
                                            control = list(eval.max=20000))
                              # --------------------------------------------------
                              .count
                              untrace(identity)
                              # --------------------------------------------------
                              list(.count, CNT)
                             }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.NLMINB, file = "Y16_result.NLMINB")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using GenSA {GenSA} Generalized Simulated Annealing 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.expSA <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.expSA <- foreach(i = 1:length(pars.INI),
                            .packages = c("GenSA"),
                            .export=c("f.HBV02"),
                            .multicombine = TRUE) %dopar% {
                             # --------------------------------------------------
                             # BLOCK: tracer
                             .count <- 0
                             trace(identity, tracer=function() .count <<- .count +1, print = F)
                             # --------------------------------------------------
                             CNT <- GenSA(par = unlist(pars.INI[i]),
                                          fn = f.HBV02, 
                                          lower = pars.LOW,
                                          upper = pars.HIG,
                                          control = list(maxit = 200,
                                                         nb.stop.improvement = 100,
                                                         max.call = 20000))
                             # --------------------------------------------------
                             .count
                             untrace(identity)
                             # --------------------------------------------------
                             list(.count, CNT)
                            }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.expSA, file = "Y16_result.expSA")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using DEoptim {DEoptim} Differential Evolution Optimization
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.expDE <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.expDE <- foreach(i = 1:length(pars.INI),
                            .packages = c("DEoptim"),
                            .export=c("f.HBV02"),
                            .multicombine = TRUE) %dopar% {
                             # --------------------------------------------------
                             # BLOCK: tracer
                             .count <- 0
                             trace(identity, tracer=function() .count <<- .count +1, print = F)
                             # --------------------------------------------------
                             CNT <- DEoptim(fn = f.HBV02, 
                                            lower = pars.LOW,
                                            upper = pars.HIG,
                                            control = DEoptim.control(itermax = 200))
                             # --------------------------------------------------
                             .count
                             untrace(identity)
                             # --------------------------------------------------
                             list(.count, CNT)
                            }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.expDE, file = "Y16_result.expDE")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using ga{GA} Genetic Algorithms
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.expGA <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.expGA <- foreach(i = 1:length(pars.INI),
                            .packages = c("GA"),
                            .export=c("f.HBV03"),
                            .multicombine = TRUE) %dopar% {
                             # --------------------------------------------------
                             # BLOCK: tracer
                             .count <- 0
                             trace(identity, tracer=function() .count <<- .count +1, print = F)
                             # --------------------------------------------------
                             CNT <- ga(type = "real-valued",
                                       fitness = f.HBV03,
                                       min = pars.LOW,
                                       max = pars.HIG,
                                       maxiter = 200,
                                       keepBest = TRUE,
                                       run = 100)
                             # --------------------------------------------------
                             .count
                             untrace(identity)
                             # --------------------------------------------------
                             list(.count, CNT)
                            }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.expGA, file = "Y16_result.expGA")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using nls.lm {minpack.lm}. Levenberg-Marquardt Algorithm
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.expMK <- NULL

# Start the clock!
ptm <- proc.time()

nls.lm.control(ftol = 1.490116e-6,
               ptol = 1.490116e-6, gtol = 0, diag = list(), epsfcn = 0,
               factor = 100, maxfev = 20000, maxiter = 500, nprint = 0)

# Optimization function is executed and an object class list is created
Y16_result.expMK.PCA01 <- foreach(i = 1:1,
                                  .packages = c("minpack.lm"),
                                  .export=c("f.HBV04"),
                                  .multicombine = TRUE) %dopar% {
                                   # --------------------------------------------------
                                   # BLOCK: tracer
                                   .count <- 0
                                   trace(identity, tracer=function() .count <<- .count +1, print = F)
                                   # --------------------------------------------------
                                   CNT <- nls.lm(par = pars.OPTI.G, # unlist(pars.INI[i]),
                                                 lower = pars.LOW.G,
                                                 upper = pars.HIG.G,
                                                 fn = f.HBV04)
                                   # --------------------------------------------------
                                   .count
                                   untrace(identity)
                                   # --------------------------------------------------
                                   list(.count, CNT)
                                  }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.expMK.PCA01, file = "Y16_result.expMK.PCA01")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package  SCEoptim {hydromad}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.optim.SCE <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.optim.SCE <- foreach(i = 1:length(pars.INI),
                                .packages = c("hydromad", "zoo", "lattice", "latticeExtra", "RColorBrewer", "polynom", "reshape"),
                                .export=c("f.HBV01"),
                                .multicombine = TRUE) %dopar% {
                                 # --------------------------------------------------
                                 # BLOCK: tracer
                                 .count <- 0
                                 trace(identity, tracer=function() .count <<- .count +1, print = F)
                                 # --------------------------------------------------
                                 CNT <- SCEoptim(FUN = f.HBV01,
                                                 par = unlist(pars.INI[i]),
                                                 lower = pars.LOW,
                                                 upper = pars.HIG,
                                                 control = list(maxit = 200,
                                                                maxeval = 20000,
                                                                elitism = 5,
                                                                initsample = "random",
                                                                maxtime = 1800))
                                 # --------------------------------------------------
                                 .count
                                 untrace(identity)
                                 # --------------------------------------------------
                                 list(.count, CNT)
                                }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.optim.SCE, file = "Y16_result.optim.SCE")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package hydroPSO {hydroPSO}
# Particle Swarm Optimisation
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.PSO <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.PSO <- foreach(i = 1:length(pars.INI),
                          .packages = c("hydroPSO"),
                          .export=c("f.HBV02"),
                          .multicombine = TRUE) %dopar% {
                           # --------------------------------------------------
                           # BLOCK: tracer
                           .count <- 0
                           trace(identity, tracer=function() .count <<- .count +1, print = F)
                           # --------------------------------------------------
                           CNT <- hydroPSO(par = unlist(pars.INI[i]),
                                           fn = f.HBV02,
                                           lower = pars.LOW,
                                           upper = pars.HIG,
                                           method=c("ipso"),
                                           control=list(maxit = 400,
                                                        maxfn = 20000,
                                                        topology = "gbest",
                                                        npart = 50))
                           # --------------------------------------------------
                           .count
                           untrace(identity)
                           # --------------------------------------------------
                           list(.count, CNT)
                          }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.PSO, file = "Y16_result.PSO")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package DIRECT {nloptr}
# DIviding RECTangles Algorithm for Global Optimization
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.DIRECT <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.DIRECT <- foreach(i = 1:length(pars.INI),
                             .packages = c("nloptr"),
                             .export=c("f.HBV02"),
                             .multicombine = TRUE) %dopar% {
                              # --------------------------------------------------
                              # BLOCK: tracer
                              .count <- 0
                              trace(identity, tracer=function() .count <<- .count +1, print = F)
                              # --------------------------------------------------
                              CNT <- direct(fn = f.HBV02,
                                            lower = pars.LOW,
                                            upper = pars.HIG,
                                            control=list(maxeval = 20000))
                              # --------------------------------------------------
                              .count
                              untrace(identity)
                              # --------------------------------------------------
                              list(.count, CNT)
                             }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.DIRECT, file = "Y16_result.DIRECT")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package csr2lm {nloptr}
# Controlled Random Search (CRS)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.CRS <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.CRS <- foreach(i = 1:length(pars.INI),
                          .packages = c("nloptr"),
                          .export=c("f.HBV02"),
                          .multicombine = TRUE) %dopar% {
                           # --------------------------------------------------
                           # BLOCK: tracer
                           .count <- 0
                           trace(identity, tracer=function() .count <<- .count +1, print = F)
                           # --------------------------------------------------
                           CNT <- crs2lm(x0 = unlist(pars.INI[i]),
                                         fn = f.HBV02,
                                         lower = pars.LOW,
                                         upper = pars.HIG,
                                         maxeval = 20000)
                           # --------------------------------------------------
                           .count
                           untrace(identity)
                           # --------------------------------------------------
                           list(.count, CNT)
                          }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.CRS, file = "Y16_result.CRS")

# OUTPUT-OK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package auglag {alabama}
# Augmented Lagrangian Minimization Algorithm
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Y16_result.ALAB <- NULL

# Start the clock!
ptm <- proc.time()

# Optimization function is executed and an object class list is created
Y16_result.ALAB <- foreach(i = 1:length(pars.INI),
                           .packages = c("nloptr"),
                           .export=c("f.HBV02"),
                           .multicombine = TRUE) %dopar% {
                            # --------------------------------------------------
                            # BLOCK: tracer
                            .count <- 0
                            trace(identity, tracer=function() .count <<- .count +1, print = F)
                            # --------------------------------------------------
                            CNT <- auglag(x0 = unlist(pars.INI[i]),
                                          fn = f.HBV02,
                                          lower = pars.LOW,
                                          upper = pars.HIG,
                                          localsolver = c("MMA"),
                                          control = list(maxeval = 20000))
                            # --------------------------------------------------
                            .count
                            untrace(identity)
                            # --------------------------------------------------
                            list(.count, CNT)
                           }

# Stop the clock
proc.time() - ptm

# Binary results are saved
save(Y16_result.ALAB, file = "Y16_result.ALAB")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Package stogo {nloptr}
# Stochastic Global Optimization
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

result.StoGo <- NULL

# Start the clock!
ptm <- proc.time()

# --------------------------------------------------
# BLOCK: tracer
.count <- 0
trace(f.HBV02, tracer=function() .count <<- .count +1)
# --------------------------------------------------

result.StoGo <- stogo(x0 = pars.OPTI,
                      fn = f.HBV02,
                      lower = pars.LOW,
                      upper = pars.HIG,
                      maxeval = 100,
                      randomized = TRUE,
                      nl.info = TRUE)


# Stop the clock
proc.time() - ptm

# --------------------------------------------------
print(.count)
untrace(f.HBV02)
# --------------------------------------------------

result.StoGo

# OUTPUT-METHODS OK;
# ============================================================================

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPTIMIZATION using Inverse modeling of ODEs with the FME package
# Inverse Modelling, Sensitivity and Monte Carlo Analysis in R Using Package FME
# Karline Soetaert -----------
# Thomas Petzoldt ------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ============================================================================

# A vector containing optimized csr2lm parameters is created, NS=0.8406343
pars.OPTI <- c(fc=800,
               lp=0.1575015,
               beta=0.194125,
               perc=13.3988,
               uzl=91.46412,
               k0=0.282918,
               k1=0.1661498,
               k2=0.02042238,
               maxbas=1.053002)

# A vector containing optimized (DEoptim !!!) R parameters is created
pars.BEST <- c(fc=798.0893,
               lp=0.3585103,
               beta=0.2928497,
               perc=9.728414,
               uzl=82.11291,
               k0=0.2642246,
               k1=0.1221602,
               k2=0.01307838,
               maxbas=1.022056)

# An observations data.frame is created based on outQOBS vector from 
# function f.HBV01. Therefore, input *.TXT must also be loaded
df.base2 <- data.frame(xx=(seq(1, (length(outQOBS)), by = 1)),yy=outQOBS)

# IMPORTANT !!!!! Total length of observations record is calculated
length(outQOBS)

# //////////////////////////////////////////////////////////////////////////
# FUNCTION: f.CostV Cost function as requested by FME
# //////////////////////////////////////////////////////////////////////////

# It dependes on named-vector "pars.OPTI" coming from ANY of the
# 16 optimization methods, either local or global
# modCost {FME} calculates the discrepancy of a model solution with observations
f.CostV <- function(par = pars.OPTI) {
 out01 <- f.HBV01(par)
 cost <- modCost(model = out01,obs = df.base2, x = "xx")
 return(cost)                
}

# Cost function is executed using the optimum PEST parameters
out.CostV <- f.CostV(par = pars.OPTI)

# Model cost is evaluated
out.CostV$model

# Variables data.frames is requested
View(out.CostV$var)

# Residuals data.frame is requested
View(out.CostV$residuals)

# A simple scatter-plot is requested
plot(out.CostV$residuals$obs, out.CostV$residuals$mod )

# A simple histogram is geneted on the residuals
hist(out.CostV$residuals$res)

# A simple boxplot is geneted on the residuals
boxplot(out.CostV$residuals$res)

# A Shapiro Test is performed on the residuals
shapiro.test(sample(out.CostV$residuals$res, 4999))

# A 1/3 Shapiro Test is performed on the residuals
shapiro.test((sample(out.CostV$residuals$res, 4999))^(1/3))

# A simple scatter-plot is requested
plot(out.CostV$residuals$x, out.CostV$residuals$res)


# ============================================================================
# 3. Local sensitivity analysis (Soetaert & Petzoldt)
# ============================================================================

# Binary format relevant data.frames are loaded 
load("df.rbind.OF")
load("df.rbind.PAR")
df.rbind.PAR$NS <- df.rbind.OF$NS

# Local Sensitivity Analysis, sensFun{FME} is applied
# Function sensFun estimates the sensitivity of the model 
# output to the parameter values in a set of so-called 
# sensitivity functions (Brun et al. 2001; Soetaert and Herman 2009)
Sfun <- sensFun(f.CostV, pars.OPTI) # pars.INI, pars.LOW, pars.HIG

# In FME, normalised, dimensionless sensitivities of model output 
# to parameters are in a sensitivity matrix. The higher the absolute
# sensitivity value, the more important the parameter
View(Sfun)

# Thus the magnitudes of the sensitivity summary values can be used
# to rank the importance of parameters on the output variables
summary(Sfun)

# As it makes no sense to finetune parameters that have little effect, this ranking serves to
# choose candidate parameters for model fitting
df.SFUN.desc <- as.data.frame(summary(Sfun))

# -------------------------
# value	   L2    ranking
# -------------------------
# perc	  0.00444	   1   Same order of magnitude
# k2	    0.00345	   2   Same order of magnitude
# maxbas	0.00340	   3   Same order of magnitude
# k1	    0.00211	   4   Same order of magnitude
# uzl	    0.00139	   5   Same order of magnitude
# k0	    0.00062	   6
# fc	    0.00048	   7
# beta	  0.00026	   8
# lp	    0.00000	   9

# The sensitivities of the modelled values to the parameter values 
# change in time, thus it makes sense to visualise 
# the sensitivity functions as they fluctuate
plot(Sfun)

# A data.frame subset for one year is created (1995)
Sfun.subset <- Sfun[367:731,]
Sfun.subset <- Sfun.subset[,-c(1:2)]

# A dummy precipitation variable is included
#Sfun.subset$RAIN01 <-outPREC[367:731]
#Sfun.subset$RAIN02 <-outPREC[367:731]
#Sfun.subset$RAIN03 <-outPREC[367:731]

# data.frame is melted
Sfun.subset.melt <- melt(Sfun.subset)

# A id repetition series is created
Sfun.subset.melt$id <- (rep(1:365,by=9))

# A ggplot object is created
ggplot() +
 geom_hline(data=Sfun.subset.melt,yintercept = 0.00,colour = '#ff0000',size = 0.75,linetype = 2) +
 geom_line(aes(x = id,y = value),data=Sfun.subset.melt,size = 0.95) +
 facet_wrap(facets = ~variable) +#, scales = 'free_y') +
 scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
 scale_y_continuous(breaks = scales::pretty_breaks(min.n = 5.0)) +
 xlab("Time (days)") + 
 ylab("Parameter dimensionless sensitivity") +
 theme_bw(base_size = 14.0)

# Correlations are calculated for the entire dataset
df.CORR <- Sfun

# Precipitation is added to data.frame
df.CORR$RAIN <- outPREC

# Correlation between precipitation and parameters is calculated
cor(df.CORR$RAIN, df.CORR$perc)    # Rank 1
cor(df.CORR$RAIN, df.CORR$k2)      # Rank 2
cor(df.CORR$RAIN, df.CORR$maxbas)  # Rank 3
cor(df.CORR$RAIN, df.CORR$k1)      # Rank 4

# Pairwise relationships are visualised with a pairs plot
pairs(Sfun, col = c("blue"))

# ============================================================================
# 4. Multivariate parameter identifiability (Soetaert & Petzoldt)
# ============================================================================

# Function collin extends the analysis to all possible parameter combinations,
# by estimating the approximate linear dependence ("collinearity") of parameter sets

# The collinearity for all parameter combinations is estimated by function collin,
# taking the previously estimated sensitivity functions as argument
df.collin <- collin(Sfun)

# A parameter set is said to be identifiable, if all parameters within the set 
# can be uniquely estimated based on (perfect) measurements. Parameters that 
# have large collinearity will NOT be identifiable.

# INF values are removed
df.collin[mapply(is.infinite, df.collin)] <- 20

# Plot  shows how the collinearity index increases as more and more parameters are included in the set
# however, collinearity values are relatively low and far below the 10-15 range 
# suggested by Brun et al. 2001.
# Therefore, the HBV parameter set is algebraically identifiable.
plot(df.collin)

# Collinearity ggplot
ggplot() +
 geom_point(aes(x = N,y = collinearity),data=df.collin,shape = 8,size = 3.5,position = position_jitter(width = 0.1)) +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 8.0,min.n = 8.0),limits = c(0.8,2)) +
 scale_x_continuous(breaks = c(c(2,3,4,5,6,7,8,9))) +
 theme(axis.line.x = element_line(),axis.ticks.x = element_line(),panel.grid.minor = element_blank()) +
 xlab("Number of parameter combinations") + 
 ylab("Collinearity Index") +
 theme_bw(base_size = 14.0)

# ============================================================================
# 6. MCMC (Soetaert & Petzoldt)
# ============================================================================

#MCMC <- modMCMC(f = f.CostV, p = pars.OPTI, niter = 25000)
#save(MCMC, file = "MCMC")
load("MCMC")

summary(MCMC)

plot(MCMC, Full = TRUE)

pairs(MCMC, nsample = 4999)

df.mcmc <- as.data.frame(MCMC$pars)

# Descriptive statistics are requested and rounded to 4 decimals
df.mcmc.desc <- round((as.data.frame(stat.desc(df.mcmc))),4)

# ============================================================================
# 7. Model prediction
# ============================================================================

sR <- sensRange(func = f.HBV01, parms = pars.BEST, parInput = MCMC$par)

summary(sR)

plot(summary(sR), xlab = "time")

# ============================================================================
