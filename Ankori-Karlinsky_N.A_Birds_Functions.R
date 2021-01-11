##################################################################
################### FUNCTIONS FOR BIRD PROJECT ###################
##################################################################

# Roi Ankori-Karlinsky, Columbia University Spring 2020

# This script: contains functions for analyses and summary table generation for figures
# 1. First functions for getting bird diversity metrics for BBA/Adj when n=BBS
# 2. Then functions for environmental variables for BBA/Adj when n=BBS
#     a) Gamma
#     b) Alpha
#     c) Beta
# 3. This can be sourced into other scripts for quick loading and use of these functions

###################################################################

# Loading packages
library(tidyverse) #tidyverse for data cleaning and manipulation
library(mosaic) # mosaic for easy calculation of summary statistics
library(vegan) # vegan pacakge for calculating beta diversity
library(iNEXT) # iNEXT package for calculating regional richness
library(data.table) # data.table package for data manipulation
library(DescTools) # DescTools package for easy calculation of CIs



###################################################################
################# 1. BIRD DIVERSITY FUNCTIONS #####################
###################################################################


## Functions!!

# BirdRand function - Regional bird richness for BBA/Adj when n=BBS
# inputs are datasets (BBA and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the observed regional species richness for the subsetted data using iNEXT package
# and puts results in a vector
# After doing that it calculates mean and 95% CIs for all 1000 runs and outputs those values

BirdRand <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  BirdDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  BBA <- BBA[14:291] # Keeping only species incidence data
  BBA <- t(BBA) #Flipping table so iNEXT could use it
  
  for(n in 1:niter){ # creating forloop
    
    x=BBA[,sample(ncol(BBA),nrow(BBS),replace=T)] # sampling from BBA data n=BBS rows
    bba=DataInfo(x,datatype = "incidence_raw") # Calculating observed richness
    gam=bba$S.obs # Saving observed regional richness
    BirdDude[n,]<-gam # storing gamma species result into a vector
  }
  BirdDude <- arrange(BirdDude,nrow) #arranging by small to large
  me = mean(BirdDude$nrow) # mean
  low = BirdDude[25,] #low CI
  up = BirdDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#***********************

# BirdLoc function - Local bird richness for BBA/Adj when n=BBS
# Inputs are datasets (BBA and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the observed mean local species richness for the subsetted data
# and puts results in a vector
# After doing that it calculates mean and 95% CIs for all 1000 runs and outputs those values

BirdLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  BirdDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Atlas2_Total) #mean local richness
    BirdDude[n,]<-alph # storing local richness result into a vector
  }
  BirdDude <- arrange(BirdDude,nrow) #arranging by small to large
  me = mean(BirdDude$nrow) # mean
  low = BirdDude[25,] #low CI
  up = BirdDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#***********************

# BirdRandProp function - Vector of difference in gamma richness for when n=BBS
# Inputs are datasets (BBA and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the observed species richness for the subsetted data using iNEXT package
# and puts results in a vector 

BirdRandProp <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  BirdDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  BBA <- BBA[14:291] # Keeping only species incidence data
  BBA <- t(BBA) #Flipping table so iNEXT could use it
  
  for(n in 1:niter){ # creating forloop
    
    x=BBA[,sample(ncol(BBA),nrow(BBS),replace=T)] # sampling from BBA data n=BBS rows
    bba=DataInfo(x,datatype = "incidence_raw") # Calculating observed richness
    gam=bba$S.obs # Saving observed regional richness
    BirdDude[n,]<-gam # storing gamma species result into a vector
  }
  BirdDude
}



#***********************


# JacRand function - Compositional turnover for BBA/Adj when n=BBS
# Inputs are datasets (BBA/Adj and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the Jaccard dissimilarity Index for the whole dataset and puts the results
# into a vector and returns means and 95% CIs

JacRand <- function(BBA,BBS){
  
  #Parameters
  
  repeats=1000 # number of iterations or repeats in loop (1000)
  JacDude <- data.frame(nrow=repeats) # creating vector with row number=iterations to store results in from for loop, which can then be averaged to give us beta score
  
  #For loop
  
  for(n in 1:repeats){ # creating forloop
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    x=x[14:291] # removing the 'orig.id' column created from the sample function so we calculate only the bird incidences
    Baccardi=vegdist(x, method="jaccard",binary=T) # Calculating the Jaccard index for the region in the subsetted data
    JackyBoy=mean(Baccardi) # getting mean jaccard index for that subsetted data
    JacDude[n,] <- JackyBoy # Storing that Jaccard score in a dataframe
  }
  
  #Return value
  JacDude <- arrange(JacDude,nrow) #arranging by small to large
  me = mean(JacDude$nrow) # mean
  low = JacDude[25,] #low CI
  up = JacDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make it easier to integrate into summary tables
  all # output
}


# *****************************************************************


###################################################################
############# 2. ENVIRONMENTAL VARIABLES FUNCTIONS ################
###################################################################


# a) GAMMA RANGES

# ******************************************************************

# EleRand function - inputs are datasets (BBA and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the elevation range for the whole dataset and puts the results
# into a vector. Then returns mean and CIs
EleRand <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  EleDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    ran=max(x$Max_Elev)-min(x$Min_Elev) # calculating elevation range for data
    EleDude[n,]<-ran # storing elevation range result into a vector
  }
  
  EleDude <- arrange(EleDude,nrow) #arranging by small to large
  me = mean(EleDude$nrow) # mean
  low = EleDude[25,] #low CI
  up = EleDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}



#**********************************************

# TempRand function
TempRand <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  TempDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    ran=max(x$Temp)-min(x$Temp) # calculating temperature range for data
    TempDude[n,]<-ran # storing temperature range result into a vector
  }
  
  TempDude <- arrange(TempDude,nrow) #arranging by small to large
  me = mean(TempDude$nrow) # mean
  low = TempDude[25,] #low CI
  up = TempDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
}



#**********************************************

# PrecRand function

PrecRand <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  PrecDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    ran=max(x$Precipitation)-min(x$Precipitation) # calculating precipitation range for data
    PrecDude[n,]<-ran # storing precipitation range result into a vector
  }
  
  PrecDude <- arrange(PrecDude,nrow) #arranging by small to large
  me = mean(PrecDude$nrow) # mean
  low = PrecDude[25,] #low CI
  up = PrecDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}



#************************************

#HetRand function (just gets most land-class richness present)
HetRand <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  HetDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    ran=max(x$Het_Count2006) # calculating land-class presence
    HetDude[n,]<-ran # storing land-class range result into a vector
  }
  
  HetDude <- arrange(HetDude,nrow) #arranging by small to large
  me = mean(HetDude$nrow) # mean
  low = HetDude[25,] #low CI
  up = HetDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


# ******************************************************************

# b) ALPHA VALUES

#Local Elevation mean for BBA when n=BBS

ElevLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  ElevDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Elev_Mean) #mean elevation in sampling unit
    ElevDude[n,]<-alph # storing result into a vector
  }
  
  ElevDude <- arrange(ElevDude,nrow) #arranging by small to large
  me = mean(ElevDude$nrow) # mean
  low = ElevDude[25,] #low CI
  up = ElevDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#*****************************************


#Local Precipitation mean for BBA when n=BBS

PrecLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  PrecDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Precipitation) #mean precipitation in sampling unit
    PrecDude [n,]<-alph # storing result into a vector
  }
  
  PrecDude  <- arrange(PrecDude ,nrow) #arranging by small to large
  me = mean(PrecDude$nrow) # mean
  low = PrecDude[25,] #low CI
  up = PrecDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#*****************************************

#Local Temp mean for BBA when n=BBS

TempLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  TempDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Temp) #mean temp in sampling unit
    TempDude[n,]<-alph # storing result into a vector
  }
  
  TempDude  <- arrange(TempDude ,nrow) #arranging by small to large
  me = mean(TempDude$nrow) # mean
  low = TempDude[25,] #low CI
  up = TempDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}

#*****************************************


#Local Elevation Range mean for BBA when n=BBS

ElevRangeLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  ElevDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Elev_Range) #mean elev range in sampling unit
    ElevDude[n,]<-alph # storing result into a vector
  }
  
  ElevDude  <- arrange(ElevDude,nrow) #arranging by small to large
  me = mean(ElevDude$nrow) # mean
  low = ElevDude[25,] #low CI
  up = ElevDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}

#*******************************************

#Local Land-Class Richness mean for BBA when n=BBS

LandRichLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  LandDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Het_Count2006) #mean land-class richness in sampling unit
    LandDude[n,]<-alph # storing result into a vector
  }
  
  LandDude  <- arrange(LandDude,nrow) #arranging by small to large
  me = mean(LandDude$nrow) # mean
  low = LandDude[25,] #low CI
  up = LandDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}

#*******************************************

#Local LandClass Diversity mean for BBA when n=BBS

LandDivLoc <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  LandDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    alph=mean(x$Shannon2006) #mean shannon index in sampling unit
    LandDude[n,]<-alph # storing result into a vector
  }
  
  LandDude  <- arrange(LandDude,nrow) #arranging by small to large
  me = mean(LandDude$nrow) # mean
  low = LandDude[25,] #low CI
  up = LandDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


# ******************************************************************

# c) BETA VALUES

##Beta elevation mean
ElevB <- function(BBA) {
  BBA <- BBA %>% filter(!is.na(Elev_Mean),Elev_Mean!=0) # filtering out 0 values 
  bet=vegdist(BBA$Elev_Mean,method="bray",na.rm=T) # getting distance metric
  bet=1-bet # turning into similarity rather than dissimilarity
  bet <- as.vector(bet) #turning into vector
  MeanCI(bet,conf.level = 0.95,na.rm=T)# return mean and CIs 
}


#Now for BBA when n=BBS

ElevBeta <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  ElevDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% filter(Elev_Mean!=0) # filtering out 0 values 
    x=sample(BBA,nrow(BBS),replace=T) # sampling from BBA data n=BBS rows
    bet=vegdist(x$Elev_Mean,method="bray",na.rm=T) # getting distance metric
    bet=1-bet #turning into similarity
    betty=mean(bet) # calculating mean
    ElevDude[n,]<-betty # storing result into a vector
  }
  
  ElevDude  <- arrange(ElevDude,nrow) #arranging by small to large
  me = mean(ElevDude$nrow) # mean
  low = ElevDude[25,] #low CI
  up = ElevDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}




#*****************************************

#Precipitation

PrecB <- function(BBA) {
  BBA <- BBA %>% filter(Precipitation!=0) # filtering out 0 values 
  bet=vegdist(BBA$Precipitation,method="bray",na.rm=T) # getting distance metric
  bet=1-bet #turning into similarity
  bet <- as.vector(bet) #turning into vector
  
  MeanCI(bet,conf.level = 0.95,na.rm=T)# return mean and CIs 
}


#Beta Precipitation mean for BBA when n=BBS

PrecBeta <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  PrecDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% filter(Precipitation!=0) # filtering out 0 values 
    x=sample(BBA,nrow(BBS),replace = T) # sampling from BBA data n=BBS rows
    bet=vegdist(x$Precipitation,method="bray",na.rm=T) # getting distance metric
    bet=1-bet # turning into similarity
    betty=mean(bet) # calculating mean
    PrecDude[n,]<-betty # storing result into a vector
  }
  
  PrecDude  <- arrange(PrecDude,nrow) #arranging by small to large
  me = mean(PrecDude$nrow) # mean
  low = PrecDude[25,] #low CI
  up = PrecDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#*****************************************

TempB <- function(BBA) {
  BBA <- BBA %>% filter(Temp!=0) # filtering out 0 values 
  bet=vegdist(BBA$Temp,method="bray",na.rm=T) # getting distance metric
  bet=1-bet #turning into similarity
  bet <- as.vector(bet) #turning into vector
  
  MeanCI(bet,conf.level = 0.95,na.rm=T)# return mean and CIs 
}


#Beta Temp mean for BBA when n=BBS

TempBeta <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  TempDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% filter(Temp!=0) # filtering out 0 values 
    x=sample(BBA,nrow(BBS),replace = T) # sampling from BBA data n=BBS rows
    bet=vegdist(x$Temp,method="bray",na.rm=T) # getting distance metric
    bet=1-bet # turning into similarity
    betty=mean(bet) # calculating mean
    TempDude[n,]<-betty # storing result into a vector
  }
  
  TempDude <- arrange(TempDude,nrow) #arranging by small to large
  me = mean(TempDude$nrow) # mean
  low = TempDude[25,] #low CI
  up = TempDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}


#*******************************************

#LandClass

LandDivB <- function(BBA) {
  bet=vegdist(BBA[318:332],method="bray",na.rm=T) # getting distance metric
  bet=1-bet # turning into similarity
  bet <- as.vector(bet) #turning into vector
  
  MeanCI(bet,conf.level = 0.95,na.rm=T)# return mean and CIs 
}


#Beta LandClass Diversity mean for BBA when n=BBS

LandDivBeta <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  LandDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    x=sample(BBA,nrow(BBS),replace = T) # sampling from BBA data n=BBS rows
    bet=vegdist(x[318:332],method="bray",na.rm=T) # getting distance metric
    bet=1-bet #turning into similarity
    bet <- as.vector(bet) # turning to vector
    betty=mean(bet) # calculating mean
    LandDude[n,]<-betty # storing result into a vector
  }
  
  LandDude <- arrange(LandDude,nrow) #arranging by small to large
  me = mean(LandDude$nrow,na.rm=T) # mean
  low = LandDude[25,] #low CI
  up = LandDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}

#*******************************************

#Multivariate Gower
GowerB <- function(BBA) {
  bet=vegdist(BBA[c(293,299,301,335)],method="gower",na.rm=T) # getting distance metric
  bet <- as.vector(bet) #turning into vector
  
  MeanCI(bet,conf.level = 0.95,na.rm=T)# return mean and CIs 
}


#Beta multivariate Gower distance mean for BBA when n=BBS

GowerBeta <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  GowerDude <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    
    x=sample(BBA,nrow(BBS),replace = T) # sampling from BBA data n=BBS rows
    bet=vegdist(x[c(293,299,301,335)],method="gower",na.rm=T) # getting distance metric
    betty=mean(bet) # calculating mean
    GowerDude[n,]<-betty # storing result into a vector
  }
  
  GowerDude <- arrange(GowerDude,nrow) #arranging by small to large
  me = mean(GowerDude$nrow) # mean
  low = GowerDude[25,] #low CI
  up = GowerDude[975,] # up CI
  
  all=data.frame(me,low,up) # putting together
  all <- t(all) # flipping to make easier to integrate into summary tables
  all # output
  
}

#****************************************************************************

#D) Proportion functions


# EleRandProp function - inputs are datasets (BBA and BBS)
# Its parameters are the amount of random samples (1000 for us) and an empty dataframe
# It samples a given BBA/Adjacent dataset n=BBS times
# It then calculates the elevation range for the whole dataset and puts the results
# into a vector. 

EleRandProp <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  test1 <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp)) # Removing NAs
    x=sample(BBA,nrow(BBS)) # sampling from BBA data n=BBS rows
    ran=max(x$Max_Elev)-min(x$Min_Elev) # calculating elevation range for data
    test1[n,]<-ran # storing elevation range result into a vector
  }
  
  test1 # return vector
}



# TempRandProp function
TempRandProp <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  test1 <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS)) # sampling from BBA data n=BBS rows
    ran=max(x$Temp)-min(x$Temp) # calculating temperature range for data
    test1[n,]<-ran # storing temperature range result into a vector
  }
  
  test1 # return vector
}


PrecRandProp <- function(BBA,BBS){
  
  niter=1000 # number of iterations
  test1 <- data.frame(nrow=niter) # creating vector with row number=iterations
  
  for(n in 1:niter){ # creating forloop
    BBA <- BBA %>% 
      filter(!is.na(Elev_Range),
             !is.na(Shannon2006),
             !is.na(Precipitation),
             !is.na(Temp))
    x=sample(BBA,nrow(BBS)) # sampling from BBA data n=BBS rows
    ran=max(x$Precipitation)-min(x$Precipitation) # calculating precipitation range for data
    test1[n,]<-ran # storing precipitation range result into a vector
  }
  
  test1 # return vector
}


