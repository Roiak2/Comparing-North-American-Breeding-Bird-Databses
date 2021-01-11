###################################################################
################ SUMMARY TABLES FOR BIRD PROJECT ##################
###################################################################

# Roi Ankori-Karlinsky, Spring 2020


#This script creates summary tables for all figures for both predictions and assumption:
# 1. Loads functions to be used
# 2. Loads raw data
# 3. Summarizes alpha, beta, and gamma bird diversity for all 48 datasets for entire region
#     and for each of the five states (MA,MI,NY,PA,VT):
#     BBA, BBA(n=BBS), BBS, Adj, Adj(n=BBS), Sampled, Sampled(n=BBS), SampledAdj
# 4. Exports those three tables to csv files
# 5. Repeats the same process for environmental variables rather than bird diversity:
#     Land-class richness & diversity (Shannon-Wiener Index), Elevation mean and range,
#     Mean annual precipitation, Mean summer temperature


###################################################################


### Loading packages

library(tidyverse) #tidyverse for data cleaning and manipulation
library(mosaic) # mosaic for easy calculation of summary statistics
library(vegan) # vegan pacakge for calculating beta diversity
library(iNEXT) # iNEXT package for calculating regional richness
library(data.table) # data.table package for data manipulation
library(DescTools) # DescTools package for easy calculation of CIs


#set working directory
setwd('C:/Users/roiak/Documents/Ronen/Final/Final Code')

### Source function code
source("Ankori-Karlinsky_N.A_Birds_Functions.R")


### Loading datasets

path_files<- 'C:/Users/roiak/Documents/Ronen/Final' # Path files for easy loading and writing

source("Ankori-Karlinsky_N.A_Birds_Datasets.R")



###################################################################
################### BIRD DIVERSITY SUMMARY TABLES #################
###################################################################


#### Regional Richness Summary Table #####

#summary for BBA and BBS
Richness_Comparison <- tidy_full %>%
  filter(Observations=="1") %>% # keeping only observations of birds 
  group_by(Dataset,State) %>% # grouping by dataset,state
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n = n_distinct(ID)) %>%  #counting only unique observations per species
  mutate(Region = State) %>%
  dplyr::select(-State)

#Now summary table with adjacent
Adjacent_Richness_Comparison <- tidy_full %>% 
  filter(Observations=="1",Dataset=="BBA",BBS_Adjacent!="No"|is.na(BBS_Adjacent)) %>% # keeping only observations of birds  that are adjacent
  group_by(State) %>% # grouping by dataset,state,
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n = n_distinct(ID)) %>% # counting only unique observations per species
  mutate(Dataset= "Adjacent",
         Region = State) %>%
  dplyr::select(-State) #Saving as Adjacent

#Now summary table with well sampled
Sampled_Richness_Comparison <- tidy_full %>% 
  filter(Observations=="1",Dataset=="BBA",!is.na(Atlas2_Effort),Well_Sampled!="No"|is.na(Well_Sampled)) %>% # keeping only observations of birds with hour info that are well sampled
  group_by(State) %>% # grouping by dataset,state, 
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n = n_distinct(ID)) %>%
  mutate(Dataset="SampledBBA",
         Region = State) %>%
  dplyr::select(-State)

#Now summary table with info on both well sampled and adjacent blocks
Sampled_Adjacent_Richness_Comparison <- tidy_full %>% 
  filter(Observations=="1",!is.na(Atlas2_Effort),
         Dataset=="BBA",
         Well_Sampled!="No"|is.na(Well_Sampled),
         BBS_Adjacent!="No"|is.na(BBS_Adjacent)) %>% # keeping only observations of birds with hour info that are well sampled and adjacent
  group_by(State) %>% # grouping by dataset,state, 
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n= n_distinct(ID)) %>%
  mutate(Dataset="SampledAdjacent",
         Region = State) %>%
  dplyr::select(-State)

#Now making tables for entire region
BBA.BBS.Region <- tidy_full %>%
  filter(Observations=="1") %>% # keeping only observations of birds 
  group_by(Dataset) %>% # grouping by dataset,state
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n= n_distinct(ID)) # counting only unique observations per species

BBA.BBS.Region$Region <- "All"

#Adjacent
Adj.Region <- tidy_full %>%
  filter(Observations=="1",Dataset=="BBA",BBS_Adjacent!="No"|is.na(BBS_Adjacent)) %>% # keeping only observations of birds 
  group_by(Dataset) %>% # grouping by dataset,state
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n= n_distinct(ID))

Adj.Region$Dataset <- "Adjacent"
Adj.Region$Region <- "All"

#Sampled
Sample.Region <-tidy_full %>% 
  filter(Observations=="1",Dataset=="BBA",!is.na(Atlas2_Effort),Well_Sampled!="No"|is.na(Well_Sampled)) %>% # keeping only observations of birds with hour info that are well sampled
  group_by(Dataset) %>% # grouping by dataset,state, 
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n = n_distinct(ID))

Sample.Region$Dataset <- "SampledBBA"
Sample.Region$Region <- "All"

#Adjacent and Sampled
SampleAdj.Region <- tidy_full %>%
  filter(Observations=="1",!is.na(Atlas2_Effort),
         Dataset=="BBA",
         Well_Sampled!="No"|is.na(Well_Sampled),
         BBS_Adjacent!="No"|is.na(BBS_Adjacent)) %>% # keeping only observations of birds 
  group_by(Dataset) %>% # grouping by dataset,state
  dplyr::summarize(Regional_Richness = n_distinct(Species),
                   n = n_distinct(ID))

SampleAdj.Region$Dataset <- "SampledAdjacent"
SampleAdj.Region$Region <- "All"

df_list3 <- list(Richness_Comparison,Adjacent_Richness_Comparison, Sampled_Richness_Comparison,
                 Sampled_Adjacent_Richness_Comparison,BBA.BBS.Region,Adj.Region,
                 Sample.Region,SampleAdj.Region)

Sample_Size_Comparison <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list3, accumulate=FALSE)

Sample_Size_Comparison <- Sample_Size_Comparison %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))

Sample_Size <- Sample_Size_Comparison %>%
  dplyr::select(-Regional_Richness)


############### Adding BBA when n=BBS ##############


#Region column
Region <- c("All","MA","MI","NY","PA","VT")

#Estimate mean regional richness based on 1000 runs when n=BBS 
Estimated_Gamma_Richness <- data.frame(BirdRand(BBA,BBS),BirdRand(Adj,BBS),BirdRand(SampledBBA,BBS),
                                 BirdRand(MABBA,MABBS),BirdRand(MA_BBA_Adj,MABBS),BirdRand(SampledMABBA,MABBS),
                                 BirdRand(MIBBA,MIBBS),BirdRand(MI_BBA_Adj,MIBBS),BirdRand(SampledMIBBA,MIBBS),
                                 BirdRand(NYBBA,NYBBS),BirdRand(NY_BBA_Adj,NYBBS),BirdRand(SampledNYBBA,NYBBS),
                                 BirdRand(PABBA,PABBS),BirdRand(PA_BBA_Adj,PABBS),BirdRand(SampledPABBA,PABBS),
                                 BirdRand(VTBBA,VTBBS),BirdRand(VT_BBA_Adj,VTBBS),BirdRand(SampledVTBBA,VTBBS))

#Labelling columns
colnames(Estimated_Gamma_Richness) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
BBA.nBBS <- Estimated_Gamma_Richness %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Regional_Richness = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Regional_Richness),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
BBA.nBBS<-BBA.nBBS %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
BBA.nBBS$Dataset <- as.character(paste(BBA.nBBS$Dataset, "(n=BBS)"))

#Combining with larger data
df_list4 <- list(Sample_Size_Comparison, BBA.nBBS)

Sample_Size_Comparison2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list4, accumulate=FALSE)

setcolorder(Sample_Size_Comparison2, c("Dataset","Region","Regional_Richness",
                                       "Low.CI","Up.CI","n"))

Sample_Size_Comparison2 <- Sample_Size_Comparison2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


#Writing datatable

write_csv(Sample_Size_Comparison2,paste0(path_files,"/Final Data/Gamma.csv"))


#removing files from environment to save memory
rm(Richness_Comparison, # summary tables for regional richness
   Adjacent_Richness_Comparison,
   Sampled_Richness_Comparison,
   Sampled_Adjacent_Richness_Comparison,
   BBA.BBS.Region,
   Adj.Region,
   Sample.Region,
   SampleAdj.Region,
   BBA.nBBS, # n=BBS estimated tables
   Estimated_Gamma_Richness,
   Sample_Size_Comparison # Final tables
   #Sample_Size_Comparison2
   )


#######################################################

#######################################################


#### Local Richness Summary

#dataframe for BBA and BBS 
Local_Richness <- data.frame(MeanCI(BBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(BBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledBBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledAdj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MABBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MABBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MA_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledMABBA$Atlas2_Total,conf.level=0.95),
                             MeanCI(SampledMA_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MIBBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MIBBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(MI_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledMIBBA$Atlas2_Total,conf.level=0.95),
                             MeanCI(SampledMI_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(NYBBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(NYBBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(NY_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledNYBBA$Atlas2_Total,conf.level=0.95),
                             MeanCI(SampledNY_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(PABBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(PABBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(PA_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledPABBA$Atlas2_Total,conf.level=0.95),
                             MeanCI(SampledPA_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(VTBBS$Atlas2_Total,conf.level = 0.95),
                             MeanCI(VTBBA$Atlas2_Total,conf.level = 0.95),
                             MeanCI(VT_BBA_Adj$Atlas2_Total,conf.level = 0.95),
                             MeanCI(SampledVTBBA$Atlas2_Total,conf.level=0.95),
                             MeanCI(SampledVT_BBA_Adj$Atlas2_Total,conf.level = 0.95)
                            )

#Column names
colnames(Local_Richness) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                             "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                             "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                             "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                             "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                             "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
                             )


#Pivoting dataframe
Local_Richness <- Local_Richness %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Local_Richness = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Local_Richness),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Local_Richness<-Local_Richness%>%
  separate(Region, into = c("Region","Dataset"))

Local_Richness <- Local_Richness %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))

#Adding sample size to local richness table

df_list <- list(Local_Richness, Sample_Size)

Local_Richness1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list, accumulate=FALSE)

Local_Richness1 <- Local_Richness1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


#### Adding data for BBA when n=BBS

#Estimate mean local richness based on 1000 runs when n=BBS 
Estimated_Alpha_Richness <- data.frame(BirdLoc(BBA,BBS),BirdLoc(Adj,BBS),BirdLoc(SampledBBA,BBS),
                                 BirdLoc(MABBA,MABBS),BirdLoc(MA_BBA_Adj,MABBS),BirdLoc(SampledMABBA,MABBS),
                                 BirdLoc(MIBBA,MIBBS),BirdLoc(MI_BBA_Adj,MIBBS),BirdLoc(SampledMIBBA,MIBBS),
                                 BirdLoc(NYBBA,NYBBS),BirdLoc(NY_BBA_Adj,NYBBS),BirdLoc(SampledNYBBA,NYBBS),
                                 BirdLoc(PABBA,PABBS),BirdLoc(PA_BBA_Adj,PABBS),BirdLoc(SampledPABBA,PABBS),
                                 BirdLoc(VTBBA,VTBBS),BirdLoc(VT_BBA_Adj,VTBBS),BirdLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Alpha_Richness) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
BBA.nBBS <- Estimated_Alpha_Richness %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Local_Richness = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Local_Richness),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
BBA.nBBS<-BBA.nBBS %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
BBA.nBBS$Dataset <- as.character(paste(BBA.nBBS$Dataset, "(n=BBS)"))

#Combining with larger data
df_list2 <- list(Local_Richness1, BBA.nBBS)

Local_Richness2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list2, accumulate=FALSE)

setcolorder(Local_Richness2, c("Dataset","Region","Local_Richness",
                                       "Low.CI","Up.CI","n"))

Local_Richness2 <- Local_Richness2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


#Writing datatable

write_csv(Local_Richness2,paste0(path_files,"/Final Data/Alpha.csv"))


#removing files from environment to save memory
rm(Local_Richness,Local_Richness1,#local richness tables
   BBA.nBBS, # estimated local richness table for n=BBS
   Estimated_Alpha_Richness # Estimated horizontal local richness
   #Local_Richness2 # final table
   )


######################################################



############ SUMMARIZING BETA DIVERSITY TABLE #################

#Calculating Jaccard Index 
bs=vegdist(BBS[14:291],method="jaccard",binary=T)
#Turning into a vector
bs <- as.vector(bs)
#Checking
#str(bs)

#now for BBA
ba=vegdist(BBA[14:291],method="jaccard",binary=T)
#Turning into a vector
ba <- as.vector(ba)
#Checking
#str(ba)

mabs=vegdist(MABBS[14:291],method="jaccard",binary=T)
mabs <- as.vector(mabs)

maba=vegdist(MABBA[14:291],method="jaccard",binary=T)
maba <- as.vector(maba)

mibs=vegdist(MIBBS[14:291],method="jaccard",binary=T)
mibs <- as.vector(mibs)

miba=vegdist(MIBBA[14:291],method="jaccard",binary=T)
miba <- as.vector(miba)

nybs=vegdist(NYBBS[14:291],method="jaccard",binary=T)
nybs <- as.vector(nybs)

nyba=vegdist(NYBBA[14:291],method="jaccard",binary=T)
nyba <- as.vector(nyba)

pabs=vegdist(PABBS[14:291],method="jaccard",binary=T)
pabs <- as.vector(pabs)

paba=vegdist(PABBA[14:291],method="jaccard",binary=T)
paba<-as.vector(paba)

vtbs=vegdist(VTBBS[14:291],method="jaccard",binary=T)
vtbs <- as.vector(vtbs)

vtba=vegdist(VTBBA[14:291],method="jaccard",binary=T)
vtba <- as.vector(vtba)


#Calculating Jaccard Index for adjacent and well sampled  blocks

#First sampled
sampled.bba=vegdist(SampledBBA[14:291],method="jaccard",binary=T)
#Turning into a vector
sampled.bba <- as.vector(sampled.bba)
#Checking
#str(ba)

sampled.ma=vegdist(SampledMABBA[14:291],method="jaccard",binary=T)
sampled.ma <- as.vector(sampled.ma)

sampled.mi=vegdist(SampledMIBBA[14:291],method="jaccard",binary=T)
sampled.mi <- as.vector(sampled.mi)

sampled.ny=vegdist(SampledNYBBA[14:291],method="jaccard",binary=T)
sampled.ny <- as.vector(sampled.ny)

sampled.pa=vegdist(SampledPABBA[14:291],method="jaccard",binary=T)
sampled.pa<-as.vector(sampled.pa)

sampled.vt=vegdist(SampledVTBBA[14:291],method="jaccard",binary=T)
sampled.vt<- as.vector(sampled.vt)


#Now adjacent

adj=vegdist(Adj[14:291],method="jaccard",binary=T)
#Turning into a vector
adj <- as.vector(adj)
#Checking
#str(adj)

maadj=vegdist(MA_BBA_Adj[14:291],method="jaccard",binary=T)
maadj <- as.vector(maadj)

miadj=vegdist(MI_BBA_Adj[14:291],method="jaccard",binary=T)
miadj <- as.vector(miadj)

nyadj=vegdist(NY_BBA_Adj[14:291],method="jaccard",binary=T)
nyadj <- as.vector(nyadj)

paadj=vegdist(PA_BBA_Adj[14:291],method="jaccard",binary=T)
paadj<-as.vector(paadj)

vtadj=vegdist(VT_BBA_Adj[14:291],method="jaccard",binary=T)
vtadj <- as.vector(vtadj)

######## Now well sampled adjacent ##########

sampled.adj=vegdist(SampledAdj[14:291],method="jaccard",binary=T)
#Turning into a vector
sampled.adj <- as.vector(sampled.adj)
#Checking
#str(ba)

sampled.maadj=vegdist(SampledMA_BBA_Adj[14:291],method="jaccard",binary=T)
sampled.maadj <- as.vector(sampled.maadj)

sampled.miadj=vegdist(SampledMI_BBA_Adj[14:291],method="jaccard",binary=T)
sampled.miadj <- as.vector(sampled.miadj)

sampled.nyadj=vegdist(SampledNY_BBA_Adj[14:291],method="jaccard",binary=T)
sampled.nyadj <- as.vector(sampled.nyadj)

sampled.paadj=vegdist(SampledPA_BBA_Adj[14:291],method="jaccard",binary=T)
sampled.paadj<-as.vector(sampled.paadj)

sampled.vtadj=vegdist(SampledVT_BBA_Adj[14:291],method="jaccard",binary=T)
sampled.vtadj <- as.vector(sampled.vtadj)


#Comparing Means by calculating means and CIs into a dataframe
Species_Dissimilarity <- data.frame(MeanCI(bs,conf.level = 0.95),MeanCI(ba,conf.level = 0.95),
                                    MeanCI(adj,conf.level = 0.95),MeanCI(sampled.bba,conf.level=0.95),
                                    MeanCI(sampled.adj,conf.level = 0.95),
                                    MeanCI(mabs,conf.level = 0.95),MeanCI(maba,conf.level = 0.95),
                                    MeanCI(maadj,conf.level = 0.95),MeanCI(sampled.ma,conf.level=0.95),
                                    MeanCI(sampled.maadj,conf.level = 0.95),
                                    MeanCI(mibs,conf.level = 0.95),MeanCI(miba,conf.level = 0.95),
                                    MeanCI(miadj,conf.level = 0.95),MeanCI(sampled.mi,conf.level=0.95),
                                    MeanCI(sampled.miadj,conf.level = 0.95),
                                    MeanCI(nybs,conf.level = 0.95),MeanCI(nyba,conf.level = 0.95),
                                    MeanCI(nyadj,conf.level = 0.95),MeanCI(sampled.ny,conf.level=0.95),
                                    MeanCI(sampled.nyadj,conf.level = 0.95),
                                    MeanCI(pabs,conf.level = 0.95),MeanCI(paba,conf.level = 0.95),
                                    MeanCI(paadj,conf.level = 0.95),MeanCI(sampled.pa,conf.level=0.95),
                                    MeanCI(sampled.paadj,conf.level = 0.95),
                                    MeanCI(vtbs,conf.level = 0.95),MeanCI(vtba,conf.level = 0.95),
                                    MeanCI(vtadj,conf.level = 0.95),MeanCI(sampled.vt,conf.level=0.95),
                                    MeanCI(sampled.vtadj,conf.level = 0.95)
                                    )

#Column names
colnames(Species_Dissimilarity) <- c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA","All_SampledAdjacent",
                                     "MA_BBS","MA_BBA", "MA_Adjacent","MA_SampledBBA","MA_SampledAdjacent",
                                     "MI_BBS","MI_BBA", "MI_Adjacent","MI_SampledBBA","MI_SampledAdjacent",
                                     "NY_BBS","NY_BBA", "NY_Adjacent","NY_SampledBBA","NY_SampledAdjacent",
                                     "PA_BBS","PA_BBA", "PA_Adjacent","PA_SampledBBA","PA_SampledAdjacent",
                                     "VT_BBS","VT_BBA", "VT_Adjacent","VT_SampledBBA","VT_SampledAdjacent")

#Turning into long format
Species_Dissimilarity_Table <- Species_Dissimilarity %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Dataset",
               values_to = "Value") %>%
  group_by(Dataset) %>%
  mutate(Mean_Dissimilarity = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean_Dissimilarity),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating region and data column
Species_Dissimilarity_Table <- Species_Dissimilarity_Table %>%
  separate(Dataset, into = c("Region","Dataset")) 


#More cleaning
setcolorder(Species_Dissimilarity_Table,c("Dataset","Region","Mean_Dissimilarity",
                                          "Low.CI","Up.CI"))

##### Adding for BBA when n=BBS

#Estimate mean regional richness based on 1000 runs when n=BBS 
Estimated_Dissimilarity <- data.frame(JacRand(BBA,BBS),JacRand(Adj,BBS),JacRand(SampledBBA,BBS),
                                 JacRand(MABBA,MABBS),JacRand(MA_BBA_Adj,MABBS),JacRand(SampledMABBA,MABBS),
                                 JacRand(MIBBA,MIBBS),JacRand(MI_BBA_Adj,MIBBS),JacRand(SampledMIBBA,MIBBS),
                                 JacRand(NYBBA,NYBBS),JacRand(NY_BBA_Adj,NYBBS),JacRand(SampledNYBBA,NYBBS),
                                 JacRand(PABBA,PABBS),JacRand(PA_BBA_Adj,PABBS),JacRand(SampledPABBA,PABBS),
                                 JacRand(VTBBA,VTBBS),JacRand(VT_BBA_Adj,VTBBS),JacRand(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Dissimilarity) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
BBA.nBBS <- Estimated_Dissimilarity %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean_Dissimilarity = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean_Dissimilarity),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
BBA.nBBS<-BBA.nBBS %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
BBA.nBBS$Dataset <- as.character(paste(BBA.nBBS$Dataset, "(n=BBS)"))

#Combining with larger data
df_list5 <- list(Species_Dissimilarity_Table, BBA.nBBS)

Beta <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list5, accumulate=FALSE)

setcolorder(Beta, c("Dataset","Region","Mean_Dissimilarity",
                                       "Low.CI","Up.CI"))

Beta <- Beta %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


#Writing datatable

write_csv(Beta,paste0(path_files,"/Final Data/Beta.csv"))


#removing objects from environment to save space for next round!
rm(
  df_list,df_list2,df_list3,df_list4,df_list5, #all the lists
  Estimated_Dissimilarity,Species_Dissimilarity,Species_Dissimilarity_Table, BBA.nBBS, #dissimilarity tables
  bs,ba,adj,sampled.bba,sampled.adj, # regional Jaccard indices
  mabs,maba,maadj,sampled.ma,sampled.maadj, # MA Jaccard indices
  mibs,miba,miadj,sampled.mi,sampled.miadj, # MI Jaccard indices
  nybs,nyba,nyadj,sampled.ny,sampled.nyadj, # NY Jaccard indices
  pabs,paba,paadj,sampled.pa,sampled.paadj, # PA Jaccard indices
  vtbs,vtba,vtadj,sampled.vt,sampled.vtadj  # VT Jaccard indices
  #Beta, # the beta final table
)


#*******************************************************************

####################################################################
###################### Environmental Variables #####################
####################################################################


###########################################################
###### Getting data for all states and datasets ###########

Environmentals <- full_d%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset,State) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = State) %>%
  dplyr::select(-State)

#for entire dataset
Environmentals1 <- full_d%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = "All")

#Joining datasets
Environmentals2 <- rbind.data.frame(Environmentals,Environmentals1) %>%
  arrange(Region,Dataset)

#setting columns in order
setcolorder(Environmentals2,c("Dataset","Region",
                                "Gamma_Elevation",
                                "Gamma_Precipitation",
                                "Gamma_Temp",
                                "Gamma_Het"))


# Now Adjacent
Adj.Enviro <- Adj%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset,State) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = State) %>%
  dplyr::select(-State)

Adj.Enviro$Dataset <- "Adjacent"

#for entire dataset
Adj.Enviro1 <- Adj%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = "All")

Adj.Enviro1$Dataset <- "Adjacent"

#Joining datasets
Adj.Enviro2 <- rbind(Adj.Enviro,Adj.Enviro1) %>%
  arrange(Region,Dataset)

#setting columns in order
setcolorder(Adj.Enviro2,c("Dataset","Region",
                              "Gamma_Elevation",
                              "Gamma_Precipitation",
                              "Gamma_Temp",
                              "Gamma_Het"))


#Now Sampled

Sampled.Enviro <- SampledBBA%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset,State) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = State) %>%
  dplyr::select(-State)

Sampled.Enviro$Dataset <- "SampledBBA"

#for entire dataset
Sampled.Enviro1 <- SampledBBA%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = "All")

Sampled.Enviro1$Dataset <- "SampledBBA"

#Joining datasets
Sampled.Enviro2 <- rbind(Sampled.Enviro,Sampled.Enviro1) %>%
  arrange(Region,Dataset)

#setting columns in order
setcolorder(Sampled.Enviro2,c("Dataset","Region",
                          "Gamma_Elevation",
                          "Gamma_Precipitation",
                          "Gamma_Temp",
                          "Gamma_Het"))



#Sampled and Adjacent

SampledAdj.Enviro <- SampledAdj%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset,State) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = State) %>%
  dplyr::select(-State)

SampledAdj.Enviro$Dataset <- "SampledAdjacent"

#for entire dataset
SampledAdj.Enviro1 <- SampledAdj%>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp)) %>%
  group_by(Dataset) %>%
  summarize(Gamma_Elevation = (max(Max_Elev)-min(Min_Elev)),
            Gamma_Precipitation = max(Precipitation)-min(Precipitation),
            Gamma_Temp = max(Temp)-min(Temp),
            Gamma_Het = max(Het_Count2006)) %>%
  mutate(Region = "All")

SampledAdj.Enviro1$Dataset <- "SampledAdjacent"

#Joining datasets
SampledAdj.Enviro2 <- rbind(SampledAdj.Enviro,SampledAdj.Enviro1) %>%
  arrange(Region,Dataset)

#setting columns in order
setcolorder(SampledAdj.Enviro2,c("Dataset","Region",
                              "Gamma_Elevation",
                              "Gamma_Precipitation",
                              "Gamma_Temp",
                              "Gamma_Het"))


##Joining all these together

df_list6 <- list(Environmentals2,Adj.Enviro2,Sampled.Enviro2,SampledAdj.Enviro2)

Enviro <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list6, accumulate=FALSE)

Enviro <- Enviro %>% 
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


############### Adding BBA when n=BBS ##############


#Region column
Region <- c("All","MA","MI","NY","PA","VT")

#Estimate regional elevation range based on 1000 runs when n=BBS 
Estimated_Elev <- data.frame(EleRand(BBA,BBS),EleRand(Adj,BBS),EleRand(SampledBBA,BBS),
                                 EleRand(MABBA,MABBS),EleRand(MA_BBA_Adj,MABBS),EleRand(SampledMABBA,MABBS),
                                 EleRand(MIBBA,MIBBS),EleRand(MI_BBA_Adj,MIBBS),EleRand(SampledMIBBA,MIBBS),
                                 EleRand(NYBBA,NYBBS),EleRand(NY_BBA_Adj,NYBBS),EleRand(SampledNYBBA,NYBBS),
                                 EleRand(PABBA,PABBS),EleRand(PA_BBA_Adj,PABBS),EleRand(SampledPABBA,PABBS),
                                 EleRand(VTBBA,VTBBS),EleRand(VT_BBA_Adj,VTBBS),EleRand(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Elev) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Elev <- Estimated_Elev %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gamma_Elevation = Value[c(NA,1,NA)],
         Elev_Low.CI = Value[c(NA,2,NA)],
         Elev_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gamma_Elevation),!is.na(Elev_Low.CI),!is.na(Elev_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Elev<-Estimated_Elev %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Elev$Dataset <- as.character(paste(Estimated_Elev$Dataset, "(n=BBS)"))

#Arranging dataset labels
Estimated_Elev <- Estimated_Elev %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA (n=BBS)",
                                   "Adjacent (n=BBS)",
                                   "SampledBBA (n=BBS)")))


#Estimate temperature range regionally based on 1000 runs when n=BBS 
Estimated_Temp <- data.frame(TempRand(BBA,BBS),TempRand(Adj,BBS),TempRand(SampledBBA,BBS),
                             TempRand(MABBA,MABBS),TempRand(MA_BBA_Adj,MABBS),TempRand(SampledMABBA,MABBS),
                             TempRand(MIBBA,MIBBS),TempRand(MI_BBA_Adj,MIBBS),TempRand(SampledMIBBA,MIBBS),
                             TempRand(NYBBA,NYBBS),TempRand(NY_BBA_Adj,NYBBS),TempRand(SampledNYBBA,NYBBS),
                             TempRand(PABBA,PABBS),TempRand(PA_BBA_Adj,PABBS),TempRand(SampledPABBA,PABBS),
                             TempRand(VTBBA,VTBBS),TempRand(VT_BBA_Adj,VTBBS),TempRand(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Temp) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Temp <- Estimated_Temp %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gamma_Temp = Value[c(NA,1,NA)],
         Temp_Low.CI = Value[c(NA,2,NA)],
         Temp_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gamma_Temp),!is.na(Temp_Low.CI),!is.na(Temp_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Temp<-Estimated_Temp %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Temp$Dataset <- as.character(paste(Estimated_Temp$Dataset, "(n=BBS)"))

#Arranging
Estimated_Temp <- Estimated_Temp %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA (n=BBS)",
                                   "Adjacent (n=BBS)",
                                   "SampledBBA (n=BBS)")))


#### Now precipitation
Estimated_Prec <- data.frame(PrecRand(BBA,BBS),PrecRand(Adj,BBS),PrecRand(SampledBBA,BBS),
                             PrecRand(MABBA,MABBS),PrecRand(MA_BBA_Adj,MABBS),PrecRand(SampledMABBA,MABBS),
                             PrecRand(MIBBA,MIBBS),PrecRand(MI_BBA_Adj,MIBBS),PrecRand(SampledMIBBA,MIBBS),
                             PrecRand(NYBBA,NYBBS),PrecRand(NY_BBA_Adj,NYBBS),PrecRand(SampledNYBBA,NYBBS),
                             PrecRand(PABBA,PABBS),PrecRand(PA_BBA_Adj,PABBS),PrecRand(SampledPABBA,PABBS),
                             PrecRand(VTBBA,VTBBS),PrecRand(VT_BBA_Adj,VTBBS),PrecRand(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Prec) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Prec <- Estimated_Prec %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gamma_Precipitation = Value[c(NA,1,NA)],
         Prec_Low.CI = Value[c(NA,2,NA)],
         Prec_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gamma_Precipitation),!is.na(Prec_Low.CI),!is.na(Prec_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Prec<-Estimated_Prec %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Prec$Dataset <- as.character(paste(Estimated_Prec$Dataset, "(n=BBS)"))

#Arranging
Estimated_Prec <- Estimated_Prec %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA (n=BBS)",
                                   "Adjacent (n=BBS)",
                                   "SampledBBA (n=BBS)")))


## Now for Land-Class Richness

Estimated_Het <- data.frame(HetRand(BBA,BBS),HetRand(Adj,BBS),HetRand(SampledBBA,BBS),
                             HetRand(MABBA,MABBS),HetRand(MA_BBA_Adj,MABBS),HetRand(SampledMABBA,MABBS),
                             HetRand(MIBBA,MIBBS),HetRand(MI_BBA_Adj,MIBBS),HetRand(SampledMIBBA,MIBBS),
                             HetRand(NYBBA,NYBBS),HetRand(NY_BBA_Adj,NYBBS),HetRand(SampledNYBBA,NYBBS),
                             HetRand(PABBA,PABBS),HetRand(PA_BBA_Adj,PABBS),HetRand(SampledPABBA,PABBS),
                             HetRand(VTBBA,VTBBS),HetRand(VT_BBA_Adj,VTBBS),HetRand(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Het) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Het <- Estimated_Het %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gamma_Het = Value[c(NA,1,NA)],
         Het_Low.CI = Value[c(NA,2,NA)],
         Het_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gamma_Het),!is.na(Het_Low.CI),!is.na(Het_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Het<-Estimated_Het %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Het$Dataset <- as.character(paste(Estimated_Het$Dataset, "(n=BBS)"))

#Arranging
Estimated_Het <- Estimated_Het %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA (n=BBS)",
                                   "Adjacent (n=BBS)",
                                   "SampledBBA (n=BBS)")))


#Combining with larger data
df_list7 <- list(Estimated_Elev,Estimated_Temp,Estimated_Prec,Estimated_Het,Enviro)

Gamma_Enviro <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list7, accumulate=FALSE)

setcolorder(Gamma_Enviro, c("Dataset","Region",
                            "Gamma_Elevation","Elev_Low.CI","Elev_Up.CI",
                            "Gamma_Precipitation","Prec_Low.CI","Prec_Up.CI",
                            "Gamma_Temp","Temp_Low.CI","Temp_Up.CI",
                            "Gamma_Het","Het_Low.CI","Het_Up.CI"))

Gamma_Enviro <- Gamma_Enviro %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


#Writing datatable

write_csv(Gamma_Enviro,paste0(path_files,"/Final Data/Gamma_Enviro.csv"))

# removing files from environment to save memory
rm(
  df_list6,df_list7, #lists
  Enviro, Environmentals, Environmentals1, Environmentals2, # joined tables
  Estimated_Elev, Estimated_Prec, Estimated_Temp, Estimated_Het, # Tables when n=BBS
  Adj.Enviro, Adj.Enviro1, Adj.Enviro2, #Observed gamma
  Sampled.Enviro, Sampled.Enviro1, Sampled.Enviro2, 
  SampledAdj.Enviro,SampledAdj.Enviro1,SampledAdj.Enviro2
  #Gamma_Enviro # Final table
  )


##############################################################################


##############################################################################
######################### ALPHA ENVIRONMENTAL CONDITIONS ##################### 
##############################################################################


# Making sure all dataframes don't have NA values 

#Reloading data without NA values
full_d <- full_d %>%
  filter(!is.na(Elev_Range),
         !is.na(Shannon2006),
         !is.na(Precipitation),
         !is.na(Temp))

#Creating dataset with only well sampled BBA blocks (minimum 30 hours, no more than 300 hours of outliers)
SampledFull <- filter(full_d,ifelse(Dataset=="BBA",
                                    Atlas2_Effort > 29 & Atlas2_Effort < 300,Atlas2_Effort>0))

#Creating dataset with only BBS-adjacent BBA blocks
Adjacent_Full <- filter(full_d,ifelse(Dataset=="BBA",
                                      BBS_Adjacent=="Yes",is.na(BBS_Adjacent)))

#Creating dataset with well sampled adjacent blocks
SampledAdj_Full <- filter(SampledFull,ifelse(Dataset=="BBA",
                                             BBS_Adjacent=="Yes",is.na(BBS_Adjacent)))

#Creating BBA only dataset
BBA <- filter(full_d, Dataset=="BBA")

#Creating BBS only dataset
BBS <- filter(full_d, Dataset=="BBS")

#Creating BBA only  dataset that is Adjacent to BBS routes (Adj)
Adj <- filter(Adjacent_Full, Dataset=="BBA")

#Creating BBA only dataset with sufficient sampling 
SampledBBA <- filter(SampledFull, Dataset=="BBA")

#Creating BBA only dataset that is well sampled and also adjacent to BBS routes
SampledAdj <- filter(SampledAdj_Full,Dataset=="BBA")


#Tidy data

#First only keeping dataset, Plot ID, State, Hours Effort and species occurrences

tidy_full <- full_d[c(1:3,7,13:291)]

#Now using pivot_longer() to turn species into a column

tidy_full<-tidy_full %>%
  pivot_longer(c("Acadian Flycatcher":"Yellow-throated Warbler"), #columns
               names_to = "Species", # name new variable
               values_to = "Observations" # name value column
  )

#Creating column that says whether BBA is sufficiently sampled
tidy_full <-tidy_full %>%
  mutate(Well_Sampled = ifelse(Dataset=="BBS",NA,ifelse(
    Dataset=="BBA" & Atlas2_Effort>29 & Atlas2_Effort<300,"Yes","No")))


############ MASSACHUSETTS ################

#All
MA <-filter(full_d,State=="MA")

#Well Sampled
SampledMA <- filter(SampledFull, State=="MA")

#Adjacent
MAAdj <- filter(Adjacent_Full, State=="MA")

#Well sampled and Adjacent
SampledMAAdj <- filter(SampledAdj_Full, State=="MA")

#BBS
MABBS <- filter(BBS, State=="MA")
#BBA
MABBA <- filter(BBA, State=="MA")
#Sampled BBA
SampledMABBA <- filter(SampledBBA, State=="MA")
#BBS-Adjacent BBA blocks
MA_BBA_Adj <- filter(Adj, State=="MA")
#Sampled adjacent
SampledMA_BBA_Adj <- filter(SampledAdj, State=="MA")


############ MICHIGAN ################

#All
MI <-filter(full_d,State=="MI")

#Well Sampled
SampledMI <- filter(SampledFull, State=="MI")

#Adjacent
MIAdj <- filter(Adjacent_Full, State=="MI")

#Well sampled and Adjacent
SampledMIAdj <- filter(SampledAdj_Full, State=="MI")

#BBS
MIBBS <- filter(BBS, State=="MI")
#BBA
MIBBA <- filter(BBA, State=="MI")
#Sampled BBA
SampledMIBBA <- filter(SampledBBA, State=="MI")
#BBS-Adjacent BBA blocks
MI_BBA_Adj <- filter(Adj, State=="MI")
#Sampled adjacent
SampledMI_BBA_Adj <- filter(SampledAdj, State=="MI")


############ NEW YORK ################

#All
NY <-filter(full_d,State=="NY")

#Well Sampled
SampledNY <- filter(SampledFull, State=="NY")

#Adjacent
NYAdj <- filter(Adjacent_Full, State=="NY")

#Well sampled and Adjacent
SampledNYAdj <- filter(SampledAdj_Full, State=="NY")

#BBS
NYBBS <- filter(BBS, State=="NY")
#BBA
NYBBA <- filter(BBA, State=="NY")
#Sampled BBA
SampledNYBBA <- filter(SampledBBA, State=="NY")
#BBS-Adjacent BBA blocks
NY_BBA_Adj <- filter(Adj, State=="NY")
#Sampled adjacent
SampledNY_BBA_Adj <- filter(SampledAdj, State=="NY")


############ PENNSYLVANIA ################

#All
PA <-filter(full_d,State=="PA")

#Well Sampled
SampledPA <- filter(SampledFull, State=="PA")

#Adjacent
PAAdj <- filter(Adjacent_Full, State=="PA")

#Well sampled and Adjacent
SampledPAAdj <- filter(SampledAdj_Full, State=="PA")

#BBS
PABBS <- filter(BBS, State=="PA")
#BBA
PABBA <- filter(BBA, State=="PA")
#Sampled BBA
SampledPABBA <- filter(SampledBBA, State=="PA")
#BBS-Adjacent BBA blocks
PA_BBA_Adj <- filter(Adj, State=="PA")
#Sampled adjacent
SampledPA_BBA_Adj <- filter(SampledAdj, State=="PA")


############ VERMONT ################

#All
VT <-filter(full_d,State=="VT")

#Well Sampled
SampledVT <- filter(SampledFull, State=="VT")

#Adjacent
VTAdj <- filter(Adjacent_Full, State=="VT")

#Well sampled and Adjacent
SampledVTAdj <- filter(SampledAdj_Full, State=="VT")

#BBS
VTBBS <- filter(BBS, State=="VT")
#BBA
VTBBA <- filter(BBA, State=="VT")
#Sampled BBA
SampledVTBBA <- filter(SampledBBA, State=="VT")
#BBS-Adjacent BBA blocks
VT_BBA_Adj <- filter(Adj, State=="VT")
#Sampled adjacent
SampledVT_BBA_Adj <- filter(SampledAdj, State=="VT")


#******************************************************************


######## Elevation Mean

Alpha.Elev <- data.frame(MeanCI(BBS$Elev_Mean,conf.level = 0.95,na.rm=T),
                                             MeanCI(BBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledBBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledAdj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MABBS$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MABBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MA_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledMABBA$Elev_Mean,conf.level=0.95),
                                             MeanCI(SampledMA_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MIBBS$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MIBBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(MI_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledMIBBA$Elev_Mean,conf.level=0.95),
                                             MeanCI(SampledMI_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(NYBBS$Elev_Mean,conf.level = 0.95),
                                             MeanCI(NYBBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(NY_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledNYBBA$Elev_Mean,conf.level=0.95),
                                             MeanCI(SampledNY_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(PABBS$Elev_Mean,conf.level = 0.95),
                                             MeanCI(PABBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(PA_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledPABBA$Elev_Mean,conf.level=0.95),
                                             MeanCI(SampledPA_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(VTBBS$Elev_Mean,conf.level = 0.95),
                                             MeanCI(VTBBA$Elev_Mean,conf.level = 0.95),
                                             MeanCI(VT_BBA_Adj$Elev_Mean,conf.level = 0.95),
                                             MeanCI(SampledVTBBA$Elev_Mean,conf.level=0.95),
                                             MeanCI(SampledVT_BBA_Adj$Elev_Mean,conf.level = 0.95)
                         )

#Column names
colnames(Alpha.Elev) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                             "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                             "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                             "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                             "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                             "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
                         )


#Creating dataframe (pivoting to long-form) 
Alpha.Elev <- Alpha.Elev %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Elevation = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Elevation),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.Elev<-Alpha.Elev%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.Elev <- Alpha.Elev %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean local elevation based on 1000 runs when n=BBS 
Estimated_Elev <- data.frame(ElevLoc(BBA,BBS),ElevLoc(Adj,BBS),ElevLoc(SampledBBA,BBS),
                                 ElevLoc(MABBA,MABBS),ElevLoc(MA_BBA_Adj,MABBS),ElevLoc(SampledMABBA,MABBS),
                                 ElevLoc(MIBBA,MIBBS),ElevLoc(MI_BBA_Adj,MIBBS),ElevLoc(SampledMIBBA,MIBBS),
                                 ElevLoc(NYBBA,NYBBS),ElevLoc(NY_BBA_Adj,NYBBS),ElevLoc(SampledNYBBA,NYBBS),
                                 ElevLoc(PABBA,PABBS),ElevLoc(PA_BBA_Adj,PABBS),ElevLoc(SampledPABBA,PABBS),
                                 ElevLoc(VTBBA,VTBBS),ElevLoc(VT_BBA_Adj,VTBBS),ElevLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Elev) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Elev <- Estimated_Elev %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Elevation = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Elevation),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Elev<-Estimated_Elev %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Elev$Dataset <- as.character(paste(Estimated_Elev$Dataset, "(n=BBS)"))

#Combining with larger data
df_list8 <- list(Alpha.Elev, Estimated_Elev)

Alpha.Elev2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list8, accumulate=FALSE)

#Column order
setcolorder(Alpha.Elev2, c("Dataset","Region","Mean.Elevation",
                               "Low.CI","Up.CI"))

#Arranging by Dataset and Region
Alpha.Elev2 <- Alpha.Elev2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Making sure column names are right
colnames(Alpha.Elev2) <- c("Dataset","Region",
                           "Mean_Elevation",
                           "Elev_LowCI",
                           "Elev_UpCI")


##### Now for Precipitation
  
Alpha.Prec <- data.frame(MeanCI(BBS$Precipitation,conf.level = 0.95,na.rm=T),
                         MeanCI(BBA$Precipitation,conf.level = 0.95),
                         MeanCI(Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledBBA$Precipitation,conf.level = 0.95),
                         MeanCI(SampledAdj$Precipitation,conf.level = 0.95),
                         MeanCI(MABBS$Precipitation,conf.level = 0.95),
                         MeanCI(MABBA$Precipitation,conf.level = 0.95),
                         MeanCI(MA_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledMABBA$Precipitation,conf.level=0.95),
                         MeanCI(SampledMA_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(MIBBS$Precipitation,conf.level = 0.95),
                         MeanCI(MIBBA$Precipitation,conf.level = 0.95),
                         MeanCI(MI_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledMIBBA$Precipitation,conf.level=0.95),
                         MeanCI(SampledMI_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(NYBBS$Precipitation,conf.level = 0.95),
                         MeanCI(NYBBA$Precipitation,conf.level = 0.95),
                         MeanCI(NY_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledNYBBA$Precipitation,conf.level=0.95),
                         MeanCI(SampledNY_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(PABBS$Precipitation,conf.level = 0.95),
                         MeanCI(PABBA$Precipitation,conf.level = 0.95),
                         MeanCI(PA_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledPABBA$Precipitation,conf.level=0.95),
                         MeanCI(SampledPA_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(VTBBS$Precipitation,conf.level = 0.95),
                         MeanCI(VTBBA$Precipitation,conf.level = 0.95),
                         MeanCI(VT_BBA_Adj$Precipitation,conf.level = 0.95),
                         MeanCI(SampledVTBBA$Precipitation,conf.level=0.95),
                         MeanCI(SampledVT_BBA_Adj$Precipitation,conf.level = 0.95)
)

#Column names
colnames(Alpha.Prec) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                         "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                         "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                         "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                         "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                         "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
                         )


#Creating dataframe 
Alpha.Prec <- Alpha.Prec %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Precipitation = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Precipitation),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.Prec<-Alpha.Prec%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.Prec <- Alpha.Prec %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean local precipitation based on 1000 runs when n=BBS 
Estimated_Prec <- data.frame(PrecLoc(BBA,BBS),PrecLoc(Adj,BBS),PrecLoc(SampledBBA,BBS),
                             PrecLoc(MABBA,MABBS),PrecLoc(MA_BBA_Adj,MABBS),PrecLoc(SampledMABBA,MABBS),
                             PrecLoc(MIBBA,MIBBS),PrecLoc(MI_BBA_Adj,MIBBS),PrecLoc(SampledMIBBA,MIBBS),
                             PrecLoc(NYBBA,NYBBS),PrecLoc(NY_BBA_Adj,NYBBS),PrecLoc(SampledNYBBA,NYBBS),
                             PrecLoc(PABBA,PABBS),PrecLoc(PA_BBA_Adj,PABBS),PrecLoc(SampledPABBA,PABBS),
                             PrecLoc(VTBBA,VTBBS),PrecLoc(VT_BBA_Adj,VTBBS),PrecLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Prec) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Prec <- Estimated_Prec %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Precipitation = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Precipitation),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Prec<-Estimated_Prec %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Prec$Dataset <- as.character(paste(Estimated_Prec$Dataset, "(n=BBS)"))

#Combining with larger data
df_list9 <- list(Alpha.Prec, Estimated_Prec)

Alpha.Prec2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list9, accumulate=FALSE)

#Column order
setcolorder(Alpha.Prec2, c("Dataset","Region","Mean.Precipitation",
                           "Low.CI","Up.CI"))

#Arranging by region and dataset
Alpha.Prec2 <- Alpha.Prec2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Making sure labels are correct
colnames(Alpha.Prec2) <- c("Dataset","Region",
                           "Mean_Prec",
                           "Prec_LowCI",
                           "Prec_UpCI")



##### Now for Temp

Alpha.Temp <- data.frame(MeanCI(BBS$Temp,conf.level = 0.95,na.rm=T),
                         MeanCI(BBA$Temp,conf.level = 0.95),
                         MeanCI(Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledBBA$Temp,conf.level = 0.95),
                         MeanCI(SampledAdj$Temp,conf.level = 0.95),
                         MeanCI(MABBS$Temp,conf.level = 0.95),
                         MeanCI(MABBA$Temp,conf.level = 0.95),
                         MeanCI(MA_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledMABBA$Temp,conf.level=0.95),
                         MeanCI(SampledMA_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(MIBBS$Temp,conf.level = 0.95),
                         MeanCI(MIBBA$Temp,conf.level = 0.95),
                         MeanCI(MI_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledMIBBA$Temp,conf.level=0.95),
                         MeanCI(SampledMI_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(NYBBS$Temp,conf.level = 0.95),
                         MeanCI(NYBBA$Temp,conf.level = 0.95),
                         MeanCI(NY_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledNYBBA$Temp,conf.level=0.95),
                         MeanCI(SampledNY_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(PABBS$Temp,conf.level = 0.95),
                         MeanCI(PABBA$Temp,conf.level = 0.95),
                         MeanCI(PA_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledPABBA$Temp,conf.level=0.95),
                         MeanCI(SampledPA_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(VTBBS$Temp,conf.level = 0.95),
                         MeanCI(VTBBA$Temp,conf.level = 0.95),
                         MeanCI(VT_BBA_Adj$Temp,conf.level = 0.95),
                         MeanCI(SampledVTBBA$Temp,conf.level=0.95),
                         MeanCI(SampledVT_BBA_Adj$Temp,conf.level = 0.95)
                         )

#Column names
colnames(Alpha.Temp) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                         "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                         "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                         "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                         "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                         "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
                         )


#Creating dataframe 
Alpha.Temp <- Alpha.Temp %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Temp = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Temp),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.Temp<-Alpha.Temp%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.Temp <- Alpha.Temp %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean local temperature based on 1000 runs when n=BBS 
Estimated_Temp <- data.frame(TempLoc(BBA,BBS),TempLoc(Adj,BBS),TempLoc(SampledBBA,BBS),
                             TempLoc(MABBA,MABBS),TempLoc(MA_BBA_Adj,MABBS),TempLoc(SampledMABBA,MABBS),
                             TempLoc(MIBBA,MIBBS),TempLoc(MI_BBA_Adj,MIBBS),TempLoc(SampledMIBBA,MIBBS),
                             TempLoc(NYBBA,NYBBS),TempLoc(NY_BBA_Adj,NYBBS),TempLoc(SampledNYBBA,NYBBS),
                             TempLoc(PABBA,PABBS),TempLoc(PA_BBA_Adj,PABBS),TempLoc(SampledPABBA,PABBS),
                             TempLoc(VTBBA,VTBBS),TempLoc(VT_BBA_Adj,VTBBS),TempLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_Temp) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_Temp <- Estimated_Temp %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Temp = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Temp),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_Temp<-Estimated_Temp %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_Temp$Dataset <- as.character(paste(Estimated_Temp$Dataset, "(n=BBS)"))

#Combining with larger data
df_list10 <- list(Alpha.Temp, Estimated_Temp)

Alpha.Temp2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list10, accumulate=FALSE)

#Column order
setcolorder(Alpha.Temp2, c("Dataset","Region","Mean.Temp",
                           "Low.CI","Up.CI"))

#Arranging by region and dataset
Alpha.Temp2 <- Alpha.Temp2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Making sure column names are right
colnames(Alpha.Temp2) <- c("Dataset","Region",
                           "Mean_Temp",
                           "Temp_LowCI",
                           "Temp_UpCI")



##### Now for Elevation Range

Alpha.Elev_Range <- data.frame(MeanCI(BBS$Elev_Range,conf.level = 0.95,na.rm=T),
                         MeanCI(BBA$Elev_Range,conf.level = 0.95),
                         MeanCI(Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledBBA$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledAdj$Elev_Range,conf.level = 0.95),
                         MeanCI(MABBS$Elev_Range,conf.level = 0.95),
                         MeanCI(MABBA$Elev_Range,conf.level = 0.95),
                         MeanCI(MA_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledMABBA$Elev_Range,conf.level=0.95),
                         MeanCI(SampledMA_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(MIBBS$Elev_Range,conf.level = 0.95),
                         MeanCI(MIBBA$Elev_Range,conf.level = 0.95),
                         MeanCI(MI_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledMIBBA$Elev_Range,conf.level=0.95),
                         MeanCI(SampledMI_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(NYBBS$Elev_Range,conf.level = 0.95),
                         MeanCI(NYBBA$Elev_Range,conf.level = 0.95),
                         MeanCI(NY_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledNYBBA$Elev_Range,conf.level=0.95),
                         MeanCI(SampledNY_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(PABBS$Elev_Range,conf.level = 0.95),
                         MeanCI(PABBA$Elev_Range,conf.level = 0.95),
                         MeanCI(PA_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledPABBA$Elev_Range,conf.level=0.95),
                         MeanCI(SampledPA_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(VTBBS$Elev_Range,conf.level = 0.95),
                         MeanCI(VTBBA$Elev_Range,conf.level = 0.95),
                         MeanCI(VT_BBA_Adj$Elev_Range,conf.level = 0.95),
                         MeanCI(SampledVTBBA$Elev_Range,conf.level=0.95),
                         MeanCI(SampledVT_BBA_Adj$Elev_Range,conf.level = 0.95)
)

#Column names
colnames(Alpha.Elev_Range) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                         "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                         "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                         "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                         "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                         "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Alpha.Elev_Range <- Alpha.Elev_Range %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Elev_Range = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Elev_Range),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.Elev_Range<-Alpha.Elev_Range%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.Elev_Range <- Alpha.Elev_Range %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean elevation range in sampling unit based on 1000 runs when n=BBS 
Estimated_ElevRange <- data.frame(ElevRangeLoc(BBA,BBS),ElevRangeLoc(Adj,BBS),ElevRangeLoc(SampledBBA,BBS),
                             ElevRangeLoc(MABBA,MABBS),ElevRangeLoc(MA_BBA_Adj,MABBS),ElevRangeLoc(SampledMABBA,MABBS),
                             ElevRangeLoc(MIBBA,MIBBS),ElevRangeLoc(MI_BBA_Adj,MIBBS),ElevRangeLoc(SampledMIBBA,MIBBS),
                             ElevRangeLoc(NYBBA,NYBBS),ElevRangeLoc(NY_BBA_Adj,NYBBS),ElevRangeLoc(SampledNYBBA,NYBBS),
                             ElevRangeLoc(PABBA,PABBS),ElevRangeLoc(PA_BBA_Adj,PABBS),ElevRangeLoc(SampledPABBA,PABBS),
                             ElevRangeLoc(VTBBA,VTBBS),ElevRangeLoc(VT_BBA_Adj,VTBBS),ElevRangeLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_ElevRange) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                             "MA_BBA","MA_Adjacent","MA_SampledBBA",
                             "MI_BBA","MI_Adjacent","MI_SampledBBA",
                             "NY_BBA","NY_Adjacent","NY_SampledBBA",
                             "PA_BBA","PA_Adjacent","PA_SampledBBA",
                             "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_ElevRange <- Estimated_ElevRange %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.Elev_Range = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.Elev_Range),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_ElevRange<-Estimated_ElevRange %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_ElevRange$Dataset <- as.character(paste(Estimated_ElevRange$Dataset, "(n=BBS)"))

#Combining with larger data
df_list11 <- list(Alpha.Elev_Range, Estimated_ElevRange)

Alpha.Elev_Range2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list11, accumulate=FALSE)

#Column order
setcolorder(Alpha.Elev_Range2, c("Dataset","Region","Mean.Elev_Range",
                           "Low.CI","Up.CI"))

#Arranging by region and dataset
Alpha.Elev_Range2 <- Alpha.Elev_Range2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Column names
colnames(Alpha.Elev_Range2) <- c("Dataset","Region",
                           "Mean_ElevRange",
                           "ElevRange_LowCI",
                           "ElevRange_UpCI")


##### Now for Land-Class Richness

Alpha.LandRich <- data.frame(MeanCI(BBS$Het_Count2006,conf.level = 0.95,na.rm=T),
                               MeanCI(BBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledBBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledAdj$Het_Count2006,conf.level = 0.95),
                               MeanCI(MABBS$Het_Count2006,conf.level = 0.95),
                               MeanCI(MABBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(MA_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledMABBA$Het_Count2006,conf.level=0.95),
                               MeanCI(SampledMA_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(MIBBS$Het_Count2006,conf.level = 0.95),
                               MeanCI(MIBBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(MI_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledMIBBA$Het_Count2006,conf.level=0.95),
                               MeanCI(SampledMI_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(NYBBS$Het_Count2006,conf.level = 0.95),
                               MeanCI(NYBBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(NY_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledNYBBA$Het_Count2006,conf.level=0.95),
                               MeanCI(SampledNY_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(PABBS$Het_Count2006,conf.level = 0.95),
                               MeanCI(PABBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(PA_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledPABBA$Het_Count2006,conf.level=0.95),
                               MeanCI(SampledPA_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(VTBBS$Het_Count2006,conf.level = 0.95),
                               MeanCI(VTBBA$Het_Count2006,conf.level = 0.95),
                               MeanCI(VT_BBA_Adj$Het_Count2006,conf.level = 0.95),
                               MeanCI(SampledVTBBA$Het_Count2006,conf.level=0.95),
                               MeanCI(SampledVT_BBA_Adj$Het_Count2006,conf.level = 0.95)
)

#Column names
colnames(Alpha.LandRich) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                               "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                               "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                               "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                               "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                               "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Alpha.LandRich <- Alpha.LandRich %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.LandRich = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.LandRich),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.LandRich<-Alpha.LandRich%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.LandRich <- Alpha.LandRich %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean local land-class richness based on 1000 runs when n=BBS 
Estimated_LandRich <- data.frame(LandRichLoc(BBA,BBS),LandRichLoc(Adj,BBS),LandRichLoc(SampledBBA,BBS),
                                  LandRichLoc(MABBA,MABBS),LandRichLoc(MA_BBA_Adj,MABBS),LandRichLoc(SampledMABBA,MABBS),
                                  LandRichLoc(MIBBA,MIBBS),LandRichLoc(MI_BBA_Adj,MIBBS),LandRichLoc(SampledMIBBA,MIBBS),
                                  LandRichLoc(NYBBA,NYBBS),LandRichLoc(NY_BBA_Adj,NYBBS),LandRichLoc(SampledNYBBA,NYBBS),
                                  LandRichLoc(PABBA,PABBS),LandRichLoc(PA_BBA_Adj,PABBS),LandRichLoc(SampledPABBA,PABBS),
                                  LandRichLoc(VTBBA,VTBBS),LandRichLoc(VT_BBA_Adj,VTBBS),LandRichLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_LandRich) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                  "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                  "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                  "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                  "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                  "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_LandRich <- Estimated_LandRich %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.LandRich = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.LandRich),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_LandRich<-Estimated_LandRich %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_LandRich$Dataset <- as.character(paste(Estimated_LandRich$Dataset, "(n=BBS)"))

#Combining with larger data
df_list12 <- list(Alpha.LandRich, Estimated_LandRich)

Alpha.LandRich2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list12, accumulate=FALSE)

#Column order
setcolorder(Alpha.LandRich2, c("Dataset","Region","Mean.LandRich",
                                 "Low.CI","Up.CI"))

#Arranging by region and dataset
Alpha.LandRich2 <- Alpha.LandRich2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Column names
colnames(Alpha.LandRich2) <- c("Dataset","Region",
                           "Mean_LandRich",
                           "LandRich_LowCI",
                           "LandRich_UpCI")


##### Now for Land-Class Diversity

Alpha.LandDiv <- data.frame(MeanCI(BBS$Shannon2006,conf.level = 0.95,na.rm=T),
                             MeanCI(BBA$Shannon2006,conf.level = 0.95),
                             MeanCI(Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledBBA$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledAdj$Shannon2006,conf.level = 0.95),
                             MeanCI(MABBS$Shannon2006,conf.level = 0.95),
                             MeanCI(MABBA$Shannon2006,conf.level = 0.95),
                             MeanCI(MA_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledMABBA$Shannon2006,conf.level=0.95),
                             MeanCI(SampledMA_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(MIBBS$Shannon2006,conf.level = 0.95),
                             MeanCI(MIBBA$Shannon2006,conf.level = 0.95),
                             MeanCI(MI_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledMIBBA$Shannon2006,conf.level=0.95),
                             MeanCI(SampledMI_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(NYBBS$Shannon2006,conf.level = 0.95),
                             MeanCI(NYBBA$Shannon2006,conf.level = 0.95),
                             MeanCI(NY_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledNYBBA$Shannon2006,conf.level=0.95),
                             MeanCI(SampledNY_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(PABBS$Shannon2006,conf.level = 0.95),
                             MeanCI(PABBA$Shannon2006,conf.level = 0.95),
                             MeanCI(PA_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledPABBA$Shannon2006,conf.level=0.95),
                             MeanCI(SampledPA_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(VTBBS$Shannon2006,conf.level = 0.95),
                             MeanCI(VTBBA$Shannon2006,conf.level = 0.95),
                             MeanCI(VT_BBA_Adj$Shannon2006,conf.level = 0.95),
                             MeanCI(SampledVTBBA$Shannon2006,conf.level=0.95),
                             MeanCI(SampledVT_BBA_Adj$Shannon2006,conf.level = 0.95)
)

#Column names
colnames(Alpha.LandDiv) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                             "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                             "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                             "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                             "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                             "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Alpha.LandDiv <- Alpha.LandDiv %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.LandDiv = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.LandDiv),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Alpha.LandDiv<-Alpha.LandDiv%>%
  separate(Region, into = c("Region","Dataset"))

#Arranging
Alpha.LandDiv <- Alpha.LandDiv %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","Adjacent","SampledBBA","SampledAdjacent")))


### Adding when BBA (n=BBS)

#Estimate mean local land-class diversity (Shannon-Wiener) based on 1000 runs when n=BBS 
Estimated_LandDiv <- data.frame(LandDivLoc(BBA,BBS),LandDivLoc(Adj,BBS),LandDivLoc(SampledBBA,BBS),
                                 LandDivLoc(MABBA,MABBS),LandDivLoc(MA_BBA_Adj,MABBS),LandDivLoc(SampledMABBA,MABBS),
                                 LandDivLoc(MIBBA,MIBBS),LandDivLoc(MI_BBA_Adj,MIBBS),LandDivLoc(SampledMIBBA,MIBBS),
                                 LandDivLoc(NYBBA,NYBBS),LandDivLoc(NY_BBA_Adj,NYBBS),LandDivLoc(SampledNYBBA,NYBBS),
                                 LandDivLoc(PABBA,PABBS),LandDivLoc(PA_BBA_Adj,PABBS),LandDivLoc(SampledPABBA,PABBS),
                                 LandDivLoc(VTBBA,VTBBS),LandDivLoc(VT_BBA_Adj,VTBBS),LandDivLoc(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_LandDiv) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_LandDiv <- Estimated_LandDiv %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Mean.LandDiv = Value[c(NA,1,NA)],
         Low.CI = Value[c(NA,2,NA)],
         Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Mean.LandDiv),!is.na(Low.CI),!is.na(Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_LandDiv<-Estimated_LandDiv %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_LandDiv$Dataset <- as.character(paste(Estimated_LandDiv$Dataset, "(n=BBS)"))

#Combining with larger data
df_list13 <- list(Alpha.LandDiv, Estimated_LandDiv)

Alpha.LandDiv2 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list13, accumulate=FALSE)

#Column order
setcolorder(Alpha.LandDiv2, c("Dataset","Region","Mean.LandDiv",
                               "Low.CI","Up.CI"))

#Arranging by region and dataset
Alpha.LandDiv2 <- Alpha.LandDiv2 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))

#Column names
colnames(Alpha.LandDiv2) <- c("Dataset","Region",
                           "Mean_LandDiv",
                           "LandDiv_LowCI",
                           "LandDiv_UpCI")



#################################

#### ADDING ALL ALPHA ENVIRONMENTAL VARIABLES INTO ONE DATASET

df_list14 <- list(Alpha.Elev2,Alpha.Elev_Range2,
                  Alpha.Prec2,Alpha.Temp2,
                  Alpha.LandRich2,Alpha.LandDiv2)

Enviro.Alpha <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list14, accumulate=FALSE) 

#Arranging by region and dataset
Enviro.Alpha <- Enviro.Alpha %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


## Writing table

write_csv(Enviro.Alpha,paste0(path_files,"/Final Data/Alpha_Enviro.csv"))

# removing files from environment to save memory
rm(
  df_list8, df_list9, df_list10, df_list11, df_list12, df_list13, df_list14, #lists
  Estimated_Elev, Estimated_ElevRange, Estimated_Prec, # Estimates when n=BBS
  Estimated_Temp, Estimated_LandRich, Estimated_LandDiv, 
  Alpha.Elev, Alpha.Elev2, Alpha.Elev_Range, Alpha.Elev_Range2, # Observed means
  Alpha.Prec, Alpha.Prec2, Alpha.Temp, Alpha.Temp2, 
  Alpha.LandRich, Alpha.LandRich2, Alpha.LandDiv, Alpha.LandDiv2
  #Enviro.Alpha # final data
)


###############################################################






###################################################################################
##################### SIMILARITY TABLE OF ENVIRONMENTAL VARIABLES #################
###################################################################################


## Elevation

Beta.Elev <- data.frame(ElevB(BBS),ElevB(BBA), ElevB(Adj),
                        ElevB(SampledBBA),ElevB(SampledAdj),
                        ElevB(MABBS),ElevB(MABBA),ElevB(MA_BBA_Adj), #MA
                        ElevB(SampledMABBA),ElevB(SampledMA_BBA_Adj),
                        ElevB(MIBBS),ElevB(MIBBA),ElevB(MI_BBA_Adj), #MI
                        ElevB(SampledMIBBA),ElevB(SampledMI_BBA_Adj),
                        ElevB(NYBBS),ElevB(NYBBA),ElevB(NY_BBA_Adj), #NY
                        ElevB(SampledNYBBA),ElevB(SampledNY_BBA_Adj),
                        ElevB(PABBS),ElevB(PABBA),ElevB(PA_BBA_Adj), #PA
                        ElevB(SampledPABBA),ElevB(SampledPA_BBA_Adj),
                        ElevB(VTBBS),ElevB(VTBBA),ElevB(VT_BBA_Adj), #VT
                        ElevB(SampledVTBBA),ElevB(SampledVT_BBA_Adj)
                        )

#Column names
colnames(Beta.Elev) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                        "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                        "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                        "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                        "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                        "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Beta.Elev <- Beta.Elev %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Elev_Bray = Value[c(NA,1,NA)],
         Elev_Low.CI = Value[c(NA,2,NA)],
         Elev_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Elev_Bray),!is.na(Elev_Low.CI),!is.na(Elev_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Beta.Elev<-Beta.Elev%>%
  separate(Region, into = c("Region","Dataset"))


#### Adding BBA when n=BBS

Estimated_ElevB <- data.frame(ElevBeta(BBA,BBS),ElevBeta(Adj,BBS),ElevBeta(SampledBBA,BBS),
                                ElevBeta(MABBA,MABBS),ElevBeta(MA_BBA_Adj,MABBS),ElevBeta(SampledMABBA,MABBS),
                                ElevBeta(MIBBA,MIBBS),ElevBeta(MI_BBA_Adj,MIBBS),ElevBeta(SampledMIBBA,MIBBS),
                                ElevBeta(NYBBA,NYBBS),ElevBeta(NY_BBA_Adj,NYBBS),ElevBeta(SampledNYBBA,NYBBS),
                                ElevBeta(PABBA,PABBS),ElevBeta(PA_BBA_Adj,PABBS),ElevBeta(SampledPABBA,PABBS),
                                ElevBeta(VTBBA,VTBBS),ElevBeta(VT_BBA_Adj,VTBBS),ElevBeta(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_ElevB) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_ElevB <- Estimated_ElevB %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Elev_Bray = Value[c(NA,1,NA)],
         Elev_Low.CI = Value[c(NA,2,NA)],
         Elev_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Elev_Bray),!is.na(Elev_Low.CI),!is.na(Elev_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_ElevB<-Estimated_ElevB %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_ElevB$Dataset <- as.character(paste(Estimated_ElevB$Dataset, "(n=BBS)"))

#Combining
df_list15 <- list(Beta.Elev,Estimated_ElevB)

Beta.Elev1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list15, accumulate=FALSE)

#Arranging by region and dataset
Beta.Elev1 <- Beta.Elev1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


####### Precipitation


Beta.Prec <- data.frame(PrecB(BBS),PrecB(BBA), PrecB(Adj),
                        PrecB(SampledBBA),PrecB(SampledAdj),
                        PrecB(MABBS),PrecB(MABBA),PrecB(MA_BBA_Adj), #MA
                        PrecB(SampledMABBA),PrecB(SampledMA_BBA_Adj),
                        PrecB(MIBBS),PrecB(MIBBA),PrecB(MI_BBA_Adj), #MI
                        PrecB(SampledMIBBA),PrecB(SampledMI_BBA_Adj),
                        PrecB(NYBBS),PrecB(NYBBA),PrecB(NY_BBA_Adj), #NY
                        PrecB(SampledNYBBA),PrecB(SampledNY_BBA_Adj),
                        PrecB(PABBS),PrecB(PABBA),PrecB(PA_BBA_Adj), #PA
                        PrecB(SampledPABBA),PrecB(SampledPA_BBA_Adj),
                        PrecB(VTBBS),PrecB(VTBBA),PrecB(VT_BBA_Adj), #VT
                        PrecB(SampledVTBBA),PrecB(SampledVT_BBA_Adj)
)

#Column names
colnames(Beta.Prec) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                        "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                        "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                        "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                        "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                        "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Beta.Prec <- Beta.Prec %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Prec_Bray = Value[c(NA,1,NA)],
         Prec_Low.CI = Value[c(NA,2,NA)],
         Prec_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Prec_Bray),!is.na(Prec_Low.CI),!is.na(Prec_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Beta.Prec<-Beta.Prec%>%
  separate(Region, into = c("Region","Dataset"))



#### Adding BBA when n=BBS

Estimated_PrecB <- data.frame(PrecBeta(BBA,BBS),PrecBeta(Adj,BBS),PrecBeta(SampledBBA,BBS),
                              PrecBeta(MABBA,MABBS),PrecBeta(MA_BBA_Adj,MABBS),PrecBeta(SampledMABBA,MABBS),
                              PrecBeta(MIBBA,MIBBS),PrecBeta(MI_BBA_Adj,MIBBS),PrecBeta(SampledMIBBA,MIBBS),
                              PrecBeta(NYBBA,NYBBS),PrecBeta(NY_BBA_Adj,NYBBS),PrecBeta(SampledNYBBA,NYBBS),
                              PrecBeta(PABBA,PABBS),PrecBeta(PA_BBA_Adj,PABBS),PrecBeta(SampledPABBA,PABBS),
                              PrecBeta(VTBBA,VTBBS),PrecBeta(VT_BBA_Adj,VTBBS),PrecBeta(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_PrecB) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                              "MA_BBA","MA_Adjacent","MA_SampledBBA",
                              "MI_BBA","MI_Adjacent","MI_SampledBBA",
                              "NY_BBA","NY_Adjacent","NY_SampledBBA",
                              "PA_BBA","PA_Adjacent","PA_SampledBBA",
                              "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_PrecB <- Estimated_PrecB %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Prec_Bray = Value[c(NA,1,NA)],
         Prec_Low.CI = Value[c(NA,2,NA)],
         Prec_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Prec_Bray),!is.na(Prec_Low.CI),!is.na(Prec_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_PrecB<-Estimated_PrecB %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_PrecB$Dataset <- as.character(paste(Estimated_PrecB$Dataset, "(n=BBS)"))


#Combining
df_list16 <- list(Beta.Prec,Estimated_PrecB)

Beta.Prec1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list16, accumulate=FALSE)

#Arranging by region and dataset
Beta.Prec1 <- Beta.Prec1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


####### Temperature

Beta.Temp <- data.frame(TempB(BBS),TempB(BBA), TempB(Adj),
                        TempB(SampledBBA),TempB(SampledAdj),
                        TempB(MABBS),TempB(MABBA),TempB(MA_BBA_Adj), #MA
                        TempB(SampledMABBA),TempB(SampledMA_BBA_Adj),
                        TempB(MIBBS),TempB(MIBBA),TempB(MI_BBA_Adj), #MI
                        TempB(SampledMIBBA),TempB(SampledMI_BBA_Adj),
                        TempB(NYBBS),TempB(NYBBA),TempB(NY_BBA_Adj), #NY
                        TempB(SampledNYBBA),TempB(SampledNY_BBA_Adj),
                        TempB(PABBS),TempB(PABBA),TempB(PA_BBA_Adj), #PA
                        TempB(SampledPABBA),TempB(SampledPA_BBA_Adj),
                        TempB(VTBBS),TempB(VTBBA),TempB(VT_BBA_Adj), #VT
                        TempB(SampledVTBBA),TempB(SampledVT_BBA_Adj)
)

#Column names
colnames(Beta.Temp) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                        "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                        "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                        "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                        "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                        "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Beta.Temp <- Beta.Temp %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Temp_Bray = Value[c(NA,1,NA)],
         Temp_Low.CI = Value[c(NA,2,NA)],
         Temp_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Temp_Bray),!is.na(Temp_Low.CI),!is.na(Temp_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Beta.Temp<-Beta.Temp%>%
  separate(Region, into = c("Region","Dataset"))



#### Adding BBA when n=BBS

Estimated_TempB <- data.frame(TempBeta(BBA,BBS),TempBeta(Adj,BBS),TempBeta(SampledBBA,BBS),
                              TempBeta(MABBA,MABBS),TempBeta(MA_BBA_Adj,MABBS),TempBeta(SampledMABBA,MABBS),
                              TempBeta(MIBBA,MIBBS),TempBeta(MI_BBA_Adj,MIBBS),TempBeta(SampledMIBBA,MIBBS),
                              TempBeta(NYBBA,NYBBS),TempBeta(NY_BBA_Adj,NYBBS),TempBeta(SampledNYBBA,NYBBS),
                              TempBeta(PABBA,PABBS),TempBeta(PA_BBA_Adj,PABBS),TempBeta(SampledPABBA,PABBS),
                              TempBeta(VTBBA,VTBBS),TempBeta(VT_BBA_Adj,VTBBS),TempBeta(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_TempB) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                              "MA_BBA","MA_Adjacent","MA_SampledBBA",
                              "MI_BBA","MI_Adjacent","MI_SampledBBA",
                              "NY_BBA","NY_Adjacent","NY_SampledBBA",
                              "PA_BBA","PA_Adjacent","PA_SampledBBA",
                              "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_TempB <- Estimated_TempB %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Temp_Bray = Value[c(NA,1,NA)],
         Temp_Low.CI = Value[c(NA,2,NA)],
         Temp_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Temp_Bray),!is.na(Temp_Low.CI),!is.na(Temp_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_TempB<-Estimated_TempB %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_TempB$Dataset <- as.character(paste(Estimated_TempB$Dataset, "(n=BBS)"))


#Combining
df_list17 <- list(Beta.Temp,Estimated_TempB)

Beta.Temp1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list17, accumulate=FALSE)

#Arranging by region and dataset
Beta.Temp1 <- Beta.Temp1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


####### Land Class Diversity

Beta.LandDiv <- data.frame(LandDivB(BBS),LandDivB(BBA), LandDivB(Adj),
                        LandDivB(SampledBBA),LandDivB(SampledAdj),
                        LandDivB(MABBS),LandDivB(MABBA),LandDivB(MA_BBA_Adj), #MA
                        LandDivB(SampledMABBA),LandDivB(SampledMA_BBA_Adj),
                        LandDivB(MIBBS),LandDivB(MIBBA),LandDivB(MI_BBA_Adj), #MI
                        LandDivB(SampledMIBBA),LandDivB(SampledMI_BBA_Adj),
                        LandDivB(NYBBS),LandDivB(NYBBA),LandDivB(NY_BBA_Adj), #NY
                        LandDivB(SampledNYBBA),LandDivB(SampledNY_BBA_Adj),
                        LandDivB(PABBS),LandDivB(PABBA),LandDivB(PA_BBA_Adj), #PA
                        LandDivB(SampledPABBA),LandDivB(SampledPA_BBA_Adj),
                        LandDivB(VTBBS),LandDivB(VTBBA),LandDivB(VT_BBA_Adj), #VT
                        LandDivB(SampledVTBBA),LandDivB(SampledVT_BBA_Adj)
)

#Column names
colnames(Beta.LandDiv) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                        "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                        "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                        "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                        "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                        "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Beta.LandDiv <- Beta.LandDiv %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(LandDiv_Bray = Value[c(NA,1,NA)],
         LandDiv_Low.CI = Value[c(NA,2,NA)],
         LandDiv_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(LandDiv_Bray),!is.na(LandDiv_Low.CI),!is.na(LandDiv_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Beta.LandDiv<-Beta.LandDiv%>%
  separate(Region, into = c("Region","Dataset"))



#### Adding BBA when n=BBS

Estimated_LandDivB <- data.frame(LandDivBeta(BBA,BBS),LandDivBeta(Adj,BBS),LandDivBeta(SampledBBA,BBS),
                              LandDivBeta(MABBA,MABBS),LandDivBeta(MA_BBA_Adj,MABBS),LandDivBeta(SampledMABBA,MABBS),
                              LandDivBeta(MIBBA,MIBBS),LandDivBeta(MI_BBA_Adj,MIBBS),LandDivBeta(SampledMIBBA,MIBBS),
                              LandDivBeta(NYBBA,NYBBS),LandDivBeta(NY_BBA_Adj,NYBBS),LandDivBeta(SampledNYBBA,NYBBS),
                              LandDivBeta(PABBA,PABBS),LandDivBeta(PA_BBA_Adj,PABBS),LandDivBeta(SampledPABBA,PABBS),
                              LandDivBeta(VTBBA,VTBBS),LandDivBeta(VT_BBA_Adj,VTBBS),LandDivBeta(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_LandDivB) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                              "MA_BBA","MA_Adjacent","MA_SampledBBA",
                              "MI_BBA","MI_Adjacent","MI_SampledBBA",
                              "NY_BBA","NY_Adjacent","NY_SampledBBA",
                              "PA_BBA","PA_Adjacent","PA_SampledBBA",
                              "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_LandDivB <- Estimated_LandDivB %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(LandDiv_Bray = Value[c(NA,1,NA)],
         LandDiv_Low.CI = Value[c(NA,2,NA)],
         LandDiv_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(LandDiv_Bray),!is.na(LandDiv_Low.CI),!is.na(LandDiv_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_LandDivB<-Estimated_LandDivB %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_LandDivB$Dataset <- as.character(paste(Estimated_LandDivB$Dataset, "(n=BBS)"))


#Combining
df_list18 <- list(Beta.LandDiv,Estimated_LandDivB)

Beta.LandDiv1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list18, accumulate=FALSE)

#arranging by region and dataset
Beta.LandDiv1 <- Beta.LandDiv1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


####### Gower multivariate


Beta.Gower <- data.frame(GowerB(BBS),GowerB(BBA), GowerB(Adj),
                           GowerB(SampledBBA),GowerB(SampledAdj),
                           GowerB(MABBS),GowerB(MABBA),GowerB(MA_BBA_Adj), #MA
                           GowerB(SampledMABBA),GowerB(SampledMA_BBA_Adj),
                           GowerB(MIBBS),GowerB(MIBBA),GowerB(MI_BBA_Adj), #MI
                           GowerB(SampledMIBBA),GowerB(SampledMI_BBA_Adj),
                           GowerB(NYBBS),GowerB(NYBBA),GowerB(NY_BBA_Adj), #NY
                           GowerB(SampledNYBBA),GowerB(SampledNY_BBA_Adj),
                           GowerB(PABBS),GowerB(PABBA),GowerB(PA_BBA_Adj), #PA
                           GowerB(SampledPABBA),GowerB(SampledPA_BBA_Adj),
                           GowerB(VTBBS),GowerB(VTBBA),GowerB(VT_BBA_Adj), #VT
                           GowerB(SampledVTBBA),GowerB(SampledVT_BBA_Adj)
)

#Column names
colnames(Beta.Gower) <-c("All_BBS","All_BBA","All_Adjacent","All_SampledBBA", "All_SampledAdjacent",
                           "MA_BBS","MA_BBA","MA_Adjacent","MA_SampledBBA", "MA_SampledAdjacent",
                           "MI_BBS","MI_BBA","MI_Adjacent","MI_SampledBBA", "MI_SampledAdjacent",
                           "NY_BBS","NY_BBA","NY_Adjacent","NY_SampledBBA", "NY_SampledAdjacent",
                           "PA_BBS","PA_BBA","PA_Adjacent","PA_SampledBBA", "PA_SampledAdjacent",
                           "VT_BBS","VT_BBA","VT_Adjacent","VT_SampledBBA", "VT_SampledAdjacent"
)


#Creating dataframe 
Beta.Gower <- Beta.Gower %>%
  pivot_longer(c("All_BBS":"VT_SampledAdjacent"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gower_Bray = Value[c(NA,1,NA)],
         Gower_Low.CI = Value[c(NA,2,NA)],
         Gower_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gower_Bray),!is.na(Gower_Low.CI),!is.na(Gower_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Beta.Gower<-Beta.Gower%>%
  separate(Region, into = c("Region","Dataset"))



#### Adding BBA when n=BBS

Estimated_GowerB <- data.frame(GowerBeta(BBA,BBS),GowerBeta(Adj,BBS),GowerBeta(SampledBBA,BBS),
                                 GowerBeta(MABBA,MABBS),GowerBeta(MA_BBA_Adj,MABBS),GowerBeta(SampledMABBA,MABBS),
                                 GowerBeta(MIBBA,MIBBS),GowerBeta(MI_BBA_Adj,MIBBS),GowerBeta(SampledMIBBA,MIBBS),
                                 GowerBeta(NYBBA,NYBBS),GowerBeta(NY_BBA_Adj,NYBBS),GowerBeta(SampledNYBBA,NYBBS),
                                 GowerBeta(PABBA,PABBS),GowerBeta(PA_BBA_Adj,PABBS),GowerBeta(SampledPABBA,PABBS),
                                 GowerBeta(VTBBA,VTBBS),GowerBeta(VT_BBA_Adj,VTBBS),GowerBeta(SampledVTBBA,VTBBS))


#Labelling columns
colnames(Estimated_GowerB) <-c("All_BBA","All_Adjacent","All_SampledBBA",
                                 "MA_BBA","MA_Adjacent","MA_SampledBBA",
                                 "MI_BBA","MI_Adjacent","MI_SampledBBA",
                                 "NY_BBA","NY_Adjacent","NY_SampledBBA",
                                 "PA_BBA","PA_Adjacent","PA_SampledBBA",
                                 "VT_BBA","VT_Adjacent","VT_SampledBBA")


#Creating dataframe for BBA (n=BBS)
Estimated_GowerB <- Estimated_GowerB %>%
  pivot_longer(c("All_BBA":"VT_SampledBBA"),
               names_to = "Region",
               values_to = "Value") %>%
  group_by(Region) %>%
  mutate(Gower_Bray = Value[c(NA,1,NA)],
         Gower_Low.CI = Value[c(NA,2,NA)],
         Gower_Up.CI = Value[c(NA,3,NA)]) %>%
  filter(!is.na(Gower_Bray),!is.na(Gower_Low.CI),!is.na(Gower_Up.CI)) %>%
  dplyr::select(-Value)

#Separating Region and Dataset Column
Estimated_GowerB<-Estimated_GowerB %>%
  separate(Region, into = c("Region","Dataset"))

#Labelling
Estimated_GowerB$Dataset <- as.character(paste(Estimated_GowerB$Dataset, "(n=BBS)"))


#Combining
df_list19 <- list(Beta.Gower,Estimated_GowerB)

Beta.Gower1 <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list19, accumulate=FALSE)

#arranging by region and dataset
Beta.Gower1 <- Beta.Gower1 %>%
  arrange(Region,
          factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                   "Adjacent","Adjacent (n=BBS)",
                                   "SampledBBA","SampledBBA (n=BBS)",
                                   "SampledAdjacent")))


####Putting it all together

df_list20 <- list(Beta.Elev1,Beta.Prec1,Beta.Temp1,Beta.LandDiv1,Beta.Gower1)

Enviro.Beta <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list20, accumulate=FALSE)

#Arranging by region and dataset
Enviro.Beta <- Enviro.Beta %>%
  arrange(Region,
        factor(Dataset, levels=c("BBA","BBS","BBA (n=BBS)",
                                 "Adjacent","Adjacent (n=BBS)",
                                 "SampledBBA","SampledBBA (n=BBS)",
                                 "SampledAdjacent")))

#Writing file

write_csv(Enviro.Beta,paste0(path_files,"/Final Data/Beta_Enviro.csv"))

#The final tables for species and environmental tables are joined in excel
#to then be called in the script "Datasets.R" to load for analyses and figures 

##############################################################################
##############################################################################


##############################################################################
####################### REGRESSIONS FOR EFFORT ON RICHNESS ################### 
##############################################################################


#Running regression on effect of hours of sampling effort on local richness
#Running this for thresholds from 0, 10, then by 10 
#until we reach a cutoff of a slope of less than 0.1
 
#Will then save the R2, slope, intercept, sample size, 
#and mean hours and richness per regression

#This summary table will be used to then decide the threshold of 
#"sufficient" sampling in BBA blocks. Will then use it to make figures.


#Regression for no cutoff

#Removing outliers
BBA <- BBA %>%
  filter(Atlas2_Effort < 300)

#Running regression
fit1 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA)

#Summarizing
p1 <- summary(fit1)

#intercept
int1 <- p1$coefficients[1]
#slope
slo1 <- p1$coefficients[2]
#R2
r2.1 <- p1$r.squared
#Sample size
n1 <- nrow(BBA)
#Average hours
eff1 <- round(mean(BBA$Atlas2_Effort),digits=2)
#Average local richness
rich1 <- round(mean(BBA$Atlas2_Total),digits=2)

#-------------------------------------

#Regression for cutoff of 10 hours
BBA1 <- BBA %>%
  filter(Atlas2_Effort>9)

#Running regression
fit2 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA1)

#Summarizing
p2 <- summary(fit2)

#intercept
int2 <- p2$coefficients[1]
#slope
slo2 <- p2$coefficients[2]
#R2
r2.2 <- p2$r.squared
#Sample size
n2 <- nrow(BBA1)
#Average hours
eff2 <- round(mean(BBA1$Atlas2_Effort),digits=2)
#Average local richness
rich2 <- round(mean(BBA1$Atlas2_Total),digits=2)


#-------------------------------------

#Regression for cutoff of 20 hours
BBA2 <- BBA %>%
  filter(Atlas2_Effort>19)

#Running regression
fit3 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA2)

#Summarizing
p3 <- summary(fit3)

#intercept
int3 <- p3$coefficients[1]
#slope
slo3 <- p3$coefficients[2]
#R2
r2.3 <- p3$r.squared
#Sample size
n3 <- nrow(BBA2)
#Average hours
eff3 <- round(mean(BBA2$Atlas2_Effort),digits=2)
#Average local richness
rich3 <- round(mean(BBA2$Atlas2_Total),digits=2)


#-------------------------------------

#Regression for cutoff of 30 hours
BBA3 <- BBA %>%
  filter(Atlas2_Effort>29)

#Running regression
fit4 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA3)

#Summarizing
p4 <- summary(fit4)

#intercept
int4 <- p4$coefficients[1]
#slope
slo4 <- p4$coefficients[2]
#R2
r2.4 <- p4$r.squared
#Sample size
n4 <- nrow(BBA3)
#Average hours
eff4 <- round(mean(BBA3$Atlas2_Effort),digits=2)
#Average local richness
rich4 <- round(mean(BBA3$Atlas2_Total),digits=2)


#-------------------------------------

#Regression for cutoff of 40 hours
BBA4 <- BBA %>%
  filter(Atlas2_Effort>39)

#Running regression
fit5 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA4)

#Summarizing
p5 <- summary(fit5)

#intercept
int5 <- p5$coefficients[1]
#slope
slo5 <- p5$coefficients[2]
#R2
r2.5 <- p5$r.squared
#Sample size
n5 <- nrow(BBA4)
#Average hours
eff5 <- round(mean(BBA4$Atlas2_Effort),digits=2)
#Average local richness
rich5 <- round(mean(BBA4$Atlas2_Total),digits=2)

#-------------------------------------

#Regression for cutoff of 50 hours
BBA5 <- BBA %>%
  filter(Atlas2_Effort>49)

#Running regression
fit6 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA5)

#Summarizing
p6 <- summary(fit6)

#intercept
int6 <- p6$coefficients[1]
#slope
slo6 <- p6$coefficients[2]
#R2
r2.6 <- p6$r.squared
#Sample size
n6 <- nrow(BBA5)
#Average hours
eff6 <- round(mean(BBA5$Atlas2_Effort),digits=2)
#Average local richness
rich6 <- round(mean(BBA5$Atlas2_Total),digits=2)


#-------------------------------------

#Regression for cutoff of 60 hours
BBA6 <- BBA %>%
  filter(Atlas2_Effort>59)

#Running regression
fit7 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA6)

#Summarizing
p7 <- summary(fit7)

#intercept
int7 <- p7$coefficients[1]
#slope
slo7 <- p7$coefficients[2]
#R2
r2.7 <- p7$r.squared
#Sample size
n7 <- nrow(BBA6)
#Average hours
eff7 <- round(mean(BBA6$Atlas2_Effort),digits=2)
#Average local richness
rich7 <- round(mean(BBA6$Atlas2_Total),digits=2)


#-------------------------------------

#Regression for cutoff of 70 hours
BBA7 <- BBA %>%
  filter(Atlas2_Effort>69)

#Running regression
fit8 <- lm(Atlas2_Total~Atlas2_Effort,
           data=BBA7)

#Summarizing
p8 <- summary(fit8)

#intercept
int8 <- p8$coefficients[1]
#slope
slo8 <- p8$coefficients[2]
#R2
r2.8 <- p8$r.squared
#Sample size
n8 <- nrow(BBA7)
#Average hours
eff8 <- round(mean(BBA7$Atlas2_Effort),digits=2)
#Average local richness
rich8 <- round(mean(BBA7$Atlas2_Total),digits=2)



#-------------------------------------

#Putting this all together into a table

#Making columns

Threshold<- c(0,10,20,30,40,50,60,70)

Sample_Size <- c(n1,n2,n3,n4,n5,n6,n7,n8)

Average_Effort <- c(eff1,eff2,eff3,eff4,eff5,eff6,eff7,eff8)

Average_Richness <- c(rich1,rich2,rich3,rich4,rich5,rich6,rich7,rich8)

R2 <- c(r2.1,r2.2,r2.3,r2.4,r2.5,r2.6,r2.7,r2.8)

Slope <- c(slo1,slo2,slo3,slo4,slo5,slo6,slo7,slo8)

Intercept <- c(int1,int2,int3,int4,int5,int6,int7,int8)

#Putting all together into a dataframe
Effort <- data.frame(Threshold,
                     Sample_Size,
                     Average_Effort,
                     Average_Richness,
                     R2,
                     Slope,
                     Intercept)

#Saving dataframe
write_csv(Effort,paste0(path_files,"/Final Data/Effort.csv"))


##############################################################################
############################## END OF SCRIPT ################################# 
##############################################################################

  