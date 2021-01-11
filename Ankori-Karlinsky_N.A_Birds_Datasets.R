##################################################################
#################### DATASETS FOR BIRD PROJECT ###################
##################################################################

# Roi Ankori-Karlinsky, Columbia University Spring 2020

# This script: loads datasets for analyses and figures for bird Dataset project
# 1. First loads dataset for entire region
# 2. Then loads datasets for each state
# 3. Then creates tidy data for certain analyses
# 4. Then loads datasets for figures

# This can be sourced into analysis and figure scripts for quick loading 

###################################################################

# Loading packages
library(tidyverse) #tidyverse for data cleaning and manipulation
library(data.table) # data.table package for data manipulation
library(readxl) #read excel


################### DATASETS FOR ENTIRE REGION ####################

#Creating file path
path_files<- 'C:/Users/roiak/Documents/Ronen/Final'

#Loading data calling it full_d as it contains all the variables
# (this was done by saving the Master_Data tab in the excel spreadsheet as a separate csv)
full_d <- read_csv(paste0(path_files,"/Final Data/Final_Data.csv"),
                   col_types = cols(.default = "n", ID="c",Dataset="c",Block_ID="c",Route_ID="c",
                                    State="c",BBS_Adjacent="c",BCR="c",BBSID="c"))

#Getting rid of blocks for which we have no effort data or no birds observed
full_d <- filter(full_d,!is.na(Atlas2_Effort),Atlas2_Total>0)

#creating log richness and log effort columns
full_d <- full_d %>%
  mutate(logRichness = log1p(Atlas2_Total),
         logEffort = log1p(Atlas2_Effort))

#Making sure dataset is a factor
full_d$Dataset <- as.factor(full_d$Dataset)

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

#Datsaet with pairs of geographically-matched sampling units of BBS and BBA
Matched <- read_csv(paste0(path_files,"/Final Data/Matched_Data.csv"),
                    col_types = cols(.default = "n", ID="c",Dataset="c",Block_ID="c",Route_ID="c",
                                     State="c",BBS_Adjacent="c",BCR="c",BBSID="c"))


#Creating Dataset with pairs of geographically-matching sampling units of BBS and BBA (based on largest overlap)

# #Keeping code of creating this dataset here for reference
# # #Filter out BBS-Adjacent BBA blocks by biggest area overlap
# Matched <- Adjacent_Full %>%
#    group_by(Dataset,State,BBSID) %>%
#    filter(BBS_Area==max(BBS_Area))
# # 
# # #check only BBA
#  tally(~Matched$Dataset)
# # 
# # #List with BBS
#  df_list <- list(Matched,BBS)
# # 
# # #Merge dataframes
#  Matched_Data <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list, accumulate=FALSE)
# # 
# # #Arrange by BBSID so that data is organized in matched sampling unit pairs
# Matched_Data <- Matched_Data %>%
#    arrange(BBSID,Dataset,State)
# # 
# # #checking that only 281 BBS routes
#  nrow(distinct(Matched_Data,BBSID)) #yes
# # 
# # #write to file
# write_csv(Matched_Data,paste0(path_files,"/Data/Matched_Data.csv"))


#Creating dataset of geographically-matched sampling units of BBS and BBA,
#where the local richness in all BBS-adjecent BBA blocks is averaged per one BBS route.
Matched_Ave <- Adj %>%
  group_by(State,BBSID) %>%
  summarize(BBA_Richness = mean(Atlas2_Total,na.rm=T),
            BBA_Richness_sd = sd(Atlas2_Total,na.rm=T),
            BBA_Effort = mean(Atlas2_Effort,na.rm=T),
            BBA_Effort_sd = sd(Atlas2_Effort,na.rm=T),
            number_blocks = n()) %>%
  mutate(BBA_Richness_CI = (qnorm(0.975)*BBA_Richness_sd)/sqrt(number_blocks))


#merging with BBS
df_list <- list(Matched_Ave,BBS)

Matched_Ave <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list, accumulate=FALSE)

#Cleaning
Matched_Ave <- Matched_Ave %>%
  dplyr::select(State,BBSID,
                Atlas2_Effort,Atlas2_Total,BBA_Richness,BBA_Richness_sd,BBA_Richness_CI,
                BBA_Effort,BBA_Effort_sd,number_blocks) %>%
  rename(BBS_Richness = Atlas2_Total,
         BBS_Effort = Atlas2_Effort) %>%
  arrange(State,BBSID)

rm(df_list)


#Creating dataset of well sampled (<30hr) geographically-matched sampling units of BBS and BBA,
#where the local richness in all BBS-adjecent BBA blocks is averaged per one BBS route. 
Sampled_Matched_Ave <- SampledAdj %>%
  group_by(State,BBSID) %>%
  summarize(BBA_Richness = mean(Atlas2_Total,na.rm=T),
            BBA_Richness_sd = sd(Atlas2_Total,na.rm=T),
            BBA_Effort = mean(Atlas2_Effort,na.rm=T),
            BBA_Effort_sd = sd(Atlas2_Effort,na.rm=T),
            number_blocks = n()) %>%
  mutate(BBA_Richness_CI = (qnorm(0.975)*BBA_Richness_sd)/sqrt(number_blocks))

#merging with BBS
df_list <- list(Sampled_Matched_Ave,BBS)

Sampled_Matched_Ave <- Reduce(function(x,y) merge(x,y,all=TRUE),df_list, accumulate=FALSE)

#Cleaning
Sampled_Matched_Ave <- Sampled_Matched_Ave %>%
  dplyr::select(State,BBSID,
                Atlas2_Effort,Atlas2_Total,BBA_Richness,BBA_Richness_sd,BBA_Richness_CI,
                BBA_Effort,BBA_Effort_sd,number_blocks) %>%
  rename(BBS_Richness = Atlas2_Total,
         BBS_Effort = Atlas2_Effort) %>%
  arrange(State,BBSID)

rm(df_list)

####################### DATASETS FOR MA ###########################

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



####################### DATASETS FOR MI ###########################

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



####################### DATASETS FOR NY ###########################

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



####################### DATASETS FOR PA ###########################

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



####################### DATASETS FOR VT ###########################

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



####################### CREATING TIDY DATA ########################


#### Creating tidy data for regional species richness calculations

#Doing this using tidyr package with pivot_longer() function.


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


####################### LOADING FIGURES DATA ######################


## These datasheets are created in the script "Summary Tables" ####

# Therefore, when running "Summary Tables.R", comment this out,
# then after creating those datatables and saving them, run this
# before doing anaylses/figures

#Loading data with regression on effect of sampling effort
Effort <- read_excel(paste0(path_files,"/Final Data/Comparing North American Bird Databases.xlsx"),
                       sheet = "Effort Regression")

#Loading data with rareness analysis
rareness <- read_excel(paste0(path_files,
                              "/Final Data/Comparing North American Bird Databases.xlsx"),
                       sheet = "Rareness Analysis",
                       col_types = c("text","text","numeric","numeric","numeric","numeric"))


#Loading figure data for assumptions (i.e. environmental variables)
Assumption_Figures <- read_excel(paste0(path_files,
                                        "/Final Data/Comparing North American Bird Databases.xlsx"),
                                 sheet = "Summary Table Environment",
                                 col_types = c("text", "text", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric",
                                               "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric",
                                               "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric",
                                               "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric",
                                               "numeric", "numeric", 
                                               "numeric", "numeric", "numeric"))


#Loading figure data for predictions (i.e. bird variables)
Prediction_Figures <- read_excel(paste0(path_files,
                                        "/Final Data/Comparing North American Bird Databases.xlsx"),
                                 sheet = "Summary Table Birds", 
                                 col_types = c("text", "text", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric"))

