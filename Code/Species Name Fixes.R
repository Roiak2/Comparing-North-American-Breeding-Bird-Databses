##################################################################
################## CURTIS FIXES FOR BIRD PROJECT #################
##################################################################

# Roi Ankori-Karlinsky, Columbia University Spring 2020

# This script: uses Curtis Flatcher's insights and
#              fixes naming errors and mismatches between BBA and BBS
# 1. First loads dataset for entire region
# 2. Resolves naming convention mismatches by lumping together same species
# 3. Removes unnecessary names (often ' errors)
# 4. Lumps together hybrids if they are only delineated in one of the datasets
# 5. Outputs fixed dataset to be used for analyses

###################################################################

# Loading packages
library(tidyverse) #tidyverse for data cleaning and manipulation
library(data.table) # data.table package for data manipulation



########################## LOADING DATA ###########################

#file path
path_files<- 'C:/Users/roiak/Documents/Ronen/Final' # Path files for easy loading and writing

#Loading data calling it full_d as it contains all the variables
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




######################### CURTIS REMARKS ##########################

#--LISTING ACTIONS NEEDED--#

#   1) Lump Alder and Willow flycatcher into Traill's Flycatcher to match BBS
#       a) Create new column called Traill's Flycatcher
#       b) Move Traills Flycatcher (Willow/Alder Flycatcher) from BBS 
#          and both Alder and Willow Flycatchers from BBA into new column
#       c) Run a check to make sure transferred correctly
#       d) Delete old columns to make sure no double counting occurs
#
#   2) Fix Brewers and Brewer's error by lumping them into the same column and delete old ones
#
#   3) Remove Brewster's and Lawrence's Warblers as they are hybrids
#   
#   4) Fix Coopers and Cooper's Hawk error by lumping them into the same column and delete old ones
#
#   5) Lump all Dark-eyed Junco regional types from BBS into BBA column that has just a single Dark-eyed Junco
#
#   6) Fix Henslows and Henslow's error
#
#   7) Fix Kirtlands and Kirtland's error
#
#   8) Fix Lecontes and LeConte's error
#
#   9) Remove Mallard x Black Duck Hybrid from BBA 
#
#   10) Lump Northern Flicker (Yellow-shafted Flicker) into just Northern Flicker
#
#   11) Lump Saltmarsh Sparrow with Saltmarsh Sharp-tailed Sparrow
#
#   12) Fix Swainsons and Swainson's error
#
#   13) Fix Wilsons and Wilson's Snipe error
#
#   14) Fix Wilsons and Wilson's Warbler error
#
#   15) Lump Yellow-rumped Warbler with Myrtle Warbler
#
#   16) Fix Lincolns and Lincoln's Sparrow error

#==================================================================

######################### FIXING DATA #############################

#### 1) TRAILL'S FLYCATCHER FIXES ####
#   a) Testing this out on a dummy dataset with only the relevant columns
#   b) After seeing it works applying to entire dataset

#Subsetting columns and getting a random 100 rows for a dummy dataset
#Then moving all values into new column (if any detection, count it), using dplyr

full_d %>%
  dplyr::select(Dataset, #select variables subset
                `Traills Flycatcher (Willow/Alder Flycatcher)`,
                `Willow Flycatcher`,
                `Alder Flycatcher`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Traill's Flycatcher`= ifelse(
    `Traills Flycatcher (Willow/Alder Flycatcher)`==1, #if this name detected,
    `Traills Flycatcher (Willow/Alder Flycatcher)`, #count it
    ifelse(
      `Willow Flycatcher`==1, #if this name detected,
      `Willow Flycatcher`, #count it
      ifelse(
        `Alder Flycatcher`==1, #if this name detected,
        `Alder Flycatcher`, #count it
        0 #otherwise count as absence
      )
    )
  )
  ) %>%
  View() #check

#Seems to work!

#Now checking that it works while removing the old columns
full_d %>%
  dplyr::select(Dataset, #select variables subset
                `Traills Flycatcher (Willow/Alder Flycatcher)`,
                `Willow Flycatcher`,
                `Alder Flycatcher`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Traill's Flycatcher`= ifelse(
    `Traills Flycatcher (Willow/Alder Flycatcher)`==1, #if this name detected,
    `Traills Flycatcher (Willow/Alder Flycatcher)`, #count it
    ifelse(
      `Willow Flycatcher`==1, #if this name detected,
      `Willow Flycatcher`, #count it
      ifelse(
        `Alder Flycatcher`==1, #if this name detected,
        `Alder Flycatcher`, #count it
        0 #otherwise count as absence
      )
    )
  )
  ) %>%
  dplyr::select(-`Traills Flycatcher (Willow/Alder Flycatcher)`,
                -`Willow Flycatcher`,
                -`Alder Flycatcher`) %>% #removing old columns
  View() #check

#Nice!

#Now doing it for the real dataset and saving it as real_data
real_data <- full_d %>%
  mutate(`Traill's Flycatcher`= ifelse(
    `Traills Flycatcher (Willow/Alder Flycatcher)`==1, #if this name detected,
    `Traills Flycatcher (Willow/Alder Flycatcher)`, #count it
    ifelse(
      `Willow Flycatcher`==1, #if this name detected,
      `Willow Flycatcher`, #count it
      ifelse(
        `Alder Flycatcher`==1, #if this name detected,
        `Alder Flycatcher`, #count it
        0 #otherwise count as absence
      )
    )
  )
  ) %>%
  dplyr::select(-`Traills Flycatcher (Willow/Alder Flycatcher)`,
                -`Willow Flycatcher`,
                -`Alder Flycatcher`) #removing old columns
#checking
View(real_data$`Traill's Flycatcher`)
mosaic::tally(~real_data$`Traill's Flycatcher`) #only 1s and 0s, which is great!

#==================================================================

#### 2) Fix Brewers and Brewer's error ####
#   a) dummy set subset check that lumping works
#   b) doing it on real data and removing old columns

full_d %>%
  #filter(Dataset=="BBA")%>%
  dplyr::select(Dataset, #select variables subset
                `Brewer's Blackbird`,
                `Brewers Blackbird`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Brewer's Blackbird`= ifelse(
    `Brewer's Blackbird`==1, #if this name detected,
    `Brewer's Blackbird`, #count it
    ifelse(
      `Brewers Blackbird`==1, #if this name detected,
      `Brewers Blackbird`, #count it
        0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check


#Works!

#Now doing it on real dataset
real_data <- real_data %>%
  mutate(`Brewer's Blackbird`= ifelse(
    `Brewer's Blackbird`==1, #if this name detected,
    `Brewer's Blackbird`, #count it
    ifelse(
      `Brewers Blackbird`==1, #if this name detected,
      `Brewers Blackbird`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(- `Brewers Blackbird`) #remove old column

#Checking
mosaic::tally(~real_data$`Brewer's Blackbird`) #232 detections
#comparing to previous one
mosaic::tally(~full_d$`Brewer's Blackbird`) #221
mosaic::tally(~full_d$`Brewers Blackbird`) #11 (adds up to 232 so success!)

#==================================================================

#### 3) Remove Brewster's and Lawrence's Warblers ####

real_data <- real_data %>%
  dplyr::select(-`Brewster's Warbler`,
                -`Lawrence's Warbler`)

#==================================================================

#### 4) Fix Coopers and Cooper's Hawk error ####
#   a) testing it on dummy set first
#   b) then doing it on real data


full_d %>%
  #filter(Dataset=="BBA")%>%
  dplyr::select(Dataset, #select variables subset
                `Cooper's Hawk`,
                `Coopers Hawk`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Cooper's Hawk Check`= ifelse(
    `Cooper's Hawk`==1, #if this name detected,
    `Cooper's Hawk`, #count it
    ifelse(
      `Coopers Hawk`==1, #if this name detected,
      `Coopers Hawk`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now on real dataset

real_data <- real_data %>%
  mutate(`Cooper's Hawk`= ifelse(
    `Cooper's Hawk`==1, #if this name detected,
    `Cooper's Hawk`, #count it
    ifelse(
      `Coopers Hawk`==1, #if this name detected,
      `Coopers Hawk`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Coopers Hawk`) # removing old column

#checking
mosaic::tally(~real_data$`Cooper's Hawk`) #4235 detections
mosaic::tally(~full_d$`Coopers Hawk`) #88
mosaic::tally(~full_d$`Cooper's Hawk`) #4147 (adds up to 4235 success!)


#==================================================================

#### 5) Lump all Dark-eyed Junco regional types from BBS into BBA column that has just a single Dark-eyed Junco ####
#   a) first testing on dummy dataset
#   b) then checking and doing on real dataset

real_data %>%
  #filter(Dataset=="BBS") %>%
  dplyr::select(Dataset, #select variables subset
                `Dark-eyed Junco`,
                `Dark-eyed Junco (Slate-colored Junco)`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Dark-eyed Junco Check`= ifelse(
    `Dark-eyed Junco`==1, #if this name detected,
    `Dark-eyed Junco`, #count it
    ifelse(
      `Dark-eyed Junco (Slate-colored Junco)`==1, #if this name detected,
      `Dark-eyed Junco (Slate-colored Junco)`, #count it
        0 #otherwise count as absence
      )
    )
  ) %>%
  View() #check

#Works!

#Now doing it on real data
real_data <- real_data %>%
  mutate(`Dark-eyed Junco`= ifelse(
    `Dark-eyed Junco`==1, #if this name detected,
    `Dark-eyed Junco`, #count it
    ifelse(
      `Dark-eyed Junco (Slate-colored Junco)`==1, #if this name detected,
      `Dark-eyed Junco (Slate-colored Junco)`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Dark-eyed Junco (Slate-colored Junco)`) #removing old column

#checking
mosaic::tally(~real_data$`Dark-eyed Junco`) #5617 detections
mosaic::tally(~full_d$`Dark-eyed Junco`) #4997
mosaic::tally(~full_d$`Dark-eyed Junco (Slate-colored Junco)`) #620 (Adds up correctly, success!)


#==================================================================

#### 6) Fix Henslows and Henslow's error ####
#   a) First testing on dummy set
#   b) checking and doing on real data

real_data %>%
  #filter(Dataset=="BBS") %>%
  dplyr::select(Dataset, #select variables subset
                `Henslow's Sparrow`,
                `Henslows Sparrow`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Henslow's Sparrow Check`= ifelse(
    `Henslow's Sparrow`==1, #if this name detected,
    `Henslow's Sparrow`, #count it
    ifelse(
      `Henslows Sparrow`==1, #if this name detected,
      `Henslows Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now doing on real data

real_data <- real_data %>%
  mutate(`Henslow's Sparrow`= ifelse(
    `Henslow's Sparrow`==1, #if this name detected,
    `Henslow's Sparrow`, #count it
    ifelse(
      `Henslows Sparrow`==1, #if this name detected,
      `Henslows Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Henslows Sparrow`) #removing old column

#checking
mosaic::tally(~real_data$`Henslow's Sparrow`) #485 detections
mosaic::tally(~full_d$`Henslow's Sparrow`) #472 detections
mosaic::tally(~full_d$`Henslows Sparrow`) #13 detections (adds up success!)

#==================================================================

#### 7) Fix Kirtlands and Kirtland's error ####
#   a) dummy test
#   b) then real data

real_data %>%
  #filter(Dataset=="BBA") %>%
  dplyr::select(Dataset, #select variables subset
                `Kirtland's Warbler`,
                `Kirtlands Warbler`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Kirtland's Warbler Check`= ifelse(
    `Kirtland's Warbler`==1, #if this name detected,
    `Kirtland's Warbler`, #count it
    ifelse(
      `Kirtlands Warbler`==1, #if this name detected,
      `Kirtlands Warbler`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now doing it on real data

real_data <- real_data %>%
  mutate(`Kirtland's Warbler`= ifelse(
    `Kirtland's Warbler`==1, #if this name detected,
    `Kirtland's Warbler`, #count it
    ifelse(
      `Kirtlands Warbler`==1, #if this name detected,
      `Kirtlands Warbler`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Kirtlands Warbler`)

#checking
mosaic::tally(~real_data$`Kirtland's Warbler`) #127 detections
mosaic::tally(~full_d$`Kirtland's Warbler`) #124 detections
mosaic::tally(~full_d$`Kirtlands Warbler`) #3 detections (adds up success!)


#==================================================================

#### 8) Fix Lecontes and LeConte's Sparrow error ####
#   a) dummy test
#   b) real data

real_data %>%
  #filter(`Le Contes Sparrow`==1) %>%
  dplyr::select(Dataset, #select variables subset
                `LeConte's Sparrow`,
                `Le Contes Sparrow`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`LeConte's Sparrow Check`= ifelse(
    `LeConte's Sparrow`==1, #if this name detected,
    `LeConte's Sparrow`, #count it
    ifelse(
      `Le Contes Sparrow`==1, #if this name detected,
      `Le Contes Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now doing it on real data

real_data <- real_data %>%
  mutate(`LeConte's Sparrow`= ifelse(
    `LeConte's Sparrow`==1, #if this name detected,
    `LeConte's Sparrow`, #count it
    ifelse(
      `Le Contes Sparrow`==1, #if this name detected,
      `Le Contes Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Le Contes Sparrow`)

#checking
mosaic::tally(~real_data$`LeConte's Sparrow`) #66 detections
mosaic::tally(~full_d$`LeConte's Sparrow`) #65 detections
mosaic::tally(~full_d$`Le Contes Sparrow`) #1 detection (adds up success!)

#==================================================================

#### 9) Remove Mallard x Black Duck Hybrid from BBA ####

real_data <- real_data %>%
  dplyr::select(-`Mallard x Black Duck Hybrid`)



#==================================================================

#### 10) Lump Northern Flicker (Yellow-shafted Flicker) into just Northern Flicker ####
#   a) checking on dummy data
#   b) doing on real data

real_data %>%
  #filter(Dataset=="BBS") %>%
  dplyr::select(Dataset, #select variables subset
                `Northern Flicker`,
                `Northern Flicker (Yellow-shafted Flicker)`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Northern Flicker Check`= ifelse(
    `Northern Flicker`==1, #if this name detected,
    `Northern Flicker`, #count it
    ifelse(
      `Northern Flicker (Yellow-shafted Flicker)`==1, #if this name detected,
      `Northern Flicker (Yellow-shafted Flicker)`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now doing it on real data

real_data <- real_data %>%
  mutate(`Northern Flicker`= ifelse(
    `Northern Flicker`==1, #if this name detected,
    `Northern Flicker`, #count it
    ifelse(
      `Northern Flicker (Yellow-shafted Flicker)`==1, #if this name detected,
      `Northern Flicker (Yellow-shafted Flicker)`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Northern Flicker (Yellow-shafted Flicker)`)

#checking
mosaic::tally(~real_data$`Northern Flicker`) #14614 detections
mosaic::tally(~full_d$`Northern Flicker`) #13079 detections
mosaic::tally(~full_d$`Northern Flicker (Yellow-shafted Flicker)`) #1535 detection (adds up success!)

#==================================================================

#### 11) Lump Saltmarsh Sparrow with Saltmarsh Sharp-tailed Sparrow ####
#   a) dummy data
#   b) real data

real_data %>%
  #filter(Dataset=="BBS") %>%
  dplyr::select(Dataset, #select variables subset
                `Saltmarsh Sparrow`,
                `Saltmarsh Sharp-tailed Sparrow`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Saltmarsh Sparrow Check`= ifelse(
    `Saltmarsh Sparrow`==1, #if this name detected,
    `Saltmarsh Sparrow`, #count it
    ifelse(
      `Saltmarsh Sharp-tailed Sparrow`==1, #if this name detected,
      `Saltmarsh Sharp-tailed Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#works!

#now doing it on real data

real_data <- real_data %>%
  mutate(`Saltmarsh Sparrow`= ifelse(
    `Saltmarsh Sparrow`==1, #if this name detected,
    `Saltmarsh Sparrow`, #count it
    ifelse(
      `Saltmarsh Sharp-tailed Sparrow`==1, #if this name detected,
      `Saltmarsh Sharp-tailed Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  dplyr::select(-`Saltmarsh Sharp-tailed Sparrow`)

#checking
mosaic::tally(~real_data$`Saltmarsh Sparrow`) #139 detections
mosaic::tally(~full_d$`Saltmarsh Sparrow`) #60 detections
mosaic::tally(~full_d$`Saltmarsh Sharp-tailed Sparrow`) #79 detection (adds up success!)


#==================================================================

#### 12) Fix Swainsons and Swainson's Thrush error ####
#   a) dummy dataset
#   b) real data

real_data %>%
  #filter(Dataset=="BBS") %>%
  select(Dataset,
         `Swainson's Thrush`,
         `Swainsons Thrush`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Swainson's Thrush Check`= ifelse(
    `Swainson's Thrush`==1, #if this name detected,
    `Swainson's Thrush`, #count it
    ifelse(
      `Swainsons Thrush`==1, #if this name detected,
      `Swainsons Thrush`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#Works!

#Now doing it on real data

real_data <- real_data %>%
  mutate(`Swainson's Thrush`= ifelse(
    `Swainson's Thrush`==1, #if this name detected,
    `Swainson's Thrush`, #count it
    ifelse(
      `Swainsons Thrush`==1, #if this name detected,
      `Swainsons Thrush`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  select(-`Swainsons Thrush`) #removing old column


#Checking
mosaic::tally(~real_data$`Swainson's Thrush`) #1223
mosaic::tally(~full_d$`Swainson's Thrush`) #1190
mosaic::tally(~full_d$`Swainsons Thrush`) #33 (adds up success!)

#==================================================================

#### 13) Fix Wilsons and Wilson's Snipe error ####
#   a) dummy datset
#   b) real dataset

real_data %>%
 # filter(Dataset=="BBS") %>%
  select(Dataset,
         `Wilson's Snipe`,
         `Wilsons Snipe`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Wilson's Snipe Check`= ifelse(
    `Wilson's Snipe`==1, #if this name detected,
    `Wilson's Snipe`, #count it
    ifelse(
      `Wilsons Snipe`==1, #if this name detected,
      `Wilsons Snipe`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check

#Works!

#Now doing it on real data

real_data <- real_data %>%
  mutate(`Wilson's Snipe`= ifelse(
    `Wilson's Snipe`==1, #if this name detected,
    `Wilson's Snipe`, #count it
    ifelse(
      `Wilsons Snipe`==1, #if this name detected,
      `Wilsons Snipe`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  select(-`Wilsons Snipe`) #removing old column


#Checking
mosaic::tally(~real_data$`Wilson's Snipe`) #1619
mosaic::tally(~full_d$`Wilson's Snipe`) #1578
mosaic::tally(~full_d$`Wilsons Snipe`) #41 (adds up success!)

#==================================================================

#### 14) Fix Wilsons and Wilson's Warbler error ####
#   a) dummy dataset
#   b) real data

real_data %>%
  #filter(`Wilsons Warbler`==1) %>%
  select(Dataset,
         `Wilson's Warbler`,
         `Wilsons Warbler`) %>%
  #.[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Wilson's Warbler Check`= ifelse(
    `Wilson's Warbler`==1, #if this name detected,
    `Wilson's Warbler`, #count it
    ifelse(
      `Wilsons Warbler`==1, #if this name detected,
      `Wilsons Warbler`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check


#Works!

#Now doing it on real data

real_data <- real_data %>%
  mutate(`Wilson's Warbler`= ifelse(
    `Wilson's Warbler`==1, #if this name detected,
    `Wilson's Warbler`, #count it
    ifelse(
      `Wilsons Warbler`==1, #if this name detected,
      `Wilsons Warbler`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  select(-`Wilsons Warbler`) #removing old column


#Checking
mosaic::tally(~real_data$`Wilson's Warbler`) #47
mosaic::tally(~full_d$`Wilson's Warbler`) #46
mosaic::tally(~full_d$`Wilsons Warbler`) #1 (adds up success!)

#==================================================================

#### 15) Lump Yellow-rumped Warbler with Myrtle Warbler ####
#   a) dummy dataset
#   b) real data

real_data %>%
  #filter(Dataset=="BBS") %>%
  select(Dataset,
         `Yellow-rumped Warbler`,
         `Yellow-rumped Warbler (Myrtle Warbler)`) %>%
  .[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Yellow-rumped Warbler Check`= ifelse(
    `Yellow-rumped Warbler`==1, #if this name detected,
    `Yellow-rumped Warbler`, #count it
    ifelse(
      `Yellow-rumped Warbler (Myrtle Warbler)`==1, #if this name detected,
      `Yellow-rumped Warbler (Myrtle Warbler)`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check


#Works!

#Now doing it on real data

real_data <- real_data %>%
  mutate(`Yellow-rumped Warbler`= ifelse(
    `Yellow-rumped Warbler`==1, #if this name detected,
    `Yellow-rumped Warbler`, #count it
    ifelse(
      `Yellow-rumped Warbler (Myrtle Warbler)`==1, #if this name detected,
      `Yellow-rumped Warbler (Myrtle Warbler)`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  select(-`Yellow-rumped Warbler (Myrtle Warbler)`) #removing old column


#Checking
mosaic::tally(~real_data$`Yellow-rumped Warbler`) #5755
mosaic::tally(~full_d$`Yellow-rumped Warbler`) #5636
mosaic::tally(~full_d$`Yellow-rumped Warbler (Myrtle Warbler)`) #119 (adds up success!)

#==================================================================


#### 16) Fix Lincolns and Lincoln's Sparrow error ####
#   a) dummy dataset
#   b) real data

real_data %>%
  #filter(`Wilsons Warbler`==1) %>%
  select(Dataset,
         `Lincolns Sparrow`,
         `Lincoln's Sparrow`) %>%
  #.[sample(nrow(.),100),] %>% #randomly sample 100 rows
  mutate(`Lincoln's Sparrow Check`= ifelse(
    `Lincoln's Sparrow`==1, #if this name detected,
    `Lincoln's Sparrow`, #count it
    ifelse(
      `Lincolns Sparrow`==1, #if this name detected,
      `Lincolns Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  View() #check


#Works!

#Now doing it on real data

real_data <- real_data %>%
  mutate(`Lincoln's Sparrow`= ifelse(
    `Lincoln's Sparrow`==1, #if this name detected,
    `Lincoln's Sparrow`, #count it
    ifelse(
      `Lincolns Sparrow`==1, #if this name detected,
      `Lincolns Sparrow`, #count it
      0 #otherwise count as absence
    )
  )
  ) %>%
  select(-`Lincolns Sparrow`) #removing old column


#Checking
mosaic::tally(~real_data$`Lincoln's Sparrow`) #587
mosaic::tally(~full_d$`Lincoln's Sparrow`) #572
mosaic::tally(~full_d$`Lincolns Sparrow`) #15 (adds up success!)


######################### WRITING TABLE ###########################

#Saving csv file with finalized clean data

write_csv(real_data,paste0(path_files,"/Data/Final_Data.csv"))

#Now saving a version without filtering for unobserved sampling units
write_csv(real_data,paste0(path_files,"/Data/Full_Correct.csv"))

