# Analyses for ApiaryDat Manuscript
# Alger & Burnham
# 7/29/18

# Clear memory of characters:
rm(list=ls())

# source my packages
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(data.table)
library(ape)
library(spdep)
library("MuMIn")
library("multcomp")

# Set working directory:
setwd("~/Documents/GitHub/AlgerProjects/VTApiaries")

#upload data
RegData <- read.csv("CSV_files/RegActiveAndDelinquent.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

FullApiaryDat <- read.csv("CSV_files/VTApiaries.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

##############################################################
# DATA PREP ###################################################
##############################################################


# Data cleaning to make two dfs- one for apiary/colony level data, one for beekeeper level data:

#Classifying beekeepers by beek type
histDat = data.table(RegData)
histDat[, `n` := .N, by = BeekeeperID]
histDat$Beektype <- ifelse(histDat$n == 1,"Hobbyist", ifelse(histDat$n <=5, "Sideliner", "Commercial"))

#Merging the registration with census data to make a full dataframe:

#select only columns we need:
histDat <- dplyr::select(histDat, LocationID, Beektype, BeeKeeperStatus, n)
FullApiaryDat <- dplyr::select(FullApiaryDat, -BeeKeeperStatus)


#Create df for apiary/colony level analyses:
Apiarydf <- merge.data.frame(FullApiaryDat,histDat, by = "LocationID", all.y = TRUE)

# Create df for beekeepers level analyses only:
Beekdf <- Apiarydf[!duplicated(Apiarydf$BeekeeperID), ]


################################## ############## ############## ############## ###############
#End Data Prep   ################ ############## ############## ############## ###############
################################ ############## ############## ############## ###############

#Analyses:

# 1. Density/Distribution of beekeepers, colonies, and apiaries… Which counties have the greatest?  ANOVA
      # for beekeepers, use Beekdf and'County' column
      # for apiaries and colonies, use Apiarydf, For colonies: column     ColonyCount and County column. For apiaries: and length of Apiary df 

# Colony Count: 
Apiarydf$CountyName <- as.factor(Apiarydf$CountyName)

mod <- glm(data = Apiarydf, ColonyCount ~ CountyName, family = poisson)
Anova(mod)
#summary(glht(mod, mcp(CountyName="Tukey")))

# Apiary Count
mod1 <- glm(data = Beekdf, n ~ CountyName, family = poisson)
Anova(mod1)
#summary(glht(mod1, mcp(CountyName="Tukey")))

# Beekeeper Count:



# 2. Distribution of colony losses-  Is there spatial clustering? If so, where are the greatest colony losses? (it appears that losses are higher in counties east of the Green Mountains) Moran's I/ANOVA
    # Use apiary df, and use PerTotLoss, which is the % colony loss for that apiary. and 'County', Lat and Long (long might be spelled wrong- longtitude)


# clean data frame by selcting needed columns 
moransDat <- select(Apiarydf, CountyName, MiteCounts, PerTotLoss, Latitude, Longtitude, Beektype)

# remove NAs
moransDat <- moransDat[complete.cases(moransDat),]

# create distnace matrix
#For DWV:
colLoss.dists <- as.matrix(dist(cbind(moransDat$Longtitude, moransDat$Latitude)))
colLoss.dists.inv <- 1/colLoss.dists
diag(colLoss.dists.inv) <- 0

# remove infinity 
colLoss.dists.inv[colLoss.dists.inv=="Inf"] <- 0


# run Moran/s I
Moran.I(moransDat$PerTotLoss, colLoss.dists.inv)

# 3. Distribution of mite monitoring- Where are highest proportions of mite monitors? (it appears that mite monitoring is more common west of the green mountains) Moran's I/ANOVA. 
# Use True/False data from 'MiteCounts' True means the beekeeper monitors mite levles

# run Moran's I
Moran.I(moransDat$MiteCounts, colLoss.dists.inv)










# 5. What factors best explain colony losses? 

 # variables to test: location (county), number of kinds of mite treatments, whether they use synthetic/organic treatments, if they monitor their mite loads, beekeeper type.. etc.? GLM

######## Some basic stats/data cleaning needed prior to GLM

# For each beekeeper, classify beekeepers as organic mitacides vs. synthetic mitacides vs. nothing, (use beekeeper ID to merge with apiary level data) 
# This is going to be complicated and probably involve 'grep' I could try and do it when I'm finished with my drafts. I have code that might be able to help with this from ApiaryDat.R

# For each beekeeper, compute the number of different mite treatments they report using (0,1,2,3...), (use beekeeper ID to merge with apiary level data)
# Also complicated and need 'grep'
  
