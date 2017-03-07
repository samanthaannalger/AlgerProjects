# Samantha A. Alger
# P. Alexander Burnham
# 07, March 2017
# qPCR Results for Plants:
# 2015 bombus survey plants
# plant transmission exp plants
# Data manipulation 

#------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
  setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# read in data:
plants2015 <- read.csv("Plants_2015Survey.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 
plantqPCR <- read.csv("PlantqPCRResults.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 

# check data:
head(plants2015)
head(plantqPCR)

library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

## Clean up for 2015 plant data___________________________________________

# take only columns that we want:
plantqPCR <- select(plantqPCR, labID, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)


# remove unwanted rows from dataframe:
plantqPCR<-plantqPCR[!(plantqPCR$labID=="No Sample"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="Gblock"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="NTC"),]

# merge 2015 eco data with qPCR results
plants2015 <- merge(plantqPCR, plants2015, by="labID")

# create a binary vector in dataframe that pulled out samples where both replicates are positive (samples that are positive will receive a 1):
plants2015$BINYprefilter <- ifelse(plants2015$quantity_sd > 0, 1, 0)

#select columns needed and make new dataframe
plants2015 <- select(plants2015, labID, ID, target_name, site, spp, dateCollected, sampleType, BINYprefilter)

#delete unique rows to show only one entry for each sample. 
plants2015 <- unique(plants2015)

write.csv(plants2015, file = "plants2015DF.csv")


##___________________________________________________________________
#Plant Transmission data

#Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# read in data:
plantTransPlants <- read.csv("plantTransPlants.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
plantqPCR <- read.csv("PlantqPCRResults.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 

# check data:
head(plantTransPlants)
head(plantqPCR)

# data cleaning whip dataframe into shape and merge with eco data:
# remove uneeded columns from DF:
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

## Clean up for 2015 plant data___________________________________________

# take only columns that we want:
plantqPCR <- select(plantqPCR, labID, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)


# remove unwanted rows from dataframe:
plantqPCR<-plantqPCR[!(plantqPCR$labID=="No Sample"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="Gblock"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="NTC"),]


# take only columns that we want:
plantqPCR <- select(plantqPCR, labID, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)


# remove unwanted rows from dataframe:
plantqPCR<-plantqPCR[!(plantqPCR$labID=="No Sample"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="Gblock"),]
plantqPCR <- plantqPCR[!(plantqPCR$labID=="NTC"),]

# merge 2015 eco data with qPCR results
plantTransPlants <- merge(plantqPCR, plantTransPlants, by="labID")

head(plantTransPlants)

# create a binary vector in dataframe that pulled out samples where both replicates are positive (samples that are positive will receive a 1):
plantTransPlants$BINYprefilter <- ifelse(plantTransPlants$quantity_sd > 0, 1, 0)

#select columns needed and make new dataframe
plantTransPlants <- select(plantTransPlants, labID, ID, target_name, date, tent, group, fieldID, spp, exp, BINYprefilter)

#delete unique rows to show only one entry for each sample. 
plantTransPlants <- unique(plantTransPlants)

write.csv(plantTransPlants, file = "plantTransPlantsDF.csv")
