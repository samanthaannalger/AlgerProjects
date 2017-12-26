###########################################################################################
# Data Analysis for Imid Experiment
# Samantha Alger
# December 11, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/ImidPilot/")

###########################################################################################
# Read in Virus Data:
ImidVirus <- read.table("csv_files/ImidqPCR_RawResults.csv", 
                       header=TRUE, 
                       sep = ",", 
                       stringsAsFactors = FALSE) 


# Read in Lab Data
ImidStat <- read.table("csv_files/ImidLabData.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 


###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# source my functions
source("Scripts/BurnhamFunctions.R")
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)

#--------------------------------------------------------------------------
# data cleaning

# preliminary data cleaning
ImidVirus <- PrelimClean(data=ImidVirus) 

# take only columns that we want:
data <- select(ImidVirus, sample_name, Treatment, dil.factor, ID)

# use dilution factors to calculate normalized virus load
ImidVirus <- VirusNorm(data=ImidVirus, number_bees = 1)

# use actin to normalize normalized viral load
ImidVirus <- actinNormal(data=ImidVirus)

# remove actin from data frame
ImidVirus <- ImidVirus[!(ImidVirus$target_name=="ACTIN"),]

# adds virus binary data and makes norm genome copy 0 if above threashold CT
ImidVirus <- CT_Threash(data=ImidVirus)

#Fixes Virus Binary column
ImidVirus$virusBINY <- ifelse(ImidVirus$NormGenomeCopy > 0, 1, 0)

# create full ID to check for inccorrect duplicates 
ImidVirus$fullID <- with(ImidVirus, paste0(sample_name, target_name))

# remove duplicates for examples with Cts above threashold (i.e. both coerced to 0s) 
ImidVirus <- ImidVirus[!duplicated(ImidVirus$fullID), ] 

# merge the virus data to the main dataset
ImidVirus <- VirusMerger2000(data1 = ImidVirus, data2 = ImidStat)

# log transform virus data:
ImidVirus$logDWV <- log(ImidVirus$DWVload + 1)
ImidVirus$logBQCV <- log(ImidVirus$BQCVload + 1)


##############################################################
#figure for DWV load:
#Select only DWV positive samples
ImidDWV <- ImidVirus[ which(ImidVirus$DWVbinary=="1"), ]

#Checking out by plant species
Imid <- ddply(ImidDWV, c("Treatment"), summarise, 
                  n = length(logDWV),
                  mean = mean(logDWV, na.rm=TRUE),
                  sd = sd(logDWV, na.rm=TRUE),
                  se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(Imid, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 20))


##############################################################
#figure for DWV load:
#Select only DWV positive samples
ImidBQCV <- ImidVirus[ which(ImidVirus$BQCVbinary=="1"), 

#Checking out by plant species
Imid <- ddply(ImidBQCV, c("Treatment"), summarise, 
              n = length(logBQCV),
              mean = mean(logBQCV, na.rm=TRUE),
              sd = sd(logBQCV, na.rm=TRUE),
              se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(Imid, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 15))

###############################################################
#Figure for Virus Prevalence

#Checking out by plant species
DWVPrev <- ddply(ImidVirus, c("Treatment"), summarise, 
                  n = length(DWVbinary),
                  mean = mean(DWVbinary, na.rm=TRUE),
                  sd = sd(DWVbinary, na.rm=TRUE),
                  se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(DWVPrev, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position= "none") + coord_cartesian(ylim = c(0, .25)) + scale_y_continuous(labels = scales::percent)
