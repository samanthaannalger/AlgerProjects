# P. Alexander Burnham
# Samantha A. Alger
# 11, October 2016
# qPCR Results
# Data manipulation 

#------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerCollaborations/2015_Bombus_Survey/CSV_Files")

# read in data:
BombSurvDF <- read.table("qPCR_results_dataframe.csv", header=TRUE, sep = ",") 
EcoSurvDF <- read.table("2015EcoSurveyData.csv", header=TRUE, sep = ",")
HbDF <- read.table("qPCR_results_dataframe_HB.csv", header=TRUE, sep =",")
EcoHbDF <- read.table("2015EcoSurveyData_HB.csv", header=TRUE, sep = ",")
NHBS <- read.table("NHBS_2015_DataSubset1.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
Diversity <- read.table("Diversity.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

length(EcoSurvDF$ID)
# check dataframe to see what the fuck is wrong with it:
head(BombSurvDF)
head(EcoSurvDF)
head(HbDF)
head(EcoHbDF)

#------------------------------------------------------------------------
# data cleaning whip dataframe into shape and merge with eco data:
# remove uneeded columns from DF:
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

# take only columns that we want:
BombSurvDF <- select(BombSurvDF, sample_name, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)


# remove duplicate rows from dataframe:
BombSurvDF<-BombSurvDF[!duplicated(BombSurvDF), ]


# remove NTC rows from dataframe:
BombSurvDF<-BombSurvDF[!(BombSurvDF$sample_name=="No Sample"),]

#create table to check for non-numeric sample names
table(BombSurvDF$sample_name)

# use grep to remove rows containing non Bombus Survey data:
#Field Experiment
BombSurvDF <- BombSurvDF[-grep("\\F",BombSurvDF$sample_name),] 
#Koppert
BombSurvDF <- BombSurvDF[-grep("\\K",BombSurvDF$sample_name),]
#Ross Conrad
BombSurvDF <- BombSurvDF[-grep("\\RC",BombSurvDF$sample_name),]
#BB- BioBest Colonies
BombSurvDF <- BombSurvDF[-grep("\\-",BombSurvDF$sample_name),] 
#gamma irradiated pollen
BombSurvDF <- BombSurvDF[-grep("\\Pollen",BombSurvDF$sample_name),]
#No template control
BombSurvDF <- BombSurvDF[-grep("\\NTC",BombSurvDF$sample_name),]
#Gblock
BombSurvDF <- BombSurvDF[-grep("\\Gblock",BombSurvDF$sample_name),]

# make sample name a numeric variable:
BombSurvDF$sample_name <- as.numeric(as.character(BombSurvDF$sample_name))

# order from low to high my sample name:
BombSurvDF <- BombSurvDF[order(BombSurvDF$sample_name),]

EcoSurvDF <- EcoSurvDF[order(EcoSurvDF$sample_name),]

# Merge datasets eco and bombus qPCR:
#Need rownames and all.x=TRUE because data frames are different sizes.
BombSurvDF$variable <- rownames(BombSurvDF)
BombSurv <- merge(BombSurvDF, EcoSurvDF, all.x =TRUE)

# remove temporary column - "variable" from DF:
BombSurv <- select(BombSurv, -variable)

# set constant values for genome copies per bee calculations:
crude_extr <- 100
eluteRNA <- 50
total_extr_vol <- 600
number_bees <- 1
cDNA_eff <- 0.1
rxn_vol <- 3

# create column for genome copies per bee:
BombSurv$genome_copbee <- ((((((BombSurv$quantity_mean / cDNA_eff) / rxn_vol) * BombSurv$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees

# line turns NAs into 0s:
#BombSurv$genome_copbee[is.na(BombSurv$genome_copbee)] <- 0


# pull only actin values out of dataframe
ActinOnly <- BombSurv[which(BombSurv$target_name=="ACTIN"),]

# create DF of ACTIN genome copies and lab ID:
ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$genome_copbee)
colnames(ActinDF) <- c("sample_name", "ACT_genome_copbee")

# merge ACTIN dataframe with main dataframe:
#Need rownames and all.x=TRUE because data frames are different sizes.

BombSurv <- merge(BombSurv, ActinDF, by="sample_name")

# find mean of all ACTIN values:
ActinMean <- mean(ActinOnly$genome_copbee, na.rm = TRUE)

# create column for normalized genome copies per bee:
BombSurv$norm_genome_copbee <- (BombSurv$genome_copbee/BombSurv$ACT_genome_copbee)*ActinMean


#----------------------------------------------------------------------
#Working with the HB data:
# take only columns that we want:
HbDF <- select(HbDF, sample_name, target_name, Ct_mean_hb, Ct_sd, quantity_mean, quantity_sd, run)

# remove duplicate rows from dataframe:
HbDF<-HbDF[!duplicated(HbDF), ]

# remove NTC rows from dataframe:
HbDF<-HbDF[!(HbDF$sample_name=="No Sample"),]

#merge dilutions with PCR data by sample name
HbDF <- merge(HbDF, EcoHbDF, by="sample_name")

# set constant values for genome copies per bee calculations:
crude_extr <- 100
eluteRNA <- 50
GITCperbee <- 200
number_bees <- 4
cDNA_eff <- 0.1
rxn_vol <- 3

#create column for total_extr_vol
HbDF$total_extr_vol <- (GITCperbee * HbDF$num_bees)

# create column for genome copies per bee:
HbDF$genome_copbeeHB <- ((((((HbDF$quantity_mean / cDNA_eff) / rxn_vol) * HbDF$dil.factor) * eluteRNA) / crude_extr) * HbDF$total_extr_vol) / HbDF$num_bees

#----------------------------------------------------------------------
# pull only actin values out of dataframe
ActinOnly <- HbDF[which(HbDF$target_name=="ACTIN"),]

# create DF of ACTIN genome copies and lab ID:
ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$genome_copbeeHB)
colnames(ActinDF) <- c("sample_name", "ACT_genome_copbeeHB")

# merge ACTIN dataframe with main dataframe:
#Need rownames and all.x=TRUE because data frames are different sizes.

HbDF <- merge(HbDF, ActinDF, by="sample_name")

# find mean of all ACTIN values:
ActinMean <- mean(ActinOnly$genome_copbeeHB, na.rm = TRUE)

# create column for normalized genome copies per bee:
HbDF$norm_genome_copbeeHB <- (HbDF$genome_copbeeHB/HbDF$ACT_genome_copbeeHB)*ActinMean


HbDF <- select(HbDF, site, target_name, norm_genome_copbeeHB,Ct_mean_hb)

#Make column for all honey bee virus data that identifies each as "honey bees collected and tested"
HbDF$HBCollected <- rep(1,length(HbDF$site))

BombSurv <- merge(BombSurv, HbDF, by=c("site","target_name"), all.x=TRUE, all.y=FALSE)

head(select(BombSurv, site, norm_genome_copbeeHB, sample_name, target_name), 200)

# make NAs 0
BombSurv$genome_copbee[is.na(BombSurv$genome_copbee)] <- 0
BombSurv$norm_genome_copbee[is.na(BombSurv$norm_genome_copbee)] <- 0
BombSurv$norm_genome_copbeeHB[is.na(BombSurv$norm_genome_copbeeHB)] <-0

# create a binary vector in dataframe
BombSurv$virusBINY_PreFilter <- ifelse(BombSurv$genome_copbee > 0, 1, 0)

#------------------------------------------------------------------------

#sam is making a dataframe to figure out which samples are problematic:
Sammy<-table(BombSurv$sample_name)

sammyscounts <-data.frame(Sammy)
head(sammyscounts)

sammy_is_so_smart <- sammyscounts[which(sammyscounts$Freq!=4),]

# rename column names! Never forget this one!
colnames(sammy_is_so_smart)<- c("ID","Freq")


#------------------------------------------------------------------------
# 2015 Bombus Survey PCR results:
# Data Cleaning Functions:
# P. Alexander Burnham
# 10, October 2016


#------------------------------------------------------------------------
# this function puts NAs in vectors where values parameters are not met then removes those columns, cleaning the data frame:

Burnhams_Fabulous_Data_Cleanser <- function(dataframe=BombSurv, variable=BombSurv$target_name){
  
  splitDF <- split(BombSurv, BombSurv$target_name)
  
  # filtering out for CT values (bombus) that are outside limit of detection, assigning either a 1 or 0 in a new column
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Ct_mean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Ct_mean > 32.525, 0, 1)
  splitDF$IAPV$virusBINY  <- ifelse(splitDF$IAPV$Ct_mean > 30.796, 0, 1)
  
  x <- rbind(splitDF$DWV,splitDF$BQCV,splitDF$IAPV)

  return(x)}


# END OF FUNCTION

# calling Burnhams_Fabulous_Data_Cleanser created above and using it on BombSurv data frame: 
keptstuff <- Burnhams_Fabulous_Data_Cleanser(BombSurv)

# make NAs 0 in virusBINY
keptstuff$virusBINY [is.na(keptstuff$virusBINY)] <- 0

#----------------------------------------------------------------------
#making genome copies 0 if CT value is below the limit of detection (Bombus)

# splitting dataframe by target name
splitkeptstuff <- split(keptstuff, keptstuff$target_name)

# make DWV norm_genome_copbee 0 if Ct value is > 32.918
splitkeptstuff$DWV$norm_genome_copbee[which(splitkeptstuff$DWV$Ct_mean > 32.918)] <- 0

# make BQCV norm_genome_copbee 0 if Ct value is > 32.525
splitkeptstuff$BQCV$norm_genome_copbee[which(splitkeptstuff$BQCV$Ct_mean > 32.525)] <- 0

# make IAPV norm_genome_copbee 0 if Ct value is > 30.796
splitkeptstuff$IAPV$norm_genome_copbee[which(splitkeptstuff$IAPV$Ct_mean > 30.796)] <- 0

# merge split dataframe back into "BombSurv" dataframe:
BombSurv <- rbind(splitkeptstuff$DWV, splitkeptstuff$BQCV, splitkeptstuff$IAPV)

#------------------------------------------------------------------------
#making genome copies 0 if CT value is below the limit of detection (apis)

# make DWV norm_genome_copbee 0 if Ct value is > 32.918
splitkeptstuff$DWV$norm_genome_copbeeHB[which(splitkeptstuff$DWV$Ct_mean_hb > 32.918)] <- 0

# make BQCV norm_genome_copbee 0 if Ct value is > 32.525
splitkeptstuff$BQCV$norm_genome_copbeeHB[which(splitkeptstuff$BQCV$Ct_mean_hb > 32.525)] <- 0

# make IAPV norm_genome_copbee 0 if Ct value is > 30.796
splitkeptstuff$IAPV$norm_genome_copbeeHB[which(splitkeptstuff$IAPV$Ct_mean_hb > 30.796)] <- 0

# merge split dataframe back into "BombSurv" dataframe:
BombSurv <- rbind(splitkeptstuff$DWV, splitkeptstuff$BQCV, splitkeptstuff$IAPV)
#-----------------------------------------------------------------------
#Change genome copy HB from 0 to NA for sites with no HB virus data (no hbs collected) using column HBCollected

BombSurv$norm_genome_copbeeHB <- ifelse(BombSurv$HBCollected == "NA", "NA", BombSurv$norm_genome_copbeeHB)

#View(BombSurv)
#-----------------------------------------------------------------------
# adding bee abundance data to BombSurv

#looking at proportion of bees that are honey bees from bee abundance data
BeeAbun <- read.csv("BeeAbun.csv", header = TRUE)

#Making new column for apis=1
BeeAbun$HB <- ifelse(BeeAbun$Morpho_species == "apis", 1, 0)

#using ddply to find proportion (mean)
BeeAbunSum<- ddply(BeeAbun, c("site"), summarise, 
                   HB_Abun = mean(HB, na.rm=TRUE))

#merge HB proportion data to bombus data
BombSurv <- merge(BombSurv, BeeAbunSum, by=c("site"), all.x=TRUE, all.y=FALSE)
#---------------------------------------------------------------------------
BombSurvSplit <- split(BombSurv, BombSurv$target_name)

BombSurvSplit$IAPV$HBSiteBin <- rep(NA, length(BombSurvSplit$IAPV$virusBINY))

BombSurvSplit$BQCV$HBSiteBin <- rep(NA, length(BombSurvSplit$BQCV$virusBINY))

BombSurvSplit$DWV$HBSiteBin <- ifelse(log(BombSurvSplit$DWV$norm_genome_copbeeHB) > 15, "High", "Low")


BombSurv <- rbind(BombSurvSplit$DWV, BombSurvSplit$BQCV, BombSurvSplit$IAPV)

BombSurv <- merge(BombSurv, Diversity, by=c("site"), all.x=TRUE, all.y=FALSE)

#---------------------------------------------------------------------------
# END OF DATA CLEANING

#Export csv

write.csv(BombSurv, file = "BombSurv.csv")

#Working with NHBS data Adding NHBS data to BombSurv
BombSurv2 <- merge(BombSurv, NHBS, by=c("site"), all.x=TRUE, all.y=FALSE)

write.csv(BombSurv2, file = "BombSurv2.csv")


#---------------------------------------------------------------------

# THE END!!!

