# Author: Samantha Alger, and P. Alexander Burnham
# Ross Conrad PCR results
# 19, October 2016

#-----------------------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/Desktop/RScripts/Ross Conrad")

# read in data:
BombSurvDF <- read.table("qPCR_results_dataframe_Ross.csv", header=TRUE, sep = ",") 

ConradDF <- BombSurvDF[grep("\\RC",BombSurvDF$sample_name),]

library(dplyr)

# take only columns that we want:
ConradDF <- select(ConradDF, sample_name, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)

#remove IAPV
ConradDF <- ConradDF[!(ConradDF$target_name=="IAPV"),]

#read in dilution csv
Conrad_dil <- read.csv("RossConrad_dilutions.csv",header = TRUE)

#merge dilutions with PCR data by sample name
ConradDF <- merge(ConradDF, Conrad_dil, by="sample_name")

# remove duplicate rows from dataframe:
ConradDF<-ConradDF[!duplicated(ConradDF), ]

#remove first row because it is an NA
ConradDF <- ConradDF[-1,]

# set constant values for genome copies per bee calculations:
crude_extr <- 100
eluteRNA <- 50
total_extr_vol <- 600
number_bees <- 4
cDNA_eff <- 0.1
rxn_vol <- 3

# create column for genome copies per bee:
ConradDF$genome_copbee_HB <- ((((((ConradDF$quantity_mean / cDNA_eff) / rxn_vol) * ConradDF$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees

# pull only actin values out of dataframe
ActinOnly <- ConradDF[which(ConradDF$target_name=="ACTIN"),]

# create DF of ACTIN genome copies and lab ID:
ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$genome_copbee)
colnames(ActinDF) <- c("sample_name", "ACT_genome_copbee")

# merge ACTIN dataframe with main dataframe:
#Need rownames and all.x=TRUE because data frames are different sizes.

ConradDF <- merge(ConradDF, ActinDF, by="sample_name")

# find mean of all ACTIN values:
ActinMean <- mean(ActinOnly$genome_copbee, na.rm = TRUE)

# create column for normalized genome copies per bee:
ConradDF$norm_genome_copbee <- (ConradDF$genome_copbee/ConradDF$ACT_genome_copbee)*ActinMean

#remove Actin column (only keep DWV column)
ConradDF <- ConradDF[which(ConradDF$target_name=="DWV"),]

#create new vector of IDs without 'RC'
ConradDF$Colony_Number <- c(1,10,11,12,2,3,4,5,7,8,9)

#sort by this new vector
ConradDF <- ConradDF[order(ConradDF$Colony_Number),]

#barplot of vial load by colony number
barplot(ConradDF$norm_genome_copbee, 
        names.arg = ConradDF$sample_name, 
        xlab="Colony ID",
        ylab="DWV Viral Load (average genome copies per bee)",
        main="Deformed Wing Virus Loads by Colony in Apiary '125'",
        col='goldenrod',
        font.lab=2
        )

#write csv file
write.csv(ConradDF, file = "ConradDF.csv")

