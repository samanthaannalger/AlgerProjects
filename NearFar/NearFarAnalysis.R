# P. Alexander Burnham
# 25, Feb. 2017
# qPCR Results
# Data cleaning and Analysis:

##################################################################################-
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/NearFar")

# read in data:
NFdata <- read.table("NearFar_qPCR_results.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

DilDF <- read.table("DilutionFactors.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 


##################################################################################-
# load packages:

# remove uneeded columns from DF:
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
##################################################################################-
# merge dil factors with main DF
NFdata <- merge(NFdata, DilDF, by="ID")



# set constant values for genome copies per bee calculations:
crude_extr <- 100
eluteRNA <- 50
total_extr_vol <- 600
number_bees <- 1
cDNA_eff <- 0.1
rxn_vol <- 3

# create column for genome copies per bee:
NFdata$genome_copbee <- ((((((NFdata$QuantityMean / cDNA_eff) / rxn_vol) * NFdata$dilFactor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees











