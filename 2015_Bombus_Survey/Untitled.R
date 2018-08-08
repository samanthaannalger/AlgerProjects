# BombSurv Map showing hot spots for DWV and BQCV in bumble bees


# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files/")

# source my packages


###########################################################
# Read in Data
SurvMap <- read.csv("spatialMerge.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)
