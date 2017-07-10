###########################################################################################
# Data Analysis for 2015 Bombus Virus Study
# Samantha Alger and P. Alexander Burnham
# July 10, 2017
###########################################################################################

#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())

# Set Working Directory 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# load in data
BombSurv <- read.table("BombSurvNHBS.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

# remove unneeded columns from the DF
BombSurv <- select(BombSurv, -X, -Ct_mean, -Ct_sd, -quantity_mean, -quantity_sd, -run, -date_processed, -dil.factor, -genome_copbee, -Ct_mean_hb, -ID, -ACT_genome_copbee, -City, -Name, -virusBINY_PreFilter, -siteNotes, -X)

# remove unwanted sites and bombus species
BombSurv<-BombSurv[!BombSurv$site==("PITH"),]
BombSurv<-BombSurv[!BombSurv$site==("STOW"),]
BombSurv<-BombSurv[!BombSurv$species==("Griseocollis"),]
BombSurv<-BombSurv[!BombSurv$species==("Sandersonii"),]

###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# Call blue color palette for graphics
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("plyr")
library("spdep")
library("ape")
library(lme4)
library(car)


###############################################################################################
# calculating Morans-I to check for spacial auto correlation using spdep

# load site level data and merge pathogen data with GIS HB colony/apiary output:
SurvData <- read.csv("MixedModelDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
SpatDat <- read.table("SpatDatBuffs.csv", header=TRUE,sep=",",stringsAsFactors=FALSE)

# merge and remove un-needed columns:
SpatialDat <- merge(SurvData, SpatDat, by = "site")
SpatialDat <- select(SpatialDat, -X, -siteNotes, -apiaryNotes)

#Split by viruses
SpatDWV <- subset(SpatialDat, target_name=="DWV")
colnames(SpatDWV)[c(5,13)] <- c("DWVload", "DWVprev")

# rename and merge with BQCV data to create a wide form data set:
SpatBQCV <- subset(SpatialDat, target_name=="BQCV")
SpatBQCV <- SpatBQCV[c(1,5,13)] 
colnames(SpatBQCV)[c(2,3)] <- c("BQCVload", "BQCVprev")

# merge data to create final APC data frame:
SpatDat <- merge(SpatDWV, SpatBQCV, by = "site")

###############################################################################################
###############################################################################################
# formatting bonbsurv to test Spatial Autocorralation on
SpatDat <- select(SpatDat, -lat, -long, -elevation, -town, -apiaryNotes, -apiary, -siteNotes)

# merge data to create final APC data frame:
tempSurv <- merge(BombSurv, SpatDat, by = "site")

# two data frames for DWV and BQCV for Morans I
DWV <- subset(tempSurv, target_name=="BQCV")
BQCV <- subset(tempSurv, target_name=="DWV")

###############################################################################################
# function to pull out AIC and pval for DWV (prev) and BQCV (prev)


###########################################################################
# function name: AICfinderPrev
# description:finds p val and AIC for glmer model 
# parameters: 
# data = data frame, yvar and xvar
# returns a list (requires library(lme4))
###########################################################################

AICfinderPrev <- function(X=Xvar, Y="virusBINY", data=DWV){
  
  data$y <- data[,Y]
  data$x <- data[,X]
  
  Fullmod <- glmer(data=data, formula = y~x + (1|site/species), 
                   family = binomial(link = "logit"))
  
  x <- summary(Fullmod)
  
  return(list(x$AICtab[1], paste("P=", x$coefficients[2,4])))
}

###########################################################################
# END OF FUNCITON
###########################################################################

# create vector of explainitory variables to test:
Xvar <- c("sumApiaries800", "sumColonies800","sumApiaries1", "sumColonies1","sumApiaries2", "sumColonies2","sumApiaries3", "sumColonies3","sumApiaries4", "sumColonies4","sumApiaries5", "sumColonies5")

# apply funciton to run though every iteration of DWV prev:
sapply(X=Xvar, FUN=AICfinderPrev, data=DWV)

# apply funciton to run though every iteration of BQCV prev:
sapply(X=Xvar, FUN=AICfinderPrev, data=BQCV)

###############################################################################################
# function to pull out AIC and pval for DWV (prev) and BQCV (prev)


###########################################################################
# function name: AICfinderLoad
# description:finds p val and AIC for glmer model 
# parameters: 
# data = data frame, yvar and xvar
# returns a list (requires library(lme4))
###########################################################################

AICfinderLoad <- function(X=Xvar, Y="norm_genome_copbee", data=DWV){
  
  data$y <- data[,Y]
  data$x <- data[,X]
  
  Fullmod <- lmer(data=data, formula = y~x + (1|site/species))
  
  z<-Anova(Fullmod)
  
  return(list(AIC(Fullmod), paste("P=", z$`Pr(>Chisq)`)))
}

###########################################################################
# END OF FUNCITON
###########################################################################

# apply funciton to run though every iteration of DWV load:
sapply(X=Xvar, FUN=AICfinderLoad, data=DWV)

# apply funciton to run though every iteration of BQCV load:
sapply(X=Xvar, FUN=AICfinderLoad, data=BQCV)















