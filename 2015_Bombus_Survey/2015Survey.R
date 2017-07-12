###########################################################################################
# Data Analysis for 2015 Bombus Virus Study
# Samantha Alger and P. Alexander Burnham
# July 10, 2017
###########################################################################################

#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())

# Call Packages
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("plyr")
library("spdep")
library("ape")
library("lme4")
library("car")
library("ape")

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


###############################################################################################

# load site level data and merge pathogen data with GIS HB colony/apiary output:
SpatDat <- read.table("SpatDatBuffs.csv", header=TRUE,sep=",",stringsAsFactors=FALSE)
SpatDat <- select(SpatDat, -elevation, -town, -apiary, -siteNotes, -apiaryNotes)

SurvData <- read.csv("MixedModelDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
SpatialDat <- merge(SurvData, SpatDat, by = "site")

###############################################################################################
###############################################################################################
# formatting bombsurv to test Spatial Autocorralation on

# merge data to create final APC data frame:
SpatDat <- select(SpatDat, -lat, -long)
BombSurv <- merge(BombSurv, SpatDat, by = "site")

# create log virus data:
BombSurv$logVirus <- log(1+BombSurv$norm_genome_copbee)
BombSurv$logHB <- log(1+BombSurv$norm_genome_copbeeHB)

# two data frames for DWV and BQCV for Morans I
DWV <- subset(BombSurv, target_name=="BQCV")
BQCV <- subset(BombSurv, target_name=="DWV")

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
# function to pull out AIC and pval for DWV (load) and BQCV (load)


###########################################################################
# function name: AICfinderLoad
# description:finds p val and AIC for glmer model 
# parameters: 
# data = data frame, yvar and xvar
# returns a list (requires library(lme4))
###########################################################################

AICfinderLoad <- function(X=Xvar, Y="logVirus", data=DWV){
  
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

# DECIDED TO USE "sumColony1" as predictor varable based on these data (AIC and P val)

###################################################################################################
# CREATING MODELS TO TEST FOR SPATIAL AUTOCORRELATION
###################################################################################################

# create data frames to test spatial AC
SpatialDatDWV <- subset(SpatialDat, target_name=="BQCV")
SpatialDatBQCV <- subset(SpatialDat, target_name=="DWV")

#----------------------------------------------------------------------------------------------------
# BQCV PREV:

BQCVprev <- lm(data=SpatialDatBQCV, BombPrev ~ sumColonies1)
BQCVprevResid <- summary(BQCVprev)

BQCVprevResid$residual

#----------------------------------------------------------------------------------------------------
# DWV PREV

DWVprev <- lm(data=SpatialDatDWV, BombPrev ~ sumColonies1)
DWVprevResid <- summary(DWVprev)

DWVprevResid$residual

#----------------------------------------------------------------------------------------------------
# DWV LOAD

DWVload <- lm(data=SpatialDatDWV, BombusViralLoad ~ sumColonies1)
DWVloadResid <- summary(DWVload)

DWVloadResid$residual

#----------------------------------------------------------------------------------------------------
# BQCV LOAD

BQCVload <- lm(data=SpatialDatBQCV, BombusViralLoad ~ sumColonies1)
BQCVloadResid <- summary(BQCVload)

BQCVloadResid$residual

#----------------------------------------------------------------------------------------------------
# DWV HB LOAD

DWVhb <- lm(data=SpatialDatDWV, HBviralLoad ~ sumColonies1)
HBbqcvResid <- summary(DWVhb)

HBbqcvResid$residual

#----------------------------------------------------------------------------------------------------
# BQCV HB LOAD

BQCVhb <- lm(data=SpatialDatBQCV, HBviralLoad ~ sumColonies1)
HBdwvResid <- summary(BQCVhb)

HBdwvResid$residual

#----------------------------------------------------------------------------------------------------

# CREATING DISTANCE MATRICES FOR MORANS.I TEST:

#For DWV:
DWV.dists <- as.matrix(dist(cbind(SpatialDatDWV$long, SpatialDatDWV$lat)))

DWV.dists.inv <- 1/DWV.dists
diag(DWV.dists.inv) <- 0


#For BQCV:
BQ.dists <- as.matrix(dist(cbind(SpatialDatBQCV$long, SpatialDatBQCV$lat)))

BQ.dists.inv <- 1/BQ.dists
diag(BQ.dists.inv) <- 0

###################################################################################################
# TESTING FOR SPATIAL AUTOCORRELATION
###################################################################################################

# BQCV PREV:
Moran.I(BQCVprevResid$residuals, BQ.dists.inv) # NO SPACIAL-AUTO COR

# DWV PREV:
Moran.I(DWVprevResid$residuals, DWV.dists.inv) # YES SPACIAL-AUTO COR (clustered)

# BQCV LOAD:
Moran.I(BQCVloadResid$residual, BQ.dists.inv) # YES SPACIAL-AUTO COR (clustered)

# DWV LOAD:
Moran.I(DWVloadResid$residual, DWV.dists.inv) # NO SPACIAL-AUTO COR

# BQCV HB LOAD
Moran.I(HBbqcvResid$residual, BQ.dists.inv) # NO SPACIAL-AUTO COR

# DWV HB LOAD:
Moran.I(HBdwvResid$residual, DWV.dists.inv) # NO SPACIAL-AUTO COR

###################################################################################################
# CREATING FULL MODELS:
###################################################################################################


#####################################################################################
# DWV PREV ##########################################################################
#####################################################################################

DWVprevMod <- glmer(data=DWV, formula = virusBINY~sumColonies1 + Density + (1|site) + (1|species) + (1|lat) + (1|long), family = binomial(link = "logit"))

summary(DWVprevMod)
Anova(DWVprevMod)

#####################################################################################
# DWV LOAD ##########################################################################
#####################################################################################

DWVloadMod <- lmer(data=DWV, formula = norm_genome_copbee ~ sumColonies1 + Density + (1|site) + (1|species) + (1|lat) + (1|long))

summary(DWVloadMod)
Anova(DWVloadMod)

#####################################################################################
# BQCV PREV #########################################################################
#####################################################################################

BQCVprevMod <- glmer(data=BQCV, formula = virusBINY~sumColonies1 + Density + (1|site) + (1|species) + (1|lat) + (1|long), family = binomial(link = "logit"))

summary(BQCVprevMod)
Anova(BQCVprevMod)

#####################################################################################
# BQCV LOAD #########################################################################
#####################################################################################

BQCVloadMod <- lmer(data=BQCV, formula = y~x + (1|site))

Summary()
Anova()

# END MODELS

###################################################################################################
# CREATING FINAL PUBLICATION GRAPHICS:
###################################################################################################

# remove unwanted target:
BombSurvNoAIPV<-BombSurv[!BombSurv$target_name==("IAPV"),]

###################################################################################################
# Load:

# using ddply to get summary of virus load by target name and sumColonies1:
VirusSum <- ddply(BombSurvNoAIPV, c("sumColonies1", "target_name"), summarise, 
                   n = length(logVirus),
                   mean = mean(logVirus, na.rm=TRUE),
                   sd = sd(logVirus, na.rm=TRUE),
                   se = sd / sqrt(n))

VirusSum$sumColonies1 <- as.character(VirusSum$sumColonies1)
levels(VirusSum$sumColonies1) <- c(0,1,2,15,19,25,27,32)


#Create plot in ggplot 
plot <- ggplot(data = VirusSum, 
               aes(x = sumColonies1, 
                   y = mean, 
                   group = target_name, 
                   colour = target_name)
) + geom_point(size=4) + scale_colour_manual(values = c("dodgerblue4", "black")) + coord_cartesian(ylim = c(0, 20)) + labs(x = "# colonies within 1km radius", y = "log(genome copies/bee)", color="Virus:") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.4)) 

# add a theme and add asterix for significance 
plot + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .85))  


###################################################################################################
# Prevalence 


# using ddply to get summary of virus load by target name and sumColonies1:
VirusSumPrev <- ddply(BombSurvNoAIPV, c("sumColonies1", "target_name"), summarise, 
                  n = length(virusBINY),
                  mean = mean(virusBINY, na.rm=TRUE),
                  sd = sd(virusBINY, na.rm=TRUE),
                  se = sd / sqrt(n))

VirusSumPrev$sumColonies1 <- as.factor(VirusSumPrev$sumColonies1)

levels(VirusSumPrev$sumColonies1) <- c(0,1,2,15,19,25,27,32)

#Create plot in ggplot 
plot <- ggplot(data = VirusSumPrev, 
               aes(x = sumColonies1, 
                   y = mean, 
                   group = target_name, 
                   colour = target_name)
) + geom_point(size=4) + scale_colour_manual(values = c("dodgerblue4", "black")) + coord_cartesian(ylim = c(0, 1)) + labs(x = "# Colonies Within 1km Radius", y = "Precent Prevalence", color="Virus:") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.4)) + scale_y_continuous(labels = scales::percent)  

# add a theme and add asterix for significance 
plot + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) 


# create a data frame of all the other fixed effects
CoFixed <- select(SpatialDatDWV, sumColonies1, sumApiaries1, ShannonDIV, Density, HBviralLoad)

# remove duplicate values for sum Colonies 
CoFixed0<-CoFixed[CoFixed$sumColonies1==0,]
CoFixed32<-CoFixed[CoFixed$sumColonies1==32,]

# and find the column means
CoFixed0 <- colMeans(CoFixed0)
CoFixed32 <- colMeans(CoFixed32)

# remove the duplicates
CoFixed<-CoFixed[!CoFixed$sumColonies1==c(0),]
CoFixed<-CoFixed[!CoFixed$sumColonies1==c(32),]

# append to create fixed dataframe
CoFixed <- rbind(CoFixed, CoFixed0, CoFixed32)

# created scaled version of Density 
CoFixed$ScaledDensity <- (CoFixed$Density - min(CoFixed$Density))/(max(CoFixed$Density)-min(CoFixed$Density))

# created scaled version of HB viral load
CoFixed$ScaledHBviralLoad <- (CoFixed$HBviralLoad - min(CoFixed$HBviralLoad))/(max(CoFixed$HBviralLoad)-min(CoFixed$HBviralLoad))

# created scaled version of Shannon Diversity
CoFixed$ScaledShannonDIV <- (CoFixed$ShannonDIV - min(CoFixed$ShannonDIV))/(max(CoFixed$ShannonDIV)-min(CoFixed$ShannonDIV))
