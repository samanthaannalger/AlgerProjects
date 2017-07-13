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
library("MuMIn")

# Set Working Directory 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# load in data
BombSurv <- read.table("BombSurvNHBS.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

# load site level data and merge pathogen data with GIS HB colony/apiary output:
SpatDat <- read.table("SpatDatBuffs.csv", header=TRUE,sep=",",stringsAsFactors=FALSE)
SpatDat <- select(SpatDat, -elevation, -town, -apiary, -siteNotes, -apiaryNotes)
SurvData <- read.csv("MixedModelDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
SpatialDat <- merge(SurvData, SpatDat, by = "site")


# merge data to create final APC data frame:
SpatDat <- select(SpatDat, -lat, -long)
BombSurv <- merge(BombSurv, SpatDat, by = "site")

# remove unneeded columns from the DF
BombSurv <- select(BombSurv, -X, -Ct_mean, -Ct_sd, -quantity_mean, -quantity_sd, -run, -date_processed, -dil.factor, -genome_copbee, -Ct_mean_hb, -ID, -ACT_genome_copbee, -City, -Name, -virusBINY_PreFilter, -siteNotes, -X)

# remove unwanted sites and bombus species
BombSurv<-BombSurv[!BombSurv$site==("PITH"),]
BombSurv<-BombSurv[!BombSurv$site==("STOW"),]
BombSurv<-BombSurv[!BombSurv$species==("Griseocollis"),]
BombSurv<-BombSurv[!BombSurv$species==("Sandersonii"),]

# create variable that bins apiaries by how many colonies are there
BombSurv$ColoniesPooled <- ifelse(BombSurv$sumColonies1 <= 0, "0", ifelse(BombSurv$sumColonies1 <= 20, "1-19","20+"))

###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################
# formatting bombsurv to test Spatial Autocorralation on

# create log virus data:
BombSurv$logVirus <- log(1+BombSurv$norm_genome_copbee)
BombSurv$logHB <- log(1+BombSurv$norm_genome_copbeeHB)

# two data frames for DWV and BQCV for Morans I
BQCV <- subset(BombSurv, target_name=="BQCV")
DWV <- subset(BombSurv, target_name=="DWV")

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
SpatialDatBQCV <- subset(SpatialDat, target_name=="BQCV")
SpatialDatDWV <- subset(SpatialDat, target_name=="DWV")

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
Moran.I(BQCVprevResid$residuals, BQ.dists.inv) # YES SPACIAL-AUTO COR (clustered)

# DWV PREV:
Moran.I(DWVprevResid$residuals, DWV.dists.inv) # NO SPACIAL-AUTO COR 

# BQCV LOAD:
Moran.I(BQCVloadResid$residual, BQ.dists.inv) # NO SPACIAL-AUTO COR 

# DWV LOAD:
Moran.I(DWVloadResid$residual, DWV.dists.inv) # YES SPACIAL-AUTO COR (clustered)

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

#r.squaredGLMM(DWVprevMod)

#####################################################################################
# DWV LOAD ##########################################################################
#####################################################################################

# remove 0s to look at viral load of infected
DWVno0 <- DWV[!DWV$virusBINY==0,]

DWVloadMod <- lmer(data=DWVno0, formula = logVirus ~ sumColonies1 + Density + (1|site) + (1|species) + (1|lat) + (1|long))

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

# remove 0s to look at viral load of infected
BQCVno0 <- BQCV[!BQCV$virusBINY==0,]

BQCVloadMod <- lmer(data=BQCVno0, formula = logVirus ~ sumColonies1 + Density + (1|site) + (1|species) + (1|lat) + (1|long))

summary(BQCVloadMod)
Anova(BQCVloadMod)

# END MODELS

###################################################################################################
# CREATING FINAL PUBLICATION GRAPHICS:
###################################################################################################

# remove unwanted target:
BombSurvNoAIPV<-BombSurv[!BombSurv$target_name==("IAPV"),]

###################################################################################################
# Load:

# remove 0s
BombSurvNoAIPVno0<-BombSurvNoAIPV[!BombSurvNoAIPV$logVirus==0,]

#Create plot in ggplot 
plot <- ggplot(data = BombSurvNoAIPVno0, 
               aes(x = ColoniesPooled, 
                   y = logVirus, 
                   fill = target_name)
) + geom_boxplot(color="black") + coord_cartesian(ylim = c(5, 20)) + labs(x = "# apis colonies within 1km radius", y = "log(genome copies/bee)", fill="Virus:") 

# add a theme 
plot + theme_bw(base_size = 17) + scale_fill_manual(values=c("white", "gray40")) 


###################################################################################################
# Prevalence 

VirusSum <- ddply(BombSurvNoAIPV, c("target_name", "ColoniesPooled"), summarise, 
                   n = length(virusBINY),
                   mean = mean(virusBINY),
                   sd = sqrt(((mean(virusBINY))*(1-mean(virusBINY)))/n))


#Create plot in ggplot 
plot1 <- ggplot(data = VirusSum, 
               aes(x = ColoniesPooled, 
                   y = mean, 
                   shape = target_name)
) + geom_point(size=4) + coord_cartesian(ylim = c(0, 1)) + labs(x = "# apis colonies within 1km radius", y = "% prevalence", shape="Virus:") + scale_y_continuous(labels = scales::percent) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2))

# add a theme 
plot1 + theme_bw(base_size = 17) + scale_shape_manual(values=c(19, 1)) + annotate(geom = "text", x = 1, y = .11, label = "n=205",cex = 4) + annotate(geom = "text", x = 2, y = .18, label = "n=71",cex = 4) + annotate(geom = "text", x = 3, y = .3, label = "n=92",cex = 4) + annotate(geom = "text", x = 1, y = .72, label = "n=188",cex = 4) + annotate(geom = "text", x = 2, y = 1, label = "n=62",cex = 4) + annotate(geom = "text", x = 3, y = .98, label = "n=88",cex = 4) 


























