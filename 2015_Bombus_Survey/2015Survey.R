###########################################################################################
# Data Analysis for 2015 Bombus Virus Study
# Samantha Alger and P. Alexander Burnham
# July 10, 2017
# Edited by Alex on June 30, 2018
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
library("lme4")
library("car")
library("ape")
library("MuMIn")

# Set Working Directory 

# FOR SAM:
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files") 

# FOR ALEX
#setwd("~/Documents/GitHub/AlgerProjects/2015_Bombus_Survey/CSV_Files") 

# load in data
BombSurv <- read.csv("BombSurvNHBS.csv", header=TRUE, stringsAsFactors=FALSE)
BQCVrun <- read.csv("NegStrandSamplesRan.csv", header=TRUE, stringsAsFactors=FALSE)

# plant virus prevalence data:
Plants <- read.table("plants2015DF.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

# load site level data and merge pathogen data with GIS HB colony/apiary output:
SpatDat <- read.table("SpatDatBuffs.csv", header=TRUE,sep=",",stringsAsFactors=FALSE)
SpatDat <- dplyr::select(SpatDat, -elevation, -town, -apiary, -siteNotes, -apiaryNotes)
SurvData <- read.csv("MixedModelDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
SpatialDat <- merge(SurvData, SpatDat, by = "site")

# merge data to create final APC data frame:
SpatDat <- dplyr::select(SpatDat, -lat, -long)
BombSurv <- merge(BombSurv, SpatDat, by = "site")

# remove unneeded columns from the DF
BombSurv <- dplyr::select(BombSurv, -X, -Ct_mean, -Ct_sd, -quantity_mean, -quantity_sd, -run, -date_processed, -dil.factor, -genome_copbee, -Ct_mean_hb, -ID, -ACT_genome_copbee, -City, -Name, -virusBINY_PreFilter, -siteNotes, -X)

# remove unwanted sites and bombus species
BombSurv<-BombSurv[!BombSurv$site==("PITH"),]
BombSurv<-BombSurv[!BombSurv$site==("STOW"),]

BombSurv<-BombSurv[!BombSurv$species==("Griseocollis"),]
BombSurv<-BombSurv[!BombSurv$species==("Sandersonii"),]

# create variable that bins apiaries by how many colonies are there
BombSurv$ColoniesPooled <- ifelse(BombSurv$sumColonies1 <= 0, "0", ifelse(BombSurv$sumColonies1 <= 20, "1-19","20+"))






###############################################################################################
################################ SPACIAL AUTO-CORRLATION ######################################
###############################################################################################


# formatting bombsurv to test Spatial Autocorralation on
BeeAbund <- read.table("BeeAbund.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# create log virus data:
BombSurv$logVirus <- log(1+BombSurv$norm_genome_copbee)
BombSurv$logHB <- log(1+BombSurv$norm_genome_copbeeHB)

BombSurv <- merge(BombSurv, BeeAbund, by = "site")
BombSurv$HBdensRatio <-  BombSurv$Density/((BombSurv$apis+0.0000000000000001)/10)

# two data frames for DWV and BQCV for Morans I
BQCV <- subset(BombSurv, target_name=="BQCV")
DWV <- subset(BombSurv, target_name=="DWV")

# create Plants dataframe:
Plants <- merge(Plants, BeeAbund, all.x=TRUE, all.y=FALSE)
Plants <- merge(Plants, SpatialDat, by=c("site","target_name"), all.x=TRUE, all.y=FALSE)




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

# END MODELS







###############################################################################################
################## MODEL SELECTION FOR COLONY AND APIARY VARIABLES FROM GIS ###################
###############################################################################################


############################################################
# function name: AICfinderPrev
# description:finds p val and AIC for glmer model 
# parameters: 
# data = data frame, yvar and xvar
# returns a list (requires library(lme4))
############################################################

AICfinderPrev <- function(X=Xvar, Y="virusBINY", data=DWV){
  
  data$y <- data[,Y]
  data$x <- data[,X]
  
  Fullmod <- glmer(data=data, formula = y~x + (1|site/species), 
                   family = binomial(link = "logit"))
  
  x <- summary(Fullmod)
  
  return(list(x$AICtab[1], paste("P=", x$coefficients[2,4])))
}

###############################################################
# END OF FUNCITON
###############################################################

# create vector of explainitory variables to test:
Xvar <- c("sumApiaries800", "sumColonies800","sumApiaries1", "sumColonies1","sumApiaries2", "sumColonies2","sumApiaries3", "sumColonies3","sumApiaries4", "sumColonies4","sumApiaries5", "sumColonies5")

# apply funciton to run though every iteration of DWV prev:
sapply(X=Xvar, FUN=AICfinderPrev, data=DWV)

# apply funciton to run though every iteration of BQCV prev:
sapply(X=Xvar, FUN=AICfinderPrev, data=BQCV)



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

# apply function to run though every iteration of DWV load:
sapply(X=Xvar, FUN=AICfinderLoad, data=DWV)

# apply function to run though every iteration of BQCV load:
sapply(X=Xvar, FUN=AICfinderLoad, data=BQCV)

# DECIDED TO USE "sumColony1" as predictor varable based on these data (AIC and P val)





##################################################################################################
################################################################################################## 
######################################### GRAPHICS!!!!! ##########################################
##################################################################################################
##################################################################################################






###################################################################################################
##################### CREATING FINAL PUBLICATION GRAPHICS FOR BOMBUS VIRUSES ######################
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



VirusSum1 <- ddply(BombSurvNoAIPV, c("target_name", "apiary_near_far"), summarise, 
                  n = length(virusBINY),
                  mean = mean(virusBINY),
                  sd = sqrt(((mean(virusBINY))*(1-mean(virusBINY)))/n))

VirusSum1$apiary_near_far <- as.character(VirusSum1$apiary_near_far)


colors <- c("white", "grey25")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(VirusSum1, aes(x=target_name, y=mean, fill=apiary_near_far)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) + labs(x="Virus", y = "% Prevalence") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_bw(base_size = 23) + scale_fill_manual(values=colors, name="Site Type:", labels=c("Apiary Absent", "Apiary Present")) + theme(legend.position=c(.8, .85),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) + annotate(geom = "text", x = 1, y = .98, label = "*",cex = 10) + annotate(geom = "text", x = 2, y = .25, label = "*",cex = 9) 




###################################################################################################
######################## CREATING PUBLICATION GRAPHICS FOR PLANT PREV #############################
###################################################################################################
# create a binary varaible for apiary or no apiary 
Plants$apiary <- ifelse(Plants$sumColonies1 <= 0, "no apiary","apiary")
Plants$HBlowHigh <- ifelse(Plants$apis <= 4, "Low HB","High HB")


#ddply summarize:
fieldPlantsSum <- ddply(Plants, c("target_name", "apiary"), summarise, 
                        n = length(BINYprefilter),
                        mean = mean(BINYprefilter, na.rm=TRUE),
                        sd = sqrt(((mean(BINYprefilter))*(1-mean(BINYprefilter)))/n))

# remove 0 (make NA) for values so they dont plot error bars
fieldPlantsSum$sd[fieldPlantsSum$sd==0] <- NA 
fieldPlantsSum$mean[fieldPlantsSum$mean==0] <- NA

#creating the figure

#choosing color pallet
colors <- c("white", "grey30")

plot1 <- ggplot(fieldPlantsSum, aes(x=apiary, y=mean, fill=target_name)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(y="% plants with virus detected", x="Site Type") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_bw(base_size = 23) + scale_fill_manual(values=colors, name="Virus", labels=c("BQCV", "DWV")) + theme(legend.position=c(.86, .8),legend.background = element_rect(color = "black", fill = "white", size = .4, linetype = "solid")) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent) 





###################################################################################################
############################## CREATING PUBLICATION GRAPHICS FOR HB ###############################
###################################################################################################

# histogram showing apis DWV load (bimodal)

# summary of viral load for by target and site
CopDist <- ddply(BombSurv, c("target_name", "site"), summarise, 
                 n = length(norm_genome_copbeeHB),
                 mean = mean(norm_genome_copbeeHB, na.rm=TRUE),
                 sd = sd(norm_genome_copbeeHB, na.rm=TRUE),
                 se = sd / sqrt(n))

# remove BQCV and IAPV:
CopDist<-CopDist[!CopDist$target_name==("BQCV"),]
CopDist<-CopDist[!CopDist$target_name==("IAPV"),]


ggplot(data=CopDist, aes(log(1 + mean))) + 
  geom_histogram(breaks=seq(5, 25, by = 1), 
                 col="black", 
                 fill="grey30") +
  labs(x="Apis DWV log(viral load)", y="Frequency") + theme_bw(base_size=23)

################################################################################################

# bar plot showing DWV level in apis by DWV prev in bombus

HBSiteSum <- ddply(DWV, c("HBSiteBin", "target_name"), summarise, 
                   n = length(virusBINY),
                   mean = mean(virusBINY, na.rm=TRUE),
                   sd = sqrt(((mean(virusBINY))*(1-mean(virusBINY)))/n))

# remove 0 (make NA) for values so they dont plot error bars
HBSiteSum$sd[HBSiteSum$sd==0] <- NA 
HBSiteSum$mean[HBSiteSum$mean==0] <- NA

colors <- c("grey30", "white", "white")

plot1 <- ggplot(HBSiteSum, aes(x=HBSiteBin, y=mean, fill=colors)) + 
  geom_bar(stat="identity", color = "black") + labs(x="Level of DWV in Apis", y = "% Prevalence in Bombus")

plot1 + theme_bw(base_size = 23) + scale_fill_manual(values=colors) + coord_cartesian(ylim = c(0, 0.2)) + scale_y_continuous(labels = scales::percent) + theme(legend.position=c(3, 3)) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))







##################################################################################################
################################################################################################## 
######################################### MODELS!!!!! ############################################
##################################################################################################
##################################################################################################


###################################################################################################
# CREATING MODELS FOR PLANT PREV:
###################################################################################################


# Full, Null and Reduced Models
PlantsFull <- glmer(data=Plants, formula = BINYprefilter ~ apis + bombus + target_name + Density + (1|site), family = binomial(link = "logit"))

PlantsNull <- glmer(data=Plants, formula = BINYprefilter ~ 1 + (1|site), family = binomial(link = "logit"))

PlantsApis <- glmer(data=Plants, formula = BINYprefilter ~ target_name + bombus + Density + (1|site), family = binomial(link = "logit"))

PlantsTarg <- glmer(data=Plants, formula = BINYprefilter ~ apis + bombus + Density + (1|site), family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa"))

PlantsBombus <- glmer(data=Plants, formula = BINYprefilter ~ apis + target_name + Density + (1|site), family = binomial(link = "logit"))

PlantsDensity <- glmer(data=Plants, formula = BINYprefilter ~ bombus + apis + target_name + (1|site), family = binomial(link = "logit"))

# liklihood ratio tests between models for significance
anova(PlantsFull, PlantsNull, test="LRT") # full model versus the null model
anova(PlantsFull, PlantsApis, test="LRT")
anova(PlantsFull,PlantsTarg, test="LRT")
anova(PlantsFull, PlantsBombus, test="LRT")
anova(PlantsFull, PlantsDensity, test="LRT")


# To view effects and std. errors of each variable:
summary(PlantsFull)



###################################################################################################
# CREATING FULL MODELS FOR HB:
###################################################################################################

# rename NAs "no apis caught"
DWV$HBSiteBin[is.na(DWV$HBSiteBin)] <- "No Apis Caught"

# Full, Null and Reduced Models
ApisFull <- glmer(data=DWV, formula = virusBINY ~ HBSiteBin + Density + apis + (1|site), family = binomial(link = "logit"))

ApisNull <- glmer(data=DWV, formula = virusBINY ~ 1 + (1|site), family = binomial(link = "logit"))

ApisNoHB <- glmer(data=DWV, formula = virusBINY ~ Density + apis + (1|site), family = binomial(link = "logit"))

ApisNoApis <- glmer(data=DWV, formula = virusBINY ~ HBSiteBin + Density + (1|site), family = binomial(link = "logit"))

ApisNoDens <- glmer(data=DWV, formula = virusBINY ~ HBSiteBin + apis + (1|site), family = binomial(link = "logit"))



# liklihood ratio tests between models for significance
anova(ApisFull, ApisNull, test="LRT")
anova(ApisFull, ApisNoHB, test="LRT")
anova(ApisFull, ApisNoApis, test="LRT")
anova(ApisFull, ApisNoDens, test="LRT")

# To view effects and std. errors of each variable:
summary(ApisFull)


###################################################################################################
# CREATING FULL MODELS FOR BOMBUS VIRUSES:
###################################################################################################




###########################################################################
# function name: TheExtractor
# description:extracts log liklihood test stats and p vals for null vs full
# and the reduced models
# parameters: 
# Full = full model (glmer or lmer)
# Null = null model
# Density = density removed
# Colonies = colonies removed
# Species = species removed
###########################################################################

TheExtractor <- function(Full, Null, Colonies, Density, Species){
  
  sumFull <- summary(Full)
  modelFit <- anova(Full, Null, test="LRT")
  Cols <- anova(Full, Colonies, test="LRT")
  Dens <- anova(Full, Density, test="LRT")
  Spec <- anova(Full, Species, test="LRT")
  
  ModFit <- list("Model Fit P"=modelFit$`Pr(>Chisq)`[2], "Model Fit Df"=modelFit$`Chi Df`[2], "Model Fit Chi2"=modelFit$Chisq[2])
  
  ColFit <- list("Colony Fit P"=Cols$`Pr(>Chisq)`[2],"Colony Fit Df"=Cols$`Chi Df`[2],"Colony Fit Chi2"=Cols$Chisq[2])
  
  DensFit <- list("Density Fit P"=Dens$`Pr(>Chisq)`[2],"Density Fit Df"=Dens$`Chi Df`[2],"Density Fit Chi2"=Dens$Chisq[2])
  
  SpecFit <- list("Species Fit P"=Spec$`Pr(>Chisq)`[2],"Species Fit Df"=Spec$`Chi Df`[2],"Species Fit Chi2"=Spec$Chisq[2])
  
  return(list(sumFull$coefficients[1:4,1:2],ModFit, ColFit, DensFit, SpecFit))
  
}

###########################################################################
# END OF FUNCITON
###########################################################################




#####################################################################################
# DWV PREV ##########################################################################
#####################################################################################

# Full, Null and Reduced Models
DWVprevModFull <- glmer(data=DWV, formula = virusBINY~apiary_near_far + Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevModNull <- glmer(data=DWV, formula = virusBINY~1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevModnoCols <- glmer(data=DWV, formula = virusBINY~ Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevModnoDens <- glmer(data=DWV, formula = virusBINY~apiary_near_far + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevModnoSpec <- glmer(data=DWV, formula = virusBINY~apiary_near_far + Density + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

# run the function to get results of models
TheExtractor(Full=DWVprevModFull, 
             Null=DWVprevModNull, 
             Colonies=DWVprevModnoCols, 
             Density=DWVprevModnoDens,
             Species = DWVprevModnoSpec)

#####################################################################################
# DWV LOAD ##########################################################################
#####################################################################################

# remove 0s to look at viral load of infected
DWVno0 <- DWV[!DWV$virusBINY==0,]

# Full, Null and Reduced Models
DWVloadModFull <- lmer(data=DWVno0, formula = logVirus ~ apiary_near_far + Density + species + (1|site) + (1|species) + (1|lat) + (1|long))

DWVloadModNull <- lmer(data=DWVno0, formula = logVirus ~ 1  + (1|site) + (1|lat) + (1|long))

DWVloadModnoCols <- lmer(data=DWVno0, formula = logVirus ~ Density + species + (1|site) + (1|species) + (1|lat) + (1|long))

DWVloadModnoDens <- lmer(data=DWVno0, formula = logVirus ~ apiary_near_far + species + (1|site) + (1|species) + (1|lat) + (1|long))

DWVloadModnoSpec <- lmer(data=DWVno0, formula = logVirus ~ apiary_near_far + Density + (1|site) + (1|species) + (1|lat) + (1|long))

# run the function to get results of models
TheExtractor(Full=DWVloadModFull, 
             Null=DWVloadModNull, 
             Colonies=DWVloadModnoCols, 
             Density=DWVloadModnoDens,
             Species = DWVloadModnoSpec )

#####################################################################################
# BQCV PREV #########################################################################
#####################################################################################

# Full, Null and Reduced Models
BQCVprevModFull <- glmer(data=BQCV, formula = virusBINY~apiary_near_far + Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

BQCVprevModNull <- glmer(data=BQCV, formula = virusBINY~1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

BQCVprevModnoCols <- glmer(data=BQCV, formula = virusBINY~Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

BQCVprevModnoDens <- glmer(data=BQCV, formula = virusBINY~apiary_near_far + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

BQCVprevModnoSpec <- glmer(data=BQCV, formula = virusBINY~apiary_near_far + Density + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

# run the function to get results of models
TheExtractor(Full=BQCVprevModFull, 
             Null=BQCVprevModNull, 
             Colonies=BQCVprevModnoCols, 
             Density=BQCVprevModnoDens,
             Species = BQCVprevModnoSpec)


#####################################################################################
# BQCV LOAD #########################################################################
#####################################################################################

# remove 0s to look at viral load of infected
BQCVno0 <- BQCV[!BQCV$virusBINY==0,]

# Full, Null and Reduced Models
BQCVloadModFull <- lmer(data=BQCVno0, formula = logVirus ~ apiary_near_far + Density + species + (1|site) + (1|lat) + (1|long))

BQCVloadModNull <- lmer(data=BQCVno0, formula = logVirus ~ 1 + (1|site) + (1|lat) + (1|long))

BQCVloadModnoCols <- lmer(data=BQCVno0, formula = logVirus ~ Density + species + (1|lat) + (1|long) + (1|site))

BQCVloadModnoDens <- lmer(data=BQCVno0, formula = logVirus ~ apiary_near_far + species + (1|lat) + (1|long) + (1|site))

BQCVloadModnoSpec <- lmer(data=BQCVno0, formula = logVirus ~ apiary_near_far + Density + (1|site) + (1|lat) + (1|long))


# run the function to get results of models
TheExtractor(Full=BQCVloadModFull, 
             Null=BQCVloadModNull, 
             Colonies=BQCVloadModnoCols, 
             Density=BQCVloadModnoDens,
             Species = BQCVloadModnoSpec)








































###############################################################################################
# REGRESSION ANALYSIS #########################################################################
###############################################################################################
# regressions run on sum colones excluding sites that do not have any colonies:

# DWV load by number of colonies 
DWVno0just_HB <- DWVno0[!DWVno0$sumColonies1==0,]

DWVloadModFullHB <- lmer(data=DWVno0just_HB, formula = logVirus ~ sumColonies1 + Density + species + (1|site) + (1|lat) + (1|long))

DWVloadSumColonies <- lmer(data=DWVno0just_HB, formula = logVirus ~ Density + species + (1|site) + (1|lat) + (1|long))

DWVloadDensity <- lmer(data=DWVno0just_HB, formula = logVirus ~ sumColonies1 + species + (1|site) + (1|lat) + (1|long))

DWVloadSpecies <- lmer(data=DWVno0just_HB, formula = logVirus ~ sumColonies1 + Density + (1|site) + (1|lat) + (1|long))

DWVloadModNullHB <- lmer(data=DWVno0just_HB, formula = logVirus ~ 1 + (1|site) + (1|lat) + (1|long))

anova(DWVloadModFullHB, DWVloadModNullHB, test="LRT")
anova(DWVloadModFullHB, DWVloadSumColonies, test="LRT")
anova(DWVloadModFullHB, DWVloadDensity, test="LRT")
anova(DWVloadModFullHB, DWVloadSpecies, test="LRT")


summary(DWVloadModFullHB)

# BQCV load by number of colonies 
BQCVno0just_HB <- BQCVno0[!BQCVno0$sumColonies1==0,]

BQCVloadModFullHB <- lmer(data=BQCVno0just_HB, formula = logVirus ~ sumColonies1 + Density + species + (1|site) + (1|lat) + (1|long))

BQCVloadSumColonies <- lmer(data=BQCVno0just_HB, formula = logVirus ~ Density + species + (1|site) + (1|lat) + (1|long))

BQCVloadDensity <- lmer(data=BQCVno0just_HB, formula = logVirus ~ sumColonies1 + species + (1|site) + (1|lat) + (1|long))

BQCVloadSpecies <- lmer(data=BQCVno0just_HB, formula = logVirus ~ sumColonies1 + Density + (1|site) + (1|lat) + (1|long))

BQCVloadModNullHB <- lmer(data=BQCVno0just_HB, formula = logVirus ~ 1 + (1|site) + (1|lat) + (1|long))

anova(BQCVloadModFullHB, BQCVloadModNullHB, test="LRT")
anova(BQCVloadModFullHB, BQCVloadSumColonies, test="LRT")
anova(BQCVloadModFullHB, BQCVloadDensity, test="LRT")
anova(BQCVloadModFullHB, BQCVloadSpecies, test="LRT")


summary(BQCVloadModFullHB)

# DWV prev by number of colonies 
DWVjust_HB <- DWV[!DWV$sumColonies1==0,]

DWVprevModFullHB <- glmer(data=DWVjust_HB, formula = virusBINY~ Density + species + sumColonies1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevSumColonies <- glmer(data=DWVjust_HB, formula = virusBINY~ Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevDensity <- glmer(data=DWVjust_HB, formula = virusBINY~ species + sumColonies1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevSpecies <- glmer(data=DWVjust_HB, formula = virusBINY~ Density + sumColonies1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

DWVprevModNullHB <- glmer(data=DWVjust_HB, formula = virusBINY~ 1 + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

anova(DWVprevModFullHB, DWVprevModNullHB, test="LRT")
anova(DWVprevModFullHB, DWVprevSumColonies, test="LRT")
anova(DWVprevModFullHB, DWVprevDensity, test="LRT")
anova(DWVprevModFullHB, DWVprevSpecies, test="LRT")
summary(DWVprevModFullHB)

# BQCV prev by number of colonies 
BQCVjust_HB <- BQCV[!BQCV$sumColonies1==0,]

BQCVprevModFullHB <- glmer(data=BQCVjust_HB, formula = virusBINY~sumColonies1 + Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

BQCVprevModNullHB <- glmer(data=BQCVjust_HB, formula = virusBINY~ Density + species + (1|site) + (1|lat) + (1|long), family = binomial(link = "logit"))

anova(BQCVprevModFullHB, BQCVprevModNullHB, test="LRT")
anova(BQCVprevModFullHB, BQCVprevSumColonies, test="LRT")
anova(BQCVprevModFullHB, BQCVprevDensity, test="LRT")
anova(BQCVprevModFullHB, BQCVprevSpecies, test="LRT")
summary(BQCVprevModFullHB)




###############################################################################################
# SCATTER PLOTS FOR SUM COLS VS PATHOGENS #####################################################
###############################################################################################

# load in data
spatDat <- read.csv("spatialMerge.csv", header=TRUE, stringsAsFactors=FALSE)
spatDat$logViralLoad <- log(spatDat$BombusViralLoad + 1) 

# split up my data frame:
splitSpat <- split(spatDat, spatDat$target_name)
spatDWV <- splitSpat$DWV
spatBQCV <- splitSpat$BQCV

# looking at floral density by presence or absence of apiaries (NOT SIG)
x <- aov(data=spatDat, Density~apiary_near_far)
summary(x)

# DWV load:
ggplot(spatDWV, aes(x=sumColonies1, y=logViralLoad)) +
  geom_point(size=4) + theme_bw(base_size = 23) + labs(x="# colonies in 1km", y = "DWV Bombus Viral Load") + coord_cartesian(ylim = c(0, 17))

# BQCV load:
ggplot(spatBQCV, aes(x=sumColonies1, y=logViralLoad)) +
  geom_point(size=4) + theme_bw(base_size = 23) + labs(x="# colonies in 1km", y = "BQCV Bombus Viral Load") + coord_cartesian(ylim = c(10, 20))

# DWV prev:
ggplot(spatDWV, aes(x=sumColonies1, y=BombPrev)) +
  geom_point(size=4) + theme_bw(base_size = 23) + labs(x="# colonies in 1km", y = "DWV Bombus Viral Prevalence") + coord_cartesian(ylim = c(0, 1))

# BQCV prev:
ggplot(spatBQCV, aes(x=sumColonies1, y=BombPrev)) +
  geom_point(size=4) + theme_bw(base_size = 23) + labs(x="# colonies in 1km", y = "BQCV Bombus Viral Prevalence") + coord_cartesian(ylim = c(0, 1))












