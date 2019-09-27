# Samantha Alger
#Negative Strand Analysis and figures
# 4/5/2018

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

library("ggplot2")
library("dplyr")
library("lme4")
library("car")
library("plyr")

# load in data
Melt <- read.csv("USDAplate1Melt.csv", header=TRUE, stringsAsFactors=FALSE)
Cq <- read.csv("USDAplate1cq.csv", header=TRUE, stringsAsFactors=FALSE)
BombSurv <- read.csv("BombSurvNHBS.csv", header=TRUE, stringsAsFactors=FALSE)
# formatting bombsurv to test Spatial Autocorralation on
BeeAbund <- read.table("BeeAbund.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
SpatDat <- read.table("SpatDatBuffs.csv", header=TRUE,sep=",",stringsAsFactors=FALSE)


# remove unwanted sites and bombus species
BombSurv<-BombSurv[!BombSurv$site==("PITH"),]
BombSurv<-BombSurv[!BombSurv$site==("STOW"),]

BombSurv<-BombSurv[!BombSurv$species==("Griseocollis"),]
BombSurv<-BombSurv[!BombSurv$species==("Sandersonii"),]

# subset BombSurv:
Bomb <- dplyr::select(BombSurv, site, Ct_mean, sample_name, species, apiary_near_far, Density, genome_copbee, norm_genome_copbeeHB, target_name, virusBINY_PreFilter, virusBINY, HBSiteBin)
names(Bomb)[3] <- "Sample"

# merge data:
#Dat <- merge(Melt, Cq, by = c("Sample", "Target"))
#str(Dat)

# Merge Dat and Bomb
Dat <- merge(Melt, Bomb, by = c("Sample","target_name"),  all.y=TRUE)

#Dat <- merge(Melt, Bomb, by = c("Sample"),  all.x=TRUE)

DatClean <- Dat
#DatClean <- DatClean[!(DatClean$Cq>33),]
#DatClean <- DatClean[!(DatClean$Melt<78),]
DatClean$BinaryNeg <- ifelse(DatClean$Melt > 0, 1, 0)
DatClean$BinaryNeg[is.na(DatClean$BinaryNeg)] <- 0


x <- merge(DatClean, BeeAbund, by="site")
DatClean <- merge(x, SpatDat, by="site") 

DatCleanPos <- DatClean[DatClean$virusBINY==1,]

DatCleanPos[DatCleanPos$apis==0,]


DatClean$isHB <- ifelse(DatClean$site=="TIRE" | 
                          DatClean$site=="CLERK" |  
                          DatClean$site=="NEK" |  
                          DatClean$site=="FLAN",
                        "noHB", "HB")

ddply(DatClean, c("target_name", "isHB"), summarise, 
      n = length(BinaryNeg),
      mean = mean(BinaryNeg),
      sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))



# Subset for the two viruses: 
# For BQCV:
BQ <- DatClean[DatClean$target_name=="BQCV",]
# For DWV:
DW <- DatClean[DatClean$target_name=="DWV",]
reducedBQ <- select(BQ, BinaryNeg, Density, apis, apiary_near_far, species, site, Sample, lat, long)
reducedDW <- select(DW, BinaryNeg, Density, apis, apiary_near_far, species, site, Sample, lat, long)



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






BQCVprevModFull <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + Density + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModNull <- glmer(data=reducedBQ, formula = BinaryNeg ~ 1 + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoApis <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + Density + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoDens <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoSpp <- glmer(data=reducedBQ, formula = BinaryNeg ~  Density + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

# run the function to get results of models
TheExtractor(Full=BQCVprevModFull, 
             Null=BQCVprevModNull, 
             Colonies=BQCVprevModnoApis, 
             Density=BQCVprevModnoDens,
             Species = BQCVprevModnoSpp)







BQCVprevModFull <- glmer(data=reducedDW, formula = BinaryNeg ~ species + Density + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModNull <- glmer(data=reducedDW, formula = BinaryNeg ~ 1 + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoApis <- glmer(data=reducedDW, formula = BinaryNeg ~ species + Density + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoDens <- glmer(data=reducedDW, formula = BinaryNeg ~ species + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

BQCVprevModnoSpp <- glmer(data=reducedDW, formula = BinaryNeg ~  Density + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

# run the function to get results of models
TheExtractor(Full=BQCVprevModFull, 
             Null=BQCVprevModNull, 
             Colonies=BQCVprevModnoApis, 
             Density=BQCVprevModnoDens,
             Species = BQCVprevModnoSpp)










#DWVprevModFull <- glmer(data=reducedDW, formula = BinaryNeg ~ species + Density + apis + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))











BQCVprevModFull2 <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + Density + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))




DWVprevModFull3 <- glmer(data=reducedDW, formula = BinaryNeg ~ species + Density + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModNull3 <- glmer(data=reducedDW, formula = BinaryNeg ~ 1 + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3noApis <- glmer(data=reducedDW, formula = BinaryNeg ~ species + Density + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3noDensity <- glmer(data=reducedDW, formula = BinaryNeg ~ species + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3nospecies <- glmer(data=reducedDW, formula = BinaryNeg ~ Density + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))



# run the function to get results of models
TheExtractor(Full=DWVprevModFull3, 
             Null=DWVprevModNull3, 
             Colonies=DWVprevModFull3noApis, 
             Density=DWVprevModFull3noDensity,
             Species =DWVprevModFull3nospecies)













DWVprevModFull3 <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + Density + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModNull3 <- glmer(data=reducedBQ, formula = BinaryNeg ~ 1 + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3noApis <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + Density + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3noDensity <- glmer(data=reducedBQ, formula = BinaryNeg ~ species + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))

DWVprevModFull3nospecies <- glmer(data=reducedBQ, formula = BinaryNeg ~ Density + apiary_near_far + (1|site) + (1|long) + (1|lat), family = binomial(link = "logit"))



# run the function to get results of models
TheExtractor(Full=DWVprevModFull3, 
             Null=DWVprevModNull3, 
             Colonies=DWVprevModFull3noApis, 
             Density=DWVprevModFull3noDensity,
             Species =DWVprevModFull3nospecies)










# Fig and stats for BQCV:
BQ <- BQ[ which(BQ$virusBINY_PreFilter=="1"), ]

#ddply summarize:
plotdat <- ddply(BQ, c("target_name", "apiary_near_far"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))

plotdat$apiary_near_far <- ifelse(plotdat$apiary_near_far==0, "No Apiary", "Apiary")

label.df <- data.frame(Group = c("S1", "S2"),
                       Value = c(6, 9))

plot1 <- ggplot(plotdat, aes(x=apiary_near_far, y=mean, fill=target_name)) +
  geom_bar(stat="identity", color="black", 
           fill =   "white",
           position=position_dodge()) + labs(y="BQCV Replication", x="Site Type") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)

DatCleanNeg <- DatClean[DatClean$target_name=="BQCV",]















#Calculate percentage of replicating infections:
mean(BQ$BinaryNeg)
# Overall 20% of BQCV positive bumble bees had replicating infections.


#ddply summarize for species:
plotdat <- ddply(BQ, c("target_name", "species"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))


plot1 <- ggplot(plotdat, aes(x=species, y=mean, fill=target_name)) + 
  geom_bar(stat="identity", color="black",fill = "white",
           position=position_dodge()) + labs(y="Prevalence", x="Species") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)

#Percentage of replication by species:
plotdat

# 28% of bimacs, 11% of Vagans


#For DWV:

# subset for virus positive bees
DW <- DW[ which(DW$virusBINY_PreFilter=="1"), ]

DW$virusBINY


#ddply summarize:

plotdat2 <- ddply(DW, c("target_name", "apiary_near_far"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))

plotdat2$apiary_near_far <- ifelse(plotdat2$apiary_near_far==0, "No Apiary", "Apiary")

label.df <- data.frame(Group = c("S1", "S2"),
                       Value = c(6, 9))

plot1 <- ggplot(plotdat2, aes(x=apiary_near_far, y=mean, fill=target_name)) +
  geom_bar(stat="identity", color="black", 
           fill =   "white",
           position=position_dodge()) + labs(y="DWV Replication", x="Site Type") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)

DatCleanNeg <- DatClean[DatClean$target_name=="DWV",]

chisq.test(DatCleanNeg$BinaryNeg, DatCleanNeg$apiary_near_far)
chisq.test(DatCleanNeg$BinaryNeg, DatCleanNeg$species)

#Calculate % of replicating infections
mean(DW$BinaryNeg)
# 16% of DWV positive bees had replicating infections.

#ddply summarize for species:
plotdat <- ddply(DW, c("target_name", "species"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))

#Replication by species
plotdat
# bimacs 22%; Vagans 12%

plot1 <- ggplot(plotdat, aes(x=species, y=mean, fill=target_name)) + 
  geom_bar(stat="identity", color="black",fill = "white",
           position=position_dodge()) + labs(y="Prevalence", x="Species") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)
plotdat

