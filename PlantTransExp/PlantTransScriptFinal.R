# Plant Transmission Experiment
# Samantha Alger
# April 24, 2018

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
#setwd("~/AlgerProjects/PlantTransExp/CSV_Files")
setwd("~/Documents/GitHub/AlgerProjects/PlantTransExp/CSV_Files")

# packages required:
library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(car)
library(multcomp)
library(lsmeans)


# read in data:
plantTrans <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 
virusLoad <- read.csv("datTransPos.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 
PlantVL <- read.csv ("CleanPlantVirus.csv", header= TRUE, sep = ",", stringsAsFactors=FALSE)


# read in data, (skip over meta data in csv file: skip = 9) :
VideoData <- read.csv("PlantTransVideoData.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE, skip = 9)
PlantDF <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
VirusDetect <- read.csv("VirusDetectionPlants.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)



##################################################################################################
################################# Data Manipulation ##############################################
##################################################################################################

#Merge video data with the detection data:
#VideoDataMerge <- merge(VideoData, VirusDetect, by=c("expID"), all.x=TRUE, all.y=TRUE)

# recode DF groups as control, treatment, pre-experiment
plantTrans$group[plantTrans$group == "C"] <- "Bombus Only"
plantTrans$group[plantTrans$group == "T"] <- "HB + Bombus"
plantTrans$group[plantTrans$group == "PRE"] <- "Pre Experiment"


# make data frame of duration data
expID <- c("RC_acute", "BFT_acute", "WC_acute", "WC_chronic-1", "WC_chronic-2", "WC_chronic-3", "BFT_diversity", "RC_diversity", "WC_diversity")
duration <- c(173.50, 210.73, 228.08, 213.02, 181.17, 177.87, 210.08, 210.08, 210.08)
durDat <- data.frame(expID, duration)


# Create table for forage time data
VidDat <- ddply(VideoData, c("expID"), summarise, 
                       visits = length(expID),
                       foragetime = sum(Forage, na.rm=TRUE))

# merge duraction data and normalize per hr:
VidDat <- merge(VidDat, durDat)
VidDat$visits <- 60*(VidDat$visits/VidDat$duration)


#merge new visitation data with the virus DF:
plantTrans <- merge(VidDat, plantTrans,by=c("expID"), all.y=TRUE)

# select columns from virus load DF and Merge with plantTrans DF:
virusLoad <- dplyr::select(virusLoad, labID, target_name, genomeCopy)
plantTrans <-merge(plantTrans, virusLoad, by=c("labID","target_name"), all.x=TRUE)

# Make all genome copyies 0 if Binyprefilter is 0, This is because the virus load data clean up did not take into account the duplicate/triplicate discrepancies...
plantTrans$genomeCopy <- ifelse(plantTrans$BINYprefilter == 0,0, plantTrans$genomeCopy)

# Add new column for 'experiment' type so that all 'comingle (1-3)' experiments are coded as 'comingle'
plantTrans$experiment <- ifelse(plantTrans$exp == "chronic-1"|plantTrans$exp =="chronic-2"| plantTrans$exp == "chronic-3", "chronic", plantTrans$exp)

# Plant Trans Exp_____________________
# using ddply to summarize data for treatment groups: 
PlantVirusSum <- ddply(plantTrans, c("target_name", "group"), summarise, 
                       n = length(BINYprefilter),
                       mean = mean(BINYprefilter, na.rm=TRUE),
                       sd = sd(BINYprefilter, na.rm=TRUE),
                       se = sd / sqrt(n))

#Re-ordering the DF for the graph
PlantVirusSum$group <- as.factor(PlantVirusSum$group)
PlantVirusSum$group <- factor(PlantVirusSum$group, levels=c("Pre Experiment", "Bombus Only", "HB + Bombus"))

# remove bombus only from the data frame
PlantVirusSumReduced <- PlantVirusSum[!(PlantVirusSum$group=="Bombus Only"), ]

# make the data frame a factor again to remove the level
# and rename levels
PlantVirusSumReduced$group <- factor(PlantVirusSumReduced$group)
levels(PlantVirusSumReduced$group) <- c("HB Foraged", "Pre Experiment")

#subset to only include plants part of the HB+Bombus group
plantTreat <- plantTrans[ which(plantTrans$group=="HB + Bombus"), ]

#Checking out by plant species
plantSpp <- ddply(plantTreat, c("target_name", "spp"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(BINYprefilter, na.rm=TRUE),
                  sd = sd(BINYprefilter, na.rm=TRUE),
                  se = sd / sqrt(n))

plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"


##################################################################################################
################################# Statistical Models #############################################
##################################################################################################

####################################################################################
##### Prevalence
####################################################################################

#Full Model:
ModDat <- plantTreat
ModDat$spp <- as.factor(ModDat$spp)
ModDat$target_name <- as.factor(ModDat$target_name)

# subset to remove diversity and comingle
ModDatAcute <- ModDat[!ModDat$experiment==c("diversity"), ]
ModDatAcute <- ModDatAcute[!ModDatAcute$experiment==c("comingle"), ]



# ++++++++++++++++++++++++++++ plant and virus +++++++++++++++++++++++++++++

##### Revision, attempt to seperate out experiments and run analysis seperately 
null <- glmer(data=ModDatAcute, BINYprefilter ~ (1|ID), family = binomial(link="logit"))
full <- glmer(data=ModDatAcute, BINYprefilter ~ spp + target_name + (spp * target_name) + (1|ID), family = binomial(link="logit"))
main <- glmer(data=ModDatAcute, BINYprefilter ~ target_name + spp + (1|ID), family = binomial(link="logit"))
virus <- glmer(data=ModDatAcute, BINYprefilter ~ spp + (1|ID), family = binomial(link="logit"))
flower <- glmer(data=ModDatAcute, BINYprefilter ~ target_name + (1|ID), family = binomial(link="logit"))


# full model test
anova(null, full)
# interaction test
anova(full, main)
# flower
anova(main, flower)
# virus
anova(main, virus)


# failed to get interaction posthoc
#ModDatAcute$INT <- interaction(ModDatAcute$spp, ModDatAcute$target_name)
#fullModINT <- glm(data=ModDatAcute, BINYprefilter ~ INT, family = binomial(link="logit"))
#ht = glht(fullModINT, mcp(INT="Tukey"))
#summary(ht)

#lsm <- lsmeans(full, ~ spp*target_name, adjust="tukey")
#contrast(lsm, alpha=0.05, method="pairwise", adjust=NULL)


# ++++++++++++++++++++++++++++ plant and experiment +++++++++++++++++++++++++++++

# remove comingle
ModDatNoCoM <- ModDat[!ModDat$experiment==c("comingle"), ]
ModDatNoCoM$expSimple <- ifelse(ModDatNoCoM$experiment=="diversity", "diversity", "acute")
ModDatNoCoM$expSimple <- as.factor(ModDatNoCoM$expSimple)


#Checking out by plant species
plantSpp <- ddply(ModDatNoCoM, c("target_name", "spp", "expSimple"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(BINYprefilter, na.rm=TRUE),
                  sd = sd(BINYprefilter, na.rm=TRUE),
                  se = sd / sqrt(n))

plantSppMon <- plantSpp[plantSpp$expSimple=="acute",]
plantSppDiv <- plantSpp[plantSpp$expSimple=="diversity",]


##### MODEL WORKS!!!! and is legal!!!!
full1 <- glmer(data=ModDatNoCoM, BINYprefilter ~ expSimple * spp + (1|ID), family = binomial(link="logit"))
null1 <- glmer(data=ModDatNoCoM, BINYprefilter ~ (1|ID), family = binomial(link="logit"))
main1 <- glmer(data=ModDatNoCoM, BINYprefilter ~ expSimple + spp + (1|ID), family = binomial(link="logit"))
exp1 <- glmer(data=ModDatNoCoM, BINYprefilter ~ spp + (1|ID), family = binomial(link="logit"))
plant1 <- glmer(data=ModDatNoCoM, BINYprefilter ~ expSimple + (1|ID), family = binomial(link="logit"))

anova(full1, null1) # null vs full model
anova(full1, main1) # interaction test
anova(main1, exp1) 
anova(main1, plant1)

# Testing vistis and forage times across plant species:
# remove duplicates:
visitation <- ModDat[c("spp", "visits", "foragetime", "experiment")]
visitation <- visitation[!duplicated(visitation),]



####################################################################################
##### Load
####################################################################################

# Remove zeros
Loadno0 <- ModDat[!ModDat$BINYprefilter==0,]

#Checking distribution:
hist(Loadno0$genomeCopy, breaks = 12)

Loadno0$loggenomeCopy <- log10(Loadno0$genomeCopy+1)

hist(Loadno0$loggenomeCopy,breaks=12)

noBFT<- Loadno0 [! (Loadno0$spp =="BFT"), ]

# We used an ANOVAs:
#Loadno0$visits <- as.factor(Loadno0$visits)
#Loadno0$foragetime <- as.factor(Loadno0$foragetime)
Loadno0$spp <- as.factor(Loadno0$spp)


# removove comingle and create binary exp var
Loadno0 <- Loadno0[!Loadno0$experiment=="comingle",]
Loadno0$expMerge <- ifelse(Loadno0$experiment=="diversity", "diversity", "acute")
Loadno0$expMerge <- as.factor(Loadno0$expMerge)
x <- split(Loadno0, Loadno0$expMerge)


# ++++++++++++++++++++++++++++ plant and virus +++++++++++++++++++++++++++++

# viral load models
null2 <- lmer(data=Loadno0, loggenomeCopy~(1|labID), REML=FALSE)
full2 <- lmer(data=Loadno0, loggenomeCopy~spp*target_name + (1|labID), REML=FALSE)
main2 <- lmer(data=Loadno0, loggenomeCopy~spp+target_name + (1|labID), REML=FALSE)
spp2 <- lmer(data=Loadno0, loggenomeCopy~target_name + (1|labID), REML=FALSE)
virus2 <- lmer(data=Loadno0, loggenomeCopy~spp + (1|labID), REML=FALSE)

# significance
anova(full2, main2)
anova(full2, null2)
anova(main2, spp2)
anova(main2, virus2)


# post hoc
ht = glht(full2, mcp(spp="Tukey"))
summary(ht)


# ++++++++++++++++++++++++++++ plant and experiment +++++++++++++++++++++++++++++

null3 <- lmer(data=Loadno0, loggenomeCopy~(1|labID), REML=FALSE)
full3 <- lmer(data=Loadno0, loggenomeCopy~spp*expMerge + (1|labID), REML=FALSE)
main3 <- lmer(data=Loadno0, loggenomeCopy~spp+expMerge + (1|labID), REML=FALSE)
spp3 <- lmer(data=Loadno0, loggenomeCopy~expMerge + (1|labID), REML=FALSE)
exp3 <- lmer(data=Loadno0, loggenomeCopy~spp + (1|labID), REML=FALSE)

anova(full3, main3)
anova(full3, null3)
anova(main3, spp3)
anova(main3, exp3)

summary(full3)
# post hoc
ht3 = glht(full3, mcp(spp="Tukey"))
summary(ht3)
Anova(main3)

boxplot(Loadno0$loggenomeCopy~Loadno0$expMerge, ylab = "log load")



####################################################################################
##### Models with visitation data (PREV and LOAD)
####################################################################################

# ++++++++++++++++++++++++++++ prev +++++++++++++++++++++++++++++

# visitation rates and duration of visits on viruses loads and prevalence
full4 <- glmer(data=ModDat,BINYprefilter~visits*foragetime + (1|labID), family = binomial(link = "logit"))
main4 <- glmer(data=ModDat,BINYprefilter~visits+foragetime + (1|labID), family = binomial(link = "logit"))
vis4 <- glmer(data=ModDat,BINYprefilter~foragetime + (1|labID), family = binomial(link = "logit"))
time4 <- glmer(data=ModDat,BINYprefilter~visits + (1|labID), family = binomial(link = "logit"))

anova(full4, main4)
anova(main4, vis4)
anova(main4, time4)

summary(full4)

fit = glm(BINYprefilter~visits, data=ModDat, family = binomial(link = "logit"))
newdat <- data.frame(visits=seq(min(ModDat$visits, na.rm=T), max(ModDat$visits,  na.rm=T), len=100))
newdat$BINYprefilter = predict(fit, newdata=newdat, type="response")
plot(BINYprefilter~visits, data=ModDat, col="red4")
lines(BINYprefilter~visits, newdat, col="green4", lwd=2)


fit1 = glm(BINYprefilter~foragetime, data=ModDat, family = binomial(link = "logit"))
newdat <- data.frame(foragetime=seq(min(ModDat$foragetime, na.rm=T), max(ModDat$foragetime,  na.rm=T), len=100))
newdat$BINYprefilter = predict(fit1, newdata=newdat, type="response")
plot(BINYprefilter~foragetime, data=ModDat, col="red4")
lines(BINYprefilter~foragetime, newdat, col="green4", lwd=2)

# ++++++++++++++++++++++++++++ load +++++++++++++++++++++++++++++

#Loadno0$visitsFactor <- as.factor(Loadno0$visits)
#Loadno0$foragetimeFactor <- as.factor(Loadno0$foragetime)


full5 <- lmer(data=Loadno0, loggenomeCopy~ visits * foragetime + (1|labID), REML = FALSE)
main5 <- lmer(data=Loadno0, loggenomeCopy~ visits+ foragetime+ (1|labID), REML = FALSE)
vis5 <- lmer(data=Loadno0, loggenomeCopy~ foragetime + (1|labID), REML = FALSE)
time5 <- lmer(data=Loadno0, loggenomeCopy~ visits + (1|labID), REML = FALSE)

anova(full5, main5)
anova(main5, vis5)
anova(main5, time5)

Anova(full5)

# post hoc
ht5 = glht(full5, mcp(visits="Tukey"))
summary(ht5)

table(ModDat$visits, ModDat$expID)


plot(y=Loadno0$loggenomeCopy, x=Loadno0$visits)
