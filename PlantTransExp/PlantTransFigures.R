# Plant Transmission Experiment
# Samantha Alger
# April 24, 2018

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
setwd("~/Documents/GitHub/AlgerProjects/PlantTransExp/CSV_Files")

# read in data:
plantTrans <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 
virusLoad <- read.csv("datTransPos.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 
PlantVL <- read.csv ("CleanPlantVirus.csv", header= TRUE, sep = ",", stringsAsFactors=FALSE)

# read in data, (skip over meta data in csv file: skip = 9) :
VideoData <- read.csv("PlantTransVideoData.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE, skip = 9)
PlantDF <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)
VirusDetect <- read.csv("VirusDetectionPlants.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

#Merge video data with the detection data:
VideoDataMerge <- merge(VideoData, VirusDetect, by=c("expID"), all.x=TRUE, all.y=TRUE)

library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(car)

# recode DF groups as control, treatment, pre-experiment
plantTrans$group[plantTrans$group == "C"] <- "Bombus Only"
plantTrans$group[plantTrans$group == "T"] <- "HB + Bombus"
plantTrans$group[plantTrans$group == "PRE"] <- "Pre Experiment"

# P12 is negative for DWV- bad meltcurve--- I manually deleted that row from the dataframe***
# Create table for forage time data
VidDat <- ddply(VideoDataMerge, c("expID"), summarise, 
                       visits = length(expID),
                       foragetime = sum(Forage, na.rm=TRUE))

#merge new visitation data with the virus DF:
plantTrans <- merge(VidDat,plantTrans,by=c("expID"),all.y=TRUE)

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


#creating the figure

#choosing color pallet
colors <- c("olivedrab", "darkolivegreen2")

# aes: aesthetics
# geom_bar = type of graph
# stat="identity.....dodge())"  part of the code required to make grouped columns

plot1 <- ggplot(PlantVirusSum, aes(x=group, y=mean, fill=target_name)) +
  geom_bar(stat="identity",
           position=position_dodge()) + labs(x=NULL, y="% flowers with virus detected")

#adding additional aesthetics to the figure:
#name....labels...= for legend (if there is a fill)
#theme(legend.position)= puts legend on plot
#coord_cartesian= set axis limits in this case, 0-1 because prevalence
#scale y continuous..= labels as percent

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Virus", labels=c("BQCV", "DWV")) + theme(legend.position=c(.8, .9)) + coord_cartesian(ylim = c(0, 0.3)) + scale_y_continuous(labels = scales::percent)

#Creating figure for presentation Eagle Hill


# remove bombus only from the data frame
PlantVirusSumReduced <- PlantVirusSum[!(PlantVirusSum$group=="Bombus Only"), ]
# make the data frame a factor again to remove the level
# and rename levels
PlantVirusSumReduced$group <- factor(PlantVirusSumReduced$group)

levels(PlantVirusSumReduced$group) <- c("HB Foraged", "Pre Experiment")

plot1 <- ggplot(PlantVirusSumReduced, aes(x=group, y=mean, fill=target_name)) +
  geom_bar(stat="identity",
           position=position_dodge()) + labs(x=NULL, y="% plants with virus detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Virus", labels=c("BQCV", "DWV")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 0.3)) + scale_y_continuous(labels = scales::percent)


#Checking out virus detected on all treatment plants by plant species (for all experiments)

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


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by plant species (aes= aesthetics):
plot1 <- ggplot(plantSpp, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pretense"), italic("T. repens")))) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

############################################
#Full Model:
ModDat <- plantTreat

# Testing vistis and forage times across plant species:
# remove duplicates:

visitation <- ModDat[c("spp", "visits", "foragetime", "experiment")]
visitation <- visitation[!duplicated(visitation),]
#visitation <- visitation[!is.na(visitation$visits), ]

hist(visitation$visits)
hist(visitation$foragetime)


#

barplot(ModDat$foragetime, ModDat$BINYprefilter)

fig1 <- ddply(ModDat, c("BINYprefilter"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(foragetime, na.rm=TRUE),
              sd = sd(foragetime, na.rm=TRUE),
              se = sd / sqrt(n))


fig1 <- ggplot(fig1, aes(x=BINYprefilter, y=mean)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="virus (yes no)", y = "sum Forage time")

fig1
# VIRUS LOAD
# Remove zeros
Loadno0 <- ModDat[!ModDat$BINYprefilter==0,]

#Checking distribution:
hist(Loadno0$genomeCopy, breaks = 12)

Loadno0$loggenomeCopy <- log(Loadno0$genomeCopy)

hist(Loadno0$loggenomeCopy,breaks=12)

noBFT<- Loadno0 [! (Loadno0$spp =="BFT"), ]

# We used an ANOVAs:
#Loadno0$visits <- as.factor(Loadno0$visits)
#Loadno0$foragetime <- as.factor(Loadno0$foragetime)
Loadno0$spp <- as.factor(Loadno0$spp)





plot(Loadno0$spp~Loadno0$loggenomeCopy)
###############################################
# VIRUS LOAD FIGURE
VirusLoadFig <- ddply(Loadno0, c("target_name", "spp"), summarise, 
                  n = length(spp),
                  mean = mean(loggenomeCopy, na.rm=TRUE),
                  sd = sd(loggenomeCopy, na.rm=TRUE),
                  se = sd / sqrt(n))


# Virus load by plant species:
#target_name spp n       mean         sd         se
#1        BQCV BFT 3 255464.421 273413.507 157855.362
#2        BQCV  RC 3  75034.743  47963.947  27691.998
#3        BQCV  WC 1   1481.667         NA         NA
#4         DWV  RC 4  63154.050  65900.381  32950.190
#5         DWV  WC 4   2642.083   4970.145   2485.073

#plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
#plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
#plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by plant species (aes= aesthetics):
plot1 <- ggplot(VirusLoadFig, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 15)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.5))



plot(Loadno0$genomeCopy~Loadno0$foragetime)


# Figure of the log virus genome copies 
LoadPlot <- ggplot(Loadno0, aes(x=target_name, y=log10(genomeCopy+1), fill=spp)) +
  labs(x=NULL, y = "virus load log(genome copies)")+
  theme_classic(base_size = 20) +  
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  geom_dotplot(binaxis='y', stackdir = 'center', dotsize = 0.5, position = position_dodge()) + 
  scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pretense"), italic("T. repens")))) +
  theme(legend.position = c(.2,.2)) +
  coord_cartesian(ylim = c(0,6))


LoadPlot



#stats for % prevalence, species differences
statsplit <- split(plantTreat, plantTreat$target_name)

# chi.sq test for BQCV detected on flowers by plant spp.
chisq.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

fisher.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

# chi.sq test for DWV detected on flowers by plant spp.
chisq.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)

fisher.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)


#Parse out each experiment and check out data

plantTreat <- plantTrans[ which(plantTrans$group=="HB + Bombus"), ]
plantComingle <- plantTreat[which(plantTreat$exp=="comingle"),]
plantDiversity <- plantTreat[which(plantTreat$exp=="diversity"),]
plantAcute <- plantTreat[ which(plantTreat$exp=="chronic-1" | plantTreat$exp =="acute" | plantTreat$exp=="chronic-2"| plantTreat$exp=="chronic-3"), ]

#Checking out by plant species, ACUTE
plantSpp <- ddply(plantAcute, c("target_name", "spp"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(BINYprefilter, na.rm=TRUE),
                  sd = sd(BINYprefilter, na.rm=TRUE),
                  se = sd / sqrt(n))



plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by plant species (aes= aesthetics):
plot1 <- ggplot(plantSpp, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

# Stats for each experiment:
#stats for % prevalence, species differences
statsplit <- split(plantAcute, plantAcute$target_name)

# chi.sq test for BQCV detected on flowers by plant spp.
chisq.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

fisher.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

# chi.sq test for DWV detected on flowers by plant spp.
chisq.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)

fisher.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)


# What's the percent detection for each plant species for each virus?
plantAcuteDWV <- plantAcute[ which(plantAcute$target_name=="DWV"), ]
plantAcuteBQCV <- plantAcute[ which(plantAcute$target_name=="BQCV"), ]

tapply(plantAcuteDWV$BINYprefilter, plantAcuteDWV$spp, mean)
#BFT        RC        WC 
#0.0000000 0.3333333 0.2307692 
tapply(plantAcuteBQCV$BINYprefilter, plantAcuteBQCV$spp, mean)
#BFT         RC         WC 
#1.00000000 0.00000000 0.08333333 


# For Diversity:
plantSpp <- ddply(plantDiversity, c("target_name", "spp"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(BINYprefilter, na.rm=TRUE),
                  sd = sd(BINYprefilter, na.rm=TRUE),
                  se = sd / sqrt(n))


plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"

# What's the percent detection for each plant species for each virus?
plantDiversityDWV <- plantDiversity[ which(plantDiversity$target_name=="DWV"), ]
plantDiversityBQCV <- plantDiversity[ which(plantDiversity$target_name=="BQCV"), ]

tapply(plantDiversityDWV$BINYprefilter, plantDiversityDWV$spp, mean)
#BFT        RC        WC 
#0            1       0 
tapply(plantDiversityBQCV$BINYprefilter, plantDiversityBQCV$spp, mean)
#BFT         RC         WC 
#0.00       0.75       0.00  


# For Comingle:
plantSpp <- ddply(plantComingle, c("target_name", "spp"), summarise, 
                  n = length(BINYprefilter),
                  mean = mean(BINYprefilter, na.rm=TRUE),
                  sd = sd(BINYprefilter, na.rm=TRUE),
                  se = sd / sqrt(n))


plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"

# What's the percent detection for each plant species for each virus?
plantComingleDWV <- plantComingle[ which(plantComingle$target_name=="DWV"), ]
plantComingleBQCV <- plantComingle[ which(plantComingle$target_name=="BQCV"), ]

tapply(plantComingleDWV$BINYprefilter, plantComingleDWV$spp, mean)
# WC 
# 0.2 
tapply(plantComingleBQCV$BINYprefilter, plantComingleBQCV$spp, mean)
# WC 
# 0




#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(plantSpp, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

#stats for % prevalence, species differences
statsplit <- split(plantTreat, plantTreat$target_name)

# chi.sq test for BQCV detected on flowers by plant spp.
chisq.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

fisher.test(statsplit$BQCV$BINYprefilter, statsplit$BQCV$spp)

# chi.sq test for DWV detected on flowers by plant spp.
chisq.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)

fisher.test(statsplit$DWV$BINYprefilter, statsplit$DWV$spp)

#removing all WC and RC
plantTreat <- plantTreat[ which(plantTreat$spp=="BFT"), ]
fisher.test(plantTreat$BINYprefilter, plantTreat$target_name)

#

plantTreat$spp <- as.factor(plantTreat$spp)
plantTreat$target_name <- as.factor(plantTreat$target_name)

library(lme4)
plantMod <- glmer(data=plantTreat, formula = BINYprefilter ~ target_name*spp + (1|labID), family=binomial(link = "logit"))

summary(plantMod)


mod3 <- glm(data=plantTreat, BINYprefilter~target_name, family = binomial)

summary(mod3)



#Looking at differences between the experiments:

#Subsetting the diversity experiment:
plantDiv <- plantTrans[ which(plantTreat$exp=="diversity"), ]
plantAcute <- plantTrans[ which(plantTrans$exp=="acute"), ]
plantComingle <- plantTrans[ which(plantTrans$exp=="comingle"), ]
plantChronic <- plantTrans[ which(plantTrans$exp=="chronic"), ]

#Checking out by plant species
plantTransDiv <- ddply(plantDiv, c("target_name", "spp", "exp"), summarise, 
                       n = length(BINYprefilter),
                       mean = mean(BINYprefilter, na.rm=TRUE),
                       sd = sd(BINYprefilter, na.rm=TRUE),
                       se = sd / sqrt(n))

plantTransDiv$spp[plantTransDiv$spp == "RC"] <- "Red Clover"
plantTransDiv$spp[plantTransDiv$spp == "BFT"] <- "BirdsFoot Trefoil"
plantTransDiv$spp[plantTransDiv$spp == "WC"] <- "White Clover"

#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses for diversity
plot1 <- ggplot(plantTransDiv, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

plantTransDiv <- ddply(plantAcute, c("target_name", "spp", "exp"), summarise, 
                       n = length(BINYprefilter),
                       mean = mean(BINYprefilter, na.rm=TRUE),
                       sd = sd(BINYprefilter, na.rm=TRUE),
                       se = sd / sqrt(n))
plantTransDiv$mean


# Testing virus load on plants for Plant Trans and BombSurv here...

# only include virus positive plant samples
PlantVL <- PlantVL[ which(PlantVL$virusBINY==1), ]

# split by experiment
PlantVLsplit <- split(PlantVL, PlantVL$Exp)

# get the Max, Min, Average VL for each experiment:

summary(PlantVLsplit$PlantTrans$genomeCopy)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   6959   20790   47830  107200  125100  554300 
# Range = 10^3 - 10^5
# Mean = 10^5
sd(PlantVLsplit$PlantTrans$genomeCopy)
# 53358.1

summary(PlantVLsplit$Survey$genomeCopy)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    17010   39470   68340  175600  181200  702400
# Range= 10^4- 10^5
# Mean = 10^4

sd(PlantVLsplit$Survey$genomeCopy)
# 246831.9

# Determine starting VL for honey bees 
#read in data:
HB <- read.csv ("CorgiTest.csv", header= TRUE, sep = ",", stringsAsFactors=FALSE)

# read in functions:
###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 50, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 100
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: actinNormal
# description: normalizes virus data with actin values 
# parameters: a data frame with actin values
# returns a dataframe with normalized virus values 
###########################################################################

actinNormal <- function(data=MigVirus){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$target_name=="ACTIN"),]
  
  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$run, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("sample_name", "run", "ACT_genomeCopy")
  
  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("sample_name", "run"), all.x=TRUE)
  
  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)
  
  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean
  
  return(data)
}


# USE FUNCTIONS:
HBNorm<- VirusNorm(data=HB, 50)
HBVL<-actinNormal(data=HBNorm)
table(HBVL$target_name, HBVL$NormGenomeCopy)

# HB: DWV- 10^9, 
#     BQCV- 10^6
#NYHB: DWV - 10^4
#     BQCV: 10^8

