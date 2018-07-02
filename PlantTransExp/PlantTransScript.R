# Plant Transmission Experiment
# Samantha Alger
# April 24, 2018

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
setwd("~/AlgerProjects/PlantTransExp/CSV_Files")

# read in data:
plantTrans <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE) 

# read in data, (skip over meta data in csv file: skip = 9) :
VideoData <- read.csv("PlantTransVideoData.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE, skip = 9)

PlantDF <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

VirusDetect <- read.csv("VirusDetectionPlants.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

#Merge video data with the detection data:
VideoDataMerge <- merge(VideoData, VirusDetect, by=c("expID"), all.x=TRUE, all.y=TRUE)

library(plyr)
library(dplyr)
library(ggplot2)

# recode DF groups as control, treatment, pre-experiment

plantTrans$group[plantTrans$group == "C"] <- "Bombus Only"
plantTrans$group[plantTrans$group == "T"] <- "HB + Bombus"
plantTrans$group[plantTrans$group == "PRE"] <- "Pre Experiment"

# P12 is negative for DWV- bad meltcurve--- I manually deleted that row from the dataframe***

# Create a table for visitation data
Visits <- table(VideoDataMerge$expID)
expID <- as.vector(names(Visits))
visits <- as.vector(Visits)

VidDat <- data.frame(expID, visits)

#merge new visitation data with the virus DF:
plantTrans <- merge(VidDat,plantTrans,by=c("expID"),all.y=TRUE) 

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
           position=position_dodge()) + labs(x=NULL, y="% plants with virus detected")

#adding additional aesthetics to the figure:
#name....labels...= for legend (if there is a fill)
#theme(legend.position)= puts legend on plot
#coord_cartesian= set axis limits in this case, 0-1 because prevalence
#scale y continuous..= labels as percent

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Virus", labels=c("BQCV", "DWV")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 0.3)) + scale_y_continuous(labels = scales::percent)

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

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)


#Full Model:
#remove NAs:
ModDat <- plantTreat[complete.cases(plantTreat), ]


Fullmod3 <- glmer(data=ModDat, BINYprefilter~ target_name * spp * visits + (1|experiment), family = binomial(link="logit"))



Fullmod3 <- glm(data=ModDat, BINYprefilter ~ target_name * spp * visits, family = binomial(link="logit"))
Anova(Fullmod3)
summary(Fullmod3)


# better than null
Null <- glmer(data=ModDat, BINYprefilter ~ 1 + (1|exp), family = binomial(link="logit"))
anova(Fullmod3, Null, test = "LRT")

Anova(Fullmod3)

summary(Fullmod3)

?convergence
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
