# Plant Transmission Experiment
# Samantha Alger
# April 24, 2018

# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
#setwd("~/AlgerProjects/PlantTransExp/CSV_Files")
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

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pratense"), italic("T. repens")))) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

############################################
library(multcomp)


#Full Model:
ModDat <- plantTreat
ModDat$spp <- as.factor(ModDat$spp)
ModDat$target_name <- as.factor(ModDat$target_name)



# subset to remove diversity and comingle
ModDatAcute <- ModDat[!ModDat$experiment==c("diversity"), ]
ModDatAcute <- ModDatAcute[!ModDatAcute$experiment==c("comingle"), ]


##### Revision, attempt to seperate out experiments and run analysis seperately 
null <- glmer(data=ModDatAcute, BINYprefilter ~ (1|ID), family = binomial(link="logit"))
full <- glmer(data=ModDatAcute, BINYprefilter ~ spp + target_name + (spp * target_name) + (1|ID), family = binomial(link="logit"))
main <- glmer(data=ModDatAcute, BINYprefilter ~ target_name + spp + (1|ID), family = binomial(link="logit"))
virus <- glmer(data=ModDatAcute, BINYprefilter ~ spp + (1|ID), family = binomial(link="logit"))
flower <- glmer(data=ModDatAcute, BINYprefilter ~ target_name + (1|ID), family = binomial(link="logit"))


<<<<<<< HEAD
chisq.test(table(ModDatAcute$BINYprefilter, ModDatAcute$target_name))

=======
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613

# full model test
anova(null, full)
# interaction test
anova(full, main)
# flower
anova(main, flower)
# virus
anova(main, virus)

# failed to get interaction posthoc
library(multcomp)
ModDatAcute$INT <- interaction(ModDatAcute$spp, ModDatAcute$target_name)
fullModINT <- glm(data=ModDatAcute, BINYprefilter ~ INT, family = binomial(link="logit"))
ht = glht(fullModINT, mcp(INT="Tukey"))
summary(ht)

library(lsmeans)
lsm <- lsmeans(full, ~ spp*target_name, adjust="tukey")
contrast(lsm, alpha=0.05, method="pairwise", adjust=NULL)


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

#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by plant species (aes= aesthetics):
plot1 <- ggplot(plantSppMon, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Single Species", y = "% of Flowers with Virus Detected")+ theme_minimal(base_size = 16) + scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pratense"), italic("T. repens")))) + theme(legend.position=c(.7, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) 


#Create a bar graph for viruses by plant species (aes= aesthetics):
plot2 <- ggplot(plantSppDiv, aes(x=target_name, y=mean, fill=spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Diversity", y = NULL) + theme_minimal(base_size = 16) + scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pratense"), italic("T. repens")))) + theme(legend.position="none", plot.margin = unit(c(0,0,0,0), "cm")) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = NULL)

library(patchwork)
plot1+plot2



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


full1 <- glm(data=ModDatNoCoM, BINYprefilter ~ expSimple * spp, family = binomial(link="logit"))
anova(full1, test="LRT")

lsm1 <- lsmeans(full1, ~ spp*expSimple)
contrast(lsm1)







# Testing vistis and forage times across plant species:
# remove duplicates:

visitation <- ModDat[c("spp", "visits", "foragetime", "experiment")]
visitation <- visitation[!duplicated(visitation),]
#visitation <- visitation[!is.na(visitation$visits), ]

hist(visitation$visits)
hist(visitation$foragetime)




head(ModDat)


# VIRUS LOAD
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

#kruskal.test(data=x$acute, loggenomeCopy ~ target_name)
#kruskal.test(data=x$acute, loggenomeCopy ~ spp)
#kruskal.test(data=Loadno0, loggenomeCopy ~ expMerge)

null2 <- lmer(data=Loadno0, loggenomeCopy~(1|labID), REML=FALSE)
full2 <- lmer(data=Loadno0, loggenomeCopy~spp*target_name + (1|labID), REML=FALSE)
main2 <- lmer(data=Loadno0, loggenomeCopy~spp+target_name + (1|labID), REML=FALSE)
spp2 <- lmer(data=Loadno0, loggenomeCopy~target_name + (1|labID), REML=FALSE)
virus2 <- lmer(data=Loadno0, loggenomeCopy~spp + (1|labID), REML=FALSE)

anova(full2, main2)
anova(full2, null2)
anova(main2, spp2)
anova(main2, virus2)







null3 <- lmer(data=Loadno0, loggenomeCopy~(1|labID), REML=FALSE)
full3 <- lmer(data=Loadno0, loggenomeCopy~spp*expMerge + (1|labID), REML=FALSE)
main3 <- lmer(data=Loadno0, loggenomeCopy~spp+expMerge + (1|labID), REML=FALSE)
spp3 <- lmer(data=Loadno0, loggenomeCopy~expMerge + (1|labID), REML=FALSE)
exp3 <- lmer(data=Loadno0, loggenomeCopy~spp + (1|labID), REML=FALSE)

anova(full3, main3)
anova(full3, null3)
anova(main3, spp3)
anova(main3, exp3)

<<<<<<< HEAD
Anova(main3)

summary(full3)




# post hoc
ht3 = glht(full3, mcp(spp="Tukey"))
summary(ht3)
Anova(main3)
=======
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613





<<<<<<< HEAD
# ++++++++++++++++++++++++++++ prev +++++++++++++++++++++++++++++
# visitation rates and duration of visits on viruses loads and prevalence
full6 <- glmer(data=ModDat, BINYprefilter~visits*spp + (1|labID), family = binomial(link = "logit"))
main6 <- glmer(data=ModDat, BINYprefilter~visits+spp + (1|labID), family = binomial(link = "logit"))
anova(full6, main6)
summary(full6)


full7 <- glmer(data=ModDat, BINYprefilter~duration*spp + (1|labID), family = binomial(link = "logit"), nAGQ = 0)
main7 <- glmer(data=ModDat, BINYprefilter~duration+spp + (1|labID), family = binomial(link = "logit"))
anova(full7, main7)
summary(full7)















=======
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613


# visitation rates and duration of visits on viruses loads and prevalence




full4 <- glmer(data=ModDat,BINYprefilter~visits*foragetime + (1|labID), family = binomial(link = "logit"))
main4 <- glmer(data=ModDat,BINYprefilter~visits+foragetime + (1|labID), family = binomial(link = "logit"))
vis4 <- glmer(data=ModDat,BINYprefilter~foragetime + (1|labID), family = binomial(link = "logit"))
time4 <- glmer(data=ModDat,BINYprefilter~visits + (1|labID), family = binomial(link = "logit"))

anova(full4, main4)
anova(main4, vis4)
anova(main4, time4)

<<<<<<< HEAD
summary(full4)


main <- glmer(data=ModDat,BINYprefilter~visits+foragetime + (1|labID), family = binomial(link = "logit"))

noForage <- glmer(data=ModDat,BINYprefilter~visits + (1|labID), family = binomial(link = "logit"))

noVis <- glmer(data=ModDat,BINYprefilter~ foragetime + (1|labID), family = binomial(link = "logit"))

anova(main, noForage)
anova(main, noVis)







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

=======
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613

full5 <- lmer(data=Loadno0, loggenomeCopy~ visits * foragetime + (1|labID), REML = FALSE)
main5 <- lmer(data=Loadno0, loggenomeCopy~ visits + foragetime + (1|labID), REML = FALSE)
vis5 <- lmer(data=Loadno0, loggenomeCopy~ foragetime + (1|labID), REML = FALSE)
time5 <- lmer(data=Loadno0, loggenomeCopy~ visits + (1|labID), REML = FALSE)

anova(full5, main5)
anova(main5, vis5)
anova(main5, time5)



fig1 <- ddply(Loadno0, c("expMerge", "target_name"), summarise, 
              n = length(loggenomeCopy),
              mean = mean(loggenomeCopy, na.rm=TRUE),
              sd = sd(loggenomeCopy, na.rm=TRUE),
              se = sd / sqrt(n))






str(Loadno0)
# effect of species: spp: p = 0.01266

max(Loadno0$genomeCopy)
# 554278.8

min(Loadno0$genomeCopy)
# 104.6

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
LoadPlot <- ggplot(Loadno0, aes(x=target_name, y=log10(genomeCopy), fill=spp)) +
  labs(x=NULL, y = "Log(virus load)")+
  theme_classic(base_size = 20) +  
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + geom_dotplot(binaxis='y', stackdir = 'center', dotsize = 0.5, position = position_dodge()) + scale_fill_manual(values=colors, name="Plant Species:", labels=c(expression(italic("L. corniculatus"),italic("T. pratense"), italic("T. repens")))) + coord_cartesian(ylim = c(0, 6)) + theme(legend.position = c(.2, .2))
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
View(HBVL)

# HB: DWV- 10^9, 
#     BQCV- 10^6
#NYHB: DWV - 10^4
#     BQCV: 10^8


########### Analysis visitation time






VD_diversity <- VideoData[VideoData$Experiment=="diversity",]
VD_diversity$ExpComp <- c(rep("Diversity", length(VD_diversity$PlantSpp)))

hist(log10(VD_diversity$Forage))

mod <- aov(log10(VD_diversity$Forage)~VD_diversity$PlantSpp)
TukeyHSD(mod)

boxplot(log10(VD_diversity$Forage)~VD_diversity$PlantSpp)







VD_acute <- VideoData[!VideoData$Experiment==c("diversity"),]
VD_acute$ExpComp <- c(rep("Single species", length(VD_acute$PlantSpp)))

hist(log10(VD_acute$Forage))
VD_acute$PlantSpp <- as.factor(VD_acute$PlantSpp)
mod <- aov(data=VD_acute, log10(Forage+1)~PlantSpp)
summary(glht(mod, mcp(PlantSpp="Tukey")))
summary(mod)
boxplot(log10(VD_acute$Forage)~VD_acute$PlantSpp)








VD <- rbind(VD_acute, VD_diversity)
VD$PlantSpp <- as.factor(VD$PlantSpp)

mod <- aov(data=VD, log10(Forage+1)~PlantSpp*ExpComp)
summary(mod)





plantSpp <- ddply(VD, c("PlantSpp", "ExpComp"), summarise, 
                  n = length(Forage),
                  mean = mean(Forage, na.rm=TRUE),
                  sd = sd(Forage, na.rm=TRUE),
                  se = sd / sqrt(n))

plantSpp$PlantSpp <- as.factor(plantSpp$PlantSpp)

RC <- expression(paste(italic("T. pratense")))
WC <- expression(paste(italic("T. repens")))
BFT <- expression(paste(italic("L. corniculatus")))


time <- c(3.1, 102, 228.6, 8.3, 36, 40.3)
ExpComp <- c(rep("Single species", 3), rep("Diversity", 3))
PlantSpp <- c(rep(c("L. corniculatus", "T. pratense", "T. repens"), 2))
DFplant <- data.frame(time, ExpComp, PlantSpp)

DFplant$trans <- DFplant$time/75


ggplot(VD, aes(x = PlantSpp, y = log10(Forage+1), color = ExpComp)) +
  geom_boxplot(outlier.shape=16,
               outlier.size=2, notch=FALSE) + theme_bw(base_size = 15) + labs(x="Plant spp.", y = "Log10(visit duration (s))") + scale_color_manual(values = c("slategrey", "black"), name = "Experiment") + coord_cartesian(ylim = c(0, 4)) + theme(axis.line = element_line(colour = "black"), legend.position=c(.17, .87)) + annotate(geom = "text", x = .8, y = 2.3, label = "AB",cex = 6) + annotate(geom = "text", x = 1.8, y = 3.3, label = "A",cex = 6) +  annotate(geom = "text", x = 2.8, y = 3.5, label = "B",cex = 6) + annotate(geom = "text", x = 1.2, y = 2.3, label = "AB",cex = 6) + annotate(geom = "text", x = 2.2, y = 3.3, label = "A",cex = 6) +  annotate(geom = "text", x = 3.2, y = 3.5, label = "B",cex = 6) + scale_x_discrete(labels=c("BFT" = BFT, "RC" = RC, "WC" = WC)) + scale_y_continuous(sec.axis = sec_axis(~.*75, name = "visitation rate (visit/hr)")) + theme(axis.line.y.right = element_line(color = "red"), axis.ticks.y.right = element_line(color = "red"), axis.text.y.right = element_text(color = "red")) + geom_segment(aes(x = .73, xend=.9, y=.111,yend=.111), color = "red", size=1) + geom_segment(aes(x = 1.1, xend=1.28, y=.041,yend=.041), color = "red", size=1) + geom_segment(aes(x = 1.73, xend=1.9, y=.48,yend=.48), color = "red", size=1) + geom_segment(aes(x = 2.1, xend=2.28, y=1.36,yend=1.36), color = "red", size=1) + geom_segment(aes(x = 2.73, xend=2.9, y=.537,yend=.537), color = "red", size=1) + geom_segment(aes(x = 3.1, xend=3.28, y=3.048,yend=3.048), color = "red", size=1)









<<<<<<< HEAD
plot(y=Loadno0$loggenomeCopy, x=Loadno0$visits)

=======
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613
