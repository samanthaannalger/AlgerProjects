---
title: "PlantsTransExp"
author: "Samantha A. Alger"
date: "2/27/2018"
---
  
# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
setwd("~/AlgerProjects/PlantTransExp/CSV_Files")

library(ggplot2)

# read in data, (skip over meta data in csv file: skip = 9) :
VideoData <- read.csv("PlantTransVideoData.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE, skip = 9)

PlantDF <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

VirusDetect <- read.csv("VirusDetectionPlants.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

#Merge video data with the detection data:
VideoDataMerge <- merge(VideoData, VirusDetect, by=c("expID"), all.x=TRUE, all.y=TRUE)

# Create a table for visitation data
Visits <- table(VideoDataMerge$expID)
expID <- as.vector(names(Visits))
visits <- as.vector(Visits)

VidDat <- data.frame(expID, visits)

#Check out visitation data
visitations <- table(VideoData$expID)

tapply(VideoData$Forage, VideoData$expID, mean)
tapply(PlantDF$BINYprefilter, PlantDF$expID, mean)
mean = mean(BINYprefilter, na.rm=TRUE)

Anova <- aov(Forage~expID, data=VideoData)
summary(Anova)
TukeyHSD(Anova)

boxplot(Forage~expID, data = VideoData, ylim=c(0,130))

hist(VideoData$Forage, breaks = 150)

#remove NAs from the data frame
VideoDataMerge <- VideoDataMerge[complete.cases(VideoDataMerge), ]

table(VideoDataMerge$detect,VideoDataMerge$spp)

VideoSum <- ddply(VideoDataMerge, c("PlantSpp", "detect", "virus"), summarise, 
                  n = length(PlantSpp),
                  mean = mean(detect, na.rm=TRUE))

VideoSum

#subset to include acute and chronic only
#AcuteMerge <- VideoDataMerge[which(VideoDataMerge$Experiment=="chronic-1" | VideoDataMerge$Experiment =="acute" | VideoDataMerge$Experiment=="chronic-2"| VideoDataMerge$Experiment=="chronic-3"), ]
AcuteMerge <- VideoDataMerge[which(VideoDataMerge$Experiment =="acute"), ]

#subset to include diversity only

DivMerge <- VideoDataMerge[which(VideoDataMerge$Experiment=="diversity"), ]

#Creating a figure of foraging time and visitation visits by plant species showing % detection.

#Separate by virus first
VideoDWV <- AcuteMerge[ which(AcuteMerge$virus=="DWV"), ]
VideoBQCV <- AcuteMerge[ which(AcuteMerge$virus=="BQCV"), ]

DWVPlot <- ggplot(VideoDWV, aes(x=Forage, y=detect, color=spp, shape=factor(spp))) +
  labs(x="Foraging time (s)", y = "% DWV Detection")+
  theme_classic() +  
  geom_point(shape=1, size = 4,
             position=position_jitter(width=0.003,height=.003)) +
  scale_color_manual(values=c("goldenrod", "violetred4", "black"), name = "Plant Species") + 
  coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

BQCVPlot <- ggplot(VideoBQCV, aes(x=Forage, y=detect, color=spp, shape=factor(spp))) +
        labs(x="Foraging time (s)", y = "% BQCV Detection")+
        theme_classic() +  
        geom_point(shape=1, size = 4,
                   position=position_jitter(width=0.003,height=.003)) +
                    scale_color_manual(values=c("goldenrod", "violetred4", "black"), name = "Plant Species") +
        coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

DWVPlot
BQCVPlot


VideoSum <- ddply(AcuteMerge, c("PlantSpp", "detect", "virus"), summarise, 
                     n = length(PlantSpp),
                     mean = mean(detect, na.rm=TRUE))

VideoSum
#plot1<- ggplot(VideoDWVSum, aes(x=n, y=mean, fill=PlantSpp)) + geom_bar(stat="identity", position=position_dodge()) + labs(x=NULL, y="% plants with virus detected")

#plot1
###########FOR DIVERSITY

#Separate by virus first
VideoDWV <- DivMerge[ which(DivMerge$virus=="DWV"), ]
VideoBQCV <- DivMerge[ which(DivMerge$virus=="BQCV"), ]

DWVPlot <- ggplot(VideoDWV, aes(x=Forage, y=detect, color=spp, shape=factor(spp))) +
  labs(x="Foraging time (s)", y = "% DWV Detection")+
  theme_classic() +  
  geom_point(shape=1, size = 4,
             position=position_jitter(width=0.3,height=.3)) +
  scale_color_manual(values=c("goldenrod", "violetred4", "black"), name = "Plant Species") + 
  coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

BQCVPlot <- ggplot(VideoBQCV, aes(x=Forage, y=detect, color=spp, shape=factor(spp))) +
  labs(x="Foraging time (s)", y = "% BQCV Detection")+
  theme_classic() +  
  geom_point(shape=1, size = 4,
             position=position_jitter(width=0.003,height=.003)) +
  scale_color_manual(values=c("goldenrod", "violetred4", "black"), name = "Plant Species") +
  coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

DWVPlot
BQCVPlot


VideoDivSum <- ddply(DivMerge, c("PlantSpp","virus", "detect"), summarise, 
                     n = length(PlantSpp),
                     mean = mean(detect, na.rm=TRUE))

# Use brewer color palettes
p+scale_color_brewer(palette="Dark2")
# Use grey scale
p + scale_color_grey()

