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

View(VideoDataMerge)

#remove NAs from the data frame so that we are only working with Acute and Chronic experiments
VideoDataMerge <- VideoDataMerge[complete.cases(VideoDataMerge), ]

table(VideoDataMerge$detect,VideoDataMerge$spp)

#Creating a figure of foraging time and visitation visits by plant species showing % detection.

#Separate by virus first
VideoDWV <- VideoDataMerge[ which(VideoDataMerge$virus=="DWV"), ]
VideoBQCV <- VideoDataMerge[ which(VideoDataMerge$virus=="BQCV"), ]

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


# Use brewer color palettes
p+scale_color_brewer(palette="Dark2")
# Use grey scale
p + scale_color_grey()

