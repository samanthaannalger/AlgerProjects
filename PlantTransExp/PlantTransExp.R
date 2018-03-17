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

# read in data, (skip over meta data in csv file: skip = 9) :
VideoData <- read.csv("PlantTransVideoData.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE, skip = 9)

PlantDF <- read.csv("plantTransPlantsDF.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

visitations <- table(VideoData$expID)
