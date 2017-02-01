# Shannon Index (plants) 2015 Survey
# P. Alexander Burnham
# 7, November 2016


ls()
rm(list=ls())


setwd("~/AlgerCollaborations/2015_BombSurv_Plants")

mydata <- read.csv("shannon.csv",header=T, row.names="Site")
mydata2 <- read.csv("shannon.csv",header=T)
dense <- read.csv("density.csv",header=T)

Density <- dense$Density
Site <- as.character(mydata2$Site)


library(vegan)


ShannonDIV <- diversity(mydata, index = "shannon", MARGIN = 1, base = exp(1))
#ShannonDIV <- cbind(data.frame(ShannonDIV), Site)

SimpsonDIV <- diversity(mydata, index = "simpson", MARGIN = 1, base = exp(1))
#SimpsonDIV  <- cbind(data.frame(SimpsonDIV), Site)

FisherDIV <- fisher.alpha(mydata, MARGIN = 1)
#FisherDIV <- cbind(data.frame(FisherDIV), Site)

Diversity <- cbind(Site, ShannonDIV, SimpsonDIV, FisherDIV, Density) 
Diversity <- data.frame(Diversity)


write.csv(Diversity, file = "Diversity.csv")
