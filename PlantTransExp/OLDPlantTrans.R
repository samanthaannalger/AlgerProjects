# Samantha A. Alger
# 21, January 2017
# Plant Transmission Experiment


#------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerCollaborations/PlantTransExp")

# read in data:
PlantVirus <- read.table("PlantTransPlantVirus.csv", header=TRUE, sep = ",") 

#------------------------------------------------------------------------
# data cleaning whip dataframe into shape and merge with eco data:
# remove uneeded columns from DF:
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

# group presence-absence by treatment and Virus:

DWV <- ddply(PlantVirus, c("C.T"), summarise, 
                 n = length(DWV),
                 mean = mean(DWV, na.rm=TRUE),
                 sd = sd(DWV, na.rm=TRUE),
                 se = sd / sqrt(n))

DWV <- DWV[-2,]

BQCV <- ddply(PlantVirus, c("C.T"), summarise,
             n = length(BQCV),
             mean = mean(BQCV, na.rm=TRUE),
             sd = sd(BQCV, na.rm=TRUE),
             se = sd / sqrt(n))

BQCV <- BQCV[-2,]

Virus <- c(rep("DWV", 2), rep("BQCV", 2))
PlantsFig <- rbind(DWV, BQCV)
PlantsFig <- cbind(PlantsFig, Virus)
PlantsFig$C.T <- ordered(PlantsFig$C.T, levels = c("Treatment", "Control", "Pre"))

#------------------------------------------------------------------------

#choosing color pallet
colors <- c("slategray3", "dodgerblue4")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(PlantsFig, aes(x=C.T, y=mean, fill=Virus)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) + labs(x="", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Virus:", labels=c("BQCV", "DWV")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 0.3)) + scale_x_discrete(labels=c("Honey Bees + Bumble Bees", "Bumble Bees Only")) + scale_y_continuous(labels = scales::percent)

##scale_x_discrete() is awesome for changing axis labels
#------------------------------------------------------------------------

# Dividing up by flower Spp. now:


DWV <- ddply(PlantVirus, c("C.T", "Spp"), summarise, 
             n = length(DWV),
             mean = mean(DWV, na.rm=TRUE),
             sd = sd(DWV, na.rm=TRUE),
             se = sd / sqrt(n))

DWV <- DWV[-c(1:6),]

BQCV <- ddply(PlantVirus, c("C.T", "Spp"), summarise,
              n = length(BQCV),
              mean = mean(BQCV, na.rm=TRUE),
              sd = sd(BQCV, na.rm=TRUE),
              se = sd / sqrt(n))

BQCV <- BQCV[-c(1:6),]

Virus <- c(rep("DWV", 3), rep("BQCV", 3))
PlantsFig1 <- rbind(DWV, BQCV)
PlantsFig1 <- cbind(PlantsFig1, Virus)
#------------------------------------------------------------------------

#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(PlantsFig1, aes(x=Virus, y=mean, fill=Spp)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position=c(.8, .85)) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

