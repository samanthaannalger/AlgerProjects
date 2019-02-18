# S. Alger, P Alexander Burnham
# Estimating Monarch population at BorderView Farm
# July 25, 2018

###### SET UP DATA SET AND VARS FOR CALCULATIONS: #######
MWm2 <- 24.25 # density of plants per 1 meter squared
plotArea <- 137088 # area of plot in meters squared
plantsSE <- 300 # number of plants sampled at each time step
eggCounts <- c(0,0,0,1,1,3,0,0,0,0,0,0,0,0) # vector of counts of eggs at each of 14 time steps

# calculate total number of plants at BV:
plantsInPlot <- MWm2 * plotArea

# calculate the number of estimated eggs in plot for each time step: 
totalEggsVec <- (eggCounts/plantsSE) * plantsInPlot 

# sum them up for total eggs at BV so far:
sumeggs<-sum(totalEggsVec)

# with 10% survivorship to adulthood:
sumeggs*.10


# Clear memory of characters:
ls()
rm(list=ls())

# source my packages
library(plyr)
library(ggplot2)
library(reshape2)

# Set Working Directory: 
setwd("~/AlgerProjects/Milkweed/Data")

# read in data:
Monarch <- read.csv("MonarchData.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

head(Monarch)

MonarchStats <- ddply(Monarch, c("Site", "Date"), summarise, 
                  egg = sum(Total_Eggs, na.rm = TRUE),
                  larvae = sum(Total_Larvae, na.rm = TRUE),
                  adults = sum(TotalAdult, na.rm = TRUE))

MonarchLong <- melt(MonarchStats, id.vars = c("Site", "Date"))

MonarchLong <- MonarchLong[order(MonarchLong$Date),]

MonarchLong$Date <- as.character(MonarchLong$Date)


p <- ggplot(data=MonarchLong, aes(x=Date, y=value, color=variable)) +
  geom_line(aes(linetype = Site), size=1.5) + 
  geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Number of Observations") + 
  theme_minimal(base_size = 15)

p + theme(legend.position="top", legend.box = "horizontal")


x <- split(MonarchLong, MonarchLong$Site)
Borderview <- x$Borderview 
Dewing <- x$Dewing 

b <- ggplot(Borderview, aes(x=Date, y=value, fill=variable)) + 
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = c(.8, .8), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title="Borderview", x ="Date", y = NULL, fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 15))
b


d <- ggplot(Dewing, aes(x=Date, y=value, fill=variable)) + 
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position =c(3,3)) +
  labs(title="Dewing", x ="Date", y = "Number of Observations", fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 15))
d

library(patchwork)
d+b
