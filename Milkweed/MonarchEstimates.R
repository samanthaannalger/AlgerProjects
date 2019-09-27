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
# setwd("~/AlgerProjects/Milkweed/Data")
setwd("~/Documents/GitHub/AlgerProjects/Milkweed/Data")


# read in data:
Monarch <- read.csv("MonarchData.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE)
Bees <- read.csv("BeeSurveyObs.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

head(Monarch)

MonarchStats <- ddply(Monarch, c("Site", "Date"), summarise, 
                  egg = sum(Total_Eggs, na.rm = TRUE),
                  larvae = sum(Total_Larvae, na.rm = TRUE),
                  adults = sum(TotalAdult, na.rm = TRUE))

MonarchLong <- melt(MonarchStats, id.vars = c("Site", "Date"))




BeeStats <- ddply(Bees, c("Site","Flower"), summarise, 
                      'Bombus-queen' = sum(BigBombus, BigOrange, na.rm = TRUE),
                      'Bombus-worker' = sum(SmallBombus, SmallOrange, na.rm = TRUE),
                      'Black-big' = sum(BigBlack, na.rm = TRUE),
                      'Black-slender' = sum(SmallBlack, na.rm = TRUE),
                      'Black-tiny' = sum(TinyBlack, na.rm = TRUE),
                      'Green' = sum(Green, na.rm = TRUE),
                      'Apis' = sum(Apis, na.rm = TRUE),
                      'Flies' = sum(Flies, na.rm = TRUE))

BeeStats


library(reshape2)
BeeStats <- melt(BeeStats, id.vars = c("Site", "Flower"))

BeeStats$logValue <- log10(BeeStats$value + 1)

# Figure for bees at Borderview:
BeeBV <- BeeStats[ which(BeeStats$Site=='Borderview'), ]



#choosing color pallet
colors <- c("steelblue", "grey30")

plot1 <- ggplot(data=BeeBV, aes(x=variable, y=logValue, fill = Flower)) + geom_bar(stat = 'identity', position = 'dodge') 
plot1 + theme_bw(base_size = 21) + scale_fill_manual(values=colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Morphotype", y = "Log10(Observations)") + ggtitle('Borderview') + ylim(0, 3)

# Figure for bees at Dewing:
BeeDew <- BeeStats[ which(BeeStats$Site=='Dewing'), ]

#choosing color pallet
colors <- c("steelblue", "grey30")

<<<<<<< HEAD
library(patchwork)
d+b



Dewing$value

Borderview$value



timePoints <- 13
blocks <- 3
Sites <- 2
plants <- 100


TotalEffortPerSite <- timePoints * blocks * plants

# Borderview
###########################################################################
Bord <- ddply(Borderview, c("variable"), summarise, 
              val = sum(value, na.rm = TRUE))

Bord$variable <- as.character(Bord$variable)
var <- c(Bord$variable, "total")
val <- c(Bord$val, 28)
Bord <- data.frame(var, val)
Bord$Normal <- Bord$val / TotalEffortPerSite
Bord$Double <- Bord$Normal*2
###########################################################################


# Dewing
###########################################################################
Dew <- ddply(Dewing, c("variable"), summarise, 
                      val = sum(value, na.rm = TRUE))

Dew$variable <- as.character(Dew$variable)
var <- c(Dew$variable, "total")
val <- c(Dew$val, 48)
Dew <- data.frame(var, val)
Dew$Normal <- Dew$val / TotalEffortPerSite
Dew$Double <- Dew$Normal*2
###########################################################################

(Dew$Normal[3]*TotalEffortPerSite)/Bord$Normal[3]





# Parameters:
BordDens <- 24.25 # density of plants at borderview per 1 meter squared
DewDens <- 12.42 # density of plants at dewing per 1 meter squared
BordArea <- 140000 # Borderview area in meters
DewArea <- 60000 # Dewing area in meters
numFlowers <- 300 # number of flowers sampled per site per time step

# vectors of egg counts at each time point for each site
DewingEggs <- Dewing[Dewing$variable=="egg",]
BorderviewEggs <- Borderview[Borderview$variable=="egg",]

# number of plants at each time point
numBord <- BordDens * BordArea
numDew <- DewDens * DewArea

# the sum of all estimated egg counts for each time step for each site
TotEggsEstimateDewing <- sum((DewingEggs$value/numFlowers) * numDew)
TotEggsEstimateBorderview <- sum((BorderviewEggs$value/numFlowers) * numBord)

TotEggsEstimateDewing 
TotEggsEstimateBorderview

# number of eggs at dewing over the course of the year: 59616
# number of eggs at borderview over the course of the year: 56583.33


=======
plot1 <- ggplot(data=BeeDew, aes(x=variable, y=logValue, fill = Flower)) + geom_bar(stat = 'identity', position = 'dodge') 
plot1 + theme_bw(base_size = 21) + scale_fill_manual(values=colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Morphotype", y = "Log10(Observations)") + ggtitle('Dewing') + ylim(0, 3)
>>>>>>> 5ee3c5086e4d4a5b59c0d6698871520473b1e613
