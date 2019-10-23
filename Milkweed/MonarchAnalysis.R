# S. Alger, P Alexander Burnham
# Estimating Monarch population at BorderView Farm
# July 25, 2018

# Clear memory of characters:
ls()
rm(list=ls())

# source my packages
library(plyr)
library(ggplot2)
library(reshape2)

# Set Working Directory: 
setwd("~/Documents/GitHub/AlgerProjects/Milkweed/Data")

# read in data:
Monarch <- read.csv("MonarchData.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 
Bees <- read.csv("1819beeObs.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 


####### MONARCH 2018 #####################################################################

MonarchStats <- ddply(Monarch, c("Site", "Date"), summarise, 
                  egg = sum(Total_Eggs, na.rm = TRUE),
                  larvae = sum(Total_Larvae, na.rm = TRUE),
                  adults = sum(TotalAdult, na.rm = TRUE))

MonarchStats$Date <- as.Date(MonarchStats$Date, "%m/%d/%y")

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = c(3, 3), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title="Borderview", x ="Date", y = NULL, fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 20))
b


d <- ggplot(Dewing, aes(x=Date, y=value, fill=variable)) + 
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position =c(.2,.8)) +
  labs(title="Dewing", x ="Date", y = "Number of Observations", fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 20))
d

library(patchwork)
d+b


####### BEES 2018 ########################################################################

Bees18 <- Bees[Bees$Year==2018,]

BeeStats <- ddply(Bees18, c("Site","Flower"), summarise, 
                  'Bombus-queen' = sum(BigBombus, BigOrange, na.rm = TRUE),
                  'Bombus-worker' = sum(SmallBombus, SmallOrange, na.rm = TRUE),
                  'Black-big' = sum(BigBlack, na.rm = TRUE),
                  'Black-slender' = sum(SmallBlack, na.rm = TRUE),
                  'Black-tiny' = sum(TinyBlack, na.rm = TRUE),
                  'Green' = sum(Green, na.rm = TRUE),
                  'Apis' = sum(Apis, na.rm = TRUE),
                  'Flies' = sum(Flies, na.rm = TRUE))



library(reshape2)
BeeStats <- melt(BeeStats, id.vars = c("Site", "Flower"))

BeeStats$logValue <- log10(BeeStats$value + 1)

# Figure for bees at Borderview:
BeeBV <- BeeStats[ which(BeeStats$Site=='Borderview'), ]

# Figure for bees at Dewing:
BeeDew <- BeeStats[ which(BeeStats$Site=='Dewing'), ]

#choosing color pallet
colors <- c("steelblue", "grey30")

plot1 <- ggplot(data=BeeBV, aes(x=variable, y=logValue, fill = Flower)) + geom_bar(stat = 'identity', position = 'dodge') 
plot1 + theme_bw(base_size = 21) + scale_fill_manual(values=colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Morphotype", y = "Log10(Observations)") + ggtitle('Borderview') + ylim(0, 3)


plot2 <- ggplot(data=BeeDew, aes(x=variable, y=logValue, fill = Flower)) + geom_bar(stat = 'identity', position = 'dodge') 
plot2 + theme_bw(base_size = 21) + scale_fill_manual(values=colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "Morphotype", y = "Log10(Observations)") + ggtitle('Dewing') + ylim(0, 3)



####### BEES 2019 ########################################################################

BeeStatsYear <- ddply(Bees, c("Site","Year"), summarise, 
                  'Bombus-queen' = sum(BigBombus, BigOrange, na.rm = TRUE),
                  'Bombus-worker' = sum(SmallBombus, SmallOrange, na.rm = TRUE),
                  'Black-big' = sum(BigBlack, na.rm = TRUE),
                  'Black-slender' = sum(SmallBlack, na.rm = TRUE),
                  'Black-tiny' = sum(TinyBlack, na.rm = TRUE),
                  'Green' = sum(Green, na.rm = TRUE),
                  'Apis' = sum(Apis, na.rm = TRUE))

BeeStatsYear <- melt(BeeStatsYear, id.vars = c("Site", "Year"))

BeeStatsYear$logValue <- log10(BeeStatsYear$value + 1)

# Figure for bees at Borderview:
BeeBVYear <- BeeStatsYear[ which(BeeStatsYear$Site=='Borderview'), ]
BeeBVYear$Year <- as.character(BeeBVYear$Year)

# Figure for bees at Dewing:
BeeDewYear <- BeeStatsYear[which(BeeStatsYear$Site=='Dewing'), ]
BeeDewYear$Year <- as.character(BeeDewYear$Year)



plot3 <- ggplot(data=BeeBVYear, aes(x=variable, y=logValue, fill = Year)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw(base_size = 18) + scale_fill_manual(values=colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position =c(.2,.8),plot.title = element_text(hjust = 0.5)) +
  labs(x = "Morphotype", y = "Log10(Observations)") + 
  ggtitle('Borderview') + 
  ylim(0, 3.5)



plot4 <- ggplot(data=BeeDewYear, aes(x=variable, y=logValue, fill = Year)) + geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw(base_size = 18) + 
  scale_fill_manual(values=colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position =c(3,3),plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Morphotype", y = "Log10(Observations)") + 
  ggtitle('Dewing') + ylim(0, 3.5)

plot3+plot4



####### MONARCH 2019 ############################################################################

Monarch19 <- read.csv("2019MonarchDat.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

Monarch19$Date <- as.Date(Monarch19$Date, "%m/%d/%y")


Monarch19$larva <- Monarch19$Instar_1 + Monarch19$Instar_2 + Monarch19$Instar_3 + Monarch19$Instar_4 + Monarch19$Instar_5

MonarchStats19 <- ddply(Monarch19, c("Site", "Date"), summarise, 
                      egg = sum(Eggs, na.rm = TRUE),
                      larvae = sum(larva, na.rm = TRUE),
                      adults = sum(Adult_Unk, na.rm = TRUE))


MonarchLong <- melt(MonarchStats19, id.vars = c("Site", "Date"))
MonarchLong <- MonarchLong[order(MonarchLong$Date),]
MonarchLong$Date <- as.character(MonarchLong$Date)


x <- split(MonarchLong, MonarchLong$Site)
Borderview <- x$Borderview 
Dewing <- x$Dewing 

b <- ggplot(Borderview, aes(x=Date, y=value, fill=variable)) + 
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = c(3, 3), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title="Borderview", x ="Date", y = NULL, fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 80))



d <- ggplot(Dewing, aes(x=Date, y=value, fill=variable)) + 
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position =c(.2,.8)) +
  labs(title="Dewing", x ="Date", y = "Number of Observations", fill = "Life Stage") + 
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 80))


d+b



Borderview







AphidDat <- ddply(Monarch19, c("Site", "Date", "Block", "PlantBloomingState"), summarise, 
                      Aphids = sum(Aphids, na.rm = TRUE))

AphidDat$Date <- as.character(AphidDat$Date)
Bord <- AphidDat[AphidDat$Site=="Borderview",]
Dew <- AphidDat[AphidDat$Site=="Dewing",]





plot5 <- ggplot(Dew[1:20,], aes(x=Date, y=log10(Aphids+1), color=Block, group=Block)) + 
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position ="top") +
  labs(title="Dewing", x ="Date", y = "log10(aphids/50 plants)", color = "Block:") + 
  geom_point(stat="identity", size=5) +
  geom_line() +
  coord_cartesian(ylim = c(0, 4)) + 
  geom_segment(mapping=aes(x=1, y=3, xend=4, yend=3), arrow=arrow(), size=1, color="black") +
  geom_segment(mapping=aes(x=4, y=3, xend=7, yend=3), arrow=arrow(), size=1, color="black") +
  annotate("text", x = 2.5, y = 3.5, label = "Pre-Bloom", size=7, color="black") +
  annotate("text", x = 5.5, y = 3.5, label = "Bloom", size=7, color="black")

plot6 <- ggplot(Bord[1:30,], aes(x=Date, y=log10(Aphids+1), color=Block, group=Block)) + 
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position ="top") +
  labs(title="Borderview", x ="Date", y = "log10(aphids/50 plants)", color = "Block:") + 
  geom_point(stat="identity", size=5) +
  geom_line() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_segment(mapping=aes(x=1, y=3, xend=4, yend=3), arrow=arrow(), size=1, color="black") +
  geom_segment(mapping=aes(x=4, y=3, xend=6, yend=3), arrow=arrow(), size=1, color="black") +
  annotate("text", x = 2.5, y = 3.5, label = "Pre-Bloom", size=7, color="black") +
  annotate("text", x = 5, y = 3.5, label = "Bloom", size=7, color="black")

plot5+plot6




####### Defoliation 2019 ############################################################################
Monarch19$PercDefoliation


Defoliation <- ddply(Monarch19, c("Site", "Date"), summarise, 
                        Defoliation = mean(PercDefoliation, na.rm = TRUE),
                        SD = sd(PercDefoliation, na.rm = TRUE),
                        N = length(PercDefoliation),
                        SE = SD/sqrt(N)
                     )

Defoliation$Date <- as.character(Defoliation$Date)

ggplot(Defoliation, aes(x=Date, y=Defoliation, color=Site, group=Site)) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_line(data=Defoliation[!is.na(Defoliation$Defoliation),], size=1.2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin = Defoliation - SE, ymax = Defoliation + SE, width = .7))+
  labs(x ="Date", y = "% Defoliation", color = "Site:") +
  scale_color_manual(values = colors)+
  theme(legend.position="top")




hist(Monarch19$PercDefoliation, breaks = 30, 
     xlab="% Defoliation", 
     main="Frequency % Defoliation",
     col="steelblue")

hist(Monarch19$Aphids, breaks = 30, 
     xlab="Number of Aphids", 
     main="Frequency of Aphid Abundance",
     col="grey")

