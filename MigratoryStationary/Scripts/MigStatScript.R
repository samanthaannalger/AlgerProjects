###########################################################################################
# Data Analysis for Migratory Stationary Study 
# P. Alexander Burnham
# May 25, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/MigratoryStationary/")

###########################################################################################

# Read in Data:
MigStat <- read.table("Data/MigratoryStationaryData.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################################

MigStat$NosemaLoad <- (MigStat$NosemaChamber1 + MigStat$NosemaChamber2)/2




aov.out <- aov(NosemaLoad ~ Treatment * SamplingEvent + Error(ID), data=MigStat)
summary(aov.out)

library(plyr)
library(ggplot2)

# summary stats for plotting purposes:
Summary <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                     n = length(NosemaLoad),
                     mean = mean(NosemaLoad, na.rm = TRUE),
                     sd = sd(NosemaLoad, na.rm = TRUE),
                     se = sd / sqrt(n))

# summary stats for plotting purposes:
Summary <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                 n = length(VarroaLoad),
                 mean = mean(VarroaLoad, na.rm = TRUE),
                 sd = sd(VarroaLoad, na.rm = TRUE),
                 se = sd / sqrt(n))

# summary stats for plotting purposes:
Summary <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                 n = length(FOB),
                 mean = mean(FOB, na.rm = TRUE),
                 sd = sd(FOB, na.rm = TRUE),
                 se = sd / sqrt(n))


Summary <- Summary[-c(1,2,3,6,9),]

#Create plot in ggplot for mass without x label



ggplot(data = Summary, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   group = Treatment)
) + geom_point(size=3) + scale_colour_manual(values = c("black", "blue")) + labs(x = "Sampling Event", y = "Nosema Load") + coord_cartesian(ylim = c(0, 30), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.1, .85), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Queen Origin")

aov.out <- aov(FOB ~ Treatment * SamplingEvent + Error(ID), data=MigStat)
summary(aov.out)
