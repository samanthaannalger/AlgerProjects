###########################################################################################
# Functions for plant trans
# Samantha Alger
# June 30, 2018

###########################################################################################
# BEGIN FUNCTIONS #########################################################################
###########################################################################################

###########################################################################
# function name: 
# description:conducts a repeated measures ANOVA on a data set and plots 
# parameters: 
# data = data frame, in this case requires columns for SamplingEvent, ID 
# and Treatment plus a column for yvar 
# column = name of continuous variable column in quotations
###########################################################################

PlantVirus <- function(data, column, key=10){
  
  data$x <- data[,column]
  
  # repeated measures on continuous variable
  aov.out <- aov(x~Treatment * SamplingEvent + Error(ID), data=data)
  
  
  #PLOT:
  library(plyr)
  
  #Checking out by plant species for each experiment:
  plantSpp <- ddply(data, c("target_name", "spp"), summarise, 
                    n = length(BINYprefilter),
                    mean = mean(BINYprefilter, na.rm=TRUE),
                    sd = sd(BINYprefilter, na.rm=TRUE),
                    se = sd / sqrt(n))
  
  
  plantSpp$spp[plantSpp$spp == "RC"] <- "Red Clover"
  plantSpp$spp[plantSpp$spp == "BFT"] <- "BirdsFoot Trefoil"
  plantSpp$spp[plantSpp$spp == "WC"] <- "White Clover"
  
  # What's the percent detection for each plant species for each virus?
  plantAcuteDWV <- data[ which(data$target_name=="DWV"), ]
  plantAcuteBQCV <- data[ which(data$target_name=="BQCV"), ]
  
  tapply(plantAcuteDWV$BINYprefilter, plantAcuteDWV$spp, mean)
  tapply(plantAcuteBQCV$BINYprefilter, plantAcuteBQCV$spp, mean)
}