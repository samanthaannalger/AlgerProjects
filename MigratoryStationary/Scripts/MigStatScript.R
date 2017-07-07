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
# Read in Virus Data:
MigVirus <- read.table("Data/MigVirus.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 


# Read in Nosema/Varroa/Eco Data:
MigStat <- read.table("Data/MigratoryStationaryData.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 


###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# source my functions
source("Scripts/BurnhamFunctions.R")
library(plyr)
library(ggplot2)
library(dplyr)

#--------------------------------------------------------------------------
# data cleaning

# preliminary data cleaning
MigVirus <- PrelimClean(data=MigVirus) 

# take only columns that we want:
data <- select(MigStat, sample_name, Treatment, dil.factor)

# merge with dilution factors
MigVirus <- merge(MigVirus, data, by=c("sample_name"), all.x =TRUE)

# use dilution factors to calcualte normalized virus load
MigVirus <- VirusNorm(data=MigVirus, number_bees = 50)

# use actin to normalize normalized viral load
MigVirus <- actinNormal(data=MigVirus)

# remove actin from data frame
MigVirus <- MigVirus[!(MigVirus$target_name=="ACTIN"),]

# adds virus binary data and makes norm genome copy 0 if above threashold CT
MigVirus <- CT_Threash(data=MigVirus)

# create full ID to check for inccorrect duplicates 
MigVirus$fullID <- with(MigVirus, paste0(sample_name, target_name))

# remove duplicates for examples with Cts above threashold (i.e. both coerced to 0s) 
MigVirus <- MigVirus[!duplicated(MigVirus$fullID), ] 

# merge the virus data to the main dataset
MigStat <- VirusMerger2000(data1 = MigVirus, data2 = MigStat)

# create average Nosema Load between chambers
MigStat$NosemaLoad <- (MigStat$NosemaChamber1 + MigStat$NosemaChamber2)/2

# normalize varroa load by number of bees sampled
MigStat$Varroa <- (MigStat$VarroaLoad / MigStat$TotBees) * 100

# log transform virus data:
MigStat$logDWV <- log(MigStat$DWVload + 1)
MigStat$logBQCV <- log(MigStat$BQCVload + 1)

# data analysis for DWV (repeated measures ANOVA)
RepANOVA(data=MigStatPositiveDWV, column="logDWV")

# data analysis for BQCV (repeated measures ANOVA)
RepANOVA(data=MigStatPositiveBQCV, column="logBQCV")

# data analysis for FOB (repeated measures ANOVA)
RepANOVA(data=MigStat, column="FOB")

# data analysis for Varroa (repeated measures ANOVA)
RepANOVA(data=MigStat, column="Varroa")

# data analysis for Varroa (repeated measures ANOVA)
RepANOVA(data=MigStat, column="NosemaLoad")

# data analysis for Varroa (repeated measures ANOVA)
RepANOVA(data=MigStat, column="BroodPattern")


# data analysis for Varroa (repeated measures ANOVA)
RepANOVA(data=MigStat, column="DWVbinary")



# substetting only positive BQCV bees
MigStatPositiveBQCV <- subset(MigStat, BQCVbinary==1)

# substetting only positive BQCV bees
MigStatPositiveDWV <- subset(MigStat, DWVbinary==1)




VirusSum1 <- ddply(MigStatPositiveDWV, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(logDWV),
                   mean = mean(logDWV, na.rm=TRUE),
                   sd = sd(logDWV, na.rm = TRUE),
                   se = sd / sqrt(n))






































ggplot(data = Summary, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   group = Treatment)
) + geom_point(size=3) + scale_colour_manual(values = c("black", "blue")) + labs(x = "Sampling Event", y = "Nosema Load") + coord_cartesian(ylim = c(0, 30), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.1, .85), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Queen Origin")


