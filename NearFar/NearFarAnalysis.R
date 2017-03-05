# P. Alexander Burnham
# 25, Feb. 2017
# qPCR Results
# Data cleaning and Analysis:

##################################################################################-
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/NearFar")

# read in data:
NFdata <- read.table("NearFar_qPCR_results.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

DilDF <- read.table("DilutionFactors.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

MeltTemp <- read.table("MeltTemp.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

##################################################################################-
# load packages:

# remove uneeded columns from DF:
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

##################################################################################-
# Calculate genome copies per bee and remove actin from data set:

# merge dil factors with main DF
NFdata <- merge(NFdata, DilDF, by="ID")

# set constant values for genome copies per bee calculations:
crude_extr <- 100
eluteRNA <- 50
total_extr_vol <- 600
number_bees <- 1
cDNA_eff <- 0.1
rxn_vol <- 3


# create column for genome copies per bee:
NFdata$genome_copbee <- ((((((NFdata$QuantityMean / cDNA_eff) / rxn_vol) * NFdata$dilFactor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees

# pull only actin values out of dataframe
ActinOnly <- NFdata[which(NFdata$Virus=="ACTIN"),]

# create DF of ACTIN genome copies and lab ID:
ActinDF <- data.frame(ActinOnly$ID, ActinOnly$genome_copbee)
colnames(ActinDF) <- c("ID", "ACT_genome_copbee")

# merge ACTIN dataframe with main dataframe:
#Need rownames and all.x=TRUE because data frames are different sizes.
NFdata <- merge(NFdata, ActinDF, by="ID")

# find mean of all ACTIN values:
ActinMean <- mean(ActinOnly$genome_copbee, na.rm = TRUE)

# create column for normalized genome copies per bee:
NFdata$norm_genome_copbee <- (NFdata$genome_copbee/NFdata$ACT_genome_copbee)*ActinMean

# Remove actin from main dataframe (NFdata)
NFdata <- NFdata[-which(NFdata$Virus=="ACTIN"),]

# merge Melt Temps with main data frame by ID and Virus
NFdata <- merge(NFdata, MeltTemp, by=c("ID","Virus"), all.x=TRUE, all.y=FALSE)
##################################################################################-

# this function puts NAs in vectors where values parameters are not met then removes those columns, cleaning the data frame:

Burnhams_Fabulous_Data_Cleanser <- function(dataframe=NFdata, variable=NFdata$Virus){
  
  splitDF <- split(NFdata, NFdata$Virus)
  
  # filtering out for CT values (bombus) that are outside limit of detection, assigning either a 1 or 0 in a new column
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$CqMean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$CqMean > 32.525, 0, 1)
  splitDF$IAPV$virusBINY  <- ifelse(splitDF$IAPV$CqMean > 30.796, 0, 1)
  
  x <- rbind(splitDF$DWV,splitDF$BQCV,splitDF$IAPV)
  
  # make NA values 0
  x$virusBINY[is.na(x$virusBINY)] <- 0
  
  return(x)
  }


# END OF FUNCTION

# calling Burnhams_Fabulous_Data_Cleanser created above and using it on NFdata data frame: 
NFdata <- Burnhams_Fabulous_Data_Cleanser(NFdata)

##################################################################################-
# create new variable called viralLoad that makes values 0 when they fail the virusBINY logical (0 or 1)
viralLoad <- NFdata$norm_genome_copbee
NFdata <- cbind(viralLoad, NFdata)

NFdata$viralLoad[NFdata$virusBINY == 0] <-0

# create log base 10 viral load column:

NFdata$LOGviralLoad <- log10(NFdata$viralLoad)
NFdata$LOGviralLoad[NFdata$LOGviralLoad == -Inf] <-0
head(NFdata)
x <- NFdata[which(NFdata$LOGviralLoad>0),]

# write out data file as .csv
write.csv(x=NFdata, file="NFdata.csv")
##################################################################################-



###########_____DATA___ANALYSIS______############

hist(NFdata$LOGviralLoad)
#----------------------------------------------------------------
# summary stats for plotting purposes:
VirusSummary <- ddply(x, c("Virus", "treatment"), summarise, 
                      n = length(LOGviralLoad),
                      mean = mean(LOGviralLoad, na.rm = TRUE),
                      sd = sd(LOGviralLoad, na.rm = TRUE),
                      se = sd / sqrt(n))
#subsetting the data to remove NA and IAPV (no viral loads)
VirusSummary <- VirusSummary[-c(3,6),]

print(VirusSummary)
#----------------------------------------------------------------

#choosing color pallet
colors <- c("slategray3", "dodgerblue4")

plot2 <- ggplot(VirusSummary, aes(x=Virus, y=mean, fill=treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Virus", 
                                                    y = "Viral Load log(genome copies)")

plot2 + theme_minimal(base_size = 17) + theme(legend.position=c(.85, .85)) + coord_cartesian(ylim = c(0, 8)) + scale_fill_manual(values=colors) 

xsplit <- split(x, x$Virus)
t.test(xsplit$BQCV$LOGviralLoad~xsplit$BQCV$treatment)
t.test(xsplit$DWV$LOGviralLoad~xsplit$DWV$treatment)

#-----------------------------------------------------------------------------------
# looking at prevalence:

#----------------------------------------------------------------
# summary stats for plotting purposes:
PrevSummary <- ddply(NFdata, c("Virus", "treatment"), summarise, 
                      n = length(virusBINY),
                      mean = mean(virusBINY, na.rm = TRUE),
                      sd = sd(virusBINY, na.rm = TRUE),
                      se = sd / sqrt(n))
#subsetting the data to remove NA and IAPV (no viral loads)
PrevSummary <- PrevSummary[-c(3,6,7,8,9),]

print(PrevSummary)
#----------------------------------------------------------------

#choosing color pallet
colors <- c("slategray3", "dodgerblue4")

plot2 <- ggplot(PrevSummary, aes(x=Virus, y=mean, fill=treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Virus", 
                                                    y = "Virus Prevalence")

plot2 + theme_minimal(base_size = 17) + theme(legend.position=c(.85, .9)) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values=colors) 

