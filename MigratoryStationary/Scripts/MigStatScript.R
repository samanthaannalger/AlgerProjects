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


###########################################################################################





###############################################################################################
################################### FUNCTIONS #################################################
###############################################################################################


###########################################################################
# function name: VirusNorm
# description: normalizes virus data with actin values and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 50, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 100
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  data$total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################




actinNormal <- function(data=MigVirus){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$target_name=="ACTIN"),]

  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("sample_name", "ACT_genomeCopy")

  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("sample_name", "run"))

  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)

  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean

return(data)
}

actinNormal(data=MigVirus)












###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

source("Scripts/BurnhamFunctions.R")


MigVirus <- PrelimClean(data=MigVirus) 


RepANOVA(data=MigStat, column="Varroa")




# create average Nosema Load between chambers
MigStat$NosemaLoad <- (MigStat$NosemaChamber1 + MigStat$NosemaChamber2)/2

MigStat$Varroa <- MigStat$VarroaLoad / MigStat$TotBees






















































library(plyr)
library(ggplot2)

ggplot(data = Summary, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   group = Treatment)
) + geom_point(size=3) + scale_colour_manual(values = c("black", "blue")) + labs(x = "Sampling Event", y = "Nosema Load") + coord_cartesian(ylim = c(0, 30), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.1, .85), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Queen Origin")


