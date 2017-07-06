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


# Read in Data:
MigStat <- read.table("Data/MigratoryStationaryData.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 

# create average Nosema Load between chambers
MigStat$NosemaLoad <- (MigStat$NosemaChamber1 + MigStat$NosemaChamber2)/2

MigStat$Varroa <- MigStat$VarroaLoad / MigStat$TotBees

###########################################################################################





###############################################################################################
################################### FUNCTIONS #################################################
###############################################################################################








###########################################################################
# function name: PrelimClean
# description: removes unneed PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=MigVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- select(data, sample_name, ID, labID, SamplingEvent, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)

  # remove duplicate rows
  data<-data[!duplicated(data), ]

  # remove NTC rows from dataframe:
  data<-data[!(data$sample_name=="No Sample"),]

  # remove Gblock rows from dataframe:
  data<-data[!(data$sample_name=="G-Block"),]

return(data)

}

###########################################################################
# END OF FUNCITON
###########################################################################








###########################################################################
# function name: VirusNorm
# description: normalizes virus data with actin values and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 4, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 100
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  num_bees <- 50
  
  #create column for total_extr_vol
  data$total_extr_vol <- (GITCperbee * num_bees)
  
  # create column for genome copies per bee:
  data$genomeCopies <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################





# pull only actin values out of dataframe
ActinOnly <- BombSurv[which(BombSurv$target_name=="ACTIN"),]

# create DF of ACTIN genome copies and lab ID:
ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$genome_copbee)
colnames(ActinDF) <- c("sample_name", "ACT_genome_copbee")

# merge ACTIN dataframe with main dataframe:
#Need rownames and all.x=TRUE because data frames are different sizes.

BombSurv <- merge(BombSurv, ActinDF, by="sample_name")

# find mean of all ACTIN values:
ActinMean <- mean(ActinOnly$genome_copbee, na.rm = TRUE)

# create column for normalized genome copies per bee:
BombSurv$norm_genome_copbee <- (BombSurv$genome_copbee/BombSurv$ACT_genome_copbee)*ActinMean









###########################################################################
# function name: RepANOVA
# description:conducts a repeated measures ANOVA on a data set and plots 
# parameters: 
# data = data frame, 
# column = name of continuous variable column in quotations
###########################################################################

RepANOVA <- function(data, column){
  
  data$x <- data[,column]
  
  # repeated measures on continuous variable
  aov.out <- aov(x~Treatment * SamplingEvent + Error(ID), data=data)
  
  
  #PLOT:
  library(plyr)
  
  sum <- ddply(data, c("Treatment", "SamplingEvent"), summarise, 
               n = length(x),
               mean = mean(x, na.rm = TRUE),
               sd = sd(x, na.rm = TRUE),
               se = sd / sqrt(n))
  
  # split data frame by origin  
  sumsplit <- split(sum, sum$Treatment)
  
  # create matrix of values and vector of time
  SamplingEvent <- c(1:3)
  Migratory <- sumsplit$Migratory$mean
  Stationary <- sumsplit$Stationary$mean
  Exposed <- sumsplit$Exposed$mean
  mat <- as.matrix(cbind(Migratory, Stationary, Exposed))
  
  # plot matrix
  matplot(x=SamplingEvent, y=mat,
          type="l",
          xlab="Sampling Event", 
          ylab=as.character(column), 
          lwd=2,
          lty=1,
          font.lab=2,
          bty="l", 
          col=2:6)
  
  grid(col="gray")
  
  # create legend for plot
  legend(x=1,y=max(Migratory),
         c("Migratory",
           "Stationary",
           "Exposed"),
         pch=19,
         col=2:6,
         bty="n",
         bg="white")
  
  # return matrix and stats for ANOVA
  return(list(mat, aov.out, summary(aov.out)))
}

###########################################################################
# END OF FUNCITON
###########################################################################










###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

MigVirus <- PrelimClean(data=MigVirus) 


RepANOVA(data=MigStat, column="Varroa")




























aov.out <- aov(FOB ~ Treatment * SamplingEvent + Error(ID), data=MigStat)
summary(aov.out)






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
