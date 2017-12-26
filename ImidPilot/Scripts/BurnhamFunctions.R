###########################################################################################
# Functions I Wrote and Use on a Daily Basis
# P. Alexander Burnham
# July 6, 2017

###########################################################################################
# BEGIN FUNCTIONS #########################################################################
###########################################################################################


###########################################################################
# function name: RepANOVA
# description:conducts a repeated measures ANOVA on a data set and plots 
# parameters: 
# data = data frame, in this case requires columns for SamplingEvent, ID 
# and Treatment plus a column for yvar 
# column = name of continuous variable column in quotations
###########################################################################

RepANOVA <- function(data, column, key=10){
  
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




###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=ImidVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- select(data, Treatment, ID, dil.factor, sample_name, target_name, Ct_mean, Ct_sd, quantity_mean, quantity_sd, run)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$sample_name=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$sample_name=="Gblock"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################



###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 100
  eluteRNA <- 50
  GITCperbee <- 600
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: actinNormal
# description: normalizes virus data with actin values 
# parameters: a data frame with actin values
# returns a dataframe with normalized virus values 
###########################################################################

actinNormal <- function(data=MigVirus){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$target_name=="ACTIN"),]
  
  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$sample_name, ActinOnly$run, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("sample_name", "run", "ACT_genomeCopy")
  
  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("sample_name", "run"), all.x=TRUE)
  
  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)
  
  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean
  
  return(data)
}


###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: CT_Threash
# description: creates binary data and makes genome copy 0 if below Ct threash
# parameters: dataframe
###########################################################################

CT_Threash <- function(data=data){
  
  splitDF <- split(data, data$target_name)

  # make DWV norm_genome_copbee 0 if Ct value is > 32.918
  splitDF$DWV$NormGenomeCopy[which(splitDF$DWV$Ct_mean > 35.918)] <- 0
  splitDF$BQCV$NormGenomeCopy[which(splitDF$BQCV$Ct_mean > 32.525)] <- 0
  splitDF$IAPV$NormGenomeCopy[which(splitDF$IAPV$Ct_mean > 30.796)] <- 0
  
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Ct_mean > 35.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Ct_mean > 32.525, 0, 1)
  splitDF$IAPV$virusBINY  <- ifelse(splitDF$IAPV$Ct_mean > 30.796, 0, 1)
  
  # merge split dataframe back into "BombSurv" dataframe:
  data <- rbind(splitDF$DWV, splitDF$BQCV, splitDF$IAPV)
  
  # norm_genomeCopy is 0 if NA
  data$virusBINY[is.na(data$virusBINY)] <- 0

return(data)

}

###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: VirusMerger2000
# description: merges binary and viral load data for 3 viruses to data frame
# parameters: needs two data frames (data1 is the main one)
###########################################################################

VirusMerger2000 <- function(data1=DF2, data2=DF2){
  
  MigVirusSplit <- split(data1, data1$target_name)
  
  BQCV <- select(MigVirusSplit$BQCV, sample_name, NormGenomeCopy, virusBINY)
  DWV <- select(MigVirusSplit$DWV, sample_name, NormGenomeCopy, virusBINY)
  
  sample_name <- BQCV$sample_name
  BQCVload <- BQCV$NormGenomeCopy
  BQCVbinary <- BQCV$virusBINY
  DWVload <- DWV$NormGenomeCopy
  DWVbinary <- DWV$virusBINY
  
  DF <- data.frame(sample_name, BQCVload, DWVload, DWVbinary, BQCVbinary)
  
  data <- merge(data2, DF, by = c("sample_name"), all.x = T)
  
  return(data)
}

###########################################################################
# END OF FUNCITON
###########################################################################

