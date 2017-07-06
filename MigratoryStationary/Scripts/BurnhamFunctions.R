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








