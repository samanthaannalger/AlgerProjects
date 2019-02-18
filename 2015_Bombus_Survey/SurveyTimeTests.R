###########################################################################################
# Data Analysis for 2015 Bombus Virus Study
# Samantha Alger and P. Alexander Burnham
# July 10, 2017
# Edited by Alex on June 30, 2018
###########################################################################################

#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())

# Call Packages
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("plyr")
library("spdep")
library("lme4")
library("car")
library("ape")
library("MuMIn")

# Set Working Directory 

setwd("~/Documents/GitHub/AlgerProjects/2015_Bombus_Survey/CSV_Files") 

dat <- read.csv("BombSurv.csv", header=TRUE, stringsAsFactors=FALSE)
ModDF <- read.csv("MixedModelDF.csv", header=TRUE, stringsAsFactors=FALSE)

datSplit <- split(dat, dat$target_name)

datSplit$DWV$Date_collected <- as.Date(datSplit$DWV$Date_collected, "%m/%d/%Y")

library(lubridate)

datSplit$DWV$Month <- month(as.POSIXlt(datSplit$DWV$Date_collected, format="%d/%m/%Y"))




smallDat <- ddply(datSplit$DWV, c("Date_collected", "site"), summarise, 
                   n = length(apiary_near_far),
                   mean = mean(apiary_near_far),
                   sd = sqrt(((mean(apiary_near_far))*(1-mean(apiary_near_far)))/n))

smallDat$Month <- month(as.POSIXlt(smallDat$Date_collected, format="%d/%m/%Y"))

chiTab <- table(smallDat$Month, smallDat$mean)

fisher.test(chiTab)












x <- split(ModDF, ModDF$target_name)
DWV <- x$DWV

df <- merge(DWV, smallDat, by="site")

plot(log(df$HBviralLoad), x=df$Date_collected)
summary(lm(log10(df$HBviralLoad+1)~df$Date_collected))

df$Month <- month(as.POSIXlt(df$Date_collected, format="%d/%m/%Y"))

table(df$Month, df$)


DWV$HBviralLoad



kruskal.test(df$HBviralLoad~df$Month)







plot(y=datSplit$DWV$apiary_near_far, x=datSplit$DWV$Date_collected)

mylogit <- glmer(data = datSplit$DWV, apiary_near_far ~ Date_collected + (1|site), family = "binomial")
summary(mylogit)


