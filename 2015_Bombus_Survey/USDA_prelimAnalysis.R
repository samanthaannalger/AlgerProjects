###########################################################################################
# Data Analysis for 2015 Bombus Virus Study Neg srand Analysis 
# P. Alexander Burnham
# Jan. 19, 2018
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
library("ape")
library("lme4")
library("car")
library("ape")
library("MuMIn")

# Set Working Directory 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# load in data
BombSurv <- read.table("BombSurvNHBS.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
USDASurv <- read.table("USDA_finalPCR_results.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)



# code to find negative bees and write out columns wanted
#NegList <- BombSurv[BombSurv$virusBINY==0,]
#NegList <- select(NegList, sample_name, target_name, Ct_mean, norm_genome_copbee)
#write.csv(NegList, file = "NegVirus_NegStrd_2015.csv")

# code to find positive bees and write out columns wanted
PosList <- BombSurv[BombSurv$virusBINY==1,]
PosList <- select(PosList, sample_name, target_name, Ct_mean, norm_genome_copbee)

# change the name of sample_name to ID
colnames(PosList)[1] <- "ID" 

#write.csv(PosList, file = "PosVirus_NegStrd_2015.csv")

# just the samples that are postive
x <- PosList$ID[!duplicated(PosList$ID)]

# read nanodrop data:
drop <- read.table("BombSurv_RNANanodropResults.csv",
                   header=TRUE,
                   sep=",",
                   stringsAsFactors=FALSE)

# find number of samples in data set that are postive:
length(PosList$ID[!duplicated(PosList$ID)])

# cross refernce entire data set
complete <- drop[(drop$ID %in% x),]

# select columns we want from DF complete
complete <- select(complete, ID, ng.ul, final_vol)

# merge data frames to get ng.ul in with main frame
mergedDF <- merge(x = PosList, y = complete, by.x = "ID")

# write out .csv file:
#write.csv(mergedDF, file = "NegStrd_PosVirus_2015.csv")

# sample for USDA 1 step checking: 87, 67, 362  
mergedDF[which(mergedDF$target_name=="DWV"),]

















# remove "None" from data frame
USDASurv$Melt.Temperature[USDASurv$Melt.Temperature=="None"] <- 0

# create wide form data set:
USDASurv <- reshape(USDASurv, idvar = "Sample", timevar = "Target", direction = "wide")

# change the names of the data frame
names(USDASurv) <- c("sample_name", "NegCt", "NegMelt", "PosCt", "PosMelt", "SelfP_Ct", "SelfP_Melt")

# change name of drop ID to sample_name
names(drop)[1] <- "sample_name"

# Merge USDA data with main data frame:
BombUSDA <- merge(x = BombSurv, y = USDASurv, by = "sample_name", all.x = TRUE)
BombUSDA <- merge(x = BombUSDA, y = drop, by = "sample_name", all.x = TRUE)


# find data for just DWV for bees we ran at USDA
BombUSDA_DWV <- BombUSDA[BombUSDA$target_name=="DWV",]
BombUSDA_DWV <- BombUSDA_DWV[!is.na(BombUSDA_DWV$PosCt),]
BombUSDA_DWV <- BombUSDA_DWV[-12,]

###############################################################################################
# data analysis on these bees:

# new variable Negative strand binary:
BombUSDA_DWV$NegStrand <- ifelse(BombUSDA_DWV$NegMelt >= 78, 1, 0)

# viral load
boxplot(BombUSDA_DWV$norm_genome_copbee~BombUSDA_DWV$NegStrand)
x <- aov(BombUSDA_DWV$norm_genome_copbee~BombUSDA_DWV$NegStrand)
summary(x)

# concentrarion 
boxplot(BombUSDA_DWV$ng.ul~BombUSDA_DWV$NegStrand)
x <- aov(BombUSDA_DWV$ng.ul~BombUSDA_DWV$NegStrand)
summary(x)






table(mergedDF$target_name)









# calculate total volume est. left in vial
complete$totVol <- 45 - complete$RNA_for_Dilution 

# calculate total RNA left in vial 
complete$totRNA <- complete$totVol * complete$ng.ul

# histgoram of USDA samples ng/ul
hist(complete$ng.ul, xlim = c(10,350))

# five number summary of the sample concentrations we brought to USDA
summary(complete$ng.ul)

# how many are above 100 ng/ul
z <- complete$ng.ul>150
length(which(z==TRUE))

# which are postive for DWV
posDWV <- PosList[PosList$target_name=="DWV",]
