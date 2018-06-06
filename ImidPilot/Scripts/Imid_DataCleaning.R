###########################################################################################
# Data Cleaning for Imid Experiment
# Samantha Alger
# December 11, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/ImidPilot/")

###########################################################################################
# Read in Virus Data:
ImidVirus <- read.table("csv_files/ImidqPCR_RawResults.csv", 
                       header=TRUE, 
                       sep = ",", 
                       stringsAsFactors = FALSE)

KoppertVirus <- read.table("csv_files/ImidPreTestingPcrResults.csv", 
                        header=TRUE, 
                        sep = ",", 
                        stringsAsFactors = FALSE)


# Read in Lab Data
ImidStat <- read.table("csv_files/ImidLabData.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 

KoppertStat <- read.table("csv_files/ImidPreTestingData.csv", 
                       header=TRUE, 
                       sep = ",", 
                       stringsAsFactors = FALSE) 

#Read in Consumption Data

Consump <- read.table("csv_files/FoodConsumpImidPil.csv", header=TRUE, sep = ",") 

###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# source my functions
source("Scripts/BurnhamFunctions.R")
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
#--------------------------------------------------------------------------
# Merge PCR and lab data for Koppert data to add dilution factors and 'treatment' (colony #) to the PCR result dataframe

KoppertVirus <- merge(x = KoppertVirus, y = KoppertStat, by.x = "sample_name")
table(KoppertVirus$sample_name)


#--------------------------------------------------------------------------
# data cleaning

# preliminary data cleaning
ImidVirus <- PrelimClean(data=ImidVirus) 
KoppertVirus <- PrelimClean(data=KoppertVirus)

# take only columns that we want:
data <- dplyr::select(ImidVirus, sample_name, Treatment, dil.factor)
dataPRE <- dplyr::select(KoppertVirus, sample_name, Treatment, dil.factor)

# Separate I1-I10 and use OLDPrelimClean function
# Use PrelimClean funciton on the rest of the data
# append two datasets

OLDImidVirus<-ImidVirus[(ImidVirus$run=="117"),]
ImidVirus<-ImidVirus[!(ImidVirus$run=="117"),]

# use dilution factors to calculate normalized virus load
OLDImidVirus <- OLDVirusNorm(data=OLDImidVirus, number_bees = 1)
ImidVirus <- VirusNorm(data=ImidVirus, number_bees = 1)
KoppertVirus <- OLDVirusNorm(data=KoppertVirus, number_bees = 1)

# Append OldImidVirus and ImidVirus as "ImidVirus"...
ImidVirus <- rbind(ImidVirus, OLDImidVirus)


# use actin to normalize normalized viral load
ImidVirus <- actinNormal(data=ImidVirus)
KoppertVirus <- actinNormal(data=KoppertVirus)

# remove actin from data frame
ImidVirus <- ImidVirus[!(ImidVirus$target_name=="ACTIN"),]
KoppertVirus <- KoppertVirus[!(KoppertVirus$target_name=="ACTIN"),]

# adds virus binary data and makes norm genome copy 0 if above threashold CT
ImidVirus <- CT_Threash(data=ImidVirus)
KoppertVirus <- CT_Threash(data=KoppertVirus)

#Fixes Virus Binary column
ImidVirus$virusBINY <- ifelse(ImidVirus$NormGenomeCopy > 0, 1, 0)
KoppertVirus$virusBINY <- ifelse(KoppertVirus$NormGenomeCopy > 0, 1, 0)

# create full ID to check for inccorrect duplicates 
ImidVirus$fullID <- with(ImidVirus, paste0(sample_name, target_name))
KoppertVirus$fullID <- with(KoppertVirus, paste0(sample_name, target_name))

# remove duplicates for examples with Cts above threashold (i.e. both coerced to 0s) 
ImidVirus <- ImidVirus[!duplicated(ImidVirus$fullID), ] 
KoppertVirus <- KoppertVirus[!duplicated(KoppertVirus$fullID), ] 

# merge the virus data to the main dataset
ImidVirus <- VirusMerger2000(data1 = ImidVirus, data2 = ImidStat)
KoppertVirus <- VirusMerger2000(data1 = KoppertVirus, data2 = KoppertStat)

# log transform virus data:
ImidVirus$logDWV <- log(ImidVirus$DWVload + 1)
ImidVirus$logBQCV <- log(ImidVirus$BQCVload + 1)

KoppertVirus$logDWV <- log(KoppertVirus$DWVload + 1)
KoppertVirus$logBQCV <- log(KoppertVirus$BQCVload + 1)
#-------------------------------------------------------------
# PreTesting


#rename column name from'treatment' to 'colony'
colnames(KoppertStat)[2] <-"colony"
colnames(KoppertVirus)[2] <-"colony"

# Subset to include only colonies 6-8
KoppertVirus <- KoppertVirus[ which(KoppertVirus$colony > 5 & KoppertVirus$colony < 9), ]
KoppertStat <- KoppertStat[ which(KoppertStat$colony > 5 & KoppertStat$colony < 9), ]

# change colony # to a categorical variable
KoppertStat$colony <- as.character(KoppertStat$colony)
KoppertVirus$colony <- as.character(KoppertVirus$colony)

# include virus positive only samples
KoppertDWV <- KoppertVirus[ which(KoppertVirus$DWVbinary=="1"), ]
KoppertBQCV <- KoppertVirus[ which(KoppertVirus$BQCVbinary=="1"), ]

KopDWV <- ddply(KoppertDWV, c("colony"), summarise, 
              n = length(logDWV),
              mean = mean(logDWV, na.rm=TRUE),
              sd = sd(logDWV, na.rm=TRUE),
              se = sd / sqrt(n))

#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black", "green")

#Create a bar graph for viruses loads per Koppert colony:
plot1 <- ggplot(KopDWV, aes(x=colony, y=mean, fill=colony)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="colony", y = "log(Virus load) ")

plot1
x <- aov(data=KoppertDWV, logDWV~colony)
summary(x)
x

x <- aov(data=KoppertBQCV, logBQCV~colony)
summary(x)
#Prevalence is 100% for both viruses in Koppert colonies


#KoppertDWV$colony <- is.factor(KoppertDWV$colony)
is.finite(KoppertDWV$logDWV)

#y <- kruskal.test(logDWV~colony, data = KoppertDWV)

##############################################################
# Get Average DWV load for each Koppert colony
# Add new column to ImidVirus for the preVL, merge by colony ID
# Compute difference between imid VL and preVL for each bee

KopDWVlogmeans <-aggregate(logDWV ~ colony, KoppertDWV, mean)
KopDWVloadmeans <-aggregate(DWVload ~ colony, KoppertDWV, mean)
KopBQCVlogmeans <-aggregate(logBQCV ~ colony, KoppertBQCV, mean)
KopBQCVloadmeans <-aggregate(BQCVload ~ colony, KoppertBQCV, mean)

colnames(KopDWVlogmeans)[2] <-"PreLogDWV"
colnames(KopDWVloadmeans)[2] <-"PreDWVLoad"
colnames(KopBQCVlogmeans)[2] <-"PreLogBQCV"
colnames(KopBQCVloadmeans)[2] <-"PreBQCVLoad"

#Add new column and merge average pre loads and logs to the Imid DF
ImidVirus <- merge(x = ImidVirus, y = KopDWVlogmeans, by = "colony", all.x=TRUE)
ImidVirus <- merge(x = ImidVirus, y = KopDWVloadmeans, by = "colony", all.x=TRUE)
ImidVirus <- merge(x = ImidVirus, y = KopBQCVlogmeans, by = "colony", all.x=TRUE)
ImidVirus <- merge(x = ImidVirus, y = KopBQCVloadmeans, by = "colony", all.x=TRUE)

# select only columns that I want from the consump df
Consump <- dplyr::select(Consump, sample_name, Consumption_g, Consumption_mL, Imid_Consump, TimeStep, Treatment, Colony)

# merge consumption df to virus df
#ImidVirus <- merge(x = ImidVirus, y = Consump, by = "sample_name", all.x = TRUE)

# determine the TOTAL amount of imidacloprid consumed by each bee (aggregate by ID) and merge this to the Imid df

ImidConsumpTotal <-aggregate(Imid_Consump ~ sample_name, Consump, sum)

colnames(ImidConsumpTotal)[2] <-"ImidConsumpTotal"

#Add new column and merge average total imid consumed
Consump <- merge(x = Consump, y = ImidConsumpTotal, by = "sample_name", all.x=TRUE)

#ImidVirus <- merge(x = ImidVirus, y = ImidConsumpTotal, by = "sample_name", all.x=TRUE)

# Write out data frames
#Sucrose/imid consumption data:
write.csv(Consump, "csv_files/ConsumpDF.csv")

#Virus data:
write.csv(ImidVirus, "csv_files/ImidDF.csv")

