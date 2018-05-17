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
data <- select(ImidVirus, sample_name, Treatment, dil.factor)
dataPRE <- select(KoppertVirus, sample_name, Treatment, dil.factor)

# use dilution factors to calculate normalized virus load
ImidVirus <- VirusNorm(data=ImidVirus, number_bees = 1)
KoppertVirus <- VirusNorm(data=KoppertVirus, number_bees = 1)
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

# Subset to include only colonies 6-9
KoppertVirus <- KoppertVirus[ which(KoppertVirus$colony > 5), ]
KoppertStat <- KoppertStat[ which(KoppertStat$colony > 5), ]

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

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(KopDWV, aes(x=colony, y=mean, fill=colony)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="colony", y = "log(Virus load) ")

plot1
x <- aov(data=KoppertDWV, logDWV~colony)
summary(x)

#KoppertDWV$colony <- is.factor(KoppertDWV$colony)
is.finite(KoppertDWV$logDWV)

y <- kruskal.test(logDWV~colony, data = KoppertDWV)

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
Consump <- select(Consump, sample_name, Consumption_g, Consumption_mL, Imid_Consump, TimeStep)

# merge consumption df to virus df
ImidVirus <- merge(x = ImidVirus, y = Consump, by = "sample_name", all.x = TRUE)

# determine the TOTAL amount of imidacloprid consumed by each bee (aggregate by ID) and merge this to the Imid df

ImidConsumpTotal <-aggregate(Imid_Consump ~ sample_name, ImidVirus, sum)

colnames(ImidConsumpTotal)[2] <-"ImidConsumpTotal"
#Add new column and merge average total imid consumed
ImidVirus <- merge(x = ImidVirus, y = ImidConsumpTotal, by = "sample_name", all.x=TRUE)

write.csv(ImidVirus, "ImidDF.csv")

ImidDF <- read.csv("csv_files/ImidDF.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)


# DATA ANALYSIS and Figures##########################
ImidDWV <- ImidDF[ which(ImidDF$DWVbinary=="1"), ]
ImidBQCV <- ImidDF[ which(ImidDF$BQCVbinary=="1"), ]

ImidBQCV$BQCVloadDif <- (ImidBQCV$BQCVload-ImidBQCV$PreBQCVLoad)
ImidBQCV$LogBQCVDif <- (ImidBQCV$logBQCV-ImidBQCV$PreLogBQCV)

ImidDWV$DWVloadDif <- (ImidDWV$DWVload-ImidDWV$PreDWVLoad)
ImidDWV$LogDWVDif <- (ImidDWV$logDWV-ImidDWV$PreLogDWV)

##############################################################
# figure for BQCV

BQCVPlot <- ggplot(ImidBQCV, aes(x=Treatment, y=BQCVloadDif)) +
  labs(x="Treatment", y = "BQCV Load Difference")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)

BQCVPlot

# figure for DWV

DWVPlot <- ggplot(ImidDWV, aes(x=Treatment, y=DWVloadDif)) +
  labs(x="Treatment", y = "DWV Load Difference")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)

##############################################################
#figure for DWV load:
#Select only DWV positive samples
ImidDWV <- ImidVirus[ which(ImidVirus$DWVbinary=="1"), ]


#Checking out by treatment
Imid <- ddply(ImidDWV, c("Treatment"), summarise, 
                  n = length(logDWV),
                  mean = mean(logDWV, na.rm=TRUE),
                  sd = sd(logDWV, na.rm=TRUE),
                  se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(Imid, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Treatment", y = "log(Virus load) ")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 15))

x <- aov(data=ImidDWV, logDWV~Treatment)
summary(x)
table(ImidDWV$Treatment, ImidDWV$DWVbinary)

##############################################################
#figure for BQCV load:
#Select only BQCV positive samples
ImidBQCV <- ImidVirus[ which(ImidVirus$BQCVbinary=="1"), ] 

#Checking out by plant species
Imid <- ddply(ImidBQCV, c("Treatment"), summarise, 
              n = length(logBQCV),
              mean = mean(logBQCV, na.rm=TRUE),
              sd = sd(logBQCV, na.rm=TRUE),
              se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(Imid, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Treatment", y = "log(virus load)")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 20))

x <- aov(data=ImidBQCV, logBQCV~Treatment)
summary(x)

###############################################################
#Figure for Virus Prevalence

#Checking out by plant species
DWVPrev <- ddply(ImidVirus, c("Treatment"), summarise, 
                  n = length(DWVbinary),
                  mean = mean(DWVbinary, na.rm=TRUE),
                  sd = sd(DWVbinary, na.rm=TRUE),
                  se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(DWVPrev, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Virus", y = "% of Flowers with Virus Detected")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position= "none") + coord_cartesian(ylim = c(0, .50)) + scale_y_continuous(labels = scales::percent)
