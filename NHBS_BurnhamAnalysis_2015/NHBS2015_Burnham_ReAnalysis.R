# Data Analysis for the 2015 Vermont NHBS
# P. Alexander Burnham
# 13 June 2016

# Metadata 
# Author: P. Alexander Burnham
# Date: 13 June 2016
#
# Data Set: These data were collected during the 2015 National Honey Bee Survey in Vermont by Samantha Alger and Alex Burnham with all molecular work being done by the Beltsville lab in Maryland.
# Data Source: 2015 Vermont National Honey Bee Survey
# Funding Source: United States Department of Agriculture (USDA), APHIS, Bee Informed Partnership 
#
# Data Collection: Collection methods are stipulated by the National Honey Bee Survey
# Columns: (from left to right) Beekeeper last name,each virus has a 3 to 4 letter abbreviation followed by (PA=presence/absence) or (CPB=genome copies per bee) VarroaTHR_PA is binary and consists of presence/absence above the threshold.  
# Rows: Data points for all columns in order from each collection event
# Missing values: NA



#-------------------------------------------------------------------
#Preliminaries:

# Clear memory of characters
ls()
rm(list=ls())

# Set Working Directory 
setwd("~/Desktop/RScripts/NHBS_BurnhamAnalysis_2015")

# read in data
NHBS_DF <- read.table("NHBS_2015_DataSubset.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

# Head my dataframe to see what's going on
head(NHBS_DF)
 
# Call blue color palette for graphics
library(RColorBrewer)
attach(NHBS_DF)


#-------------------------------------------------------------------
# calculate each viruses mean viral load and assign to abbreviation 

SBPV <- mean(SBPV_CPB)
ABPV <- mean(ABPV_CPB)
IAPV <- mean(IAPV_CPB)
DWV <- mean(DWV_CPB)
LSV2 <- mean(LSV.2_CPB) 
CBPV <- mean(CBPV_CPB)
KBV <- mean(KBV_CPB)

# concatenate mean viral loads for each virus and create XLabs
Viruses <- c(SBPV,	ABPV,	IAPV,	DWV,	LSV2,	CBPV, KBV)
VirusLab <- c("SBPV",	"ABPV",	"IAPV",	"DWV",	"LSV-2",	"CBPV", "KBV")

# take the log of all viruses that are present 
LogViruses <- c(SBPV,	ABPV,	log(IAPV),	log(DWV),	log(LSV2),	log(CBPV), KBV) 

# colors for each virus in plot
colors<-colorRampPalette(brewer.pal(9,"Blues"))(7)
colors<-rev(colors)

# Create barplot showing log of mean viral load for each virus 
plot1 <- barplot(height=LogViruses,
                 names.arg=VirusLab,
                 xlab="Viruses",
                 ylab="log(Genome Copies per Bee)",
                 ylim=c(0,30),
                 #main=paste("Log of Mean Viral Load"),
                 col.main="black",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors,
                 border="black"
)


#-------------------------------------------------------------------
# barplot showing mite load for each of 24 anonamous apiaries
# colors for each apiary in plot
colors<-colorRampPalette(brewer.pal(9,"Blues"))(24)
colors<-rev(colors)

yards <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X")

#OrdVarroa <- sort(x=Varroa100,decreasing = FALSE)

plot2 <- barplot(height=Varroa100,
                 names.arg=yards,
                 ylab="Varroa Mites per 100 Bees",
                 xlab="Apiaries",
                 ylim=c(0,10),
                 #main=paste("Mite load by Apiary"),
                 col.main="black",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors,
                 border="black"
)

abline(a = 3, b = 0, col="red", lwd=3)


#-------------------------------------------------------------------
# barplot showing nosema load for each of 24 anonamous apiaries
# colors for each apiary in plot
colors<-colorRampPalette(brewer.pal(9,"Blues"))(24)
colors<-rev(colors)

yards <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X")


plot3 <- barplot(height=NosMil,
                 names.arg=yards,
                 ylab="Nosema spores per 100 bees (in millions)",
                 xlab="Apiaries",
                 ylim=c(0,1.5),
                 #main=paste("Nosema load by Apiary"),
                 col.main="black",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors,
                 border="black"
)

abline(a = 1, b = 0, col="red", lwd=3)

#-------------------------------------------------------------------
# Looking at viral prevalence across 7 RNA viruses

SBPV <- mean(SBPV_PA)
ABPV <- mean(ABPV_PA)
IAPV <- mean(IAPV_PA)
DWV <- mean(DWV_PA)
LSV2 <- mean(LSV.2_PA) 
CBPV <- mean(CBPV_PA)
KBV <- mean(KBV_PA)

Virus <- c("SBPV",	"ABPV",	"IAPV",	"DWV",	"LSV2",	"CBPV", "KBV")
VirusPrev <- c(SBPV,	ABPV,	IAPV,	DWV,	LSV2,	CBPV, KBV)
virus_data <- data.frame(Virus, VirusPrev)
virus_data

# colors for each virus in plot
colors<-"cadetblue3"


# plot in ggplot using minamal theme with some graphical additions
plot1 <- ggplot(virus_data, aes(x = virus_data$Virus, 
                                   y = virus_data$VirusPrev, fill=colors)) + geom_bar(stat="identity", position=position_dodge()) + labs(x="Virus", y = "% of Infected Apiaries")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values=colors) + theme(legend.position=c(-1, -1)) + scale_y_continuous(labels=percent)



#-------------------------------------------------------------------
# Looking at viral prevalence across 7 RNA viruses
Nosema <- mean(NHBS_DF$NosemaPA)
NosemaThreshold <- mean(NHBS_DF$NosemaTHR)
Varroa <- mean(NHBS_DF$VarroaPA)
VarroaThreshold <- mean(NHBS_DF$VarroaTHR_PA)


Parasites <- c(Nosema, NosemaThreshold, Varroa, VarroaThreshold)

ParasiteNames <- c("Nosema", "Nosema Threshold", "Varroa", "Varroa Threshold")

# colors for each virus in plot
colors<-colorRampPalette(brewer.pal(9,"Blues"))(2)
colors<-rev(colors)

# Create barplot showing prevalence for each virus

plot5 <- barplot(height=Parasites,
                 names.arg=ParasiteNames,
                 #xlab="Virus",
                 ylab="Parasite Prevalence",
                 ylim=c(0,1),
                 #main=paste("Prevalence of Varroa and Nosema with Thresholds"),
                 col.main="black",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = c("blue", "grey"),
                 border="black"
)

legend(x=0.3,y=1,
       c("Total Prevalence",
         "Exceeded Threshold"),
       fill=c("blue", "grey"),
       bty="n",
       bg="white")


#-------------------------------------------------------------------

# Create vectors for Migratory, Genome copies and virus type
MigBinary <- as.character(rep(Migratory,4))
GenCopies <- c(IAPV_CPB, DWV_CPB, LSV.2_CPB, CBPV_CPB)
RNA_Vlab <- c(rep("IAPV", 24), rep("DWV", 24), rep("LSV-2", 24), rep("CBPV", 24))

# Create a data frame for it
VirusDF <- data.frame(RNA_Vlab, MigBinary, GenCopies)

# create presence/abscence vector (binary) for viruses:
VirusPA <- ifelse(VirusDF$GenCopies > 0, 1, 0)

VirusDF <- data.frame(VirusDF, VirusPA)

library(plyr)
library(ggplot2)
library(scales)

VirusSum <- ddply(VirusDF, c("MigBinary", "RNA_Vlab"), summarise, 
                  n = length(VirusPA),
                  mean = mean(VirusPA),
                  sd = sd(VirusPA),
                  se = sd / sqrt(n))


colors <- c("slategray3", "dodgerblue4")

plot1 <- ggplot(VirusSum, aes(x=RNA_Vlab, y=mean, fill=MigBinary)) + geom_bar(stat="identity", position=position_dodge()) + labs(x="Virus", y = "% of Infected Apiaries")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values=colors, name="Operation:", labels=c("Stationary", "Migratory")) + theme(legend.position=c(.8, .85)) + scale_y_continuous(labels=percent)

# stats - chi-sq 

# split DF:
splitVDF <- split(VirusDF, RNA_Vlab)

#test for DWV
chisq.test(splitVDF$DWV$VirusPA, splitVDF$DWV$MigBinary)

#test for CBPV
chisq.test(splitVDF$CBPV$VirusPA, splitVDF$CBPV$MigBinary)

#testfor LSV-2
chisq.test(splitVDF$`LSV-2`$VirusPA, splitVDF$`LSV-2`$MigBinary)
#-------------------------------------------------------------------
# aggregate mean viral load by virus and apiary type (migratory/nonmigratory) and calculate se and sd

# GeekGasm = log(GenCopies) 
VirusDF$GeekGasm <- log(VirusDF$GenCopies)
VirusDF$GeekGasm[VirusDF$GeekGasm=="-Inf"] <- NA
VirusDF

library(plyr)
VirusSum1 <- ddply(VirusDF, c("MigBinary", "RNA_Vlab"), summarise, 
              n = length(GeekGasm),
              mean = mean(GeekGasm, na.rm=TRUE),
              sd = sd(GeekGasm, na.rm = TRUE),
              se = sd / sqrt(n))

#-------------------------------------------------------------------

plot1 <- ggplot(VirusSum1, aes(x=RNA_Vlab, y=mean, fill=MigBinary)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Virus", y = "log(Viral Load (genome copies/bee))")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 30)) + scale_fill_manual(values=colors, name="Operation:", labels=c("Stationary", "Migratory")) + theme(legend.position=c(.8, .85)) + annotate(geom = "text", x = 2, y = 28, label = "*",cex = 8) 

#-------------------------------------------------------------------
# Stats -> t.test of DWV p<0.05

VirusSplit <- split(VirusDF, RNA_Vlab)

t.test(VirusSplit$DWV$GeekGasm~VirusSplit$DWV$MigBinary, alternative="less")

#-------------------------------------------------------------------
# varroa by yard type:

VarroaSum <- ddply(NHBS_DF, c("MigBinary"), summarise, 
                   n = length(Varroa100),
                   mean = mean(Varroa100, na.rm=TRUE),
                   sd = sd(Varroa100, na.rm = TRUE),
                   se = sd / sqrt(n))

plot2 <- ggplot(VarroaSum, aes(y=mean, x="", fill=MigBinary)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Varroa", y = "Varroa Load (per 100 bees)")

plot2 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 8)) + scale_fill_manual(values=colors, name="Operation:", labels=c("Stationary", "Migratory")) + annotate(geom = "text", x = 1, y = 7, label = "**",cex = 10) 

#-------------------------------------------------------------------
# varroa by yard type:

x <- aov(NHBS_DF$Varroa100~NHBS_DF$Migratory)
summary(x)


#-------------------------------------------------------------------
# varroa by yard type:

NosemaSum <- ddply(NHBS_DF, c("MigBinary"), summarise, 
                   n = length(NosemaCount),
                   mean = mean(NosemaCount, na.rm=TRUE),
                   sd = sd(NosemaCount, na.rm = TRUE),
                   se = sd / sqrt(n))

plot2 <- ggplot(NosemaSum, aes(y=mean, x="", fill=MigBinary)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Nosema", y = "Nosema Load (per 50 bees)")

plot2 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 6)) + scale_fill_manual(values=colors, name="Operation:", labels=c("Stationary", "Migratory")) + annotate(geom = "text", x = 1, y = 7, label = "**",cex = 10) 

#-------------------------------------------------------------------
#  N0sema by yard type:

x <- aov(NHBS_DF$NosemaCount~NHBS_DF$Migratory)
summary(x)


