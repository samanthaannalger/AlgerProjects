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

#attatch column names to column data
attach(NHBS_DF)

# Head my dataframe to see what's going on
head(NHBS_DF)
 
# Call blue color palette for graphics
library(RColorBrewer)


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
# histogram of nosema count

hist(NosemaCount,
     breaks=35, 
     col="steelblue",
     #ylim=c(0,15),
     #xlim=c(0,1300),
     xlab="Nosema Count per 100 Bees",
     main = "Histogram of Nosema Count"
)

#-------------------------------------------------------------------
# histogram of varroa counts

hist(Varroa100,
     breaks=15, 
     col="steelblue",
     ylim=c(0,5),
     xlim=c(0,10),
     xlab="Varroa Mites per 100 Bees",
     main = "Histogram of Mite Load"
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
                 ylab="Spores per 100 bees (in millions)",
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

Viruses1 <- c(SBPV,	ABPV,	IAPV,	DWV,	LSV2,	CBPV, KBV)

# colors for each virus in plot
colors<-colorRampPalette(brewer.pal(9,"Blues"))(7)
colors<-rev(colors)

# Create barplot showing prevalence for each virus

plot4 <- barplot(height=Viruses1,
                 names.arg=VirusLab,
                 xlab="Viruses",
                 ylab="Viral Prevalence",
                 ylim=c(0,1),
                 #main=paste("Prevalence of 7 RNA Viruses"),
                 col.main="black",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors,
                 border="black"
)



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
MigBinary <- rep(Migratory,4)
GenCopies <- c(IAPV_CPB, DWV_CPB, LSV.2_CPB, CBPV_CPB)
RNA_Vlab <- c(rep("IAPV", 24), rep("DWV", 24), rep("LSV-2", 24), rep("CBPV", 24))

# Create a data frame for it
VirusDF <- data.frame(RNA_Vlab, MigBinary, GenCopies)

# remove rows that include 0 genome copies
VirusDF <- VirusDF[apply(VirusDF[c(3:3)],1,function(z) !any(z==0)),]

# take the log of genome copies and repackage into a dataframe
logGenCopies <- log(VirusDF$GenCopies)
VirusDF <- data.frame(VirusDF,logGenCopies)


#-------------------------------------------------------------------

# aggregate mean viral load by virus and apiary type (migratory/nonmigratory) and calculate se and sd

library(plyr)
VirusDF1 <- ddply(VirusDF, c("MigBinary", "RNA_Vlab"), summarise, 
              n = length(GenCopies),
              mean = mean(GenCopies),
              sd = sd(GenCopies),
              se = sd / sqrt(n))

VirusDF1$ulim <- VirusDF1$mean + VirusDF1$se
VirusDF1$llim <- VirusDF1$mean - VirusDF1$se

VirusDF1$err <- ifelse(VirusDF1$MigBinary == 0, -0.175, 0.175)
VirusDF1


#-------------------------------------------------------------------

# aggregate mean log(viral load) by virus and apiary type (migratory/nonmigratory) and calculate se and sd

library(plyr)
VirusDF2 <- ddply(VirusDF, c("MigBinary", "RNA_Vlab"), summarise, 
                  n = length(logGenCopies),
                  mean = mean(logGenCopies),
                  sd = sd(logGenCopies),
                  se = sd / sqrt(n))

VirusDF2$ulim <- VirusDF2$mean + VirusDF2$se
VirusDF2$llim <- VirusDF2$mean - VirusDF2$se

VirusDF2$err <- ifelse(VirusDF2$MigBinary == 0, -0.175, 0.175)
VirusDF2



#-------------------------------------------------------------------

# call package lattice and use barchart instead of base-R barplot
library(lattice)


barchart(VirusDF2$mean ~ VirusDF2$RNA_Vlab, 
         groups = VirusDF2$MigBinary, #Be sure to add the groups argument to specify how to group the bars
         ylim=c(0,30),
         par.settings=list(superpose.polygon=list(col=c("grey","blue"))),
         main=paste("Viral Load by Apiary Function"),
         auto.key = list(corner=c(0.9,0.9), rectangles = TRUE, points = FALSE, text = c("Non-Migratory","Migratory")), #add legend
         xlab = list(label = "Viruses", fontsize = 20),
         ylab = list(label = "log(genome copies per bee)", fontsize = 20), 
         scales = list(alternating = FALSE, tck = c(1,0), cex=1.25),
         panel = function(x, y, ..., subscripts) 
         
         {panel.barchart(x, y, subscripts = subscripts, ...)
           ll = VirusDF2$llim[subscripts]
           ul = VirusDF2$ulim[subscripts]
           #vertical error bars
           panel.segments(as.numeric(x) + VirusDF2$err[subscripts], ll, #added VirusDF2$err[subscripts]
                          as.numeric(x) + VirusDF2$err[subscripts], ul, #to specify offset of error
                          col = 'black', lwd = 1)                    #bars
           #lower horizontal cap
           panel.text(1.95, 27, "*", cex = 1.8, font = 3)
           panel.segments(as.numeric(x) + VirusDF2$err[subscripts] - 0.1, ll, #same as above 
                          as.numeric(x) + VirusDF2$err[subscripts] + 0.1, ll,
                          col = 'black', lwd = 1)
           #upper horizontal cap
           panel.segments(as.numeric(x) + VirusDF2$err[subscripts] - 0.1, ul, #same as above 
                          as.numeric(x) + VirusDF2$err[subscripts] + 0.1, ul,
                          col = 'black', lwd = 1)  
         })
         


#-------------------------------------------------------------------


# Create vectors for Migratory, Genome copies and virus type
MigBinary <- rep(Migratory,4)
GenCopies <- c(IAPV_CPB, DWV_CPB, LSV.2_CPB, CBPV_CPB)
RNA_Vlab <- c(rep("IAPV", 24), rep("DWV", 24), rep("LSV-2", 24), rep("CBPV", 24))

# Create a data frame for it
VirusDF <- data.frame(RNA_Vlab, MigBinary, GenCopies)

# remove rows that include 0 genome copies
VirusDF <- VirusDF[apply(VirusDF[c(3:3)],1,function(z) !any(z==0)),]

# take the log of genome copies and repackage into a dataframe
logGenCopies <- log(VirusDF$GenCopies)
VirusDF <- data.frame(VirusDF,logGenCopies)
VirusDF

# Run ANOVA


MyModel <- aov(logGenCopies~RNA_Vlab+MigBinary, data=VirusDF) 
z <- summary(MyModel)
z


#-------------------------------------------------------------------


# Create vectors for Migratory, Genome copies and virus type
MigBinary <- Migratory
GenCopies <- DWV_CPB
RNA_Vlab <- c(rep("DWV", 24))
MiteLoad <- Varroa100
NosemaL <- NosMil

# Create a data frame for it
VirusDF <- data.frame(RNA_Vlab, MigBinary, GenCopies, MiteLoad,NosemaL)

# remove rows that include 0 genome copies
VirusDF <- VirusDF[apply(VirusDF[c(3:3)],1,function(z) !any(z==0)),]

# take the log of genome copies and repackage into a dataframe
logGenCopies <- log(VirusDF$GenCopies)
VirusDF <- data.frame(VirusDF,logGenCopies)
VirusDF

# Run ANOVA


MyModel <- aov(logGenCopies~MigBinary, data=VirusDF) 
z <- summary(MyModel)
z


#-------------------------------------------------------------------
# run a linear regression on defmed wing viral load by mite count

VirusDF

plot(x=VirusDF$logGenCopies, 
     y=VirusDF$MiteLoad,
     xlim=c(15,30),
     ylim=c(0,10),
     font.lab=2,
     pch=19,
     ylab="Varroa Mites per 100 bees",
     xlab="DWV log(genome copies per bee)"
     #main=paste("log(viral load) by Varroa Count")
)

grid(col="grey")
# fit these data to a linear model and plot line of best fit (generate equation:
LineBF <- lm(MiteLoad~logGenCopies, data=VirusDF)
line<-abline(LineBF, col = "blue", lwd=3)

# Extract slope and intercept from linear model 

z <- summary(LineBF)
x <- names(z)
z$coefficients[c(1,2)]











