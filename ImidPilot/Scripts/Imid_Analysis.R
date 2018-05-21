# P. Alexander Burnham
# 29, December 2016
# 17, May 2018
# Data Analysis (Figures and Stats) 

#------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/ImidPilot/")

# Read in Data:
ImidDF <- read.csv("csv_files/ImidDF.csv", 
                   header=TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)

ConsumpDF <- read.csv("csv_files/ConsumpDF.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE)

EvapTrials <- read.table("csv_files/EvaporationTrial.csv",
                         header=TRUE, 
                         sep = ",")

# load plyr for data manipulation and ggplot plotting package
library(plyr)
library(ggplot2)
library(grid)
library(dplyr)
library(scales)

# DATA ANALYSIS and Figures##########################
ImidDWV <- ImidDF[ which(ImidDF$DWVbinary=="1"), ]
ImidBQCV <- ImidDF[ which(ImidDF$BQCVbinary=="1"), ]

ImidBQCV$BQCVloadDif <- (ImidBQCV$BQCVload-ImidBQCV$PreBQCVLoad)
ImidBQCV$LogBQCVDif <- log(ImidBQCV$BQCVload/ImidBQCV$PreBQCVLoad)
#ImidBQCV$LogBQCVDif <- (ImidBQCV$logBQCV-ImidBQCV$PreLogBQCV)

ImidDWV$DWVloadDif <- (ImidDWV$DWVload-ImidDWV$PreDWVLoad)
ImidDWV$LogDWVDif <- log(ImidDWV$DWVload/ImidDWV$PreDWVLoad)
#ImidDWV$LogDWVDif <- (ImidDWV$logDWV-ImidDWV$PreLogDWV)

##############################################################
# figure for BQCV

BQCVPlot <- ggplot(ImidBQCV, aes(x=Treatment, y=LogBQCVDif)) +
  labs(x="Treatment", y = "Log BQCV change")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)

BQCVPlot

#checking for normality:
hist(ImidBQCV$BQCVload, data = ImidBQCV)
hist(ImidBQCV$logBQCV, data = ImidBQCV)

# ANOVA Using log transformed data to improve normality:
mod1 <- aov(logBQCV~Treatment, data = ImidBQCV)
summary(mod1)
TukeyHSD(mod1)


# figure for DWV load:

# remove 0.1 and 1 treatment groups
ImidDWV <- ImidDWV[!(ImidDWV$Treatment=="0.1"), ]
ImidDWV <- ImidDWV[!(ImidDWV$Treatment=="1"), ]


ConsumpDF<-ConsumpDF[!(ConsumpDF$sample_name=="I-4" & ConsumpDF$TimeStep==2),]

DWVPlot <- ggplot(ImidDWV, aes(x=Treatment, y=LogDWVDif)) +
  labs(x="Treatment", y = "Log DWV change")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)
DWVPlot

#ANOVA testing for difference, using log transformed value to improve normality
mod1 <- aov(logDWV~Treatment, data = ImidDWV)
summary(mod1)

TukeyHSD(mod1)


##############################################################
###############################################################
#Figure for Virus Prevalence

# DWV prevalece:
DWVPrev <- ddply(ImidDF, c("Treatment"), summarise, 
                 n = length(DWVbinary),
                 mean = mean(DWVbinary, na.rm=TRUE),
                 sd = sd(DWVbinary, na.rm=TRUE),
                 se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black", "green")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(DWVPrev, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Treatment", y = "DWV prevalence")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position= "none") + coord_cartesian(ylim = c(0, .50)) + scale_y_continuous(labels = scales::percent)

##################################################################
# BQCV Prev, 100% of all treatments are infected.

BQCVPrev <- ddply(ImidDF, c("Treatment"), summarise, 
                 n = length(BQCVbinary),
               mean = mean(BQCVbinary, na.rm=TRUE),
                 sd = sd(BQCVbinary, na.rm=TRUE),
                 se = sd / sqrt(n))


#choosing color pallet
colors <- c("goldenrod", "violetred4", "snow1", "black")

#Create a bar graph for viruses by bombus species (aes= aesthetics):
plot1 <- ggplot(BQCVPrev, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Treatment", y = "BQCV prevalence")

plot1 + theme_minimal(base_size = 17) + scale_fill_manual(values=colors, name="Plant Species:", labels=c("Birdsfoot Trefoil", "Red Clover", "White Clover")) + theme(legend.position= "none") + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent)

##############################################################
##############################################################
##############################################################
# CONSUMPTION STATS

# Check for differences between evaporation rates of groups
Evap <- aov(evap~Treatment, data=EvapTrials)
summary(Evap)
# p = 0.46, no difference in amount of sucrose lost due to evaporation.

#############################################
# Time Series for Sucrose Consumption FIGURE:

# I manually made blank sample_names "NA", these are bees that died and were not part of the experiment 

#Remove two outliers- measuring error- values are more than the sucrose and vial can possibly weigh: sample_name: I-33, Treatment:20, time step 3: & sample I-4, Treatment: C, time step 2

ConsumpDF<-ConsumpDF[!(ConsumpDF$sample_name=="I-33" & ConsumpDF$TimeStep==3),]
ConsumpDF<-ConsumpDF[!(ConsumpDF$sample_name=="I-4" & ConsumpDF$TimeStep==2),]
ConsumpDF<-ConsumpDF[!is.na(ConsumpDF$sample_name),]
ConsumpDF<-ConsumpDF[!is.na(ConsumpDF$Consumption_g),]


# summary stats for plotting purposes:
ConsumpSummary <- ddply(ConsumpDF, c("Treatment", "TimeStep"), summarise, 
                     n = length(Consumption_mL),
                     mean = mean(Consumption_mL, na.rm = TRUE),
                     sd = sd(Consumption_mL, na.rm = TRUE),
                     se = sd / sqrt(n))


#Create plot in ggplot 
plot <- ggplot(data = ConsumpSummary, 
               aes(x = TimeStep, 
                   y = mean, 
                   group = Treatment, 
                   colour = Treatment) 
) + geom_line(size=1) + geom_point(size=3) + scale_colour_manual(values = c("dodgerblue4", "black", "darkgreen", "purple", "brown")) + labs(x = "Time (days)", y = "Consumption (mL)") + coord_cartesian(ylim = c(0, 0.5)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.1)) 
plot
# add a theme and add asterix for significance 
#plot + scale_fill_brewer(palette = "Paired") + theme_minimal(base_size = 17) + annotate(geom = "text", x = 1, y = 0.37, label = "**",cex = 6) + annotate(geom = "text", x = 2, y = 0.5, label = "***",cex = 6) + annotate(geom = "text", x = 3, y = 0.41, label = "***",cex = 6) + annotate(geom = "text", x = 4, y = 0.45, label = "***",cex = 6) + annotate(geom = "text", x = 5, y = 0.4, label = "***",cex = 6) 



# Repeated Measures ANOVA to look for differences between treatment groups:

ConsumpDF$TimeStep <- as.factor(ConsumpDF$TimeStep)
ConsumpDF$Treatment <- as.factor(ConsumpDF$Treatment)

mod <- glm(data=ConsumpDF, Consumption_mL~Treatment*TimeStep, family=Gamma)

summary(mod)

library("multcomp")
summary(glht(mod, mcp(Treatment="Tukey")))

# Results:
# ONLY the 20 ppb is significantly different from the control 
# p < 0.001
#------------------------------------------------------------------------
# Figure showing the total amount of Imid consumed by treatment group
# determine the TOTAL amount of imidacloprid consumed by each bee (aggregate by ID) and merge this to the Imid df


## NEED TO FIX IMID TOTALS in the DF, Right now, if there isn't a sample_name (lab ID), then the IMID totals given are incorrect- so all 10 ppb group is wrong, and any other bees that died part way through the experiement. Need to finish lab work for the 10 ppb and remove all rows for bees that died part way through the experiment.
ImidConsumpTotal <-aggregate(Imid_Consump ~ sample_name, ConsumpDF, sum)

colnames(ImidConsumpTotal)[2] <-"ImidConsumpTotal"

ImidConsumpTotal<-ConsumpDF[!(ConsumpDF$sample_name=="")]



#------------------------------------------------------------------------
# statistics looking for differnces between treaments at each of 5 time points using ANOVAs:

# split dataframe by time points 
TimeSplit <- split(ImidDF, ImidDF$TimeStep)

# anova for time point 1
mod1 <- aov(TimeSplit$`1`$Consumption~TimeSplit$`1`$Treatment)
summary(mod1)

TukeyHSD(mod1)


# anova for time point 2
mod2 <- aov(TimeSplit$`2`$Consumption~TimeSplit$`2`$Treatment)
summary(mod2)

TukeyHSD(mod2)

# anova for time point 3
mod3 <- aov(TimeSplit$`3`$Consumption~TimeSplit$`3`$Treatment)
summary(mod3)

TukeyHSD(mod3)


# anova for time point 4
mod4 <- aov(TimeSplit$`4`$Consumption~TimeSplit$`4`$Treatment)
summary(mod4)

TukeyHSD(mod4)


# anova for time point 5
mod5 <- aov(TimeSplit$`5`$Consumption~TimeSplit$`5`$Treatment)
summary(mod5)

TukeyHSD(mod5)


#Food Consumption Analysis:
ImidDF$Treatment <- factor(ImidDF$Treatment, levels = c("Control", "0.1 ppb", "1 ppb", "10 ppb", "20 ppb"))

# summary stats for plotting purposes:
MassSummary <- ddply(ImidDF, c("Treatment", "TimeStep"), summarise, 
                     n = length(Consumption),
                     mean = mean(Consumption, na.rm = TRUE),
                     sd = sd(Consumption, na.rm = TRUE),
                     se = sd / sqrt(n))

# remove time steps with missing data:
print(MassSummary)

colors <- c("dodgerblue4", "slategray3", "blue", "darkgreen", "gray")

plot2 <- ggplot(MassSummary, aes(x=TimeStep, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Time (days)", y = "Consumption (grams)")

plot2 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0.1, 0.5)) + scale_fill_manual(values=colors) + annotate(geom = "text", x = 0.6, y = 0.34, label = "A",cex = 4) 


