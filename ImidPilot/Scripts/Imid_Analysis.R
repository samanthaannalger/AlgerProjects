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
library(lme4)
library(car)

# DATA ANALYSIS and Figures##########################

# Only include samples positive for viruses:
ImidDWV <- ImidDF[ which(ImidDF$DWVbinary=="1"), ]
ImidBQCV <- ImidDF[ which(ImidDF$BQCVbinary=="1"), ]

#Calculate the differences between the starting virus loads (taken directly from Koppert colonies) and the final virus loads, create new columns:
ImidBQCV$BQCVloadDif <- (ImidBQCV$BQCVload-ImidBQCV$PreBQCVLoad)
ImidBQCV$LogBQCVDif <- log(ImidBQCV$BQCVload/ImidBQCV$PreBQCVLoad)

ImidDWV$DWVloadDif <- (ImidDWV$DWVload-ImidDWV$PreDWVLoad)
ImidDWV$LogDWVDif <- log(ImidDWV$DWVload/ImidDWV$PreDWVLoad)

##############################################################
# figure for BQCV

#Change order of factors
ImidBQCV$Treatment <- factor(ImidBQCV$Treatment, levels = c("C", "0.1", "1", "10", "20"))

# Figure of the log virus genome copies 
BQCVPlot <- ggplot(ImidBQCV, aes(x=Treatment, y=logBQCV, fill=Treatment)) +
  labs(x="Dose (ppb)", y = "Log(BQCV titer)")+
  theme_classic() +  
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + scale_x_discrete(labels=c("Control", "0.1","1","10","20")) + scale_fill_grey(start = 1, end = .4, guide=FALSE)

BQCVPlot

# STATISTICS FOR BQCV:##############################################################

#checking for normality:
hist(ImidBQCV$BQCVload, data = ImidBQCV)
hist(ImidBQCV$logBQCV, data = ImidBQCV)

# ANOVA Using log transformed data to improve normality:
mod1 <- aov(logBQCV~Treatment, data = ImidBQCV)
summary(mod1)
TukeyHSD(mod1)

# Mixed Model for BQCV load:
BQCVMod <- lmer(logBQCV~Treatment + (1|colony), data=ImidBQCV)
summary(BQCVMod)
Anova(BQCVMod)

library("multcomp")
summary(glht(BQCVMod, mcp(Treatment="Tukey")))

############################################################
# figure for DWV load:

#Change order of factors
ImidDWV$Treatment <- factor(ImidDWV$Treatment, levels = c("C", "0.1", "1", "10", "20"))

# remove 0.1 and 1 treatment groups
ImidDWV <- ImidDWV[!(ImidDWV$Treatment=="0.1"), ]
ImidDWV <- ImidDWV[!(ImidDWV$Treatment=="1"), ]

# Figure of the log virus genome copies 
DWVPlot <- ggplot(ImidDWV, aes(x=Treatment, y=logDWV, fill=Treatment)) +
  labs(x="Dose (ppb)", y = "Log(DWV titer)")+
  theme_classic() +  
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + scale_x_discrete(labels=c("Control","10","20")) + scale_fill_grey(start = 1, end = .4, guide=FALSE)

DWVPlot

# STATISTICS FOR DWV:##############################################################

#ANOVA testing for difference, using log transformed value to improve normality
mod1 <- aov(logDWV~Treatment, data = ImidDWV)
summary(mod1)
TukeyHSD(mod1)

# Mixed model for DWV load:
DWVMod <- lmer(logDWV~Treatment + (1|colony), data=ImidDWV)
summary(DWVMod)
Anova(DWVMod)

library("multcomp")
summary(glht(DWVMod, mcp(Treatment="Tukey")))


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

# Remove two outliers- measuring error- values are more than the sucrose and vial can possibly weigh: sample_name: I-33, Treatment:20, time step 3: & sample I-4, Treatment: C, time step 2

ConsumpDF<-ConsumpDF[!(ConsumpDF$sample_name=="I-33" & ConsumpDF$TimeStep==3),]
ConsumpDF<-ConsumpDF[!(ConsumpDF$sample_name=="I-4" & ConsumpDF$TimeStep==2),]
ConsumpDF<-ConsumpDF[!is.na(ConsumpDF$sample_name),]
ConsumpDF<-ConsumpDF[!is.na(ConsumpDF$Consumption_g),]

# reorder factors for plotting
ConsumpDF$Treatment <- factor(ConsumpDF$Treatment, levels = c("C", "0.1", "1", "10", "20"))

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
) + geom_line(size=1) + geom_point(size=3) + scale_colour_manual(values = c("dodgerblue4", "black", "darkgreen", "purple", "brown"), name = "Dose(ppb)", labels=c("Control","0.1", "1", "10","20")) + labs(x = "Time (days)", y = "Sucrose consumption/bee/day (mL)") + coord_cartesian(ylim = c(0, 0.5)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.1)) + theme_bw()

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
# the 20 ppb is significantly different from the control, 1, 0.1 (p < 0.001) ; and 10 ppb is different from 0.1 (p=0.006)


################################################################ TOTAL Sucrose CONSUMED:

# Figure showing the total amount of sucrose consumed by treatment group, Did bees consume different amounts total?

#Only take one timestep bee to plot
ConsumpDF<-ConsumpDF[(ConsumpDF$TimeStep=="1"),]

# Plot figure showing the total amount of sucrose consumed by group (ng)

ConsumpDF$logImidConsumpTotal <- log(ConsumpDF$ImidConsumpTotal)

TotalImid <- ggplot(ConsumpDF, aes(x=Treatment, y=ImidConsumpTotal, fill=Treatment)) +
  labs(x="Dose (ppb)", y = "Total Imidacloprid Consumed (ng)")+
  theme_classic() +  
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + scale_x_discrete(labels=c("0.1","1", "10","20")) + scale_fill_grey(start = 1, end = .4, guide=FALSE)

TotalImid




################################################################ TOTAL IMID CONSUMED:

# Figure showing the total amount of Imid consumed by treatment group, Did bees consume different amounts of Imid?

#Only take one timestep bee to plot
ConsumpDF<-ConsumpDF[(ConsumpDF$TimeStep=="1"),]
#Remove Control group from dataset
ConsumpDF<-ConsumpDF[!(ConsumpDF$Treatment=="C"),]
ConsumpDF$Treatment <- as.character(ConsumpDF$Treatment)
ConsumpDF$Treatment <- as.factor(ConsumpDF$Treatment)


# Check distribution, data are gamma distribution
hist(ConsumpDF$ImidConsumpTotal)

# Plot figure showing the total amount of Imid consumed by group (ng)

ConsumpDF$logImidConsumpTotal <- log(ConsumpDF$ImidConsumpTotal)

TotalImid <- ggplot(ConsumpDF, aes(x=Treatment, y=ImidConsumpTotal, fill=Treatment)) +
  labs(x="Dose (ppb)", y = "Total Imidacloprid Consumed (ng)")+
  theme_classic() +  
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + scale_x_discrete(labels=c("0.1","1", "10","20")) + scale_fill_grey(start = 1, end = .4, guide=FALSE)

TotalImid

#Test for differences with ANOVA:
mod3 <- aov(ImidConsumpTotal~Treatment, data = ConsumpDF)
summary(mod3)
TukeyHSD(mod3)

# Mixed model for Imid Consumption, colony as a random effect:
ConsumpFull <- glmer(ImidConsumpTotal~Treatment + (1|Colony), data=ConsumpDF, family = Gamma)

Mod <- glm(ImidConsumpTotal~Treatment, data=ConsumpDF, family = Gamma)

summary(Mod)
########################################################
#SURVIVAL ANALYSIS
########################################################
########################################################

# Read in Data:
BeeMort <- read.csv("csv_files/SurvivDF.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

#install.packages("survival")
# Loading the package
library("survival")

#Fitting the survival model
mod <- survdiff(Surv(time, status) ~ Treatment, data=BeeMort, rho = 0)
mod
# NO difference in mortality (df = 4, chisq = 4.3, p = 0.4), "Kaplan-Meier estimate of survival, Tests if there is a difference between two or more survival curves using the Gρ family of tests, or for a single curve against a known alternative." Using 'survival' package, function = "survdiff"

#Creating the survival curve:
survival_func=survfit(Surv(BeeMort$time,BeeMort$status)~BeeMort$Treatment)

plot(survival_func, ylab="Survival", xlab="Days", col=c("red", "blue", "green", "black", "purple"))
legend(.6, .6, legend=c("0.1", "1", "10","20","Control"), title = "Dose (ppb)",
       col=c("red", "blue", "green", "black", "purple"), lty=1, cex=0.8)








# END OF FIGURES AND STATS##################################
################################################################################################################################################################################################################################################




#Old Figures:

#BQCV:
BQCVPlot <- ggplot(ImidBQCV, aes(x=Treatment, y=LogBQCVDif)) +
  labs(x="Treatment", y = "Log BQCV change")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)

BQCVPlot

#DWV:
DWVPlot <- ggplot(ImidDWV, aes(x=Treatment, y=LogDWVDif)) +
  labs(x="Treatment", y = "Log DWV change")+
  theme_classic() +  
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=1) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red") + geom_hline(yintercept=0, linetype="solid", color = "blue", size=1.5)
DWVPlot

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

Summary()

x<- glmer(data = ImidDF, DWVbinary~Treatment + (1|colony), family = binomial(link = "logit"))

Anova(x)

#           Chisq Df Pr(>Chisq)
#Treatment 3.836  4     0.4287

# Number of bees infected in each group:
# C 0.1   1  10  20 
# 5   2   1   4   7 


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
  



#------------------------------------------------------------------------
# OLD CONSUMPTION STATS:
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


