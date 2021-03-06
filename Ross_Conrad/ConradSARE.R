###########################################################################################
# Data Analysis for Ross Conrad SARE project
# Samantha Alger and P. Alexander Burnham
# September 12, 2017
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
library("MASS")

# Set Working Directory 
setwd("~/Documents/GitHub/AlgerProjects/Ross_Conrad")

# load in data
#Conrad <- read.table("ConradSARE.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
#Conrad <- read.table("conradSAREtwo.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
Conrad <- read.table("SareFinalData.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

###############################################################
# Varroa Figure
# repeated measures anova for varroa load

# Subset the data to exclude Varroa sampling conducting "pre" quick strips
ConradSub<- Conrad[-which(Conrad$PrePostTreatVarroa=="Pre"), ]

# Summary of Varroa
VarSum <- ddply(ConradSub, c("Treatment", "SamplingEvent"), summarise, 
                n = length(Varroa),
                mean = mean(Varroa, na.rm=TRUE),
                sd = sd(Varroa, na.rm = TRUE),
                se = sd / sqrt(n))

#VarSum <- VarSum[!VarSum$SamplingEvent==7,]

# plotting DWV prev. for experiment 1
colors <- c("slategray3", "dodgerblue4", "black", "blue")

ggplot(data = VarSum, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/300 bees)") + coord_cartesian(ylim = c(0, 40), xlim = c(1,7)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color="Treatment:") + scale_x_continuous(breaks=c(1,2,3,4,5,6,7)) + scale_color_manual(values=colors)

ConradSub1 <- ConradSub[!ConradSub$SamplingEvent==7,]

mod <- glmer(data=ConradSub1, formula = Varroa~Treatment * SamplingEvent + (1+SamplingEvent|ID), family = poisson)

Anova(mod)

###############################################################

# Honey Figure
HoneySum <- ddply(Conrad, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(Honey),
                   mean = mean(Honey, na.rm=TRUE),
                   sd = sd(Honey, na.rm = TRUE),
                   se = sd / sqrt(n))

HoneySum <-HoneySum[complete.cases(HoneySum),]
#HoneySum$SamplingEvent <- as.character(HoneySum$SamplingEvent)
#HoneySum$se[5] <- NA


# plotting DWV prev. for experiment 1
colors <- c("slategray3", "dodgerblue4", "black", "blue")

ggplot(data = HoneySum, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Treatment)
) + geom_point(size=4) + labs(x = "Year", y = "Honey Harvested (# Supers)") + coord_cartesian(ylim = c(0, 4), xlim = c(2,6)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color="Treatment:") + scale_x_continuous(breaks=c(2,4,6), labels =c("2016","2017","2018")) + scale_color_manual(values=colors)




colors <- c("slategray3", "blue", "grey")

plot1 <- ggplot(HoneySum, aes(x=Treatment, y=mean, fill=SamplingEvent)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Treatment", y = "Honey Harvested (# Supers)")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 4)) + scale_fill_manual(values=colors, name="", labels=c("2016", "2017", "2018")) + theme(legend.position=c(.8, .88))

#ANOVA testing honey
HoneyModel <- aov(data=Conrad, Honey~Treatment*SamplingEvent)
summary(HoneyModel)


# Significant difference in honey between treatments. 
# C and QS=A 
# TF and TFQ=B
###############################################################





Conrad$Date <- as.Date(Conrad$DateDied, "%m/%d/%y")
TimeCourse <- Conrad[1:60,]

TimeCourse$startDate <- rep("2016-04-01", 60)

TimeCourse$SurvDays <- difftime(TimeCourse$Date ,TimeCourse$startDate , units = c("days"))

library(ggfortify)
library(survival)

fit <- survfit(Surv(SurvDays, status) ~ Treatment, data = TimeCourse)

autoplot(fit, surv.linetype = 'dashed', conf.int = FALSE,
         censor.shape = '*', censor.size = 10, 
         xlab = "Time (days)", ylab = "Percent Survival")




# Survival Figure
Conrad$SurvivalBINY <- ifelse(Conrad$Survival=="Yes",1,0)
ConradTone <- Conrad[1:60,]


SurvivalSum <- ddply(ConradTone, c("Treatment"), summarise, 
                  n = length(SurvivalBINY),
                  mean = mean(SurvivalBINY, na.rm=TRUE),
                  sd = sd(SurvivalBINY, na.rm = TRUE),
                  se = sd / sqrt(n))

SurvivalSum$se[1] <- NA

colors <- c("slategray3", "dodgerblue4", "black", "blue")

plot1 <- ggplot(SurvivalSum, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Treatment", y = "% Survival")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 1)) + scale_fill_manual(values=colors, name="", labels=c("Stationary", "Migratory")) + theme(legend.position=c(2, 2)) + scale_y_continuous(labels = scales::percent)

# chi square testing Survival

chisq.test(x=ConradTone$Treatment, y=ConradTone$Survival)
# significant difference between treatments
library("lme4")

ConradNoC <- Conrad[!(Conrad$Treatment == "C"),]

Fullmod4 <- glmer(data=ConradNoC, formula = SurvivalBINY~Treatment + (1|ID), family = binomial(link = "logit"))
summary(Fullmod4)


###############################################################
# Strength Figure

StrengthSum <- ddply(Conrad, c("Treatment", "SamplingEvent"), summarise, 
                     n = length(Strength),
                     mean = mean(Strength, na.rm=TRUE),
                     sd = sd(Strength, na.rm = TRUE),
                     se = sd / sqrt(n))


strength <- StrengthSum[!is.na(StrengthSum$mean),]

ggplot(data = strength, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/300 bees)") + coord_cartesian(ylim = c(300, 600), xlim = c(1,6)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color="Treatment:") + scale_x_continuous(breaks=c(1,2,3)) + scale_color_manual(values=colors)










colors <- c("slategray3", "dodgerblue4", "black", "blue")

plot1 <- ggplot(StrengthSum, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Treatment", y = "Brood Area (sq. inches)")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 600)) + scale_fill_manual(values=colors, name="", labels=c("Stationary", "Migratory")) + theme(legend.position=c(2, 2))

#ANOVA on Strength
StrengthModel <- aov(data=Conrad, Strength~Treatment)
summary(StrengthModel)

#No significant difference in brood area between treatments












#############################################################################
#NEW code for survorship figure:

ls()
rm(list=ls())

Conrad <- read.table("ConradSARE_Survivorship.csv",header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings=c("","NA"))

# find death binary 
Conrad$death <- ifelse(Conrad$survivorship==1,0,1)

Conrad2016<- Conrad[(Conrad$year=="2016"),]
Conrad2017<- Conrad[(Conrad$year=="2017"),]



# summary of death for 2016 conrad 
death <- ddply(Conrad2016, c("treatment", "season"), summarise, 
               n = length(death),
               mean = n / 15,
               sd = sd(death, na.rm = TRUE),
               se = sd / sqrt(n))

death <- death[! is.na(death$season=="NA"), ]


library(scales)
library(wesanderson)

z <- (wes_palettes)

ggplot(data = death, 
       aes(x=treatment, y=mean, 
           fill = season)
) + geom_bar(stat = "identity") + labs(x = "Treatment", y = "2016 Percent Losses") + coord_cartesian(ylim = c(0, 1)) + theme_minimal(base_size = 17) + theme(legend.position="top") + scale_y_continuous(labels = percent) + scale_fill_manual(values=c(wes_palettes$Moonrise2))

# summary of death for 2017 conrad 
death <- ddply(Conrad2017, c("treatment", "season"), summarise, 
               n = length(death),
               mean = n / 15,
               sd = sd(death, na.rm = TRUE),
               se = sd / sqrt(n))

death <- death[! is.na(death$season=="NA"), ]


ggplot(data = death, 
       aes(x=treatment, y=mean, 
           fill = season)
) + geom_bar(stat = "identity") + labs(x = "Treatment", y = "2017 Percent Losses") + coord_cartesian(ylim = c(0, 1)) + theme_minimal(base_size = 17) + theme(legend.position="top") + scale_y_continuous(labels = percent) + scale_fill_manual(values=c(wes_palettes$Moonrise2))



