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
setwd("~/AlgerProjects/Ross_Conrad/")

# load in data
Conrad <- read.table("ConradSARE.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
#View(Conrad)

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

# plotting DWV prev. for experiment 1
colors <- c("slategray3", "dodgerblue4", "black", "blue")

ggplot(data = VarSum, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/300 bees)") + coord_cartesian(ylim = c(0, 40), xlim = c(1,4)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color="Treatment:") + scale_x_continuous(breaks=c(1,2,3,4)) + scale_color_manual(values=colors)


mod <- lmer(data=ConradSub, formula = Varroa~Treatment * SamplingEvent + (1|ID))

Anova(mod)

###############################################################
# Honey Figure

HoneySum <- ddply(Conrad, c("Treatment"), summarise, 
                   n = length(Honey),
                   mean = mean(Honey, na.rm=TRUE),
                   sd = sd(Honey, na.rm = TRUE),
                   se = sd / sqrt(n))

HoneySum[3,5] <- NA

colors <- c("slategray3", "dodgerblue4", "black", "blue")

plot1 <- ggplot(HoneySum, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Treatment", y = "Honey Harvested (# Supers)")

plot1 + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 3)) + scale_fill_manual(values=colors, name="", labels=c("Stationary", "Migratory")) + theme(legend.position=c(2, 2))

#ANOVA testing honey

HoneyModel <- aov(data=Conrad, Honey~Treatment)
summary(HoneyModel)
TukeyHSD(HoneyModel)

# Significant difference in honey between treatments. 
# C and QS=A 
# TF and TFQ=B
###############################################################
# Survival Figure
Conrad$SurvivalBINY <- ifelse(Conrad$Survival=="Yes",1,0)

ConradTone <- Conrad[1:60,]

SurvivalSum <- ddply(ConradTone, c("Treatment"), summarise, 
                  n = length(SurvivalBINY),
                  mean = mean(SurvivalBINY, na.rm=TRUE),
                  sd = sd(SurvivalBINY, na.rm = TRUE),
                  se = sd / sqrt(n))

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

Fullmod4 <- glmer(data=Conrad, formula = SurvivalBINY~Treatment + (1|ID), family = binomial(link = "logit"))
summary(Fullmod4)


###############################################################
# Strength Figure

StrengthSum <- ddply(Conrad, c("Treatment"), summarise, 
                     n = length(Strength),
                     mean = mean(Strength, na.rm=TRUE),
                     sd = sd(Strength, na.rm = TRUE),
                     se = sd / sqrt(n))

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



