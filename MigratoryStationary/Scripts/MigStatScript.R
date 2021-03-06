###########################################################################################
# Data Analysis for Migratory Stationary Study 
# P. Alexander Burnham
# May 25, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/MigratoryStationary/")

###########################################################################################
# Read in Virus Data:
MigVirus <- read.table("Data/MigVirus.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 


# Read in Nosema/Varroa/Eco Data:
MigStat <- read.table("Data/MigratoryStationaryData.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 


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
library(lmerTest)
library(lsmeans)

#--------------------------------------------------------------------------
# data cleaning

# preliminary data cleaning
MigVirus <- PrelimClean(data=MigVirus) 

# take only columns that we want:
data <- select(MigStat, sample_name, Treatment, dil.factor)

# merge with dilution factors
MigVirus <- merge(MigVirus, data, by=c("sample_name"), all.x =TRUE)

# use dilution factors to calcualte normalized virus load
MigVirus <- VirusNorm(data=MigVirus, number_bees = 50)

# use actin to normalize normalized viral load
MigVirus <- actinNormal(data=MigVirus)

# remove actin from data frame
MigVirus <- MigVirus[!(MigVirus$target_name=="ACTIN"),]

# adds virus binary data and makes norm genome copy 0 if above threashold CT
MigVirus <- CT_Threash(data=MigVirus)

# create full ID to check for inccorrect duplicates 
MigVirus$fullID <- with(MigVirus, paste0(sample_name, target_name))

# remove duplicates for examples with Cts above threashold (i.e. both coerced to 0s) 
MigVirus <- MigVirus[!duplicated(MigVirus$fullID), ] 

# merge the virus data to the main dataset
MigStat <- VirusMerger2000(data1 = MigVirus, data2 = MigStat)

# create average Nosema Load between chambers
MigStat$NosemaLoad <- (MigStat$NosemaChamber1 + MigStat$NosemaChamber2)/2

# create average Nosema Load between chambers for recounted
MigStat$NosemaLoadRecount <- (MigStat$NosemaChamber1Recount + MigStat$NosemaChamber2Recount)/2

# normalize varroa load by number of bees sampled
MigStat$Varroa <- (MigStat$VarroaLoad / MigStat$TotBees) * 100

# log transform virus data:
MigStat$logDWV <- log(MigStat$DWVload + 1)
MigStat$logBQCV <- log(MigStat$BQCVload + 1)





# create binary variable for Nosema:
MigStat$NosemaBinary <- ifelse(MigStat$NosemaLoadRecount == 0, 0, 1)

# create binary variable for Varroa:
MigStat$VarroaBinary <- ifelse(MigStat$VarroaLoad == 0, 0, 1)

# calculate how many pathogens are in each sample (pathgogen richness)
# PathRich <- rowSums(select(MigStat, ChalkBrood, AFB, EFB, PMS, SBV, Snot, BSB, DWV, WaxMoth, # SmallBeetle, DWVbinary, IAPVbinary, BQCVbinary, NosemaBinary, VarroaBinary), na.rm = TRUE)

# write out the clean .csv file to directory 
# write.csv(MigStat, file = "MigStatClean.csv")


#################################################################################################
# DATA ANLAYSIS AND GRAPHICS FOR EXPERIMENT 1:
#################################################################################################

# Setting up dataframe without Exposed for Experiment 1 data (Mig vs Stationary)
MigStatExp_1<-MigStat[!(MigStat$Treatment=="Exposed"),]

# create data frame for running preliminary T1 tests
MigStatExp_1_T1<-MigStatExp_1[(MigStatExp_1$SamplingEvent=="1"),]


#-----------------------------------------------------------------------------------
# DWV PREV:

# Initial T1 test:
x <- table(MigStatExp_1_T1$Treatment, MigStatExp_1_T1$DWVbinary)
chisq.test(x)


# DWV prevalence using glmer
Fullmod <- glmer(data=MigStatExp_1, formula = DWVbinary~Treatment * SamplingEvent + (1+SamplingEvent|ID), family = binomial(link = "logit"))
Anova(Fullmod)

# Summary of DWV prev. for experiment 1
VirusSum <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(DWVbinary),
                   mean = mean(DWVbinary, na.rm=TRUE),
                   sd = sd(DWVbinary, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
ggplot(data = VirusSum, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Prevalance of DWV") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# DWV VL:

# repeated measures anova for DWV
aov.DWV <- aov(logDWV~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.DWV)

# Summary of DWV prev. for experiment 1
VirusSum2 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(logDWV),
                   mean = mean(logDWV, na.rm=TRUE),
                   sd = sd(logDWV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
ggplot(data = VirusSum2, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "DWV log(genome copies/bee)") + coord_cartesian(ylim = c(0, 15), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# BQCV VL:

# repeated measures anova for BQCV
aov.BQCV <- aov(logBQCV~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.BQCV)

# Summary of BQCV prev. for experiment 1
VirusSum1 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                  n = length(logBQCV),
                  mean = mean(logBQCV, na.rm=TRUE),
                  sd = sd(logBQCV, na.rm = TRUE),
                  se = sd / sqrt(n))

# plotting BQCV prev. for experiment 1
ggplot(data = VirusSum1, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "BQCV log(genome copies/bee)") + coord_cartesian(ylim = c(10, 25), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# BQCV PREV:

#BQCV prevalence using glmer
Fullmod2 <- glmer(data=MigStatExp_1, formula = BQCVbinary~Treatment * SamplingEvent + (1+SamplingEvent|ID), family = binomial(link = "logit"))

Anova(Fullmod2)

# Summary of BQCV prev. for experiment 1
VirusSum3 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                  n = length(BQCVbinary),
                  mean = mean(BQCVbinary, na.rm=TRUE),
                  sd = sd(BQCVbinary, na.rm = TRUE),
                  se = sd / sqrt(n))

# plotting BQCV prev. for experiment 1
ggplot(data = VirusSum3, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Prevalance of BQCV") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .2),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# Varroa Load:

# repeated measures anova for Varroa
aov.Var <- aov(Varroa~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.Var)

# Summary of Varroa prev. for experiment 1
VarSum <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(Varroa),
                   mean = mean(Varroa, na.rm=TRUE),
                   sd = sd(Varroa, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting Varroa prev. for experiment 1
ggplot(data = VarSum, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/100 bees))") + coord_cartesian(ylim = c(0, 3), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# Varroa PREV:

#Varroa prevalence using glmer
Fullmod2 <- glmer(data=MigStatExp_1, formula = VarroaBinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

Anova(Fullmod2)

# Summary of Varroa prev. for experiment 1
VarroaBin <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(VarroaBinary),
                   mean = mean(VarroaBinary, na.rm=TRUE),
                   sd = sd(VarroaBinary, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting Varroa prev. for experiment 1
ggplot(data = VarroaBin, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Prevalance of Varroa") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.8, .2),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))



#-----------------------------------------------------------------------------------
# Frames of Bees:

# repeated measures anova for FOB
aov.FOB <- aov(FOB~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.FOB)

# Summary of FOB for experiment 1
FOB <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                n = length(FOB),
                mean = mean(FOB, na.rm=TRUE),
                sd = sd(FOB, na.rm = TRUE),
                se = sd / sqrt(n))

# plotting FOB for experiment 1
ggplot(data = FOB, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Frames of Bees") + coord_cartesian(ylim = c(5, 25), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))


#-----------------------------------------------------------------------------------
# Brood Pattern:

# repeated measures anova for Brood Pattern
aov.BroodPat <- aov(BroodPattern~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.BroodPat)

# Summary of Brood pat for experiment 1
BroodPat <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
             n = length(BroodPattern),
             mean = mean(BroodPattern, na.rm=TRUE),
             sd = sd(BroodPattern, na.rm = TRUE),
             se = sd / sqrt(n))

# plotting brood pat for experiment 1
ggplot(data = BroodPat, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Brood Pattern") + coord_cartesian(ylim = c(2, 8), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))




#################################################################################################
# DATA ANLAYSIS AND GRAPHICS FOR EXPERIMENT 2:
#################################################################################################

# Setting up dataframe without time point 1 for Experiment 2 data (Stationary vs Exposed)
MigStatExp_2_plot<-MigStat[!(MigStat$SamplingEvent=='1'),]

# same data frame but without migratory for analysis purposes 
MigStatExp_2_analysis<-MigStatExp_2_plot[!(MigStatExp_2_plot$Treatment=="Migratory"),]

#-----------------------------------------------------------------------------------
# DWV PREV:

#DWV prevalence using glmer
Fullmod3 <- glmer(data=MigStatExp_2_analysis, formula = DWVbinary~Treatment * SamplingEvent + (1+SamplingEvent|ID), family = binomial(link = "logit"))

Anova(Fullmod3)

# Summary of DWV prev. for experiment 2
VirusSum5 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                  n = length(DWVbinary),
                  mean = mean(DWVbinary, na.rm=TRUE),
                  sd = sd(DWVbinary, na.rm = TRUE),
                  se = sd / sqrt(n))

# plotting DWV prev. for experiment 2
ggplot(data = VirusSum5, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Prevalance of DWV") + coord_cartesian(ylim = c(0, 1), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 


#-----------------------------------------------------------------------------------
# DWV VL:

# repeated measures anova for DWV
aov.DWV2 <- aov(logDWV~Treatment * SamplingEvent + Error(ID), data=MigStatExp_2_analysis)
summary(aov.DWV2)


DWV2full <- lmer(logDWV~Treatment * SamplingEvent + (1|ID), data=MigStatExp_2_analysis)

DWVnull <- lmer(logDWV~SamplingEvent + (1|ID), data=MigStatExp_2_analysis)

anova(DWV2full, DWVnull)
Anova(DWV2full)


# Summary of DWV VL for experiment 2
VirusSum6 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(logDWV),
                   mean = mean(logDWV, na.rm=TRUE),
                   sd = sd(logDWV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV VL for experiment 2
ggplot(data = VirusSum6, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "DWV log(genome copies/bee)") + coord_cartesian(ylim = c(0, 17), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 


#-----------------------------------------------------------------------------------
# Varroa Load:

# repeated measures anova for Varroa
aov.Var2 <- aov(Varroa~Treatment * SamplingEvent + Error(ID), data=MigStatExp_2_analysis)
summary(aov.Var2)

# Summary of varroa load for experiment 2
VarSum2 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(Varroa),
                   mean = mean(Varroa, na.rm=TRUE),
                   sd = sd(Varroa, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting varroa load for experiment 2
ggplot(data = VarSum2, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Varroa (mites/100 bees)") + coord_cartesian(ylim = c(0, 5), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 



#-----------------------------------------------------------------------------------

# Varroa PREV:

#Varroa prevalence using glmer
Fullmod10 <- glmer(data=MigStatExp_2_analysis, formula = VarroaBinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

Anova(Fullmod10)

# Summary of Varroa prev. for experiment 2
VirusSum5 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(DWVbinary),
                   mean = mean(DWVbinary, na.rm=TRUE),
                   sd = sd(DWVbinary, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 2
ggplot(data = VirusSum5, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Prevalance of DWV") + coord_cartesian(ylim = c(0, 1), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 


#-----------------------------------------------------------------------------------

# BQCV PREV:

#BQCV prevalence using glmer
Fullmod4 <- glmer(data=MigStatExp_2_analysis, formula = BQCVbinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

Anova(Fullmod4)

# Summary of BQCV prev. for experiment 2
VirusSum7 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(BQCVbinary),
                   mean = mean(BQCVbinary, na.rm=TRUE),
                   sd = sd(BQCVbinary, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting BQCV prev. for experiment 2
ggplot(data = VirusSum7, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Prevalance of BQCV") + coord_cartesian(ylim = c(0, 1), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .2), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 


#-----------------------------------------------------------------------------------
# BQCV VL:

# repeated measures anova for BQCV
aov.BQCV2 <- aov(logBQCV~Treatment * SamplingEvent + Error(ID), data=MigStatExp_2_analysis)
summary(aov.BQCV2)

# Summary of BQCV VL for experiment 2
VirusSum8 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(logBQCV),
                   mean = mean(logBQCV, na.rm=TRUE),
                   sd = sd(logBQCV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting BQCV VL for experiment 2
ggplot(data = VirusSum8, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "BQCV log(genome copies/bee)") + coord_cartesian(ylim = c(10, 25), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 

#-----------------------------------------------------------------------------------
# Brood Pattern:

# repeated measures anova for BroodPat
aov.BroodPat2 <- aov(BroodPattern~Treatment * SamplingEvent + Error(ID), data=MigStatExp_2_analysis)
summary(aov.BroodPat2)

# Summary of BroodPat for experiment 2
BroodPat2 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                 n = length(BroodPattern),
                 mean = mean(BroodPattern, na.rm=TRUE),
                 sd = sd(BroodPattern, na.rm = TRUE),
                 se = sd / sqrt(n))

# plotting BroodPat for experiment 2
ggplot(data = BroodPat2, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Brood Pattern") + coord_cartesian(ylim = c(3, 5), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 

#-----------------------------------------------------------------------------------
# Frames of Bees:


MigStat


# repeated measures anova for FOB
aov.FOB2 <- aov(FOB~Treatment * SamplingEvent + Error(ID), data=MigStatExp_2_analysis)
summary(aov.FOB2)

# Summary of FOB for experiment 2
FOB2 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
              n = length(FOB),
              mean = mean(FOB, na.rm=TRUE),
              sd = sd(FOB, na.rm = TRUE),
              se = sd / sqrt(n))

# plotting FOB for experiment 2
ggplot(data = FOB2, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Frames of Bees") + coord_cartesian(ylim = c(10, 30), xlim = c(2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(2,3)) 

fit <- lm(DWVload~FOB, data=MigStat)
summary(fit)
plot(fit)


##################

NosMod <- lm(MigStat$NosemaLoadRecount~MigStat$NosemaLoad)
summary(NosMod)

plot(MigStat$NosemaLoad, MigStat$NosemaLoadRecount)

mean(MigStat$NosemaLoad, na.rm=TRUE)

mean(MigStat$NosemaLoadRecount, na.rm=TRUE)






#-----------------------------------------------------------------------------------
# DWV VL:


MigStatT1<-MigStat[(MigStat$SamplingEvent=="1"),]
MigStatT1 <- MigStatT1[!(MigStatT1$Treatment=="Exposed"),]

T1anonva <- aov(NosemaLoadRecount~Treatment, data=MigStatT1)
summary(T1anonva)

# repeated measures anova for DWV
aov.Nos <- aov(NosemaLoadRecount~Treatment * SamplingEvent + Error(ID), data=MigStat)
summary(aov.Nos)

# Summary of DWV prev. for experiment 1
NosSum <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(NosemaLoadRecount),
                   mean = mean(NosemaLoadRecount, na.rm=TRUE),
                   sd = sd(NosemaLoadRecount, na.rm = TRUE),
                   se = sd / sqrt(n))


ggplot(data = NosSum, 
       aes(x = SamplingEvent, 
           y = mean, 
           col = Treatment,
           linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("black", "darkgrey", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Nosema Load") + coord_cartesian(ylim = c(0, 20), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.8, .85), legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(breaks=c(1,2,3)) 


########################################################################################################
#### ANALYSIS OF MITES VS FOB ##########################################################################
########################################################################################################

# clean data set for analysis, take varroa from time point three and merge with data set from 
# time step 2 for only migratory and stationary:
MiteFob <- MigStat[!MigStat$SamplingEvent==1, ]
MiteFob <- MiteFob[!MiteFob$Treatment=="Exposed", ]
MiteFob <- split(MiteFob, MiteFob$SamplingEvent)
MiteFob2 <- MiteFob$`2`
MiteFob3 <- MiteFob$`3`
MiteFob2$Varroa3 <- MiteFob3$Varroa

# main plot parameters
p <- ggplot(MiteFob3, aes(Varroa, FOB)) + geom_point(aes(color=Treatment), size = 3.5) + coord_cartesian(xlim = c(0,10), ylim=c(0,25))  

# aestetic parameters 
p + theme_classic(base_size = 17) + theme(legend.position=c(.85, .8), legend.background = element_rect(color = "black", size = .5)) + labs(color="Treatment", y="Varroa at Time Point 3", x="Frames of Bees at Time Point 2") + scale_color_manual(values = c("slategrey", "black")) + geom_smooth(aes(color=Treatment), method  = lm, se = FALSE)

x <- split(MiteFob3, MiteFob3$Treatment)
Mig <- x$Migratory
Stat<- x$Stationary

modMig <- lm(data=Mig, Varroa~FOB)
summary(modMig)

modStat <- lm(data=Stat, Varroa~FOB)
summary(modStat)



MiteFobFull <- MigStat[!MigStat$Treatment=="Exposed", ]


mod <- lmer(data = MiteFobFull, formula = Varroa ~  SamplingEvent * (FOB * Treatment) + (1|ID) + (1|Yard))

library(lmtest)
acf(residuals(mod))



resid(mod)
Anova(mod)
summary(mod)
coef(summary(mod))

lsmeans(mod, list(pairwise ~ Treatment), adjust = "tukey")

plot(MiteFobFull$FOB, MiteFobFull$Varroa)

