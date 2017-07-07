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

# normalize varroa load by number of bees sampled
MigStat$Varroa <- (MigStat$VarroaLoad / MigStat$TotBees) * 100

# log transform virus data:
MigStat$logDWV <- log(MigStat$DWVload + 1)
MigStat$logBQCV <- log(MigStat$BQCVload + 1)

#################################################################################################
# DATA ANLAYSIS AND GRAPHICS FOR EXPERIMENT 1:
#################################################################################################

# Setting up dataframe without Exposed for Experiment 1 data (Mig vs Stationary)
MigStatExp_1<-MigStat[!(MigStat$Treatment=="Exposed"),]

#-----------------------------------------------------------------------------------
# DWV PREV:

#DWV prevalence using glmer
Fullmod <- glmer(data=MigStatExp_1, formula = DWVbinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

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

# Summary of DWV prev. for experiment 1
VirusSum1 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                  n = length(logBQCV),
                  mean = mean(logBQCV, na.rm=TRUE),
                  sd = sd(logBQCV, na.rm = TRUE),
                  se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
ggplot(data = VirusSum1, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "BQCV log(genome copies/bee)") + coord_cartesian(ylim = c(10, 25), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# BQCV PREV:

#DWV prevalence using glmer
Fullmod2 <- glmer(data=MigStatExp_1, formula = BQCVbinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

Anova(Fullmod2)

# Summary of DWV prev. for experiment 1
VirusSum3 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                  n = length(BQCVbinary),
                  mean = mean(BQCVbinary, na.rm=TRUE),
                  sd = sd(BQCVbinary, na.rm = TRUE),
                  se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
ggplot(data = VirusSum3, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Prevalance of BQCV") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .2),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))

#-----------------------------------------------------------------------------------
# Varroa Load:

# repeated measures anova for BQCV
aov.Var <- aov(Varroa~Treatment * SamplingEvent + Error(ID), data=MigStatExp_1)
summary(aov.Var)

# Summary of DWV prev. for experiment 1
VarSum <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(Varroa),
                   mean = mean(Varroa, na.rm=TRUE),
                   sd = sd(Varroa, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
ggplot(data = VarSum, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/100 bees))") + coord_cartesian(ylim = c(0, 3), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.2, .85),legend.key.width=unit(5,"line"), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))


#-----------------------------------------------------------------------------------
# Frames of Bees:

# repeated measures anova for BQCV
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
Fullmod3 <- glmer(data=MigStatExp_2_analysis, formula = DWVbinary~Treatment * SamplingEvent + (1|ID), family = binomial(link = "logit"))

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

# Summary of DWV prev. for experiment 2
VirusSum6 <- ddply(MigStatExp_2_plot, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(logDWV),
                   mean = mean(logDWV, na.rm=TRUE),
                   sd = sd(logDWV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 2
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



