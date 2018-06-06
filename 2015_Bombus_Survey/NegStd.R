# Samantha Alge
#Negative Strand Analysis and figures
# 4/5/2018

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

library("ggplot2")
library("dplyr")
library("lme4")
library("car")
library("plyr")

# load in data
Melt <- read.csv("USDAplate1Melt.csv", header=TRUE, stringsAsFactors=FALSE)
Cq <- read.csv("USDAplate1cq.csv", header=TRUE, stringsAsFactors=FALSE)
BombSurv <- read.csv("BombSurvNHBS.csv", header=TRUE, stringsAsFactors=FALSE)

# remove unwanted sites and bombus species
BombSurv<-BombSurv[!BombSurv$site==("PITH"),]
BombSurv<-BombSurv[!BombSurv$site==("STOW"),]

BombSurv<-BombSurv[!BombSurv$species==("Griseocollis"),]
BombSurv<-BombSurv[!BombSurv$species==("Sandersonii"),]

# subset BombSurv:
Bomb <- select(BombSurv, site, Ct_mean, sample_name, species, apiary_near_far, Density, genome_copbee, norm_genome_copbeeHB, target_name)
names(Bomb)[3] <- "Sample"
Bomb <- Bomb[Bomb$target_name=="BQCV",]

# merge data:
#Dat <- merge(Melt, Cq, by = c("Sample", "Target"))
#str(Dat)

# Merge Dat and Bomb
Dat <- merge(Melt, Bomb, by = c("Sample"),  all.y=TRUE)

#Dat <- merge(Melt, Bomb, by = c("Sample"),  all.x=TRUE)

DatClean <- Dat
#DatClean <- DatClean[!(DatClean$Cq>33),]
#DatClean <- DatClean[!(DatClean$Melt<78),]
DatClean$BinaryNeg <- ifelse(DatClean$Melt > 0, 1, 0)
DatClean$BinaryNeg[is.na(DatClean$BinaryNeg)] <- 0

#ddply summarize:
plotdat <- ddply(DatClean, c("target_name", "apiary_near_far"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))

plotdat$apiary_near_far <- ifelse(plotdat$apiary_near_far==0, "No Apiary", "Apiary")

label.df <- data.frame(Group = c("S1", "S2"),
                       Value = c(6, 9))

plot1 <- ggplot(plotdat, aes(x=apiary_near_far, y=mean, fill=target_name)) +
  geom_bar(stat="identity", color="black", 
           fill =   "white",
           position=position_dodge()) + labs(y="Prevalence", x="Site Type") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)

DatCleanNeg <- DatClean[DatClean$target_name=="BQCV",]

chisq.test(DatCleanNeg$BinaryNeg, DatCleanNeg$apiary_near_far)

xFull <- glmer(data=DatCleanNeg,BinaryNeg~apiary_near_far + species + (1|Sample), family = binomial(link = "logit"))

#for apiary near/far
xNull <- glmer(data=DatCleanNeg,BinaryNeg ~ species + (1|Sample), family = binomial(link = "logit"))

#for species
xNullSP <- glmer(data=DatCleanNeg,BinaryNeg ~ apiary_near_far + (1|Sample), family = binomial(link = "logit"))

anova(xFull, xNull, type="LRT")

anova(xFull,xNullSP, type = "LRT")

summary(xFull)

xNullSP

#ddply summarize for species:
plotdat <- ddply(DatClean, c("target_name", "species"), summarise, 
                 n = length(BinaryNeg),
                 mean = mean(BinaryNeg, na.rm=TRUE),
                 sd = sqrt(((mean(BinaryNeg))*(1-mean(BinaryNeg)))/n))

#plotdat <- plotdat[-3,]


plot1 <- ggplot(plotdat, aes(x=species, y=mean, fill=target_name)) + 
  geom_bar(stat="identity", color="black",fill = "white",
           position=position_dodge()) + labs(y="Prevalence", x="Species") + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 18) + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) + guides(fill=FALSE)
plotdat
