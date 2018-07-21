#####################################################################################
# P. Alexander Burnham
# Final Publication Graphics
# October 25, 2017
# Migratory Stationary Study
# Intended Journal: PLoS ONE
#####################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/MigratoryStationary/")

# Read in Nosema/Varroa/Eco Data:
MigStat <- read.table("Data/MigStatClean.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE) 

# create log base 10
MigStat$log10BQCV <- log10(MigStat$BQCVload+1)
MigStat$log10DWV <- log10(MigStat$DWVload + 1)

MigStat$Treatment <- ifelse(MigStat$Treatment=="Stationary","Stationary\n(Isolated)", ifelse(MigStat$Treatment=="Migratory", "Migratory", "Exposed"))


MigStat$NosemaLoadRecount1 <- ((MigStat$NosemaLoadRecount*4000000)/80)/100000


# required packages:
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(MASS)
library(cowplot)

#####################################################################################
###### EXPERIMENT 1 #################################################################
#####################################################################################

# Setting up dataframe without Exposed for Experiment 1 data (Mig vs Stationary)
MigStatExp_1<-MigStat[!(MigStat$Treatment=="Exposed"),]


# DWV LOAD
#####################################################################################

# Summary of DWV prev. for experiment 1
VirusSum2 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(log10DWV),
                   mean = mean(log10DWV, na.rm=TRUE),
                   sd = sd(log10DWV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV prev. for experiment 1
DWVload1 <- ggplot(data = VirusSum2, 
                   aes(x = SamplingEvent, 
                       y = mean, 
                       group = Treatment)
) + geom_point(size=4) + labs(x = NULL, y = "DWV log(genome copies/bee)") + coord_cartesian(ylim = c(0, 8), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position = "none") + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))
DWVload1







# BQCV LOAD
#####################################################################################

# Summary of BQCV load. for experiment 1
VirusSum1 <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(log10BQCV),
                   mean = mean(log10BQCV, na.rm=TRUE),
                   sd = sd(log10BQCV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting BQCV prev. for experiment 1
BQCVload1 <- ggplot(data = VirusSum1, 
                    aes(x = SamplingEvent, 
                        y = mean, 
                        group = Treatment)
) + geom_point(size=4) + labs(x = NULL, y = "BQCV log(genome copies/bee)") + coord_cartesian(ylim = c(3, 12), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.3, .85),legend.key.width=unit(5,"line")) + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))
BQCVload1




# Varroa LOAD
#####################################################################################

# Summary of Varroa prev. for experiment 1
VarSum <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
                n = length(Varroa),
                mean = mean(Varroa, na.rm=TRUE),
                sd = sd(Varroa, na.rm = TRUE),
                se = sd / sqrt(n))

# plotting Varroa prev. for experiment 1
VARload1 <- ggplot(data = VarSum, 
                   aes(x = SamplingEvent, 
                       y = mean, 
                       group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Varroa (mites/100 bees)") + coord_cartesian(ylim = c(0, 3), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position="none") + labs(linetype="Operation Type:") + scale_x_continuous(breaks=c(1,2,3))
VARload1

nbVAR <- lmer(data=MigStatExp_1, formula = Varroa~Treatment * SamplingEvent + (1|ID))
Anova(nbVAR)


hist(log(MigStatExp_1$Varroa+1))

# FOB
#####################################################################################

# Summary of FOB for experiment 1
FOB <- ddply(MigStatExp_1, c("Treatment", "SamplingEvent"), summarise, 
             n = length(FOB),
             mean = mean(FOB, na.rm=TRUE),
             sd = sd(FOB, na.rm = TRUE),
             se = sd / sqrt(n))

# plotting FOB for experiment 1
FOB1 <- ggplot(data = FOB, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   group = Treatment)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Frames of Bees") + coord_cartesian(ylim = c(5, 25), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(aes(linetype=Treatment), size=1.5) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position = "none") + scale_x_continuous(breaks=c(1,2,3))
FOB1



#####################################################################################
###### EXPERIMENT 2 #################################################################
#####################################################################################



# DWV LOAD
#####################################################################################


# Summary of DWV VL for experiment 2
VirusSum6 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(log10DWV),
                   mean = mean(log10DWV, na.rm=TRUE),
                   sd = sd(log10DWV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting DWV VL for experiment 2
DWVload2 <- ggplot(data = VirusSum6, 
                   aes(x = SamplingEvent, 
                       y = mean, 
                       col = Treatment,
                       linetype= Treatment)
) + geom_point(size=8) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "DWV log(genome copies/bee)") + coord_cartesian(ylim = c(0, 8), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=3) + theme_classic(base_size = 30) + theme(legend.position="none")+ scale_x_continuous(breaks=c(1,2,3)) 
DWVload2


# BQCV LOAD
#####################################################################################

# Summary of BQCV VL for experiment 2
VirusSum8 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                   n = length(log10BQCV),
                   mean = mean(log10BQCV, na.rm=TRUE),
                   sd = sd(log10BQCV, na.rm = TRUE),
                   se = sd / sqrt(n))

# plotting BQCV VL for experiment 2
BQCVload2 <- ggplot(data = VirusSum8, 
                    aes(x = SamplingEvent, 
                        y = mean, 
                        col = Treatment,
                        linetype= Treatment)
) + geom_point(size=8) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "BQCV log(genome copies/bee)") + coord_cartesian(ylim = c(3, 12), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=3) + theme_classic(base_size = 30) + theme(legend.position=c(.3, .85), legend.key.width=unit(10,"line"), legend.key.height = unit(3, "line")) + scale_x_continuous(breaks=c(1,2,3))
BQCVload2

# Varroa LOAD
#####################################################################################




# Summary of varroa load for experiment 2
VarSum2 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
                 n = length(Varroa),
                 mean = mean(Varroa, na.rm=TRUE),
                 sd = sd(Varroa, na.rm = TRUE),
                 se = sd / sqrt(n))

# plotting varroa load for experiment 2
VARload2 <- ggplot(data = VarSum2, 
                   aes(x = SamplingEvent, 
                       y = mean, 
                       col = Treatment,
                       linetype= Treatment)
) + geom_point(size=8) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Varroa (mites/100 bees)") + coord_cartesian(ylim = c(0, 5), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=3) + theme_classic(base_size = 30) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
VARload2


# FOB
#####################################################################################


# Summary of FOB for experiment 2
FOB2 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
              n = length(FOB),
              mean = mean(FOB, na.rm=TRUE),
              sd = sd(FOB, na.rm = TRUE),
              se = sd / sqrt(n))

# plotting FOB for experiment 2
FOB3 <- ggplot(data = FOB2, 
               aes(x = SamplingEvent, 
                   y = mean, 
                   col = Treatment,
                   linetype= Treatment)
) + geom_point(size=8) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Frames of Bees") + coord_cartesian(ylim = c(6, 30), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=3) + theme_classic(base_size = 30) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
FOB3








# EXP2 LDA
#######################################################

# create seperate data frames for each time point:
MigStatSplit <- split(MigStat, MigStat$SamplingEvent)
Mig1 <- MigStatSplit$`1`
Mig2 <- MigStatSplit$`2`
Mig3 <- MigStatSplit$`3`



# run LDA
time2 <- lda(Treatment~Varroa + NosemaLoadRecount + logBQCV + logDWV + DWVbinary + FOB + BroodPattern + NosemaBinary + VarroaBinary, data=Mig2, na.action="na.omit")

# LDA prep function 
ggplotLDAPrep <- function(x){
  if (!is.null(Terms <- x$terms)) {
    data <- model.frame(x)
    X <- model.matrix(delete.response(Terms), data)
    g <- model.response(data)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
      X <- X[, -xint, drop = FALSE]
  }
  means <- colMeans(x$means)
  X <- scale(X, center = means, scale = FALSE) %*% x$scaling
  rtrn <- as.data.frame(cbind(X,labels=as.character(g)))
  rtrn <- data.frame(X,labels=as.character(g))
  return(rtrn)
}

# fit graph for LDA
fitGraph <- ggplotLDAPrep(time2)

# graph for LDA
LDAone <- ggplot(fitGraph, aes(LD1,LD2, color=labels))+geom_point(size=6) + theme_minimal(base_size = 30) + theme(legend.key.width = unit(3, "line"), legend.key.height = unit(3, "line"), legend.position = c(x=.82, y=.85), legend.background = element_rect(fill="white", size=.5), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + stat_ellipse(show.legend = FALSE, level=.7) + labs(x="LD1 (88.36%)", y="LD2 (11.64%)") + scale_color_manual(values=c("black", "red", "blue"), name="Treatment",breaks=c("Exposed", "Migratory", "Stationary"), labels=c("Exposed", "Migratory", "Stationary\n(Isolated)"))
LDAone






# EXP3 LDA
#######################################################

# run LDA
time3 <- lda(Treatment~Varroa + NosemaLoadRecount + logBQCV + logDWV + DWVbinary + FOB + BroodPattern + NosemaBinary + VarroaBinary, data=Mig3, na.action="na.omit")

# graph for LDA
fitGraph1 <- ggplotLDAPrep(time3)

# graph for LDA
LDAtwo <- ggplot(fitGraph1, aes(LD1,LD2, color=labels))+geom_point(size=6) + theme_minimal(base_size = 30) + scale_colour_manual(values = c("black", "red","blue" )) + theme(legend.position= "none", axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) + stat_ellipse(show.legend = FALSE, level=.7) + labs(color="Treatment:", x="LD1 (68.73%)", y="LD2 (31.27%)")
LDAtwo




# merge into one panel
#####################################################################################
#####################################################################################


# experiment 1:
plot_grid(BQCVload1, DWVload1, VARload1, FOB1, labels = c("A", "B", "C", "D"), ncol = 2,align="hv")

# experiment 2:
plot_grid(BQCVload2, DWVload2, VARload2, FOB3, labels = c("A", "B", "C", "D"), ncol = 2,  align="hv", label_size = 30)


# LDA 1 and 2:
plot_grid(LDAone, LDAtwo, labels = c("A", "B"), align="hv",  label_size = 30)













#####################################################################################
# Supplementary figure:
#####################################################################################




# Summary of BP for experiment 2
A <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
             n = length(BroodPattern),
             mean = mean(BroodPattern, na.rm=TRUE),
             sd = sd(BroodPattern, na.rm = TRUE),
             se = sd / sqrt(n))

# plotting BP for experiment 2
#A1 <- ggplot(data = A, 
              # aes(x = SamplingEvent, 
               #    y = mean, 
                #   col = Treatment,
                 #  linetype= Treatment)
#) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "Brood Pattern") + coord_cartesian(ylim = c(3, 5), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
#A1






# Summary of BP for experiment 2
B <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
           n = length(DWVbinary),
           mean = mean(DWVbinary, na.rm=TRUE),
           sd = sd(DWVbinary, na.rm = TRUE),
           se = sd / sqrt(n))

# plotting BP for experiment 2
B1 <- ggplot(data = B, 
             aes(x = SamplingEvent, 
                 y = mean, 
                 col = Treatment,
                 linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "DWV Prevalence") + coord_cartesian(ylim = c(0, 1), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
B1



















# Summary of BP for experiment 2
C <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
           n = length(VarroaBinary),
           mean = mean(VarroaBinary, na.rm=TRUE),
           sd = sd(VarroaBinary, na.rm = TRUE),
           se = sd / sqrt(n))

# plotting BP for experiment 2
C1 <- ggplot(data = C, 
             aes(x = SamplingEvent, 
                 y = mean, 
                 col = Treatment,
                 linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "Varroa Prevalence") + coord_cartesian(ylim = c(0, 1), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
C1
















# Summary of BP for experiment 2
D <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
           n = length(NosemaBinary),
           mean = mean(NosemaBinary, na.rm=TRUE),
           sd = sd(NosemaBinary, na.rm = TRUE),
           se = sd / sqrt(n))

# plotting BP for experiment 2
D1 <- ggplot(data = D, 
             aes(x = SamplingEvent, 
                 y = mean, 
                 col = Treatment,
                 linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Nosema Prevalence") + coord_cartesian(ylim = c(0, 1), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
D1











# Summary of BP for experiment 2
E <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
           n = length(NosemaLoadRecount),
           mean = mean(NosemaLoadRecount, na.rm=TRUE),
           sd = sd(NosemaLoadRecount, na.rm = TRUE),
           se = sd / sqrt(n))

# plotting BP for experiment 2
E1 <- ggplot(data = E, 
             aes(x = SamplingEvent, 
                 y = mean, 
                 col = Treatment,
                 linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = "Sampling Event", y = "Nosema (spores/bee) x100000") + coord_cartesian(ylim = c(0, 25), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position="none") + scale_x_continuous(breaks=c(1,2,3)) 
E1







# Summary of BP for experiment 2
F2 <- ddply(MigStat, c("Treatment", "SamplingEvent"), summarise, 
           n = length(BQCVbinary),
           mean = mean(BQCVbinary, na.rm=TRUE),
           sd = sd(BQCVbinary, na.rm = TRUE),
           se = sd / sqrt(n))

# plotting BP for experiment 2
F15 <- ggplot(data = F2, 
             aes(x = SamplingEvent, 
                 y = mean, 
                 col = Treatment,
                 linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "BQCV Prevalence") + coord_cartesian(ylim = c(0, 1), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(.6, .5), legend.key.width=unit(5,"line")) + scale_x_continuous(breaks=c(1,2,3)) 



# plotting BP for experiment 2
F1 <- ggplot(data = F2, 
              aes(x = SamplingEvent, 
                  y = mean, 
                  col = Treatment,
                  linetype= Treatment)
) + geom_point(size=4) + scale_colour_manual(values = c("darkgrey","black", "black")) + scale_linetype_manual(values = c(1, 1, 2)) + labs(x = NULL, y = "BQCV Prevalence") + coord_cartesian(ylim = c(0, 1), xlim = c(1,2,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05), linetype=1, show.legend=FALSE) + geom_line(size=1.5) + theme_classic(base_size = 17) + theme(legend.position=c(3, 3), legend.key.width=unit(5,"line")) + scale_x_continuous(breaks=c(1,2,3)) 
F1







legend <- get_legend(F15)


A1 <- ggplot(A) + geom_blank()

A1 <- (plot_grid(plot_grid(A1, ncol=1, align='v'),
                 plot_grid(NULL, legend, ncol=1),
                 rel_widths=c(.5, 20)))

# experiment 2:
plot_grid(A1,F1,B1,C1,D1,E1, labels = c(" ","A", "B", "C", "D", "E"), ncol = 2,  align="hv")



plot_grid(DWVload1, BQCVload1, labels = c("A", "B"), ncol = 2,  align="hv")


