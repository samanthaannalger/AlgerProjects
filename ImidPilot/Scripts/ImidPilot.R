# P. Alexander Burnham
# 29, December 2016
# Imid. Pilot Food Consumption
# Data Analysis (Figures and Stats) 

#------------------------------------------------------------------------
# Preliminaries:

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/Desktop/ImidPilot")

# Read in Data:
ImidDF <- read.table("FoodConsumpImidPil.csv", header=TRUE, sep = ",") 
head(ImidDF)


# load plyr for data manipulation and ggplot plotting package
library(plyr)
library(ggplot2)
library(grid)
library(dplyr)
library(scales)

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






















#Create plot in ggplot 
plot <- ggplot(data = MassSummary, 
               aes(x = TimeStep, 
                   y = mean, 
                   group = Treatment, 
                   colour = Treatment) 
) + geom_line(size=1) + geom_point(size=3) + scale_colour_manual(values = c("dodgerblue4", "black", "darkgreen", "purple", "brown")) + labs(x = "Time (days)", y = "Consumption (grams)") + coord_cartesian(ylim = c(0.1, 0.5)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.1)) 

# add a theme and add asterix for significance 
plot + scale_fill_brewer(palette = "Paired") + theme_minimal(base_size = 17) + annotate(geom = "text", x = 1, y = 0.37, label = "**",cex = 6) + annotate(geom = "text", x = 2, y = 0.5, label = "***",cex = 6) + annotate(geom = "text", x = 3, y = 0.41, label = "***",cex = 6) + annotate(geom = "text", x = 4, y = 0.45, label = "***",cex = 6) + annotate(geom = "text", x = 5, y = 0.4, label = "***",cex = 6) 





plot(ImidDF$TreatmentNum, 
     ImidDF$Consumption,
     ylim=c(0,0.8)
     )



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

