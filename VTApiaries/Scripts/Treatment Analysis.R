#Data Analysis for Treatments
# S Alger
#Aug. 9, 2018


#Preliminaries
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/VTApiaries/CSV_files/")

Shinydf <- read.csv("Shinydf.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)




# What 'meds' are used at the apiary level? (To be used in big analysis to predict colony loss)
MedUse <- Shinydf[! is.na(Shinydf$TreatmentsForColonyHealthMitigation), ]

# USE FOR BEEKEEPER LEVEL DATA:
#MedUse <- MedUse[!duplicated(MedUse$BeekeeperID), ]



df<- c("None","Fluvalinate","Coumaphos","Amitraz","Apiguard","Api-life VAR","Sucrocide","Powdered sugar","Oxalic acid","Formic acid","Menthol","Hopguard","Honey B Healthy","Fumagillin-B","Tylosin","Lincomycin","Oxytetracycline","Herbal antibiotics")

# preallocate matrix space for for loop operation
mat <- matrix(nrow = length(MedUse$ID), ncol = length(df))
colnames(mat) <- df

for (i in 1:length(df)) {
  mat[,i] <-  as.integer(grepl(df[i], MedUse$TreatmentsForColonyHealthMitigation))
  
}
# convert matrix to a df
mat.DF <- as.data.frame(mat)

#Creating columns:
# sum organic mitacides used:
mat.DF$oSum <- rowSums(mat.DF[,c("Formic acid", "Oxalic acid", "Apiguard", "Api-life VAR", "Menthol", "Hopguard", "Sucrocide")])

# sum synthetic mitacides used:
mat.DF$sSum <- rowSums(mat.DF[,c("Fluvalinate","Coumaphos","Amitraz")])

# synthetic mitacide used? T/F?
mat.DF$sTF <- ifelse(mat.DF$sSum >= 1, 1,0)

# organic mitacide used? T/F?
mat.DF$oTF <- ifelse(mat.DF$oSum >= 1, 1,0)

# Sum of synthetic T/F and organic TF to determine next column: 
mat.DF$soSum <- rowSums(mat.DF[,c("sTF","oTF")])

# Both organic and synthetic mitacides used (T/F)
mat.DF$soTF <- ifelse(mat.DF$soSum >= 2, 1,0)

# Sum of all mitacides used:, if 0, beek used no mitacides
mat.DF$sum <- rowSums(mat.DF[,c("oSum","sSum")])

# Mitacide used?
mat.DF$mitTF <- ifelse(mat.DF$sum >= 1, 1,0)

# How many beeks did not use any treatments:
mean(mat.DF$None)
# 0.2309198

# How many beeks reported to using a mitacide treatment?
mean(mat.DF$mitTF)
# 0.6790607

# of those, how many used synthetic vs. organic treatments?
# include only beeks who use mitacide:
treat <- mat.DF[which(mat.DF$mitTF == 1), ]

mean(treat$sTF)
# 0.03458213
mean(treat$oTF)
# 0.9855908

mean(treat$`Formic acid`)
# 0.6772334
Formic <- treat[which(treat$`Formic acid` == 1), ]
table(Formic$sum)
# 79.1% use formic acid only
mean(treat$Apiguard)
mean(treat$`Api-life VAR`)
mean(treat$Sucrocide)
mean(treat$`Oxalic acid`)
mean(treat$Menthol)
mean(treat$Hopguard)

#
table(mat.DF$sum)
164+293+51+3
293/511
51/511
3/511



# Bind to original data set to get ready for Big model:
MedUse<-cbind(MedUse, mat.DF)
SubMedUse<- dplyr::select(MedUse, LocationID, oSum, sSum, sTF, oTF, soTF, soSum, mitTF, sum, None)

# Merge climate data with treatment data (had to read in climate data from other script, Need to compile the entire script here instead)
Modeldf<- merge(SubMedUse, ClimDiv, by = 'LocationID')

# Write out CSV to be used in the model:
#write.csv(Modeldf, file = "Modeldf.csv")

mod<-lm(Modeldf$PerTotLoss~Modeldf$sum)
summary(mod)



# Full, Null and Reduced Models
LossModFull <- lmer(data=Modeldf, formula = PerTotLoss ~ MiteCounts + mitTF + SupplementalFeed + Div + Beektype + (1|BeekeeperID))

LossInteraction <- lmer(data=Modeldf, formula = PerTotLoss ~mitTF * Div + (1|BeekeeperID))

LossModNull <- lmer(data=Modeldf, formula = PerTotLoss ~ 1  + (1|BeekeeperID))

# Analysis of deviance, Type II Wald Chi Square test:
library(car)
Anova(LossModFull, test="Chisq")
Anova(LossInteraction, test="Chisq")


library("lmtest")
acf(residuals(LossModFull))

resid(LossModFull)
Anova(LossModFull)
summary(LossModFull)
coef(summary(LossModFull))

#Multiple comparisons tests each factor:
library("multcomp")
summary(glht(LossModFull, mcp(Div="Tukey")))

boxplot(Modeldf$PerTotLoss~Modeldf$Div)


#anova(LossModFull, LossModBeektype, test="LRT")
Modeldf$mitTF<-as.factor(Modeldf$mitTF)
n <- aov(Modeldf$PerTotLoss~Modeldf$mitTF)
summary(n)
TukeyHSD(n)

plot(Modeldf$PerTotLoss~Modeldf$mitTF)

# Determining the number of days since the last inspection:
Modeldf$Date <- as.Date(Modeldf$LastInspectionDate, "%d-%b-%y")
Modeldf$DaysInspection <- as.Date("2018-08-10")-Modeldf$Date 
plot(Modeldf$PerTotLoss~Modeldf$DaysInspection)

Modeldf$DaysInspection<-as.numeric(Modeldf$DaysInspection)

Modeldf$DaysInspection
