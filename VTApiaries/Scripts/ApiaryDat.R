###########################################################################################
# Vermont Apiary Inspeciton Data
# Samantha Alger
# June 10, 2018
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/VTApiaries/")

# source my packages
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(tidyr)
library(slam)
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

###########################################################
# Read in Data
FullApiaryDat <- read.csv("CSV_files/VTApiaries.csv", 
                        header=TRUE, 
                        sep = ",", 
                        stringsAsFactors = FALSE)

ApiaryDat <- read.csv("CSV_files/singles.csv", 
                      header=TRUE, 
                      sep = ",", 
                      stringsAsFactors = FALSE)

states <- read.csv("CSV_files/USstates.csv", 
                   header=TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)

#Remove duplicate data so that there is only a single row for each beekeeper (to be used in all analyses except apiary-level analyses and vendor information)

Singles <- FullApiaryDat[!duplicated(FullApiaryDat$BeekeeperID), ]
table(Singles$BeekeeperID)

write.csv(Singles, file = "singles.csv")

#########################################################
#########################################################
# Creating a word cloud from the 'opinion' column:
#########################################################
#########################################################

text <- readLines("CSV_files/WordCloud.txt")
docs <- Corpus(VectorSource(text))
#docs <- tm_map(docs, content_transformer(tolower))
inspect(docs)
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word

# specify your stopwords as a character vector
docs <- tm_map(docs, removeWords, c("bees", "year", "keeping","will", "use", "enough", "getting", "one", "many", "just", "last", "keep", "good", "one", "also", "sure", "bee", "lost","hives","hive","left", "hard", "proper", "long", "someone", "need", "every", "make", "low", "late", "right", "finding", "two", "area", "best", "put", "well", "using", "trying", "due", "lot", "start", "find", "may", "going", "even", "get", "past", "biggest","check", "another", "much", "purchased","now", "can"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
# docs <- tm_map(docs, stemDocument)
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)

#Generating word cloud
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 5,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

#########################################################
#########################################################
# Mite Monitoring
#########################################################
#########################################################
##########################################################
# Percentage of beekeepers that count mites:

#Create a Pie chart for mite counts, (yes/no)
table(ApiaryDat$MiteCounts)

mitedf <- data.frame(
  group = c("Counted mites", "Did not count mites"),
  value = c(232, 403)
)
head(mitedf)

#     group value
# 1 Counted   232....36.5%
# 2 Did not   403...63.5%

library(ggplot2)
library(scales)
# Barplot
bp<- ggplot(mitedf, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")

#Create a blank theme
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#Piechart
pie <- bp + coord_polar("y", start=0) +
  theme_minimal() +  blank_theme +
  theme(axis.text.x=element_blank()) + guides(fill=guide_legend(title="Percentage of beekeepers"))
pie
#######################################################
#How are beekeepers monitoring for mites?

#Import data:
MiteMon<- Singles

#Create two separate dfs for analysis with specified columns:
MonMethods <- dplyr::select(MiteMon, BeekeeperID,SugarShakeYN, AlcoholWashYN, BottomBoardYN, DroneSurveyYN, OtherMiteCountYN)

MonFreq <- dplyr::select(MiteMon, BeekeeperID, SugarShake, AlcoholWash,BottomBoard,DroneSurvey, OtherMiteCount)

OtherMiteCountSpecify <- dplyr::select(CensusDat,OtherMiteCountSpecify)

# Change dataframes from wide to long format
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

MonMethods <- gather(MonMethods, question, response, SugarShakeYN:OtherMiteCountYN, factor_key=TRUE)

MonFreq <- gather(MonFreq, question, response, SugarShake:OtherMiteCount, factor_key=TRUE)

#Convert all True and False to '0' and '1' in the 'ReasonLoss' dataframe
MonMethods$response<-as.integer(as.logical(MonMethods$response))

# Preparing Data for Bar Plot
MiteMon <- ddply(MonMethods, c("question"), summarise, 
                   n = length(response),
                   mean = mean(response, na.rm=TRUE),
                   sd = sd(response, na.rm=TRUE),
                   se = sd / sqrt(n))

plot1 <- ggplot(MiteMon, aes(x=question, y=mean, fill=question)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + 
  labs(x="Mite Monitoring Method", y = "Frequency") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none") + 
  scale_y_continuous(labels = scales::percent) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + 
                      se, width = 0.2)) + scale_fill_brewer() +
  scale_x_discrete(labels=c("SugarShakeYN" = "Sugar Shake", "AlcoholWashYN" = "Alcohol Wash", "BottomBoardYN" = "Bottom Board", "DroneSurveyYN"= "Drone Survey", "OtherMiteCountYN"= "Other"))

plot1

################################################################
# State of Origin for bees purchased
# DataCleaning....
OriginData <- FullApiaryDat

#Get libraries
library(usmap)
library(maps)
library(ggplot2)

#download fips codes dataset
data(state.fips)


# Remove all rows without vendor listed:
OriginData <- OriginData[! is.na(OriginData$BeePurchaseVendorName), ]

# Remove duplicates based on beekeeper ID and vendor, there are some beekeeper IDS associated with multiple vendors if beekeepers purchased bees from multiple places!
OriginData <- OriginData[!duplicated(OriginData[,c("BeekeeperID","BeePurchaseVendorName")]),]

# Change column name state of origin -> StateID
names(OriginData)[names(OriginData) == "BeePurchaseStateOfOrigin"] <- "StateID"


# merge dataframes so that states are listed with their 'fips IDs' (needed to map)
OriginData <- merge(OriginData, states, by= "StateID")
OriginData <- merge(OriginData, state.fips, by = "abb")


#subset specific columns for analysis:
OriginData <- dplyr::select(OriginData, StateID,BeekeeperID,BeePurchase,BeePurchaseVendorName, abb, fips, BeePurchasePackages,BeePurchaseNucleusColonies,BeePurchaseFullColonies,BeePurchaseQueens,BeePurchaseOther,OtherBeePurchaseSpecify)

# Subset data for mapping
Origin <- ddply(OriginData, c("fips","abb"), summarise, 
                     n = length(abb))

#Create Map:
usmap::plot_usmap(data = Origin, values = "n", lines = "black") + 
  scale_fill_continuous(low = "white", high = "blue", name = "Honey Bee Purchases", label = scales::comma) + 
  theme(legend.position = "right")

# Histogram figure for Vendors, filled by state of origin
Vendor <- ddply(OriginData, c("BeePurchaseVendorName","abb"), summarise, 
                n = length(BeePurchaseVendorName))

VendorPlot <- ggplot(data=Vendor, aes(x=BeePurchaseVendorName, y=n, fill = abb)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Medications", y = "Colony Loss")


VendorPlot
#+ scale_x_discrete(labels=c("Formic acid (Mite Away Quick Strips)" = "MAQS", "Amitraz, Formic acid (Mite Away Quick Strips)" = "Amitraz, MAQS", "Apiguard, Formic acid (Mite Away Quick Strips)" = "Apiguard, MAQS", "Formic acid (Mite Away Quick Strips), Oxalic acid" = "Oxalic acid, MAQS", "Formic acid (Mite Away Quick Strips), Honey B Healthy" = "Honey B Healthy, MAQS")) + scale_y_continuous(labels = scales::percent) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2)) + scale_fill_brewer()
####################################################################
# MEDICATIONS:


###########################################################

#Remove rows with NA in percent loss:
FullApiaryDat <- FullApiaryDat[! is.na(FullApiaryDat$PerWintLoss), ]


# reorder factors for plotting
FullApiaryDat$TreatmentsForColonyHealthMitigation <- factor(FullApiaryDat$TreatmentsForColonyHealthMitigation, levels = c("None", "Powdered sugar", "Honey B Healthy", "Apiguard", "Api-life VAR", "Formic acid (Mite Away Quick Strips)", "Oxalic acid", "Formic acid (Mite Away Quick Strips) Honey B Healthy", "Apiguard, Formic acid (Mite Away Quick Strips)", "Formic acid (Mite Away Quick Strips), Oxalic acid", "Amitraz, Formic acid (Mite Away Quick Strips)"))

Medications <- ddply(FullApiaryDat, c("TreatmentsForColonyHealthMitigation"), summarise, 
                n = length(PerTotLoss),
                mean = mean(PerTotLoss, na.rm=TRUE),
                sd = sd(PerTotLoss, na.rm=TRUE),
                se = sd / sqrt(n))

#View the table to see average colony losses by all treatment combinations
#View(Medications)

# Reduce Medication combinations to only include those that appear more than 10 times:
Medications <- Medications[ which(Medications$n > 9), ]
Medications <- Medications[! is.na(Medications$TreatmentsForColonyHealthMitigation), ]

#choosing color pallet
#colors <- c("goldenrod", "violetred4", "snow1", "black", "green")


#Create a bar graph for viruses loads per Koppert colony:
plot1 <- ggplot(Medications, aes(x=TreatmentsForColonyHealthMitigation, y=mean, fill=TreatmentsForColonyHealthMitigation)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x="Medications", y = "Colony Loss") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none") + scale_x_discrete(labels=c("Formic acid (Mite Away Quick Strips)" = "MAQS", "Amitraz, Formic acid (Mite Away Quick Strips)" = "Amitraz, MAQS", "Apiguard, Formic acid (Mite Away Quick Strips)" = "Apiguard, MAQS", "Formic acid (Mite Away Quick Strips), Oxalic acid" = "Oxalic acid, MAQS", "Formic acid (Mite Away Quick Strips), Honey B Healthy" = "Honey B Healthy, MAQS")) + scale_y_continuous(labels = scales::percent) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2)) + scale_fill_brewer()

plot1

#Prepare raw data set for ANOVA:

# Subset data in R

MedANOVA <- subset(FullApiaryDat, TreatmentsForColonyHealthMitigation == "None"| TreatmentsForColonyHealthMitigation =="Powdered sugar" | TreatmentsForColonyHealthMitigation == "Honey B Healthy"| TreatmentsForColonyHealthMitigation == "Apiguard" | TreatmentsForColonyHealthMitigation == "Api-life VAR" |  TreatmentsForColonyHealthMitigation == "Formic acid (Mite Away Quick Strips)" | TreatmentsForColonyHealthMitigation == "Oxalic acid" | TreatmentsForColonyHealthMitigation == "Formic acid (Mite Away Quick Strips) Honey B Healthy" | TreatmentsForColonyHealthMitigation == "Apiguard Formic acid (Mite Away Quick Strips)" | TreatmentsForColonyHealthMitigation == "Formic acid (Mite Away Quick Strips), Oxalic acid" | TreatmentsForColonyHealthMitigation == "Amitraz, Formic acid (Mite Away Quick Strips)")

#remove 'NA's from percent total loss column
MedANOVA <- MedANOVA[! is.na(MedANOVA$PerTotLoss), ]

#Perform ANOVA and post hoc test:
Meds <- aov(MedANOVA$PerTotLoss~MedANOVA$TreatmentsForColonyHealthMitigation)
summary(Meds)
TukeyHSD(Meds)
########################################################
# Reasons for Colony Losses:

#Import data:
CensusDat<- Singles

#Create two separate dfs for analysis with specified columns:
ReasonLoss <- dplyr::select(CensusDat, BeekeeperID, ColonyLossVarroaMiteYN, ColonyLossStarvationYN, ColonyLossBearsYN, ColonyLossAmericanFoulbroodYN, ColonyLossSwarmingYN, ColonyLossPesticidesYN, ColonyLossMitacidesYN, OtherColonyLossYN)

NumLoss <- dplyr::select(CensusDat, BeekeeperID, ColonyLossVarroaMite, ColonyLossStarvation, ColonyLossBears, ColonyLossAmericanFoulbrood, ColonyLossSwarming, ColonyLossPesticides, ColonyLossMitacides, OtherColonyLoss)

OtherLoss <- dplyr::select(CensusDat,OtherColonyLossCausesSpecify)

# Change dataframes from wide to long format
# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

ReasonLoss <- gather(ReasonLoss, question, response, ColonyLossVarroaMiteYN:OtherColonyLossYN, factor_key=TRUE)

NumLoss <- gather(NumLoss, question, response, ColonyLossVarroaMite:OtherColonyLoss, factor_key=TRUE)

#Convert all True and False to '0' and '1' in the 'ReasonLoss' dataframe
ReasonLoss$response<-as.integer(as.logical(ReasonLoss$response))

# Preparing Data for Bar Plot
LossCause <- ddply(ReasonLoss, c("question"), summarise, 
                     n = length(response),
                     mean = mean(response, na.rm=TRUE),
                     sd = sd(response, na.rm=TRUE),
                     se = sd / sqrt(n))

plot1 <- ggplot(LossCause, aes(x=question, y=mean, fill=question)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + 
    labs(x="Causes", y = "Colony Loss") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none") + 
  scale_y_continuous(labels = scales::percent) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + 
                      se, width = 0.2)) + scale_fill_brewer() +
  scale_x_discrete(labels=c("ColonyLossVarroaMiteYN" = "Varroa", "ColonyLossStarvationYN"= "Starvation", "ColonyLossBearsYN" = "Bears", "ColonyLossAmericanFoulbroodYN" = "AFB", "ColonyLossSwarmingYN" = "Swarming", "ColonyLossPesticidesYN" = "Pesticides", "ColonyLossMitacidesYN" = "Mitacides", "OtherColonyLossYN"= "Other"))

plot1

######################################################
# Create Word Cloud for "OtherLoss"

#Combine all Text:
OtherLossText <- paste(unlist(OtherLoss), collapse =" ")
# I then copied and pasted the above text to a .txt file and saved it in .csv folder. Is there a way to do this in R?

text <- readLines("CSV_files/WordCloud2.txt")
docs <- Corpus(VectorSource(text))
#docs <- tm_map(docs, content_transformer(tolower))
inspect(docs)
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word
docs <- tm_map(docs, removeWords, c("they", "due", "left","sure","know", "other", "had", "were", "into", "with", "all", "going", "was","very","over","from"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
# docs <- tm_map(docs, stemDocument)
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)

#Generating word cloud
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 5,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

##################################################
# Antibiotic Use
Antibx <- Singles
AntibxFull <- FullApiaryDat

table(Antibx$AntibioticTreatments)
#antibiotic use: 53 of 582 beekeepers (~10 %)
table(AntibxFull$AntibioticTreatments)
# of apiaries = 66
SumAntibx <- AntibxFull[which(AntibxFull$AntibioticTreatments == "TRUE"),]
sum(SumAntibx$ColonyCount)
# 356 hives
