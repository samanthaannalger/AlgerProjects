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

HistDat <- read.csv("CSV_files/histDat.csv", 
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

RegData <- read.csv("CSV_files/RegActiveAndDelinquent.csv", 
                   header=TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE)
ApiaryRegs2017 <- read.csv("CSV_files/2017ApiaryRegs.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)


#Remove duplicate data so that there is only a single row for each beekeeper (to be used in all analyses except apiary-level analyses and vendor information)

Singles <- FullApiaryDat[!duplicated(FullApiaryDat$BeekeeperID), ]

#write.csv(Singles, file = "singles.csv")

##################################################
#########################################################
#########################################################
# Basic Statistics:
#########################################################
#########################################################


# Calculating state level losses:
AnnLoss<- RegData[! is.na(RegData$PerTotLoss), ]
mean(AnnLoss$PerTotLoss)
# 38.8%: Total Annual

WinterLoss<- RegData[! is.na(RegData$PerWinLoss), ]
mean(WinterLoss$PerWinLoss)
# 35.8% Winter Loss
# Summer loss??

#Calculating apiary and colony totals:

# Subset to only include status="Active" beekeepers (remove NAs too)
ActiveOnly<- RegData[! is.na(RegData$BeeKeeperStatus), ]
ActiveOnly<- ActiveOnly[ which(ActiveOnly$BeeKeeperStatus == "Active"), ]
nrow(ActiveOnly)
# 1091 Registered Apiaries

#How many apiaries are delinquent?
DelOnly<- RegData[ which(RegData$BeeKeeperStatus == "Delinquent"), ]
nrow(DelOnly)
# 465 delinquent apiaries

#How many registered Colonies?
ActiveOnly<- ActiveOnly[! is.na(ActiveOnly$ColonyCount), ]
sum(ActiveOnly$ColonyCount)
# 8,450 registered colonies

#Calculate the number of registered beekeepers:
ActiveBeeks <- ActiveOnly[!duplicated(ActiveOnly$BeekeeperID), ]
nrow(ActiveBeeks)
# 743 Active registered beekeepers
DelBeeks <- DelOnly[!duplicated(DelOnly$BeekeeperID), ]
nrow(DelBeeks)
#254 Delinquent beekeepers


#Determine the number of apiaries each beekeeper has registered and create a histogram....

#For registered beekeepers:
library(data.table)
histDat = data.table(ActiveOnly)
histDat[, `n` := .N, by = BeekeeperID]
#histDat <-  histDat[!duplicated(histDat$BeekeeperID),]
# Histogram showing the number of apiaries belonging to beekeepers with 1 yard, two yards...etc...
hist(histDat$n, freq=TRUE, breaks=25)

###############################################

#Code apiaries as owned by either hobbists, sideliners, or commericial apiaries

histDat$Beektype <- ifelse(histDat$n == 1,"Hobbyist", ifelse(histDat$n <=5, "Sideliner", "Commercial"))

# Write.csv hist dat for shiny app demo
#write.csv(histDat,"histDat.csv")

BeekTypeStats <- ddply(histDat, c("Beektype"), summarise, 
                     apiaries = length(n),
                     colonies = sum(ColonyCount),
                     loss = mean(PerTotLoss, na.rm = TRUE),
                     sd = sd(PerTotLoss, na.rm=TRUE),
                     se = sd / sqrt(apiaries))

BeekTypeStats$perApiaries <- BeekTypeStats$apiaries/sum(BeekTypeStats$apiaries)

BeekTypeStats$perColonies <- BeekTypeStats$colonies/sum(BeekTypeStats$colonies)

BeekTypeStats

BeekTypeDF <- rbind(BeekTypeStats, BeekTypeStats)
MeasureType <- c(rep("Apiary", 3), rep("Colony", 3))
BeekTypeDF <- cbind(BeekTypeDF, MeasureType)
BeekTypeDF$Percent <- c(BeekTypeDF$perApiaries[1:3], BeekTypeDF$perColonies[4:6])

BeekTypeDF
#Are there significant differences in hive losses among beekeeper types?
#Currenlty no significant differences! p = 0.09
BeekTypeLoss <- aov(histDat$PerTotLoss~histDat$Beektype)
summary(BeekTypeLoss)
BeekTypeLoss
# Check to see which beekeepers told us about hive losses... #MISSING LOTS OF COLONY LOSS DATA FROM COMMERICIAL BEEKEEPERS!!
NAcheck<-histDat[is.na(histDat$PerTotLoss),]
table(NAcheck$Beektype)


#Create figure showing number of colonies and number of apiaries as BeekType.
#Create labels.....

#Reorder factors:
BeekTypeDF$Beektype <- factor(BeekTypeDF$Beektype, levels = c("Sideliner", "Hobbyist", "Commercial"))

BeekTypeDF <- ddply(BeekTypeDF, .(MeasureType), transform, pos = cumsum(Percent) - (0.5 * Percent))

BeekTypeFig <- ggplot(BeekTypeDF, aes(y = Percent, x = MeasureType , fill = Beektype)) + geom_bar(stat="identity") + theme_classic() + labs(x=NULL, y = "Percent of total in Vermont")+ scale_y_continuous(labels = scales::percent) + 
scale_fill_brewer() +
  scale_x_discrete(labels=c("Apiary" = "Apiaries", "Colony" = "Colonies")) + guides(fill=guide_legend(title="Beekeeper Type")) + geom_text(data=BeekTypeDF, aes(x = MeasureType, y = pos ,label = paste0(round(Percent*100, digits = 1),"%")), size=4)

BeekTypeFig


#########################################################
#########################################################
# Creating a word cloud from the 'opinion' column:
#########################################################
#########################################################

#Combine all Text:
OpinionText <- paste(unlist(ApiaryDat$Opinion), collapse =" ")
# I then copied and pasted the above text to a .txt file and saved it in .csv folder. Is there a way to do this in R?
#text <- readLines("CSV_files/WordCloud2.txt")

docs <- Corpus(VectorSource(OpinionText))

text <- readLines("CSV_files/WordCloud.txt")
#docs <- Corpus(VectorSource(text))
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
wordcloud(words = d$word, freq = d$freq, min.freq = 10,
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

mitedf <- data.frame(
  group = c("Did not count mites", "Counted mites"),
  value = c(table(ApiaryDat$MiteCounts))
)
head(mitedf)

#     group value
# 1 Counted   221....34.8%
# 2 Did not   414...65.2%

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
  labs(x="Mite Monitoring Method", y = "% Reported Use") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none", axis.text=element_text(size=15), axis.title=element_text(size=18,face="bold")) + 
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


####################################################################
# MEDICATIONS:


###########################################################
# What 'meds' do beekeepers use in their hives?
MedUse <- Singles[! is.na(Singles$TreatmentsForColonyHealthMitigation), ]

# What 'meds' do beekeepers use in their hives?
MedUse <- Singles[! is.na(Singles$TreatmentsForColonyHealthMitigation), ]

vec <- c("None","Fluvalinate","Coumaphos","Amitraz","Apiguard","Api-life VAR","Sucrocide","Powdered sugar","Oxalic acid","Formic acid (Mite Away Quick Strips)","Menthol","Hopguard","Honey B Healthy","Fumagillin-B","Tylosin (Tylan, Tylosin, Tylovet)","Lincomycin (Lincomix)","Oxytetracycline (TM, OXTC, Pennox, Terramycin)","Herbal antibiotics")


mat <- matrix(nrow = length(MedUse$TreatmentsForColonyHealthMitigation), ncol = length(vec))

for (i in 1:length(vec)){
  
  mat[,i] <- grepl(vec[i], MedUse$TreatmentsForColonyHealthMitigation, fixed = TRUE)
  
}

# turn mat into a data frame, give it column names and merge with MedUse
mat <- data.frame(mat)
names(mat) <- vec
MedUse <- cbind(MedUse, mat)

TreatPercents <- colMeans(mat, na.rm=TRUE)
Treatments <- vec

MedUse <- data.frame(Treatments, TreatPercents)


#Reorder for plotting:
MedUse$Treatments <- factor(MedUse$Treatments, levels = c("None", "Powdered sugar", "Honey B Healthy", "Apiguard", "Api-life VAR", "Menthol","Sucrocide","Hopguard", "Formic acid (Mite Away Quick Strips)", "Oxalic acid", "Amitraz","Fluvalinate","Coumaphos","Tylosin (Tylan, Tylosin, Tylovet)","Lincomycin (Lincomix)","Oxytetracycline (TM, OXTC, Pennox, Terramycin)","Herbal antibiotics", "Fumagillin-B"))


plot2 <- ggplot(MedUse, aes(x=Treatments, y=TreatPercents)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x=NULL, y = "% Reported Use") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")  + scale_x_discrete(labels=c("None","Fluvalinate"="Fluvalinate (Apistan)","Coumaphos","Amitraz"="Amitraz(Apivar)","Apiguard","Api-life VAR","Sucrocide","Powdered sugar","Oxalic acid","Formic acid (Mite Away Quick Strips)"="Formic acid (MAQS)","Menthol","Hopguard","Honey B Healthy","Fumagillin-B","Tylosin (Tylan, Tylosin, Tylovet)"= "Tylosin","Lincomycin (Lincomix)"= "Lincomycin","Oxytetracycline (TM, OXTC, Pennox, Terramycin)"="Oxytetracycline","Herbal antibiotics")) + scale_fill_brewer() + scale_y_continuous(labels = scales::percent) + theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"))

plot2


################################################
# How Does med use correlate with colony losses....



#Remove rows with NA in percent loss:
MedData <- FullApiaryDat[! is.na(FullApiaryDat$PerTotLoss), ]


# reorder factors for plotting
MedData$TreatmentsForColonyHealthMitigation <- factor(MedData$TreatmentsForColonyHealthMitigation, levels = c("None", "Powdered sugar", "Honey B Healthy", "Apiguard", "Api-life VAR", "Formic acid (Mite Away Quick Strips)", "Oxalic acid", "Formic acid (Mite Away Quick Strips) Honey B Healthy", "Apiguard, Formic acid (Mite Away Quick Strips)", "Formic acid (Mite Away Quick Strips), Oxalic acid", "Amitraz, Formic acid (Mite Away Quick Strips)"))


Medications <- ddply(MedData, c("TreatmentsForColonyHealthMitigation"), summarise, 
                n = length(PerTotLoss),
                mean = mean(PerTotLoss, na.rm=TRUE),
                sd = sd(PerTotLoss, na.rm=TRUE),
                se = sd / sqrt(n))

#View the table to see average colony losses by all treatment combinations
#View(Medications)

# Reduce Medication combinations to only include those that appear more than 10 times:
Medications <- Medications[ which(Medications$n > 9), ]
Medications <- Medications[! is.na(Medications$TreatmentsForColonyHealthMitigation), ]


#Create a bar graph
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
    labs(x="Causes", y = "% Reported Causes") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none", axis.text=element_text(size=15), axis.title=element_text(size=18,face="bold")) + 
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
#text <- readLines("CSV_files/WordCloud2.txt")

docs <- Corpus(VectorSource(OtherLossText))
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


#######################################################
# How do you make up for colony losses?
#######################################################

ColLoss <- Singles[! is.na(Singles$ColonyLossMakeup), ]

vec <- c("Make Splits or Divides","Purchase Colonies in Hives","Purchase Packages","Purchase Nucleus Colonies")


mat <- matrix(nrow = length(ColLoss$ColonyLossMakeup), ncol = length(vec))

for (i in 1:length(vec)){
  
  mat[,i] <- grepl(vec[i], ColLoss$ColonyLossMakeup, fixed = TRUE)
  
}

# turn mat into a data frame, give it column names and merge with MedUse
mat <- data.frame(mat)
names(mat) <- vec
ColLoss <- cbind(ColLoss, mat)

Percents <- colMeans(mat, na.rm=TRUE)
Treatments <- vec

ColLoss <- data.frame(Treatments, Percents)


#Reorder for plotting:
ColLoss$Treatments <- factor(ColLoss$Treatments, levels = c("Make Splits or Divides","Purchase Colonies in Hives","Purchase Packages","Purchase Nucleus Colonies"))


plot3 <- ggplot(ColLoss, aes(x=Treatments, y=Percents)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(x=NULL, y = "% Responses") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position="none")  + scale_x_discrete(labels=c("Splits/Divides","Purchase Full Hives","Purchase Packages","Purchase Nucs")) + scale_fill_brewer() + scale_y_continuous(labels = scales::percent) + theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"))

plot3

##############################################################
# SUBSET Spacial data by lat long Euclidean Distance for mapping purposes:

# Data cleaning to make final df:
histDat = data.table(RegData)
histDat[, `n` := .N, by = BeekeeperID]
#histDat <-  histDat[!duplicated(histDat$BeekeeperID),]

histDat$Beektype <- ifelse(histDat$n == 1,"Hobbyist", ifelse(histDat$n <=5, "Sideliner", "Commercial"))

####### 
#Merging the two dataframes for shiny app:

#select only columns we need:
histDat <- dplyr::select(histDat, LocationID, Beektype, AccountName, BeeKeeperStatus, n)

FullApiaryDat <- dplyr::select(FullApiaryDat, -AccountName, -BeeKeeperStatus)

Shinydf <- merge.data.frame(FullApiaryDat,histDat, by = "LocationID", all.y = TRUE)

###########################################################
# Begin Functions
###########################################################


####################################################################
# function name: LatLongMat
# description: Creates matrix with lat and long columns
# parameters: 
# data = data frame that includes "Latitude" and "Longitude"
# returns: a matrix with two columns (lat and long)
####################################################################

LatLongMat <- function(data=data){
# create matrix of lats and longs from data set
Latitude <- data$Latitude
Longtitude <- data$Longtitude
x <- cbind(Longtitude, Latitude)
x <- x[complete.cases(x), ]
x <- as.matrix(x)

return(x)

}

#####################################################################
# END OF FUNCTION
####################################################################

####################################################################
# function name: SubSetMap
# description: uses lat long matrix from "LatLongMat" to calculate distance in miles from center point (lat long) and find all apiaries within that radius
# parameters: 
# data = data frame that includes "Latitude" and "Longitude"
# rad = radius to query (numeric (miles))
# lat = center point latitude
# long = center point longtitude
# matrix = matrix of lats and longs in two columns
# returns: a list with 4 elements: data frame where all rows are apairies within rad, rad, lat and long
####################################################################

SubSetMap <- function(data = data, 
                      rad = 100, 
                      lat = 42, 
                      long = -72,
                      matrix = x){

library(geosphere)
# use the sp package to determine Euc. Dist between points in matrix "y" and central point "x"
m <- distm(x = c(lat, long), y = matrix, fun = distHaversine)
m <- as.vector(m)
distance <- m/1609.334

# merge data back to original data frame:
matrix <- as.data.frame(matrix)
temp <- as.data.frame(cbind(matrix$Latitude, matrix$Longtitude, distance))
names(temp) <- c("Latitude", "Longtitude", "distance")

t <- as.data.frame(merge(x=data, y=temp, by = c("Latitude", "Longtitude"), all.y=TRUE, sort=FALSE))

# which are within radius
queryDF <- t[t$distance<=rad,]

return(list(queryDF, rad, lat, long))

}

#####################################################################
# END OF FUNCTION
####################################################################

####################################################################
# function name: Mapfunc
# description:plots all apiaries from subsetted list within a certain mile radious from a given point, provides drop down information for apiaries when point is clicked radius
# parameters: 
# data = data frame that includes "Latitude" and "Longitude"
# rad = radius to query (numeric (miles))
# lat = center point latitude
# long = center point longtitude
# returns: a map with center point marked and all apiaries surrounding that point. interactive.
####################################################################

Mapfunc <- function(data=data, rad, lat, long) {
  library(leaflet)
  
  content <- paste("Account Name:", data$AccountName, "<br/>", 
                   "BeekeeperID:", data$BeekeeperID, "<br/>",
                   "Status:", data$BeeKeeperStatus, "<br/>",
                   "# Colonies:", data$ColonyCount, "<br/>",
                   "Annual Loss:", data$PerTotLoss*100, "%","<br/>",
                   "Beekeeper Type:", data$Beektype, "<br/>",
                   "Last Inspection Date:", data$LastInspectionDate)
m <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group="background 1") %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=long , lat=lat,
             popup="Outbreak") %>% setView(long, lat, zoom = 8) %>% addCircles(data$Longtitude,data$Latitude, popup=content, weight = 3, radius=40, 
                 color="#ffa500", stroke = TRUE, fillOpacity = 0.8) 

return(m)
}
#####################################################################
# END OF FUNCTION
####################################################################


############################################################
############ PROGRAM BODY ##################################
############################################################

LLmat <- LatLongMat(data = Shinydf)

SSdat <- SubSetMap(data = Shinydf, rad = 20, lat = -72.746286, long = 44.278876, matrix = LLmat)

Mapfunc(data=SSdat[[1]], rad= SSdat[[2]], lat = SSdat[[4]], long= SSdat[[3]])










# SHINY CODE:



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$map <- renderLeaflet({
    SSdat <- SubSetMap(data = Shinydf, rad = 20, lat = -72.746286, long = 44.278876, matrix = LLmat)
    
    Mapfunc(data=SSdat[[1]], rad= SSdat[[2]], lat = SSdat[[4]], long= SSdat[[3]]) 
  })
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Apiary Locator"),
  
  # Sidebar with a slider input for distance to central point 
  sidebarLayout(
    sidebarPanel(
      sliderInput("distance",
                  "Radius (Miles):",
                  min = 1,
                  max = 4000,
                  value = 30)
    ),
    
    # Show a map plot
    mainPanel(
      plotOutput("map")
    )
    
    
    
    
    br(),
    leafletOutput("map", height="600px"),
    absolutePanel(top=20, left=70, textInput("target_zone", "" , "Ex: Burlington, Vermont")),
    br()
    )
  )
)



# Run the application 
shinyApp(ui = ui, server = server)



