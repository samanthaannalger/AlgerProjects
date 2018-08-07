# Data Analysis for Apiary Inspection Programs
# Samantha Alger, Faith Novella
#

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/VTApiaries/CSV_files/")

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
InspectSurvey <- read.csv("InspectorSurvey.csv", 
                         header=TRUE, 
                         sep = ",", 
                         stringsAsFactors = FALSE)
InspectCloud <- read.csv("InspectWordCloud.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)
#########################################################
#########################################################
# Creating a word cloud from the 'opinion' column:
#########################################################
#########################################################

#Combine all Text:
OpinionText <- paste(unlist(InspectCloud$BeeChallenges), collapse =" ")
# I then copied and pasted the above text to a .txt file and saved it in .csv folder. Is there a way to do this in R?
#text <- readLines("CSV_files/WordCloud2.txt")

docs <- Corpus(VectorSource(OpinionText))

#text <- readLines("CSV_files/WordCloud.txt")
#docs <- Corpus(VectorSource(text))
docs <- tm_map(docs, content_transformer(tolower))
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
docs <- tm_map(docs, removeWords, c("michigan"))
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
wordcloud(words = d$word, freq = d$freq, min.freq = 2,
          max.words=200, random.order=FALSE, rot.per=0.4, 
          colors=brewer.pal(8, "Dark2"))

##################################################

table(InspectSurvey$Registration)

6/16
4/16
6/16
3/16

0.42+0.42+0.28
