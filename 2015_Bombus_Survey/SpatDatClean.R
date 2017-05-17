# Spatial Data Cleaning for 2015 Survey
# 2017-04-05
# SAA

#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())

#
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# Read in data
SpatDatBuffs<-read.csv("SpatDatBuffs.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
MixedModelDF <- read.csv("MixedModelDF.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

#replace Nulls with 0
SpatDatBuffs[SpatDatBuffs == "<Null>"] <- 0

#write file into the CSV folder to be used in mixed model analyses:
write.csv(SpatDatBuffs, file = "ApiaryBuffers.csv")
