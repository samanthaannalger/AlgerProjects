# ggmap, point density layer
# 2017, March, 9
# SAA


# Clear memory of characters:
ls()
rm(list=ls())

#Set Working Directory: 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files")

# read in data:
SpacDF <-read.table("2015SurveySpatial.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
apiaries <-read.csv("ApiaryData.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

#get required libraries
library(ggmap)
library(ggplot2)
library(maps)

vtMap<-get_map(location=c(-72.6954,43.8),zoom=8, scale="auto", source="google",maptype="satellite")
ggmap(vtMap) 

ggmap(vtMap) +
  stat_density2d(
    aes(x=long, y=lat, fill=..level..),
    data=apiaries,
    geom="polygon",
    bins=30) +
    scale_fill_gradient(low="gray", high="red")
  
  
