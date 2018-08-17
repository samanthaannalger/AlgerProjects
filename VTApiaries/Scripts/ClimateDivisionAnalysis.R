# Analyses for ApiaryDat Manuscript
# Alger & Burnham
# 7/29/18

# Clear memory of characters:
rm(list=ls())

# source my packages
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(data.table)
library(rgeos)
library(sf)
library(maps)
library(maptools)
library(sp)
library(rgdal)
# Set working directory:
setwd("~/AlgerProjects/VTApiaries/CSV_files/")

Shinydf <- read.csv("Shinydf.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

Clim <- read.csv("DivDat.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE,
                    skip = 8)

# GeoDataConverter:
# https://mygeodata.cloud/converter/

#read in data:
division <- read_sf("~/AlgerProjects/VTApiaries/CSV_files/CONUS_CLIMATE_DIVISIONS.shp")


# https://stackoverflow.com/questions/13316185/r-convert-zipcode-or-lat-long-to-county

# Pull the lat and long coordinates out of Shinydf
coords<- dplyr::select(Shinydf, Longtitude,Latitude)

#remove NAs
coords <- coords[complete.cases(coords),]

# Begin spatial data analysis...
points <- SpatialPoints(coords)
#SpatialPolygonDataFrame - I'm using a shapefile of climatological divisions
counties <- rgdal::readOGR(dsn="CONUS_CLIMATE_DIVISIONS.shp")
#select Vermont climatological divisions only based on state 'FIPS' code.
counties <- counties[counties$STATE_FIPS == 50, ]
#assume same proj as shapefile!
proj4string(points) <- proj4string(counties)
#get division polygon point is in
Div <- as.character(over(points, counties)$NAME)

#Bind it back with the coordinates
Div<-cbind(coords, Div)

#merge it back with the original df
ClimDiv <- merge(Div, Shinydf, by = c("Latitude","Longtitude"))

# Now need to add climatological data:

# Testing to see if colony loss differed by climatological divsision.
mod<-aov(ClimDiv$PerTotLoss~ClimDiv$Div)
summary(mod)
mod
TukeyHSD(mod)

