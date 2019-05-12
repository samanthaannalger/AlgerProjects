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
install.packages("rgdal")
install.packages("raster")

vtMap<-get_map(location=c(-72.6954,43.8),zoom=8, scale="auto", source="google",maptype="satellite")
ggmap(vtMap) 

ggmap(vtMap) +
  stat_density2d(
    aes(x=long, y=lat, fill=..level..),
    data=apiaries,
    geom="polygon",
    bins=20) +
  scale_fill_gradient(low = "green", high = "red", guide=FALSE) + 
  #Adding points
  geom_point(aes(x = long, y = lat), data = apiaries, alpha = .5, color = "blue", size = 1) +
  theme( #removing grid lines:
    panel.grid.major = element_blank (), # remove major grid
    panel.grid.minor = element_blank (),  # remove minor grid
    axis.text = element_blank (), 
    axis.title = element_blank (),
    axis.ticks = element_blank ()
  )

library(rgdal)
library(raster)
library(ggmap)

  vt <- map("state", "VERMONT")
  vt <- ggmap_rast(map = vt)
  vt
  

  
rastervt <- crop(rasterVt, vt)
  