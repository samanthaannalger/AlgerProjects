# BombSurv Map showing hot spots for DWV and BQCV in bumble bees
# Salger
# 8-5-2018


# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/AlgerProjects/2015_Bombus_Survey/CSV_Files/")

# source my packages
library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(rgdal)


###########################################################
# Read in Data
SurvMap <- read.csv("spatialMerge.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)
SurvMap$x <- SurvMap$long
SurvMap$y <- SurvMap$lat

MapSplit<- split(SurvMap, SurvMap$target_name)
coordinates(MapSplit$BQCV) = ~ x+ y
plot(MapSplit$BQCV)

x.range <- as.numeric(c(-74.7, -70.5))
y.range <- as.numeric(c(44, 45))

grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.05), y = seq(from = y.range[1], 
to = y.range[2], by = 0.05))  # expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE

plot(grd,cex=1.5, col="grey")
points(MapSplit$BQCV, pch = 1, col = "red", cex =1)

summary(MapSplit$BQCV$lat)
summary(MapSplit$BQCV$long)

idw <- krige(formula = MapSplit$BQCV$BombusViralLoad~1, locations = MapSplit$BQCV, newdata = grd) 

idw.output = as.data.frame(idw)  # output is defined as a data table
names(idw.output)[1:3] <- c("long", "lat", "var1.pred")  # give names to the modelled variables
ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = log(var1.pred))) +   geom_point(data = SurvMap, aes(x = long, y = lat), shape = 21,  colour = "red")

VT_contour <- rgdal::readOGR(dsn="~/AlgerProjects/2015_Bombus_Survey/CSV_Files/VT.shp")
est_contour <- fortify(est_contour, region = "name")


idw <- krige(formula = MapSplit$DWV$BombusViralLoad~1, locations = MapSplit$BQCV, newdata = grd) 

idw.output = as.data.frame(idw)  # output is defined as a data table
names(idw.output)[1:3] <- c("long", "lat", "var1.pred")  # give names to the modelled variables
ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = log(var1.pred))) +   geom_point(data = SurvMap, aes(x = long, y = lat), shape = 21,  colour = "red") + scale_fill_gradient(low = "cyan", high = "orange")
