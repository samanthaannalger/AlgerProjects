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
library(rgdal)
library(rgeos)
library(leaflet)
# Set working directory:
setwd("~/AlgerProjects/VTApiaries/CSV_files/")


# read in geojson data (county basemaps)
division <- rgdal::readOGR(dsn="AK_divisions_NAD83.geojson")

(division@data)

pal <- colorNumeric("viridis", NULL)

leaflet(division) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal)
            #  label = ~paste0(county, ": ", formatC(pop, big.mark = ","))) %>%
 # addLegend(pal = pal, values = ~log10(pop), opacity = 1.0,
  #          labFormat = labelFormat(transform = function(x) round(10^x)))