#Making Maps

# Preliminaries:

# Clear memory of characters:
rm(list=ls())


library(rgdal)
library(rgeos)
library(plyr)
library(data.table)
library(leaflet)
library(rgdal)

#set working director
setwd("~/AlgerProjects/VTApiaries/CSV_files/")


#upload data
RegData <- read.csv("RegActiveAndDelinquent.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

FullApiaryDat <- read.csv("VTApiaries.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

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

write.csv(Shinydf, file ="Shinydf.csv")
table(Shinydf$CountyName)


# Mapping points: https://andrewbtran.github.io/NICAR/2017/maps/mapping-census-data.html

# Use for making a density map!

# Chloropeth mapping:

# https://rstudio.github.io/leaflet/json.html
# https://stackoverflow.com/questions/43446802/how-to-download-ny-state-all-county-data-in-r-for-leaflet-map

#Create Summary data for Chloropleth mapping

CountySummary <- ddply(Shinydf,c("CountyName"), summarise,
                     Loss = 100* mean(PerTotLoss, na.rm=TRUE),
                     apiaries = length(n),
                     colonies = sum(ColonyCount),
                     nLoss = length(PerTotLoss [!is.na(PerTotLoss)]),
                     MiteTrue = length(MiteCounts[MiteCounts==TRUE]),
                     MiteFalse = length(MiteCounts[MiteCounts==FALSE]),
                     PerMiteTrue = 100*(MiteTrue/length(MiteCounts [!is.na(MiteCounts)])))

#rename countyname column
names(CountySummary)[1] <- "NAME"

# read in geojson data (county basemaps)
vtcounties <- rgdal::readOGR(dsn="cb_2017_us_county_5m.geojson")

#subset to only include vermont
vermont <- vtcounties[vtcounties$STATEFP == 50, ]

#merge the summary data with the base map data
vermont <- merge(vermont, CountySummary)

#set color pallet
pal <- colorNumeric("YlOrRd", NULL)

#Map for # of apiaries by county
leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(apiaries),
              label = ~paste0(NAME, ": ", formatC(apiaries, big.mark = ","))) %>%
   addLegend(pal = pal, values = ~apiaries, opacity = 1.0, title = "# Apiarires")

#Map for # colonies by county
leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(colonies),
              label = ~paste0(NAME, ": ", formatC(colonies, big.mark = ","))) %>%
  addLegend(pal = pal, values = ~colonies, opacity = 1.0, title = "# Colonies")


#Map for % colony loss by county
pal <- colorNumeric("YlOrRd", NULL)

leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(Loss),
              label = ~paste0(NAME, ": ", formatC(Loss, big.mark = ","),"%, n=",nLoss)) %>%
  addLegend(pal = pal, values = ~Loss, opacity = 1.0, title = "% Annual Loss")

# Map for mite monitoring by county

pal <- colorNumeric("Blues", NULL)
leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(PerMiteTrue),
              label = ~paste0(NAME, ": ", formatC(PerMiteTrue, big.mark = ","), "%, n=",apiaries)) %>%
  addLegend(pal = pal, values = ~PerMiteTrue, opacity = 1.0, title = "% Monitor Mites")




  