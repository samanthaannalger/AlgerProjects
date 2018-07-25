##SHINY CODE:

# S. Alger
# July 16, 2018
# Shiny App Demno for ApiaryMap

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

setwd("~/AlgerProjects/VTApiaries/shinyappApiaries")
RegData <- read.csv("RegActiveAndDelinquent.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

FullApiaryDat <- read.csv("VTApiaries.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

library(shiny)
library(leaflet)
library(data.table)
library(geosphere)
library(ggmap)
library(markdown)
library(DT)
library(expss)

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

# Select only columns we need:
Shinydf<- dplyr::select(Shinydf, Latitude, Longtitude, LocationID, BeekeeperID, PaPlantsID, AccountName, City, BeePurchaseVendorName, MiteCounts,ColonyCount, PerTotLoss, BeeKeeperStatus, n, Beektype)
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

####################################################################
# function name: sumtable
# description: creates summary table for subsetted df of apiaries
# parameters: 
# data = subsetted dataframe to be used
# returns: a formatted table
####################################################################

sumtable<- function(data=data) {
  
  dat$PerTotLoss2 <- dat$PerTotLoss*100
  
  dat = apply_labels(dat,
                     ColonyCount = "Summary",
                     PerTotLoss2 = "Annual Colony Loss",
                     Beektype = "Beekeeper Type")
  
  
  sumtable<- dat %>%
    tab_cells(ColonyCount) %>%
    tab_cols(Beektype, total()) %>%
    tab_stat_sum("Colonies") %>%
    tab_stat_valid_n("Apiaries") %>%
    tab_cells(PerTotLoss2) %>%
    tab_stat_mean("% Average") %>%
    tab_pivot()
  
  return(sumtable)
  
}

#####################################################################
# END OF FUNCTION
####################################################################


# Pre-Data Cleaning:

LLmat <- LatLongMat(data = Shinydf)
#SSdat <- SubSetMap(data = Shinydf, rad = 20, lat = -72.746286, long = 44.278876, matrix = LLmat)


#####################################################################
# BEGIN SHINY CODE
####################################################################


# User interface, creating tabs
ui <- navbarPage("Apiary Locator",
           tabPanel("Map",
                    sidebarLayout(
                      sidebarPanel(
                        sliderInput("distance", # Slider bar for distance
                                    "Radius (Miles):",
                                    min = 1,
                                    max = 40,
                                    value = 30)
                      ),
                      
                      mainPanel(
                        leafletOutput("mymap", # creating the main map output
                                      height="300px"),
                                absolutePanel(top=0, 
                                      left=55, 
                                      textInput("target_zone", 
                                                "" ,
                                                "Ex: Burlington, Vermont")),
                #Create summary table
                        h5("Summary Table"),
                        DT::dataTableOutput("sum")
                        )
                      )
            ),
           # creating an ouput for the table
           tabPanel("Table",
                    DT::dataTableOutput("table")
           ),
           navbarMenu("More",
                      tabPanel("Summary",
                               DT::dataTableOutput("summary")
                      ),
                      tabPanel("About",
                               fluidRow(
                                 column(6),
                                 column(3)
         )
      )
   )
)


# Server code for Shiny App:
#Starts with the code for mymap
server <- function(input, output, session) {
                 output$mymap<- renderLeaflet({
                  rad <- input$distance # Uses the slider bar input for the distance radius
                   # Uses latitude and longitude entered into the bar
                              if(input$target_zone=="Ex: Burlington"){
                              ZOOM=2
                              lat=0
                              long=0
                              }else{
                              target_pos=geocode(input$target_zone)
                              lat=target_pos$lat
                              long=target_pos$lon
                              ZOOM=12
                        }
    # Program body that uses the functions defined at the beginning of the strip
    # Create a subsetted DF that uses the input radious, long, and lat
    SSdat <<- SubSetMap(data = Shinydf, rad = rad, lat = long ,long =, lat, matrix = LLmat)
    # create an output table of this dataframe (Full table on tab)
    output$table <- DT::renderDataTable({
      DT::datatable(SSdat[[1]])
    })
  #create summary table on main panel under map
    x <- sumtable(SSdat[[1]])
    output$sum <- DT::renderDataTable({
    as.datatable_widget(sumtable(x))
    })
    # Map function that maps the new df subsetted above
    Mapfunc(data=SSdat[[1]], rad= SSdat[[2]], lat = SSdat[[4]], long= SSdat[[3]]) 
    
    
    
  })
  
  
  output$summary <- renderPrint({
    summary(cars)
  })
  
}

shinyApp(ui, server)
