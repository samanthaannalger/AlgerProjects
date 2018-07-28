# Shiny App for Beeks!
# July 26, 2018
# Samantha Alger

# Clear memory of characters:
rm(list=ls())

library(rgdal)
library(rgeos)
library(plyr)
library(data.table)
library(leaflet)
library(rgdal)
library(shiny)
library(shinythemes)
library(markdown)
library(knitr)
library(DT)
library(expss)
#library(lemon)
#library(kableExtra)

#set working director
setwd("~/AlgerProjects/VTApiaries/shinyappApiaries")


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
# DATA PREP:
##############################################################

# Data cleaning to make final df:
histDat = data.table(RegData)
histDat[, `n` := .N, by = BeekeeperID]
#histDat <-  histDat[!duplicated(histDat$BeekeeperID),]

histDat$Beektype <- ifelse(histDat$n == 1,"Hobbyist", ifelse(histDat$n <=5, "Sideliner", "Commercial"))

################################################
#Merging the two dataframes for shiny app:

#select only columns we need:
histDat <- dplyr::select(histDat, LocationID, Beektype, AccountName, BeeKeeperStatus, n)

FullApiaryDat <- dplyr::select(FullApiaryDat, -AccountName, -BeeKeeperStatus)

Shinydf <- merge.data.frame(FullApiaryDat,histDat, by = "LocationID", all.y = TRUE)

#Create Summary data for Chloropleth mapping

CountySummary <- ddply(Shinydf,c("CountyName"), summarise,
                       Loss = round(100* mean(PerTotLoss, na.rm=TRUE), digits = 2),
                       apiaries = length(n),
                       colonies = sum(ColonyCount),
                       beeks = length(unique(unlist(BeekeeperID[!is.na(BeekeeperID)]))),
                       n = length(n),
                       nLoss = length(PerTotLoss [!is.na(PerTotLoss)]),
                       MiteTrue = length(MiteCounts[MiteCounts==TRUE]),
                       MiteFalse = length(MiteCounts[MiteCounts==FALSE]),
                       PerMiteTrue = round(100*(MiteTrue/length(MiteCounts [!is.na(MiteCounts)])),digits = 2))

#rename countyname column
names(CountySummary)[1] <- "NAME"

# read in geojson data (county basemaps)
vtcounties <- rgdal::readOGR(dsn="cb_2017_us_county_5m.geojson")

#subset to only include vermont
vermont <- vtcounties[vtcounties@data$STATEFP == 50, ]

#merge the summary data with the base map data
vermont <- merge(vermont, CountySummary)

# Preparing a df for the resultant summary tables:

#rename the 'data' slot of the vermont spatial dataframe. Use the '@' symbol to refer to the slot name 'data', and drop unused levels to reduce the df to the VT counties only:
vermontdf <- droplevels(vermont@data)


#select only columns we will need for the summary tables
vermontdf<- dplyr::select(vermontdf, NAME, Loss, apiaries, colonies, nLoss, n, PerMiteTrue, beeks)


#############################################################################

####################################################################
# function name: sumtable
# description: creates summary table for subsetted
# parameters: 
# data = subsetted dataframe to be used
# returns: a formatted table
####################################################################

sumtable<- function(data=data, x) {
  
  data2 = dplyr::select(data, x, col.names = c(y))
  p <- kable(data2)
  
  return(p)
  
}



#vermontdf = dplyr::select(vermontdf, NAME, apiaries)
#sumtab <- kable(vermontdf, col.names = c("County", "# Apiaries"))

#####################################################################
# END OF FUNCTION
####################################################################

# Set default color for maps:
pal <- colorNumeric("Blues", NULL)

ui <- fluidPage(
  theme = shinytheme("cerulean"),
        navbarPage("BeekApp",
           tabPanel("Home",
                h4("Home of Vermont's Registered Apiary Data")),
           navbarMenu("Maps", 
                    tabPanel("Apiary Density",
                            mainPanel(
                              h3("Apiary Density"),
                              br(),
                                tabsetPanel(type = "tabs",
                                  tabPanel("Map", leafletOutput("mymap")),
                                  tabPanel("Table",
                                           DT::dataTableOutput("sum")),
                                  tabPanel("Description",
                                           h5("Apiary Density"),
                                           p("Map and table display the total number of registered apiaries in Vermont by county"))))),
                    tabPanel("Beekeeper Density",
                             mainPanel(
                               h3("Beekeeper Density"),
                               br(),
                               tabsetPanel(type = "tabs",
                                  tabPanel("Map", leafletOutput("mymap2")),
                                  tabPanel("Table",
                                            DT::dataTableOutput("sum2")),
                                  tabPanel("Description",
                                           h5("Beekeeper Density"),
                                           p("Map and table display the total number of registered beekeepers in Vermont by county"))))),
                    tabPanel("Colony Loss",
                             mainPanel(
                              h3("% Annual Colony Loss (2017)"),
                              br(),
                                tabsetPanel(type = "tabs",
                                  tabPanel("Map", leafletOutput("mymap3")),
                                  tabPanel("Table",
                                           DT::dataTableOutput("sum3")),
                                  tabPanel("Description",
                                           h5("% Annual Colony Loss (2017)"),
                                           p("Map and table display county averages for % colony loss during the 2016-2017 season"))))),
                    tabPanel("Mite Monitoring",
                            mainPanel(
                              h3("% Mite Monitoring"),
                              br(),
                                tabsetPanel(type = "tabs",
                                  tabPanel("Map", leafletOutput("mymap4")),
                                  tabPanel("Table",
                                           DT::dataTableOutput("sum4")),
                                  tabPanel("Description",
                                           h5("% Mite Monitoring"),
                                           p("Map and table display county level data for the % of apiaries managed by beekeepers who reported using mite monitoring practices")))))),

                        # creating an ouput for the table
          navbarMenu("Analyses",
                      tabPanel("Colony Loss"),
                      tabPanel("Management"),
                      tabPanel("Basic Statistics")
          )
))

server <- function(input,output, session){
    output$mymap <- renderLeaflet({
    m <- leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 0.90,
              fillColor = ~pal(apiaries),
              label = ~paste0(NAME, ": ", formatC(apiaries, big.mark = ","))) %>%
              addLegend(pal = pal, values = ~apiaries, opacity = 1, title = "# Apiaries")
    m })
    
apiaryTab = dplyr::select(vermontdf, NAME, apiaries)
  output$sum <- DT::renderDataTable ({
                  DT::datatable(apiaryTab, 
                              rownames = FALSE, 
                              colnames = c("County", "# Apiaries"), 
                              options = list(pageLength=17, dom = 'ft'))
  })
  
  output$mymap2 <- renderLeaflet({
    m <- leaflet(vermont) %>%
      addTiles() %>%
      addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 0.9,
                  fillColor = ~pal(beeks),
                  label = ~paste0(NAME, ": ", formatC(beeks, big.mark=","))) %>%
                  addLegend(pal = pal, values = ~beeks, opacity = 1, title = "# Beekeepers")
    m })
  
  beekTab = dplyr::select(vermontdf, NAME, beeks)
                output$sum2 <- DT::renderDataTable ({
                  DT::datatable(beekTab, 
                                rownames = FALSE, 
                                colnames = c("County", "# Beekeepers"), 
                                options = list(pageLength=17, dom = 'ft'))
  })


output$mymap3 <- renderLeaflet({
  pal <- colorNumeric("YlOrRd", NULL)
  m <- leaflet(vermont) %>%
    addTiles() %>%
    addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 0.9,
                fillColor = ~pal(Loss),
                label = ~paste0(NAME, ": ", formatC(Loss, big.mark = ","),"%, n=",nLoss)) %>%
                addLegend(pal = pal, values = ~Loss, opacity = 1, title = "% Annual Loss")
  m })

lossTab = dplyr::select(vermontdf, NAME, Loss, nLoss)
          output$sum3 <- DT::renderDataTable ({
                  DT::datatable(lossTab, 
                                      rownames = FALSE, 
                                      colnames = c("County", "% Annual Colony Loss", "N"), 
                                      options = list(pageLength=17, dom = 'ft'))
  
})
output$mymap4 <- renderLeaflet({
  m <- leaflet(vermont) %>%
  addTiles() %>%
  addPolygons(stroke = FALSE, smoothFactor = 0.3, fillOpacity = 1,
              fillColor = ~pal(PerMiteTrue),
              label = ~paste0(NAME, ": ", formatC(PerMiteTrue, big.mark = ","),"%, n=",n)) %>%
              addLegend(pal = pal, values = ~PerMiteTrue, opacity = 1.0, title = "% mite monitoring")
  m })

MiteTab = dplyr::select(vermontdf, NAME, PerMiteTrue, n)
                  output$sum4 <- DT::renderDataTable ({
                        DT::datatable(MiteTab, 
                                      rownames = FALSE, 
                                      colnames = c("County", "% Mite Monitoring", "N"), 
                                      options = list(pageLength=17, dom = 'ft'))
  
})
}



shinyApp(ui, server)






#Code for table using kable:

# ui:
#tableOutput("sum")),

#server:
#vermontdf = dplyr::select(vermontdf, NAME, apiaries)
#    output$sum <- function() {
#        vermontdf %>%
#        mutate(County = rownames(.)) %>%
#        select(NAME, apiaries) %>%
#        knitr::kable("html") %>%
#        kable_styling(latex_options = c("striped", "hold_position"))
#    }
