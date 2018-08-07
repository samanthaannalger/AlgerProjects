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
library(shiny)
library(shinythemes)
library(markdown)
library(knitr)
library(DT)
library(expss)
library(ggplot2)
library(plotly)
library(tidyr)
library(mgcv)
library(Matrix)

#set working director
setwd("~/AlgerProjects/VTApiaries/shinyappApiaries/BeekApp/")


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
histDat <- dplyr::select(histDat, LocationID, Beektype, BeeKeeperStatus, n)

FullApiaryDat <- dplyr::select(FullApiaryDat, -BeeKeeperStatus)

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

# Preparing df for map summary tables:

#rename the 'data' slot of the vermont spatial dataframe. Use the '@' symbol to refer to the slot name 'data', and drop unused levels to reduce the df to the VT counties only:
vermontdf <- droplevels(vermont@data)


#select only columns we will need for the summary tables
vermontdf<- dplyr::select(vermontdf, NAME, Loss, apiaries, colonies, nLoss, n, PerMiteTrue, beeks)

####################################################################
# Data prep for Analyses Tab
###################################################################
# Apiary and Colony number by beekeeper Type:

BeekTypeStats <- ddply(Shinydf, c("Beektype"), summarise, 
                       apiaries = length(n),
                       colonies = sum(ColonyCount, na.rm = TRUE),
                       loss = round(100*mean(PerTotLoss, na.rm = TRUE),digits=2),
                       sd = sd(100*PerTotLoss, na.rm=TRUE),
                       se = sd / sqrt(apiaries))
BeekTypeStats$perApiaries <- round(100*BeekTypeStats$apiaries/sum(BeekTypeStats$apiaries), digits=2)

BeekTypeStats$perColonies <- round(100*BeekTypeStats$colonies/sum(BeekTypeStats$colonies), digits=2)

# Code for stacked bar plot
BeekTypeDF <- rbind(BeekTypeStats, BeekTypeStats)
MeasureType <- c(rep("Apiary", 3), rep("Colony", 3))
BeekTypeDF <- cbind(BeekTypeDF, MeasureType)
BeekTypeDF$Percent <- c(BeekTypeDF$perApiaries[1:3], BeekTypeDF$perColonies[4:6])

#Reorder factors for stacked bar plot:
BeekTypeDF$Beektype <- factor(BeekTypeDF$Beektype, levels = c("Hobbyist","Sideliner", "Commercial"))

#rename factors for plotting:
levels(BeekTypeDF$Beektype) <- c("Hobbyist (1 apiary)", "Sideliner (2-5 apiaries)", "Commercial (5+ apiaries)")

# Reorder factors for colony loss bar plot:
BeekTypeStats$Beektype <- factor(BeekTypeStats$Beektype, levels = c("Hobbyist","Sideliner", "Commercial"))

#rename factors for plotting:
levels(BeekTypeStats$Beektype) <- c("Hobbyist (1 apiary)", "Sideliner (2-5 apiaries)", "Commercial (5+ apiaries)")

# Create df for Pie Chart:
Singles <- Shinydf[!duplicated(Shinydf$BeekeeperID), ]

mitedf <- data.frame(group = c("Did not count mites", "Counted mites"), value = c(table(Singles$MiteCounts)))

####################################################
#Figure for Mite count methods:

#Import data:
MiteMon<- Singles

#Create  separate dfs for analysis with specified columns:
MonMethods <- dplyr::select(MiteMon, BeekeeperID,SugarShakeYN, AlcoholWashYN, BottomBoardYN, DroneSurveyYN, OtherMiteCountYN)

#Switch df format to 'long' format:
MonMethods <- gather(MonMethods, question, response, SugarShakeYN:OtherMiteCountYN, factor_key=TRUE)

#Convert all True and False to '0' and '1' in the 'ReasonLoss' dataframe
MonMethods$response<-as.integer(as.logical(MonMethods$response))

# Preparing Data for Bar Plot
MiteMon <- ddply(MonMethods, c("question"), summarise, 
                 n = length(response),
                 mean =round(100*mean(response, na.rm=TRUE), digits=2),
                 sd = round(100*sd(response, na.rm=TRUE), digits=2),
                 se = sd / sqrt(n))

#########################################################
# Reasons for Losses

#Import data:
LossDat<- Singles

#Create two separate dfs for analysis with specified columns:
ReasonLoss <- dplyr::select(LossDat, BeekeeperID, ColonyLossVarroaMiteYN, ColonyLossStarvationYN, ColonyLossBearsYN, ColonyLossAmericanFoulbroodYN, ColonyLossSwarmingYN, ColonyLossPesticidesYN, ColonyLossMitacidesYN, OtherColonyLossYN)

#Convert to long format
ReasonLoss <- gather(ReasonLoss, question, response, ColonyLossVarroaMiteYN:OtherColonyLossYN, factor_key=TRUE)

#Convert all True and False to '0' and '1' in the 'ReasonLoss' dataframe
ReasonLoss$response<-as.integer(as.logical(ReasonLoss$response))

# Preparing Data for Bar Plot
LossCause <- ddply(ReasonLoss, c("question"), summarise, 
                   n = length(response),
                   mean = round(100*mean(response, na.rm=TRUE),digits=2),
                   sd = round(100*sd(response, na.rm=TRUE),digits=2),
                   se = sd / sqrt(n))


####################################################################
# END Data prep for Analyses Tab
###################################################################
####################################################################
# function name: sumtable
# description: creates summary table for subsetted df of apiaries
# parameters: 
# data = subsetted dataframe to be used
# returns: a formatted table
####################################################################

sumtable<- function(data=data) {
  
  data$PerTotLoss2 <- data$PerTotLoss*100
  
  data = apply_labels(data,
                      ColonyCount = "Sum",
                      PerTotLoss2 = "Colony Loss",
                      Beektype = "Beekeeper Type")
  
  
  sumtable<- data %>%
    tab_cells(ColonyCount) %>%
    tab_cols(Beektype, total()) %>%
    tab_stat_sum("Colonies") %>%
    tab_stat_valid_n("Apiaries") %>%
    tab_cells(PerTotLoss2) %>%
    tab_stat_mean("% Ave") %>%
    tab_pivot()
  
  return(sumtable)
  
}

####################################################################
# BEGIN SHINY CODE
###################################################################
# Set default color for maps:
pal <- colorNumeric("Blues",NULL)

ui <- fluidPage(
  theme = shinytheme("cerulean"),
        navbarPage("VT BeekApp",
           tabPanel("About",
                h4("Welcome to BeekApp"),
                   h6("Home of Vermont's Registered Apiary Data"),
                p("In 2017, Vermont's Apiary Inspection Program began collecting data on beekeeping practices and colony loss. Now, you can explore and visualize results through BeekApp!"),
                p("BeekApp is open to the public and intended for use by beekeepers, beekeeping clubs, and researchers. It allows users to explore Vermont state trends and identify opportunities for education and new research."), 
                  p("All data were collected through Vermont Apiary Registrations and Beekeeper Censuses. State wide trends are shown. Beekeeper personal information and apiary locations are kept confidential.")),
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
          navbarMenu("Data",
                      tabPanel("Registrations",
                            mainPanel(
                              br(),
                              h3("Registrations by beekeeper type"),
                              br(),
                              plotlyOutput('BeekPlot', height = "350px"),
                              DT::dataTableOutput("BeekTable"))),
                      tabPanel("Colony Loss",
                               mainPanel(
                              h3("Colony losses"),
                              br(),
                              tabsetPanel(type = "tabs",
                                          tabPanel("Loss Summary",
                                              br(),
                                              p("Annual colony loss by beekeepers type (2017)."),
                                              plotlyOutput("BeekLoss", height = "350px"),
                                              DT::dataTableOutput("BeekTable2")),
                                          tabPanel("Losses Explained",
                                              br(),
                                              p("Explanations provided by beekeepers for their colony losses."),
                                              plotlyOutput("LossExp", height = "500px"),
                                              plotlyOutput("Other")))
                              )),
                      tabPanel("Pest Management",
                               mainPanel(
                                 h3("Pest Management"),
                                 br(),
                                 tabsetPanel(type = "tabs",
                                      tabPanel("Mite Monitoring",
                                               br(),
                                               p("Figure shows the percentage of beekeepers who reported monitoring Varroa mites in their colonies."),
                                             plotlyOutput("pie", height = "350px")),
                                      tabPanel("Monitoring Methods",
                                               br(),
                                            p("The methods beekeepers report using to monitor Varroa mites in their colonies."),
                                            plotlyOutput("MiteMethods", height = "500px")),
                                      tabPanel("Treatments ",
                                               br(),
                                               br(),
                                               p("under construction"),
                                               plotlyOutput("TreatPlot", height = "500px"))))),
                      tabPanel("Challenges",
                              mainPanel(
                                h3("Challenges"),
                                p("We asked VT beekeepers to tell us about the biggest challenges they face as beekeepers and here are the results"),
                                br(),
                                br(),
                                p("under construction")
          )))))

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
                  
# Plot for Beekeeper type and apiary/colonies:
# Use ggplot sntax and then use ggplotly to make the figure "interactive"
output$BeekPlot <- renderPlotly({
                    BeekTypeFig <- ggplot(BeekTypeDF, 
                                          aes(y = Percent, 
                                              x = MeasureType , 
                                              fill = Beektype)) + 
                      geom_bar(stat="identity") + theme_classic() + 
                      labs(x=NULL, y = "% Registered in VT") +
                      scale_fill_brewer() +
                      scale_x_discrete(labels=c("Apiary" = "Apiaries", "Colony" = "Colonies")) + 
                      guides(fill=guide_legend(title="Beekeeper Type"))
  
#Use ggplotly to make the figure:
  ggplotly(BeekTypeFig)  %>% 
    layout(height = input$plotHeight, autosize=TRUE) %>% # set the size, specified in the ui
    config(displayModeBar = F) %>% # Removes the hover bar
    layout(xaxis=list(fixedrange=TRUE)) %>%  # disables the zoom option
    layout(yaxis=list(fixedrange=TRUE)) #disables the zoom option

})

# Plot for Beekeeper type and apiary/colonies:
# Use ggplot sntax and then use ggplotly to make the figure "interactive"
output$BeekLoss <- renderPlotly({
                    BeekPlot <- ggplot(BeekTypeStats, 
                                  aes(y = loss, x = Beektype, fill = Beektype)) + 
                                  geom_bar(stat="identity") + 
                                  theme_classic() + 
                                  labs(x=NULL, y = "% of total in Vermont") +
                                  scale_fill_brewer() +
                                  scale_x_discrete(labels=c("Hobbyist", "Sideliner", "Commercial")) +
                                  geom_errorbar(aes(ymin = loss - se, ymax = loss + se, width = 0.2)) +
                                  guides(fill=guide_legend(title="Beekeeper Type"))
  
  #Use ggplotly to make the figure:
  ggplotly(BeekPlot)  %>% 
    layout(height = input$plotHeight, autosize=TRUE) %>% # set the size, specified in the ui
    config(displayModeBar = F) %>% # Removes the hover bar
    layout(xaxis=list(fixedrange=TRUE)) %>%  # disables the zoom option
    layout(yaxis=list(fixedrange=TRUE)) #disables the zoom option
  
})

#create summary tables below figures
output$BeekTable <- DT::renderDataTable({
  as.datatable_widget(sumtable(Shinydf))
})

output$BeekTable2 <- DT::renderDataTable({
  as.datatable_widget(sumtable(Shinydf))
})

output$pie <- renderPlotly({
  pie<-plot_ly(mitedf, 
               labels = ~group, 
               values = ~value, 
               type = 'pie',
               marker = list(colors=c("dodgerblue", "lightskyblue"))) %>%
    layout( xaxis = list(showgrid = FALSE, 
                         zeroline = FALSE, showticklabels = FALSE), 
            yaxis = list(showgrid = FALSE, zeroline = FALSE, 
                         showticklabels = FALSE))
  ggplotly(pie) %>% 
    layout(height = input$plotHeight, autosize=TRUE) %>% # set the size, specified in the ui
    config(displayModeBar = F) %>% # Removes the hover bar
    layout(xaxis=list(fixedrange=TRUE)) %>%  # disables the zoom option
    layout(yaxis=list(fixedrange=TRUE)) #disables the zoom option
  
  })

output$MiteMethods <- renderPlotly ({
          MitePlot <- ggplot(MiteMon,
                      aes(x=question, y=mean, fill=question)) + 
                      geom_bar(stat="identity", position=position_dodge()) + 
                      theme_classic() +
                      labs(x="Mite Monitoring Method", y = "% Reported Use") + 
                      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), 
                      legend.position="none", axis.text=element_text(size=15), 
                      axis.title=element_text(size=18,face="bold")) + 
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2)) +                       scale_fill_brewer() + 
                      scale_x_discrete(labels=c("SugarShakeYN" = "Sugar Shake", "AlcoholWashYN" = "Alcohol Wash", "BottomBoardYN" = "Bottom Board", "DroneSurveyYN"= "Drone Survey", "OtherMiteCountYN"= "Other"))
  
  
  ggplotly(MitePlot)  %>% 
    layout(height = input$plotHeight, autosize=TRUE) %>% # set the size, specified in the ui
    config(displayModeBar = F) %>% # Removes the hover bar
    layout(xaxis=list(fixedrange=TRUE)) %>%  # disables the zoom option
    layout(yaxis=list(fixedrange=TRUE)) #disables the zoom option
  
  
})

output$LossExp <- renderPlotly ({
                    LossExpPlot <- ggplot(LossCause, aes(x=question, y=mean, fill=question)) + 
                    geom_bar(stat="identity", color="black",
                    position=position_dodge()) + 
                    labs(x="Causes", y = "% Reported Causes") + 
                    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), 
                          legend.position="none", axis.text=element_text(size=9), 
                          axis.title=element_text(size=18,face="bold")) + 
                    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2)) + 
                       scale_fill_brewer() +
                      scale_x_discrete(labels=c("ColonyLossVarroaMiteYN" = "Varroa", "ColonyLossStarvationYN"= "Starvation", "ColonyLossBearsYN" = "Bears", "ColonyLossAmericanFoulbroodYN" = "AFB", "ColonyLossSwarmingYN" = "Swarming", "ColonyLossPesticidesYN" = "Pesticides", "ColonyLossMitacidesYN" = "Mitacides", "OtherColonyLossYN"= "Other"))

ggplotly(LossExpPlot)  %>% 
        layout(height = input$plotHeight, autosize=TRUE) %>% # set the size, specified in the ui
        config(displayModeBar = F) %>% # Removes the hover bar
        layout(xaxis=list(fixedrange=TRUE)) %>%  # disables the zoom option
        layout(yaxis=list(fixedrange=TRUE)) #disables the zoom option
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


