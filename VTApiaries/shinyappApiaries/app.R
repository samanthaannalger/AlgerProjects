# Samantha Alger
# July 14,2018
# Shiny App for Apiary Data DEMO
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())
#setwd("~/AlgerProjects/VTApiaries/shinyappApiaries/")

library(shiny)
library(data.table)
library(rsconnect)

# Read in Data
FullApiaryDat <- read.csv("VTApiaries.csv", 
                          header=TRUE, 
                          sep = ",", 
                          stringsAsFactors = FALSE)

RegData <- read.csv("RegActiveAndDelinquent.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

histDat = data.table(RegData)
histDat[, `n` := .N, by = BeekeeperID]
#histDat <-  histDat[!duplicated(histDat$BeekeeperID),]
# Histogram showing the number of apiaries belonging to beekeepers with 1 yard, two yards...etc...
hist(histDat$n, freq=TRUE, breaks=25)

histDat$Beektype <- ifelse(histDat$n == 1,"Hobbyist", ifelse(histDat$n <=5, "Sideliner", "Commercial"))

####### 
#Merging the two dataframes for shiny app:

#select only columns we need:
histDat <- dplyr::select(histDat,LocationID, Beektype, AccountName, n)

Shinydf <- merge.data.frame(FullApiaryDat,histDat, by = "LocationID", all.y = TRUE)

# Define UI for colonyloss app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Colony Loss"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for variable to plot against loss ----
      selectInput("variable", "Variable:",
                  c("Beekeeper Type" = "Beektype",
                    "Mite Counts" = "MiteCounts",
                    "Supplemental Feed" = "SupplementalFeed"
                    )),
      
      # Input: Checkbox for whether outliers should be included ----
      checkboxInput("outliers", "Show outliers", TRUE)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Formatted text for caption ----
      h3(textOutput("caption")),
      
      # Output: Plot of the requested variable against mpg ----
      plotOutput("LossPlot")
      
    )
  )
)



# Data Preprocessing:
LossDat <- Shinydf

#rename TRUE and FALSE as YEs and NO, ** need to fix so it's not 'hard coded' for future analyses.
LossDat$MiteCounts <- factor(LossDat$MiteCounts, labels = c("No", "Yes"))
LossDat$SupplementalFeed <- factor(LossDat$SupplementalFeed, labels = c("No", "Yes"))

# Remove all rows without vendor listed:
LossDat <- LossDat[! is.na(LossDat$Beektype), ]
LossDat <- LossDat[! is.na(LossDat$MiteCounts), ]



# Define server logic to plot various variables against colony loss ----
server <- function(input, output) {
  
  # Compute the formula text ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$LossPlot functions
  formulaText <- reactive({
    paste("PerTotLoss ~", input$variable)
  })
  
  # Return the formula text for printing as a caption ----
  output$caption <- renderText({
    formulaText()
  })
  
  # Generate a plot of the requested variable against loss ----
  # and only exclude outliers if requested
  output$LossPlot <- renderPlot({
    boxplot(as.formula(formulaText()),
            data = LossDat,
            #outline = input$outliers,
            col = "#75AADB", pch = 19)
  })
  
}

shinyApp(ui, server)

