

SSdat[[1]] <- dplyr::select(SSdat[[1]], Latitude, Longtitude, LocationID, BeekeeperID, PaPlantsID, BeePurchaseVendorName, MiteCounts,ColonyCount, PerTotLoss, BeeKeeperStatus, n, Beektype, City)
colnames(SSdat[[1]])

library(markdown)

navbarPage("Navbar!",
           tabPanel("Plot",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("plotType", "Plot type",
                                     c("Scatter"="p", "Line"="l")
                        )
                      ),
                      mainPanel(
                        plotOutput("plot")
                      )
                    )
           ),
           tabPanel("Summary",
                    verbatimTextOutput("summary")
           )))


## Only run this example in interactive R sessions
if (interactive()) {
  # table example
  shinyApp(
    ui = fluidPage(
      fluidRow(
        column(12,
               tableOutput('table')
        )
      )
    ),
    server = function(input, output) {
      output$table <- renderTable(SSdat[[1]])
    }
  )
  
  
  # DataTables example
  shinyApp(
    ui = fluidPage(
      fluidRow(
        column(12,
               dataTableOutput('table')
        )
      )
    ),
    server = function(input, output) {
      output$table <- renderDataTable(SSdat[[1]])
    }
  )
}

