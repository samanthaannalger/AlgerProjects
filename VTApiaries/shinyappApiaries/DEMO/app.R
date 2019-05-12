library(shiny)

ui <- fluidPage(
  
  # Application title
  titlePanel("mtcars"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("mpg", "mpg Limit",
                  min = 11, max = 33, value = 20)
    ),
    
    mainPanel(
      tableOutput("mtcars_kable")
    )
  )
)

server <- function(input, output) {
  library(dplyr)
  library(kableExtra)
  output$mtcars_kable <- function() {
    mtcars %>%
     # mutate(car = rownames(.)) %>%
     # select(car, everything()) %>%
      knitr::kable("html") %>%
     kable_styling(latex_options = c("striped", "hold_position"))
     # add_header_above(c(" ", "Group 1" = 5, "Group 2" = 6))
  }
}

# Run the application
shinyApp(ui = ui, server = server)