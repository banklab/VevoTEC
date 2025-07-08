library(tidyverse)
library(circlize)
library(igraph)
library(shiny)
library(bslib)
source("main.R")

ui <- page_sidebar( # create a sidebar
  # App title ----
  title = "VevoTEC: Visualizing Evolutionary Transformations of Ecological Communities",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    accordion(
      accordion_panel(
        "Upload Data",
        fileInput(
          inputId = "file",
          label = h3("Select a File:"), 
          multiple = FALSE
        )
      ),
      accordion_panel(
        "Highlight",
        checkboxGroupInput(
          inputId = "selection",
          label = "Select nodes to visualize",
          choices = list("000" = 0, "001" = 1)
        )
      ),
      accordion_panel(
        "Calculate Mutational Distance",
        selectInput(
          inputId = "source", 
          label = h3("Select Source Node"), 
          choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
          selected = 1,
          multiple = FALSE
        ),
        selectInput(
          inputId = "target", 
          label = h3("Select Target Node"),
          choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3),
          selected = 2
        )
      ),
      position = "right"
    )
  ),
  # Output: Histogram ----
  navset_card_underline(
    nav_panel(
      "Plot",
      plotOutput(outputId = "circosPlot"),
      card_footer("plot description here")
    ),
    nav_panel(
      "Summary",
      card_footer("summary stats description here")
    ),
    nav_panel(
      "Raw Data",
      tableOutput("rawData")
    )
  )
)

server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$rawData <- renderTable({
    file1 <- input$file
    read.table(file1$datapath, sep = " ", header = T, check.names = F) %>% 
    as.matrix()
  })
  output$circosPlot <- renderPlot({
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#424242", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
}

shinyApp(ui = ui, server = server)

