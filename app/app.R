library(tidyverse)
library(circlize)
library(igraph)
library(shiny)
library(bslib)
source("utils.R")

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
      tableOutput("labels"),
      card_footer("summary stats description here")
    ),
    nav_panel(
      "Raw Data",
      tableOutput("rawData")
    )
  )
)

server <- function(input, output) {
  
  processed_data <- reactive({
    dataset <- input$file
    read.table(dataset$datapath, sep = " ", header = T, check.names = F) %>% 
      as.matrix()
  })
  labels <- reactive({
    input_labels <- colnames(processed_data()) # creating the default label vector for later use
    names(input_labels) <- input_labels
  })
  output$rawData <- renderTable({
    processed_data()
  })
  output$labels <- renderTable({
    labels()
  })
  output$circosPlot <- renderPlot({
    dataset <- processed_data()
    custom_order <- sort_labels(labels())
    rownames(dataset) <- colnames(dataset) # for now, we will only accept a matrix with identically ordered labels
    
    eco_transition_plot(dataset, sector_order = custom_order, plot_labels = labels())

  })
  
}

shinyApp(ui = ui, server = server)
