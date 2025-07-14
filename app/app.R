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
        uiOutput("peak_choice")
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
      )
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
      #tableOutput("peak_choice"),
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
    req(input$file) # require a file so the whole thing doesn't break on init
    read.table(input$file$datapath, sep = " ", header = T, check.names = F) %>% 
      as.matrix()
  })
  labels <- reactive({
    req(processed_data())
    input_labels <- colnames(processed_data()) # creating the default label vector for later use
    names(input_labels) <- input_labels
    input_labels
  })
  output$rawData <- renderTable({
    req(processed_data())
    processed_data()
  })
  output$labels <- renderTable({
    labels()
  })
  output$peak_choice <- renderUI({
    req(labels())    
    options <- label_convert(list_extrema(processed_data(), labels())$peaks)
    options <- setNames(names(options), options)
    checkboxGroupInput(inputId = "peak_choice", 
                       label = "Select states to visualize:",
                       choices = options,
                       selected = options)
  })
  output$circosPlot <- renderPlot({
    req(processed_data(), labels())
    dataset <- processed_data()
    custom_order <- sort_labels(labels())
    rownames(dataset) <- colnames(dataset) # for now, we will only accept a matrix with identically ordered labels
    eco_transition_plot(dataset, sector_order = custom_order, plot_labels = labels())
  })
  
  observe({
    print(label_convert(list_extrema(processed_data(), labels())$peaks))
  })
  
}

shinyApp(ui = ui, server = server)
