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
        title = h5("Upload Data"),
        value = "Upload Data",
        fileInput(
          inputId = "file",
          label = "Select a File:", 
          multiple = FALSE
        )
      ),
      accordion_panel(
        title = h5("Highlight adaptive basin(s)"),
        value = "Highlight adaptive basin(s)",
        uiOutput("peak_choice"),
        actionButton("submit_highlight", "Plot")
      ),
      accordion_panel(
        title = h5("Display shortest path(s)"),
        value = "Display shortest path(s)",
        uiOutput("source_choice"),
        uiOutput("target_choice"),
        actionButton("submit_path", "Display Paths"),
        uiOutput("path_choice")
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
                       label = "Select sink(s) to visualize:",
                       choices = options,
                       selected = options)
  })
  output$source_choice <- renderUI({
    options_valleys <- label_convert(labels()[!(labels() %in% list_extrema(processed_data(), labels())$peaks)])
    options_valleys <- setNames(names(options_valleys), options_valleys)
    selectInput(
      inputId = "source_choice", 
      label = "Select source state:",
      choices = options_valleys,
      multiple = FALSE
    )
  })
  output$target_choice <- renderUI({
    req(labels())
    options_peaks <- label_convert(list_extrema(processed_data(), labels())$peaks)
    options_peaks <- setNames(names(options_peaks), options_peaks)
    selectInput(
      inputId = "target_choice", 
      label = "Select target sink state:", 
      choices = options_peaks, 
      multiple = FALSE
    )
  })
  selected_source <- reactiveVal(NULL)
  observeEvent(input$submit_path, {
    selected_source(input$source_choice)
  })
  selected_target <- reactiveVal(NULL)
  observeEvent(input$submit_path, {
    selected_target(input$target_choice)
  })
  output$path_choice <- renderUI({
    req(input$submit_path)
    path_options <- list_shortest_paths(processed_data(), selected_source(), selected_target())
    radioButtons(
      inputId = "path",
      label = "Select path:",
      choices = path_options
    )
  })
 
  selected_peaks <- reactiveVal(NULL) # set to null by default
  observeEvent(input$submit_highlight, { # changes only when event is observed (button press)
    if (is.null(selected_path()) == FALSE){
      selected_path(NULL)
    }
    selected_peaks(input$peak_choice)
  })
  selected_path <- reactiveVal(NULL)
  observeEvent(input$path, {
    if (is.null(selected_peaks()) == FALSE){
      selected_peaks(NULL)
    }
    if (input$path == "No possible paths"){
      selected_path(NULL)
    } else {
      selected_path(input$path)
    }
  })
    
  output$circosPlot <- renderPlot({
    req(processed_data(), labels())
    dataset <- processed_data()
    custom_order <- sort_labels(labels())
    plot_title <- NULL
    highlighting <- NULL
    highlight_mode <- "all"
    if (is.null(selected_peaks()) == FALSE){
      highlighting <- list()
      for (peak in selected_peaks()){
        highlight_mode <- "basin"
        highlighting <- c(highlighting, list_basin(dataset, peak))
      }
      cleaned_title <- str_replace(paste(label_convert(selected_peaks()), collapse = " & "), "\n", " ")
      plot_title <- str_glue("Adaptive basin(s) for:\n{cleaned_title}")
    } else if (is.null(selected_path()) == FALSE){
      highlight_mode <- "path"
      highlighting = strsplit(selected_path(), ", ")[[1]] # here we split the string back into a list of nodes for plotting
      plot_title <- str_glue("Shortest path between:\n{label_convert(selected_source())} and {label_convert(selected_target())}")
    }
    rownames(dataset) <- colnames(dataset) # for now, we will only accept a matrix with identically ordered labels
    eco_transition_plot(dataset = dataset, highlighting = highlighting, sector_order = custom_order, plot_labels = labels(), highlight_mode = highlight_mode)
    title(plot_title)
  })
  
  observe({
    #print(selected_path())
  })
  
}

shinyApp(ui = ui, server = server)
