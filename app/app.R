library(tidyverse)
library(circlize)
library(igraph)
library(shiny)
library(bslib)
source("utils.R")

ui <- page_sidebar(
  title = "VevoTEC: Visualizing Evolutionary Transformations of Ecological Communities",
  
  ##############################################################################
  # Sidebar panel for uploading data and highlighting functions
  ##############################################################################
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
  
  ##############################################################################
  # Output panels: circos plot, summary tables, raw data
  ##############################################################################
  navset_card_underline(
    nav_panel(
      "Plot",
      plotOutput(outputId = "circosPlot"),
      card_footer("plot description here")
    ),
    nav_panel(
      "Summary",
      tableOutput("peaks"),
      textOutput("transitions"),
      card_footer("Network summary descriptions for: \nAdaptive basins & their sizes,
                  number of transitions between states with differing numbers of interaction partners")
    ),
    nav_panel(
      "Raw Data",
      tableOutput("rawData"),
      card_footer("Rendering of imported data as an adjacency matrix. If you encounter issues with plot 
                  generation or have unexpected results, check here to make sure the data are being parsed correctly.")
    )
  )
)

server <- function(input, output){
  ##############################################################################
  # Code pertaining to the import and initialization of data + variables
  ##############################################################################
  processed_data <- reactive({ # Importing the uploaded data (fitness table or adjacency matrix)
    req(input$file) # Require a file so the whole thing doesn't break on init
    if (tools::file_ext(input$file$datapath) == "csv"){
      parse_fitness_table(input$file$datapath) %>% 
        fitnesses_to_adjacency()
    } else { # else assume it is already an adjacency matrix
      read.table(input$file$datapath, sep = " ", header = T, check.names = F) %>% 
        as.matrix()
    }
  })
  labels <- reactive({ # Labels are propagated as a named vector in base 10
    req(processed_data())
    input_labels <- colnames(processed_data()) # Creating the default label vector for later use
    names(input_labels) <- input_labels
    input_labels
  })
  output$rawData <- renderTable({ # For the displaying of the raw adjacency matrix to the user
    req(processed_data())
    processed_data()
  })
  
  ##############################################################################
  # Code pertaining to stats displayed in the Summary pane
  ##############################################################################
  output$transitions <- renderText({ 
    req(processed_data())
    str_glue("Number of inter-order transitions: {interaction_transitions(processed_data())}")
  })
  output$peaks <- renderTable({ # Table listing basins and sizes
    peaks <- list_extrema(processed_data(), labels())$peaks
    basins <- lapply(peaks, function(x) list_basin(processed_data(), x)) %>% sapply(function(x) toString(label_convert(x[-1])))
    size <- sapply(basins, function(x) length(str_split(x, ",")[[1]]))
    data.frame("Peak" = label_convert(peaks), "Basin" = basins, "Size" = size)
  })
  
  ##############################################################################
  # Code pertaining to the Adaptive Basins widget
  ##############################################################################
  output$peak_choice <- renderUI({ # Reactive UI element to calculate and display possible peaks to select
    req(labels())    
    options <- label_convert(list_extrema(processed_data(), labels())$peaks)
    options <- setNames(names(options), options)
    checkboxGroupInput(inputId = "peak_choice", 
                       label = "Select sink(s) to visualize:",
                       choices = options,
                       selected = options)
  })
  selected_peaks <- reactiveVal(NULL) # Reactive value recording the selected peaks in the UI
  observeEvent(input$submit_highlight, { # Changes only when event is observed (button press)
    if (is.null(selected_path()) == FALSE){
      selected_path(NULL)
    }
    selected_peaks(input$peak_choice)
  })
  
  ##############################################################################
  # Code pertaining to the Shortest Paths widget
  ##############################################################################
  output$source_choice <- renderUI({ # reactive UI element to display a list of all source states (non-peaks) 
    options_valleys <- label_convert(labels()[!(labels() %in% list_extrema(processed_data(), labels())$peaks)])
    options_valleys <- setNames(names(options_valleys), options_valleys)
    selectInput(
      inputId = "source_choice", 
      label = "Select source state:",
      choices = options_valleys,
      multiple = FALSE
    )
  })
  output$target_choice <- renderUI({ # Reactive UI element to display all possible sinks (peaks)
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
  selected_source <- reactiveVal(NULL) # Reactive value recording the selected source state in the UI
  observeEvent(input$submit_path, {
    selected_source(input$source_choice)
  })
  selected_target <- reactiveVal(NULL) # Reactive value recording the selected target state in the UI
  observeEvent(input$submit_path, {
    selected_target(input$target_choice)
  })
  output$path_choice <- renderUI({ # Reactive UI element displaying choices for shortest path to display, if there are multiple
    req(input$submit_path)
    path_options <- list_shortest_paths(processed_data(), selected_source(), selected_target())
    radioButtons(
      inputId = "path",
      label = "Select path:",
      choices = path_options
    )
  })
  selected_path <- reactiveVal(NULL) # Reactive value recording the selected shortest path
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
 
  ##############################################################################
  # Rendering of the main Circos plot, based on which widgets are used
  ##############################################################################
  output$circosPlot <- renderPlot({
    req(processed_data(), labels())
    dataset <- processed_data()
    custom_order <- sort_labels(labels())
    plot_title <- NULL
    highlighting <- NULL
    highlight_mode <- "all" # By default, we visualize everything in the plot
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
      highlighting = strsplit(selected_path(), ", ")[[1]] # Here we split the string back into a list of nodes for plotting
      plot_title <- str_glue("Shortest path between:\n{label_convert(selected_source())} and {label_convert(selected_target())}")
    }
    rownames(dataset) <- colnames(dataset) # For now, we will only accept a matrix with identically ordered labels
    eco_transition_plot(dataset = dataset, highlighting = highlighting, sector_order = custom_order, plot_labels = labels(), highlight_mode = highlight_mode)
    title(plot_title)
  })
}

shinyApp(ui = ui, server = server)
