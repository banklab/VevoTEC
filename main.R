library(tidyverse)
library(shiny)
#library(plotly) #we actually don't need this since chord diagrams aren't compatible with ggplot
library(circlize)

# ok now trying a dataset with higher order interactions (two resources)
z <- read.table("example_data/transitions_matrix.dat", sep = " ", header = T, check.names = F) %>% 
  as.matrix()
rownames(z) <- colnames(z) # for now, we will only accept a matrix with identically ordered labels
input_labels <- colnames(z) # creating the default label vector for later use
names(input_labels) <- input_labels

to_binary <- function(x, width = ceiling(log2(7))) { # should later handle whatever the largest number is
  paste0(rev(as.integer(intToBits(as.integer(x)))[1:width]), collapse = "")
}

# Extract numbers from each column name
extract_numbers <- function(label) {
  # Remove parentheses
  nums <- gsub("[()]", "", label)
  # Split by comma if needed
  strsplit(nums, ",")[[1]]
}

# Convert and format
binary_labels <- sapply(input_labels, function(name) {
  nums <- extract_numbers(name)
  bin <- sapply(nums, function(n) to_binary(n))
  paste(bin, collapse = "\n")
})

custom_order <- c("(0)", "(1)", "(2)", "(4)",
                  "(3)", "(5)", "(6)", "(7)",
                  "(5,7)","(2,7)","(2,3)","(1,5)",
                  "(1,2)","(0,7)","(0,3)")
names(custom_order) = custom_order

eco_transition_plot <- function(dataset, 
                                sorting = TRUE, 
                                sector_order = NULL, 
                                plot_labels = NULL, ...){
  
  if(is.null(sector_order) && sorting == FALSE){ # if no custom order is defined
    sector_order <- input_labels # Defaults to input label order
  } else if(is.null(sector_order) && sorting == TRUE){ # If no custom order but sorting desired
    sector_order <- sort_labels(input_labels)
  } # else will provide custom sector order
  if (is.null(plot_labels)){ # if custom labels are not provided, copy from sector_order
    plot_labels <- sector_order
  } 
  # Yes, I have ensured they are in the correct order. ignore warning
  sector_colors <- c(rep("#FF0000", times = 3),"#B20000",rep("#FF0000", times = 4), rep("#c8a772", times = 2), "#B08133", rep("#c8a772", times = 4))
  
  # create gaps between different interaction orders
  circos.par$gap.degree <- c(rep(1, times = 7), 15, rep(1, times = 6), 15)
  
  #calculate the starting offset such that 0 is at the bottom
  zero_num_connections <- sum(dataset[1,]) + sum(dataset[,1])
  total_connections <- 2 * sum(z)
  degrees_per_connection <- (360 - sum(circos.par$gap.degree[circos.par$gap.degree > 1])) / 
    total_connections
  circos.par$start.degree <- -90 + (zero_num_connections * degrees_per_connection / 2)
  
  values = runif(length(sector_order))
  chordDiagram(dataset, order = sector_order,
               annotationTrack = "grid",
               directional = 1, 
               grid.col = sector_colors,
               direction.type = "arrows",
               preAllocateTracks = list(
                 list(track.height = 0.05),
                 list(track.height = 0.025),
                 list(track.height = 0.05)
               ))
#  circos.track(
#    track.index = 3, 
#    panel.fun = function(x, y) {
#      xlim <- CELL_META$xlim
#      ylim <- CELL_META$ylim
#      sector.name <- CELL_META$sector.index
#      xcenter <- CELL_META$xcenter
#      circos.rect(
#        xleft = xlim[1], ybottom = ylim[1],
#        xright = xlim[2], ytop = ylim[2],
#        col = rand_color(1), # need to figure out an opacity parameter
#        border = NA)
#      }, 
#    bg.border = NA
#  )
  circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      xcenter <- CELL_META$xcenter
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], plot_labels[sector.name],
                  facing = "downward", niceFacing = TRUE)
      },
    bg.border = NA
  )
  circos.clear()
}
eco_transition_plot(z, sector_order = custom_order, plot_labels = binary_labels)
