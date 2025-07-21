library(tidyverse)
library(circlize)
library(igraph)

# functions and code to convert to create a vector of binary labels
to_binary <- function(x, width = ceiling(log2(7))) { # should later handle whatever the largest number is
  paste0(rev(as.integer(intToBits(as.integer(x)))[1:width]), collapse = "")
}

extract_numbers <- function(label) {
  # extract numbers from each column name
  # remove parentheses
  nums <- gsub("[()]", "", label)
  # split by comma
  strsplit(nums, ",")[[1]]
}

label_convert <- function(input_labels){
  # Convert and format
  binary_labels <- sapply(input_labels, function(name) {
    nums <- extract_numbers(name)
    bin <- sapply(nums, function(n) to_binary(n))
    paste(bin, collapse = "\n")
  })
  return(binary_labels)
}

sort_labels <- function(input_labels){
  # function to sort labels by increasing distance from 0. Otherwise the provided order is used.
  binary_labels <- label_convert(input_labels)
  distances <- sapply(binary_labels, function(x) min_distance("000", x))
  names(distances) <- input_labels
  
  # determine which indices belong to which interaction orders
  interaction_order <- c()
  sector = 1
  sector_begin <- 1
  i <- 1
  for (label in binary_labels){
    if (nchar(label) > nchar(binary_labels[[sector_begin]])){ # track where the transition between interaction orders is
      sector <- sector + 1
      sector_begin <- i
    }
    interaction_order <- c(interaction_order, sector)
    i <- i + 1
  }
  
  # split 
  sorted_labels <- c(names(sort(distances[which(interaction_order == 1)], decreasing = F)), 
                     names(sort(distances[which(interaction_order == 2)], decreasing = T)))
  names(sorted_labels) <- sorted_labels
  return(sorted_labels)
}

generate_sector_colors <- function(dataset, input_labels, highlighting = NULL){
  # function to detect the number of sectors in a binary label vector, and assign a color for plotting
  color_options <- c("#f0554a", "#ba2014", "#f0554a","#e0110d11",  
                     "#cfc197", "#a19574", "#cfc197","#e0810d11",
                     "#42993c", "#267021", "#75c76f","#42993c11",
                     "#8d52a8", "#5d2e73", "#a577ba","#8d52a811")
  sector_colors <- c()
  sector_order <- 1
  sector_index <- 1
  binary_labels <- label_convert(input_labels)
  if (is.null(highlighting) == FALSE){
    binary_highlights <- label_convert(highlighting)
  }
  i <- 1
  for(label in binary_labels){
    if (nchar(label) > nchar(binary_labels[[sector_index]])){ # track where the transition between interaction orders is
      sector_index <- i
      sector_order <- sector_order + 1
    }
    if(is.null(highlighting) == FALSE && !(label %in% binary_highlights)){ # if we pass the highlighting, set everything else to opaque
      sector_colors <- c(sector_colors, color_options[[(sector_order * 4)]])
    } else if(nchar(label) == nchar(binary_labels[[sector_index]]) && sum(dataset[i,]) != 0 && sum(dataset[,i]) != 0){ # else if is not a peak
      sector_colors <- c(sector_colors, color_options[[(sector_order * 4) - 3]])
    } else if(nchar(label) == nchar(binary_labels[[sector_index]]) && sum(dataset[,i]) == 0){ # else if is a valley
     sector_colors <- c(sector_colors, color_options[[(sector_order * 4) - 1]])
    } else if(nchar(label) == nchar(binary_labels[[sector_index]]) && sum(dataset[i,]) == 0){ # if it is within the current interaction order and is a peak:
      sector_colors <- c(sector_colors, color_options[[(sector_order * 4) - 2]])
    }
    i <- i + 1
  }
  names(sector_colors) = input_labels
  return(sector_colors)
}

generate_link_colors <- function(dataset, highlight_mode = "all", highlighting = NULL){
  color_matrix <- matrix(nrow = nrow(dataset), ncol = ncol(dataset), "#00000000")
  dimnames(color_matrix) <- dimnames(dataset)
  color_options <- c("#f0554abb", "#ba2014bb", "#f0554abb","#e0110d11",  
                     "#cfc197bb", "#a19574bb", "#cfc197bb","#e0810d11",
                     "#42993cbb", "#267021bb", "#42993cbb","#42993c11",
                     "#8d52a8bb", "#5d2e73bb", "#8d52a8bb","#8d52a811")
  if (highlight_mode == "all"){
    for (node in colnames(dataset)){ #iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1){ # if link exists
          order <- length(strsplit(link, ",")[[1]]) # find which interaction order (of the source) based on commas in label
          if (sum(dataset[link,]) == 0){ # if is a peak (no links out over the row, only has links in over the column)
            color_matrix[link, node] <- color_options[[(order * 4) - 2]] # grabs peak color for the respective order
          } else if (sum(dataset[link,]) != 0){
            color_matrix[link, node] <- color_options[[(order * 4) - 3]] # else it must not be a peak
          }
          #color_matrix[link,node] <- "#000000BB"
        }
      }  
    }
  }
  if (highlight_mode == "basin" && is.null(highlighting) == FALSE){ # double check for highlighting anyway
    for (node in colnames(dataset)){ #iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1 & node %in% highlighting){ # if link exists and target node is in list of those to visualize:
          order <- length(strsplit(link, ",")[[1]]) # find which interaction order based on commas in label
          if (sum(dataset[link,]) == 0 && node %in% highlighting){ # if is a peak (no links out over the row, only has links in over the column)
            color_matrix[link, node] <- color_options[[(order * 4) - 2]] # grabs peak color for the respective order
          } else if (sum(dataset[link, ]) != 0 && node %in% highlighting){
            color_matrix[link, node] <- color_options[[(order * 4) - 3]] # else it must not be a peak
          }
        }
      }  
    }
  } else if (highlight_mode == "path"){
    print("yay")
  }
  return(color_matrix)
}

generate_arrow_colors <- function(dataset, highlight_mode = "all", highlighting){
  color_matrix <- matrix(nrow = nrow(dataset), ncol = ncol(dataset), "#00000000")
  dimnames(color_matrix) <- dimnames(dataset)
  if (highlight_mode == "all"){
    for (node in colnames(dataset)){ #iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1){ # if link exists
          color_matrix[link,node] <- "#000000BB"
        }
      }  
    }
  }
  if (highlight_mode == "basin" && is.null(highlighting) == FALSE){ # double check for highlighting anyway
    for (node in colnames(dataset)){ #iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1 & node %in% highlighting){ # if link exists and target node is in list of those to visualize:
          color_matrix[link,node] <- "#000000BB"
        }
      }  
    }
  } else if (highlight_mode == "path"){
    print("yay")
  }
  return(color_matrix)
}

list_extrema <- function(dataset, input_labels){
  peaks <- c()
  valleys <- c()
  for (i in 1:ncol(dataset)){
    if (sum(dataset[i,]) == 0){ # if there are no outgoing connections
      peaks <- c(peaks, input_labels[i])
    } else if (sum(dataset[,i]) == 0){ # if there are no incoming connections
      valleys <- c(valleys, input_labels[i])
    }
  }
  names(peaks) = peaks
  names(valleys) = valleys
  return(list("peaks" = peaks, "valleys" = valleys))
}

list_basin <- function(dataset, root){
  basin <- bfs(graph_from_adjacency_matrix(dataset), root, mode = "in", unreachable = FALSE, dist = TRUE)
  return(c(root, names(basin$dist[basin$dist > 0]))) # return the root and all nodes part of the basin (nodes not part of basin return -1)
}

list_shortest_paths <- function(dataset, source, target){
  data_matrix <- graph_from_adjacency_matrix(dataset)
  V(data_matrix)$name <- colnames(dataset)
  paths <- all_shortest_paths(data_matrix, from = source, to = target, mode = "out")$vpaths
  if (length(paths) == 0){
    return("No possible paths")
  }
  paths <- lapply(paths, function(x) as_ids(x))
  names(paths) <- lapply(paths, function(x) paste0(label_convert(x), collapse = "→"))
  paths <- lapply(paths, toString) # we will have to coerce to a string here bc shiny cannot handle lists as a selectable value. we split these later
  return(paths)
}

hamming_matrix <- function(n) {
  # Generate an adjacency matrix for a hamming space with dimension n
  nodes <- as.matrix(expand.grid(rep(list(c(0,1)), n)))
  num_nodes <- nrow(nodes)
  
  adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
      if (sum(nodes[i, ] != nodes[j, ]) == 1) {
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1  # symmetric since undirected
      }
    }
  }
  
  rownames(adj_matrix) <- apply(nodes, 1, paste0, collapse = "")
  colnames(adj_matrix) <- rownames(adj_matrix)
  
  return(adj_matrix)
}

hamming_dist <- function(g1, g2) {
  # takes two binary strings and returns the hamming distance between them
  return(sum(strsplit(g1, "")[[1]] != strsplit(g2, "")[[1]]))
}

eco_transition_plot <- function(dataset, 
                                sorting = TRUE, 
                                highlighting = NULL,
                                sector_order = NULL, 
                                plot_labels = NULL,
                                highlight_mode,
                                ...){
  if(is.null(highlighting) == FALSE){
    sector_colors <- generate_sector_colors(dataset, plot_labels, highlighting)
    link_colors <- generate_link_colors(dataset, highlight_mode, highlighting)
    arrow_colors <- generate_arrow_colors(dataset, highlight_mode, highlighting)
  } else if (is.null(highlighting) == TRUE){
    sector_colors <- generate_sector_colors(dataset, plot_labels)
    link_colors <- generate_link_colors(dataset)
    arrow_colors <- generate_arrow_colors(dataset)
  }
  if(is.null(sector_order) && sorting == FALSE){ # if no custom order is defined
    sector_order <- input_labels # Defaults to input label order
  } else if(is.null(sector_order) && sorting == TRUE){ # If no custom order but sorting desired
    sector_order <- sort_labels(input_labels)
  } # else will provide custom sector order
  if (is.null(plot_labels)){ # if custom labels are not provided, copy from sector_order
    plot_labels <- sector_order
  } else {
    plot_labels <- label_convert(plot_labels) # assuming the provided labels are in base 10, we convert to binary
  }
  
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
               col = link_colors,
               direction.type = "arrows",
               link.arr.col = arrow_colors,
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

min_distance <- function(source, targets){
  # Function to determine the min number of mutations required to produce a state. this may not actually be achievable on a given landscape 
  bitwidth <- nchar(source)
  full_network_matrix <- hamming_matrix(bitwidth) # make a general (full landscape) matrix
  full_network_graph <- graph_from_adjacency_matrix(full_network_matrix) # convert matrix to graph
  # find the shortest path between the source and target: if target is monotypic, take hamming dist
  targets <- str_split_1(targets, "\n")
  names(targets) <- targets
  if (length(targets) == 1){
    return(hamming_dist(source, targets))
  } else {
  # if target has two or more, then find the shortest combined path by:
  # find the min hamming dist between all genotypes in state
    min_dist_target <- "0"
    hamdist <- Inf
    for (target in targets){
      if (hamming_dist(source, target) < hamdist){
        hamdist = hamming_dist(source, target)
        min_dist_target = target
      }
    }
    targets <- targets[names(targets) != min_dist_target]
    # calc shortest path from source to this target
    second_source_indices <- as.vector(shortest_paths(full_network_graph, source, min_dist_target)[[1]][[1]])
    # nodes are stored in the vector V(g)$name, which is 1 indexed. You cannot reference by names
    second_sources <- V(full_network_graph)$name[second_source_indices]
    # find min hamming dist between each genotype in path and the second genotype
    second_min_dist_target <- "0"
    second_hamdist <- Inf
    for (second_source in second_sources){
      if (hamming_dist(second_source, targets) < second_hamdist){
        second_hamdist <- hamming_dist(second_source, targets)
        second_min_dist_source <- second_source
      }
    }
    # add hamming distances together
    hamdist <- hamdist + second_hamdist
  }
  return(hamdist)
}
