library(tidyverse)
library(circlize)
library(igraph)

# Function to determine the bit width (number of loci) for a dataset from its labels
# Labels must be in binary format e,g. "(00,00)", and we assume homogenous genome size.
# Input: a column of labels from a dataset (adjacency matrix)
set_bitwidth <- function(input_column){
  return(nchar(extract_numbers(input_column[1])[1]))
}

# Parses an input fitness table with column format:
# ("species1 = "00", species2 = "00",..., speciesN, fitness = <dbl>)
# Column names are unimportant, except for the final fitness column.
# Concatenates all provided (potentially variable numbers of) species into one ID label e.g. "(sp1,sp2,...)"
# Returns a cleaned dataframe with columns (fitness = <dbl>, ID = <chr>)
parse_fitness_table <- function(input_data){
  data <- read.csv(input_data)
  n <- ncol(data)
  alleles <- list("0", "1") # If the provided data are not already, we will map them to [0,1]
  data <- data %>% 
    mutate("ID" = do.call(paste, c(data[, 1:(n-1)], sep = ","))) %>% # Paste together all but the last column of the dataframe
    select(`fitness`, `ID`)
  data$ID <- sapply(data$ID, function(x) gsub(",NA", "", x)) # Remove any NAs introduced (they are just handled as strings?)
  names(alleles) <- sort(unique(strsplit(gsub(",", "", paste0(data$ID, collapse = "")), "")[[1]])) # ID the 2 alleles present and assign to [0,1]
  alleles <- c(alleles, "," = ",") # Add the comma back in after mapping (unreliable sorting otherwise)
  data$ID <- sapply(data$ID, function(x) { # Update ID column with binary alphabet
    chars <- strsplit(x, "")[[1]]
    paste0("(", paste0(alleles[chars], collapse = ""), ")") # Coerce IDs to a single string in proper format: "(g1,g2,...)"
  })
  return(data)
}

fitnesses_to_adjacency <- function(input_data){
  # we first determine which transitions are possible (have hamming distance of 1)
    # we do this by comparing the min sum of hamming distances between the two comparison possibilities e.g. a1b1+a2b2 vs a1b2 + a2b1 
  n <- length(input_data$ID)
  transition_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(input_data$ID, input_data$ID))
  for (i in 1:(n-1)){ # to make pairwise comparisons between all IDs
    for (j in 2:n){
      state1 <- strsplit(input_data$ID[i], ",")[[1]]
      state2 <- strsplit(input_data$ID[j], ",")[[1]]
      # This section only applies if we consider comparing multiple genotypes from the different species 
      # (as opposed to assuming all interaction partners are the same species, in which case we'd compare all permutations)
      hamdist <- 0
      for (l in 1:length(state1)){
        hamdist <- hamdist + hamming_dist(state1[l], state2[l])
      }
      if (hamdist == 1 && input_data$fitness[i] < input_data$fitness[j]){
        transition_matrix[i,j] <- as.integer(1)
      } else if (hamdist == 1 && input_data$fitness[i] > input_data$fitness[j]){
        transition_matrix[j,i] <- as.integer(1)
      }
    }
  }
  return(transition_matrix)
}

# Converts a single numeric string to binary given a bit width e.g. "7" & width = 3 -> "111" 
to_binary <- function(number, width) { # should later handle whatever the largest number is
  paste0(rev(as.integer(intToBits(as.integer(number)))[1:width]), collapse = "")
}

# Extract numbers from a label in format "(00,00)"; remove parentheses, split by comma.
extract_numbers <- function(label){
  nums <- gsub("[()]", "", label)
  strsplit(nums, ",")[[1]]
}

# Convert and format labels from base 10 to base 2. 
# Input: named list of format: list(`(x)` = "(x)", `(x,y)` = "(x,y)", ...)
# Output: named list of identical format with numbers in base 2. Commas in IDs are replaced with newlines (\n)
label_convert <- function(input_labels){
  width = set_bitwidth(input_labels)
  binary_labels <- sapply(input_labels, function(name){
    nums <- extract_numbers(name)
    bin <- sapply(nums, function(n) to_binary(n, width = width))
    paste(bin, collapse = "\n")
  })
  return(binary_labels)
}

# takes two binary strings (e.g "000" & "010") and returns the hamming distance between them
hamming_dist <- function(g1, g2){
  return(sum(strsplit(g1, "")[[1]] != strsplit(g2, "")[[1]]))
}

# Generate an adjacency matrix for a hamming space with dimension n
hamming_matrix <- function(n){
  nodes <- as.matrix(expand.grid(rep(list(c(0,1)), n)))
  num_nodes <- nrow(nodes)
  adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  for (i in 1:(num_nodes - 1)){
    for (j in (i + 1):num_nodes) {
      if (sum(nodes[i, ] != nodes[j, ]) == 1){
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1  # symmetric since undirected
      }
    }
  }
  rownames(adj_matrix) <- apply(nodes, 1, paste0, collapse = "")
  colnames(adj_matrix) <- rownames(adj_matrix)
  return(adj_matrix)
}

# Function to determine the min number of mutations required between two states. This may not actually be achievable on a given landscape. 
# For now, we assume the source is monotypic and the targets can be mono or polytypic
# Input: a single binary string (e.g. "000") and newline separated string (e.g. "010\n111"))
# Output: a number representing the minimum number of mutations from source to the targets
# i.e. the sum of all distances on a min span tree between the source and targets.
# this cannot handle more than two targets
min_distance <- function(source, targets){
  # This is a form of the k-min span tree problem, which is NP-hard. Probably can be generalized to 3+ species with a recursive function?
  bitwidth <- nchar(source)
  full_network_matrix <- hamming_matrix(bitwidth) # make a general (full landscape) matrix
  full_network_graph <- graph_from_adjacency_matrix(full_network_matrix) # convert matrix to graph
  targets <- str_split_1(targets, "\n")
  names(targets) <- targets
  # find the shortest path between the source and target: if target is monotypic, take hamming dist
  if (length(unique(targets)) == 1){
    return(hamming_dist(source, unique(targets)))
  } else {
    # if target has two or more, then find the shortest combined path by:
    # find the min hamming dist between all genotypes in state
    min_dist_target <- "0" # Dummy value
    hamdist <- Inf
    for (target in targets){
      if (hamming_dist(source, target) < hamdist){
        hamdist = hamming_dist(source, target)
        min_dist_target = target
      }
    }
    targets <- targets[targets != min_dist_target]
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

# Function to sort labels by increasing distance from 0.
sort_labels <- function(input_labels){
  binary_labels <- label_convert(input_labels)
  zero <- paste0(rep("0", times = set_bitwidth(input_labels)), collapse = "") # create a 0 state with correct bit width
  distances <- sapply(binary_labels, function(x) min_distance(zero, x)) # Figuring out distance of everything from 0
  names(distances) <- input_labels
  # Determine which indices belong to which interaction orders
  interaction_order <- c()
  sector = 1
  sector_begin <- 1
  i <- 1
  for (label in binary_labels){
    if (nchar(label) > nchar(binary_labels[[sector_begin]])){ # Track where the transition between interaction orders is
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

# Function to detect the number of sectors in a binary label vector, and create a vector of colors for the chord diagram
generate_sector_colors <- function(dataset, input_labels, highlighting = NULL){
  # works with up to a four-species community
  color_options <- c("#f0554a", "#ba2014", "#f0554a","#e0110d11",  
                     "#cfc197", "#a19574", "#cfc197","#e0810d11",
                     "#42993c", "#267021", "#42993c","#42993c11",
                     "#8d52a8", "#5d2e73", "#8d52a8","#8d52a811")
  #sector_colors <- sapply(binary_labels, function(x) )
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

# Function for generating color palettes of links to use with the chord diagram
generate_link_colors <- function(dataset, highlight_mode = "all", highlighting = NULL){
  color_matrix <- matrix(nrow = nrow(dataset), ncol = ncol(dataset), "#00000000")
  dimnames(color_matrix) <- dimnames(dataset)
  color_options <- c("#f0554abb", "#ba2014bb", "#f0554abb","#e0110d11",  
                     "#cfc197bb", "#a19574bb", "#cfc197bb","#e0810d11",
                     "#42993cbb", "#267021bb", "#42993cbb","#42993c11",
                     "#8d52a8bb", "#5d2e73bb", "#8d52a8bb","#8d52a811")
  if (highlight_mode == "all"){
    for (node in colnames(dataset)){ # Iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1){ # If link exists
          order <- length(strsplit(link, ",")[[1]]) # Find which interaction order (of the source) based on commas in label
          if (sum(dataset[link,]) == 0){ # If is a peak (no links out over the row, only has links in over the column)
            color_matrix[link, node] <- color_options[[(order * 4) - 2]] # Grabs peak color for the respective order
          } else if (sum(dataset[link,]) != 0){
            color_matrix[link, node] <- color_options[[(order * 4) - 3]] # Else it must not be a peak
          }
        }
      }  
    }
  }
  if (highlight_mode == "basin" && is.null(highlighting) == FALSE){ # Double check for highlighting anyway
    for (node in colnames(dataset)){ # Iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1 & node %in% highlighting){ # If link exists and target node is in list of those to visualize:
          node_order <- length(strsplit(link, ",")[[1]]) # Find which interaction order based on commas in label
          if (sum(dataset[link,]) == 0 && node %in% highlighting){ # If is a peak (no links out over the row, only has links in over the column)
            color_matrix[link, node] <- color_options[[(node_order * 4) - 2]] # Grabs peak color for the respective order
          } else if (sum(dataset[link, ]) != 0 && node %in% highlighting){
            color_matrix[link, node] <- color_options[[(node_order * 4) - 3]] # Else it must not be a peak
          }
        }
      }  
    }
  } else if (highlight_mode == "path" && is.null(highlighting) == FALSE){
      for (i in 1:(length(highlighting)-1)){ # Since this is a directed path, we simply highlight the edges between nodes in the path
        node_order <- length(strsplit(highlighting[[i]], ",")[[1]])
        color_matrix[highlighting[[i]], highlighting[[i+1]]] <- color_options[[(node_order * 4) - 3]] # None of them will be peaks, since only the last value in the path is   
      }
  }
  return(color_matrix)
}

# Function for generating color palettes of arrows to use with the chord diagram. Opacity should correspond to which links are also highlighted
generate_arrow_colors <- function(dataset, highlight_mode = "all", highlighting){
  color_matrix <- matrix(nrow = nrow(dataset), ncol = ncol(dataset), "#00000000")
  dimnames(color_matrix) <- dimnames(dataset)
  if (highlight_mode == "all"){
    for (node in colnames(dataset)){ # Iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1){ # If link exists
          color_matrix[link,node] <- "#000000BB"
        }
      }  
    }
  }
  if (highlight_mode == "basin" && is.null(highlighting) == FALSE){ # Double check for highlighting anyway
    for (node in colnames(dataset)){ # Iterate over matrix
      for (link in rownames(dataset)){
        if (dataset[link,node] == 1 & node %in% highlighting){ # If link exists and target node is in list of those to visualize:
          color_matrix[link,node] <- "#000000BB"
        }
      }  
    }
  } else if (highlight_mode == "path" && is.null(highlighting) == FALSE){
      for (i in 1:(length(highlighting)-1)){ # Since this is a directed path, we simply highlight the edges between nodes in the path
        node_order <- length(strsplit(highlighting[[i]], ",")[[1]])
        color_matrix[highlighting[[i]], highlighting[[i+1]]] <- "#000000BB" 
      }
  }
  return(color_matrix)
}

list_extrema <- function(dataset, input_labels){
  peaks <- c()
  valleys <- c()
  for (i in 1:ncol(dataset)){
    if (sum(dataset[i,]) == 0){ # If there are no outgoing connections
      peaks <- c(peaks, input_labels[i])
    } else if (sum(dataset[,i]) == 0){ # If there are no incoming connections
      valleys <- c(valleys, input_labels[i])
    }
  }
  names(peaks) = peaks
  names(valleys) = valleys
  return(list("peaks" = peaks, "valleys" = valleys))
}

list_basin <- function(dataset, root){
  basin <- bfs(graph_from_adjacency_matrix(dataset), root, mode = "in", unreachable = FALSE, dist = TRUE)
  return(c(root, names(basin$dist[basin$dist > 0]))) # Return the root and all nodes part of the basin (nodes not part of basin return -1)
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
  paths <- lapply(paths, toString) # We will have to coerce to a string here bc shiny cannot handle lists as a selectable value. We split these later
  return(paths)
}

interaction_transitions <- function(dataset){
# function to return the number of transitions between interaction orders in a dataset
  binary_labels <- label_convert(colnames(dataset))
  total_transitions <- sum(dataset)
  num_transitions_plus <- 0
  num_transitions_minus <- 0
  for (i in 1:ncol(dataset)){
    for (j in 1:nrow(dataset)){
      if (dataset[i,j] == 1 && length(str_split_1(binary_labels[i], "\n")) > length(str_split_1(binary_labels[j], "\n"))){
        num_transitions_plus <- num_transitions_plus + 1
      } else if (dataset[i,j] == 1 && length(str_split_1(binary_labels[i], "\n")) < length(str_split_1(binary_labels[j], "\n"))){
        num_transitions_minus <- num_transitions_minus + 1
      }
    }
  }
  return(list(`total` = total_transitions,
              `increasing` = num_transitions_plus,
              `decreasing` = num_transitions_minus))
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
    sector_order <- plot_labels # Defaults to input label order
  } else if(is.null(sector_order) && sorting == TRUE){ # If no custom order but sorting desired
    sector_order <- sort_labels(input_labels)
  } # else will provide custom sector order
  if (is.null(plot_labels)){ # if custom labels are not provided, copy from sector_order
    plot_labels <- sector_order
  } else {
    plot_labels <- label_convert(plot_labels)
  }
  
  # create gaps between different interaction orders
  #circos.par$gap.degree <- c(rep(1, times = 7), 15, rep(1, times = 6), 15)
  
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
