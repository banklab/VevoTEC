library(tidyverse)
library(circlize)
library(igraph)

#' Determine the bit width of a dataset
#'
#' Determines the bit width (number of loci) from a dataset’s labels.
#' Labels must be in binary format, e.g. "(00,00)", and genome size is assumed homogeneous.
#'
#' @param input_column A character vector of binary labels from the dataset (e.g. from an adjacency matrix).
#'
#' @return An integer representing the bit width (number of bits per locus).
#'
#' @examples
#' labels <- c("(00,00)", "(11,11)")
#' set_bitwidth(labels)
set_bitwidth <- function(input_column){
  return(nchar(extract_numbers(input_column[1])[1]))
}

#' Parse an input fitness table
#'
#' Reads and cleans a CSV file of fitness values and species genotype combinations.
#' Concatenates all genotype columns into a single binary ID column.
#' 
#' Input columns: ("species1 = "00", species2 = "00",..., speciesN, fitness = <dbl>)
#' Column names for species are unimportant, except for the final 'fitness' column.
#' Concatenates all provided (potentially variable numbers of) species into one ID label e.g. "(sp1,sp2,...)"
#' 
#' @param input_data File path to the input CSV with columns for species and a final column `fitness`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{fitness}{Numeric fitness values.}
#'   \item{ID}{Binary string label of concatenated genotypes, e.g. "(00,11)".}
#' }
#'
#' @examples
#' \dontrun{
#' df <- parse_fitness_table("data/fitness.csv")
#' head(df)
#' }
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

#' Create an adjacency matrix from a dataframe of fitnesses
#'
#' Builds a directed adjacency matrix of fitness-increasing transitions between ecological states.
#' Two states are connected if they differ by one locus (Hamming distance = 1),
#' with edges directed from lower to higher fitness.
#'
#' @param input_data A data frame with columns `ID` and `fitness`, e.g. from `parse_fitness_table()`.
#'
#' @return A square adjacency matrix (`matrix`) with row and column names corresponding to state IDs.
#' Entries are `1` for fitness-increasing transitions, `0` otherwise.
#'
#' @examples
#' data <- data.frame(
#'   ID = c("(00)", "(01)", "(11)"),
#'   fitness = c(0.1, 0.2, 0.3)
#' )
#' fitnesses_to_adjacency(data)
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

#' Convert decimal to binary string
#'
#' Converts an integer to a binary string with fixed bit width.
#'
#' @param number Numeric or character scalar representing the decimal number.
#' @param width Integer specifying the number of bits.
#'
#' @return A binary string of length `width`.
#'
#' @examples
#' to_binary(7, width = 3)
#' to_binary(2, width = 4)
to_binary <- function(number, width) { # should later handle whatever the largest number is
  paste0(rev(as.integer(intToBits(as.integer(number)))[1:width]), collapse = "")
}

#' Extract numeric components from a label
#'
#' Parses a label in the form "(00,01)" into its component genotypes.
#'
#' @param label Character string label in format "(g1,g2,...)"
#'
#' @return A character vector of genotype strings, e.g. c("00", "01").
#'
#' @examples
#' extract_numbers("(01,10)")
extract_numbers <- function(label){
  nums <- gsub("[()]", "", label)
  strsplit(nums, ",")[[1]]
}

#' Convert and format a set of labels from decimal to binary
#'
#' Converts all numeric values in a label list to binary form with the appropriate bit width.
#' Commas separating genotypes are replaced by newlines.
#'
#' @param input_labels Named list of character labels, e.g. list("(1)"="(1)", "(1,2)"="(1,2)").
#'
#' @return A named list with identical structure but containing binary representations.
#'
#' @examples
#' labels <- list("(1)"="(1)", "(1,2)"="(1,2)")
#' label_convert(labels)
label_convert <- function(input_labels){
  width = set_bitwidth(input_labels)
  binary_labels <- sapply(input_labels, function(name){
    nums <- extract_numbers(name)
    bin <- sapply(nums, function(n) to_binary(n, width = width))
    paste(bin, collapse = "\n")
  })
  return(binary_labels)
}

#' Calculate Hamming distance between two binary strings
#'
#' @param g1 Character string, binary (e.g. "010").
#' @param g2 Character string, binary (e.g. "000").
#'
#' @return Integer number of differing bits between `g1` and `g2`.
#'
#' @examples
#' hamming_dist("000", "010")
hamming_dist <- function(g1, g2){
  return(sum(strsplit(g1, "")[[1]] != strsplit(g2, "")[[1]]))
}

#' Generate adjacency matrix for n-dimensional Hamming space
#'
#' Creates a symmetric adjacency matrix where nodes represent binary strings,
#' and edges connect nodes differing by one bit.
#'
#' @param n Integer, number of bits (Hamming space dimension).
#'
#' @return A binary adjacency matrix with row/column names as binary strings.
#'
#' @examples
#' hamming_matrix(3)
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

#' Minimum number of mutations between two states
#'
#' Calculates the minimum number of mutations required to connect a source genotype
#' to one or more target genotypes, via the shortest path in Hamming space.
#'
#' @param source Character string, binary source genotype (e.g. "000").
#' @param targets Character string of one or more target genotypes separated by newlines.
#'
#' @return Integer representing the minimal number of mutational steps.
#'
#' @examples
#' min_distance("000", "010")       # 1
#' min_distance("000", "010\n111")  # >1
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

#' Sort binary labels by genetic distance
#'
#' Sorts a set of labels by increasing Hamming distance from the all-zero state.
#'
#' @param input_labels Character vector or named list of binary labels.
#'
#' @return A character vector of labels sorted by increasing genetic distance.
#'
#' @examples
#' labels <- c("(00)", "(01)", "(11)")
#' sort_labels(labels)
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

#' Generate colors for circos plot sectors
#'
#' Assigns consistent color palettes to sectors in a circos plot based on interaction order,
#' with optional highlighting.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param input_labels Character vector of labels for sectors.
#' @param highlighting Optional subset of labels to highlight.
#'
#' @return Named character vector of color hex codes for each sector.
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' generate_sector_colors(dataset, colnames(dataset))
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
    order <- length(strsplit(label, "\n")[[1]]) # Find which interaction order (of the source) based on commas in label
    if(is.null(highlighting) == FALSE && !(label %in% binary_highlights)){ # if we pass the highlighting, set everything else to opaque
      sector_colors <- c(sector_colors, color_options[[(order * 4)]])
    } else if(sum(dataset[i,]) != 0 && sum(dataset[,i]) != 0){ # else if is not a peak
      sector_colors <- c(sector_colors, color_options[[(order * 4) - 3]])
    } else if(sum(dataset[,i]) == 0){ # else if is a valley
     sector_colors <- c(sector_colors, color_options[[(order * 4) - 1]])
    } else if(sum(dataset[i,]) == 0){ # if it is a peak:
      sector_colors <- c(sector_colors, color_options[[(order * 4) - 2]])
    }
    i <- i + 1
  }
  names(sector_colors) = input_labels
  return(sector_colors)
}

#' Generate link colors for circos plot
#'
#' Creates a color matrix mapping transition edges to colors for a chord diagram.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param highlight_mode Character; one of `"all"`, `"basin"`, or `"path"`.
#' @param highlighting Optional character vector of node labels to emphasize.
#'
#' @return A character matrix of color hex codes corresponding to adjacency matrix positions.
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' generate_link_colors(dataset, "all")
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

#' Generate arrow colors for circos plot
#'
#' Produces an opacity-coded color matrix for directional arrows in chord diagrams.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param highlight_mode Character; one of `"all"`, `"basin"`, or `"path"`.
#' @param highlighting Vector of labels to emphasize.
#'
#' @return Character matrix of color codes (e.g. "#000000BB" for visible arrows).
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' generate_arrow_colors(dataset, "all")
#' generate_arrow_colors(dataset, "basin", highlighting = "(01)")
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

#' Identify peaks and valleys in fitness landscape
#'
#' Returns lists of local maxima (peaks) and minima (valleys) based on an adjacency matrix.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param input_labels Vector of corresponding node labels.
#'
#' @return A list with two named elements:
#' - `peaks`: labels with no outgoing edges.
#' - `valleys`: labels with no incoming edges.
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' list_extrema(dataset, colnames(dataset))
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

#' List all states within a basin of accessibility
#'
#' Traverses graph to find all nodes reachable from a given root.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param root Character string; label of the root node.
#'
#' @return Character vector of labels in the same basin as the root.
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' list_basin(dataset, "(11)")
list_basin <- function(dataset, root){
  basin <- bfs(graph_from_adjacency_matrix(dataset), root, mode = "in", unreachable = FALSE, dist = TRUE)
  return(c(root, names(basin$dist[basin$dist > 0]))) # Return the root and all nodes part of the basin (nodes not part of basin return -1)
}

#' Find shortest paths between two states
#'
#' Computes all shortest directed paths between two nodes in the transition network.
#'
#' @param dataset Adjacency matrix.
#' @param source Character string giving starting node label.
#' @param target Character string giving target node label.
#'
#' @return A list of shortest paths (as strings with arrows between labels), or
#' `"No possible paths"` if none exist.
#'
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' list_shortest_paths(dataset, "(00)", "(11)")
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

#' Calculate the number of transitions that change community diversity
#'
#' Different ecological states may have different numbers of species in them. 
#' This function sums the number of mutations in the system that cause gain or loss of community diversity.
#'
#' @param dataset Adjacency matrix
#' 
#' @returns A list with three named elements: 
#' - `total`: Number of transitions in the dataset
#' - `increasing`: Of which the number of species increases
#' - `decreasing`: Of which the number of species decreases
#' 
#' @examples
#' dataset <- data.frame(ID = c("(00)", "(01)", "(11)"), fitness = c(0.1, 0.2, 0.3)) %>% 
#'   fitnesses_to_adjacency()
#' interaction_transitions(dataset)
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

#' Generate ecological transition circos plot
#'
#' Visualizes fitness-increasing transitions between ecological states using a circos chord diagram.
#'
#' @details
#' Combines adjacency data, sector sorting, and color schemes to produce a
#' directed network visualization of evolutionary and ecological transitions.
#'
#' @param dataset Adjacency matrix of transitions.
#' @param sorting Logical; whether to automatically sort sectors.
#' @param highlighting Optional vector of nodes to highlight.
#' @param sector_order Optional custom order of sectors.
#' @param plot_labels Optional custom label vector.
#' @param highlight_mode Character; one of `"all"`, `"basin"`, or `"path"`.
#' @param ... Additional arguments passed to `circlize::chordDiagram()`.
#'
#' @return A circos plot rendered to the active graphical device.
#'
#' @examples
#' \dontrun{
#' eco_transition_plot(dataset, highlight_mode = "all")
#' }
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
  total_connections <- 2 * sum(dataset)
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
