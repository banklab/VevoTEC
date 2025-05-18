library(tidyverse)
library(shiny)
#library(plotly) #we actually don't need this since chord diagrams aren't compatible with ggplot
library(circlize)

# data may take the form of an adjacency matrix or dataframe
# i am starting with a simple csv of first-order transitions
x <- read.table("example_data/transitions_data.d", header = T, sep = ",")

# we have two arrow types to indicate directionality
chordDiagram(x, directional = 1, 
             direction.type = "arrows", 
             link.arr.type = "big.arrow", 
             diffHeight = -mm_h(2), 
             # color ramp can be set to, for instance, the indegree
             grid.col = colorRampPalette(c("lightblue", "red"))(8),
             annotationTrack = c("name", "grid"))
circos.clear()
      
#colors of the arrows can be changed
chordDiagram(x, directional = 1, 
             direction.type = "arrows", 
             grid.col = colorRampPalette(c("lightblue", "red"))(8),
             annotationTrack = c("name", "grid"))
circos.clear()



# ok now trying a dataset with higher order interactions (two resources)
y <- read.table("example_data/transitions_data_full.d", header = T, sep = ";") %>% 
  drop_na()
sector_order <- c("(0)", "(1)", "(2)", "(3)",
           "(4)", "(5)", "(6)", "(7)",
           "(0,3)","(0,7)","(1,2)","(1,5)",
           "(2,3)","(2,7)","(5,7)")
ratio <- c()
for (state in unique(c(y$From, y$To))){
  ratio <- c(ratio, sum(str_count(state, y$From)) / sum(c(str_count(state, y$From), str_count(state,y$To))))
}

# there's some strange problem with (6);(7), it detects the destination 7 as different than
# the source 7? I cut it out of the data for the time being

circos.par$gap.degree <- c(rep(1, times = 7), 15, rep(1, times = 6), 15)
circos.par$start.degree <- -7.5 
values = runif(length(sector_order))
chordDiagram(y, order = sector_order,
             annotationTrack = c("name", "grid"),
             directional = 1, 
             direction.type = "arrows",
             preAllocateTracks = list(
               track.height = 0.05  # reserves one custom track (you can adjust height)
             ))
xlim <- get.cell.meta.data("xlim")
ylim <- get.cell.meta.data("ylim")
circos.track(
  track.index = 2,  # This is the index of the preallocated track (1-based)
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xcenter <- get.cell.meta.data("xcenter")
    circos.rect(
      xleft = xlim[1], ybottom = ylim[1],
      xright = xlim[2], ytop = ylim[2],
      col = rand_color(1),
      border = NA)
    }, 
  bg.border = NA
)
circos.clear()

