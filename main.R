library(tidyverse)
library(shiny)
# library(plotly) we actually don't need this since chord diagrams aren't compatible with ggplot
library(chorddiag)
g <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") + 
  xlim(1, 6) + ylim(40, 100)
ggplotly(g)
