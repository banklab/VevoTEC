# VevoTEC: a tool for <ins>V</ins>isualizing <ins>evo</ins>lutionary <ins>T</ins>ransformations of <ins>E</ins>cological <ins>C</ins>ommunities

# Introduction

Interactions between species have traditionally been studied by <i>Ecology</i>, whereas genetic changes occurring within species have been the focus of <i>Evolutionary Biology</i>. However, it is increasingly appreciated that evolution is mediated by complex interactions between species - for example, genomes of parasites and their hosts evolve in response to each other. Their study thus requires a single formalism combining ecological and evolutionary principles. The genetic bases of such co-evolving systems are increasingly better understood, and there is a growing need for better data visualization for such systems that aid understanding from academics and the public. To address this need, we present VevoTEC, an Rshiny tool to transform eco-evolutionary data into clear visualizations of how interacting organisms coevolve at the genetic level.

# Installation

VevoTEC can be installed from the command line with `git clone git@github.com:banklab/VevoTEC.git`.

# Usage

## Prerequisites

`R` version 4.5.1 or later and libraries `circlize`, `igraph`, `shiny`, `bslib`, & `tidyverse`.

## Run

VevoTEC is an Rshiny dashboard that allows for the visualization of multi-species eco-evo dynamics using circos plots. It can be run through an R IDE of choice, or on the command line via `R -e "shiny::runApp('VevoTEC/app')"`, and visiting the corresponding localhost port.   

Contained are three widgets that allow for the upload of datasets, and different visualizations of data:
Widget | Functionality
---|---
Upload Data | 
Highlight Adaptive Basin(s) | 
Display Shortest Path(s) | 
