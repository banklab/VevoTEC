# VevoTEC: a tool for <ins>V</ins>isualizing <ins>evo</ins>lutionary <ins>T</ins>ransformations of <ins>E</ins>cological <ins>C</ins>ommunities

## Introduction

Interactions between species have traditionally been studied by <i>Ecology</i>, whereas genetic changes occurring within species have been the focus of <i>Evolutionary Biology</i>. However, it is increasingly appreciated that evolution is mediated by complex interactions between species - for example, genomes of parasites and their hosts evolve in response to each other. Their study thus requires a single formalism combining ecological and evolutionary principles. The genetic bases of such co-evolving systems are increasingly better understood, and there is a growing need for better data visualization for such systems that aid understanding from academics and the public. To address this need, we present VevoTEC, an Rshiny tool to transform eco-evolutionary data into clear visualizations of how interacting organisms coevolve at the genetic level.

## Installation

VevoTEC can be installed from the command line with `git clone git@github.com:banklab/VevoTEC.git`.

## Usage

### Prerequisites

`R` version 4.5.1 or later and libraries `circlize`, `igraph`, `shiny`, `bslib`, & `tidyverse`.

### Run

VevoTEC is an Rshiny dashboard that allows for the visualization of multi-species eco-evo dynamics using circos plots. It can be run through an R IDE of choice, or on the command line via `R -e "shiny::runApp('VevoTEC/app')"`, and visiting the corresponding localhost port.   

Contained are three widgets that allow for the upload of datasets, and different visualizations of data:
Widget | Functionality
---|---
Upload Data | Allows for the upload of a data file and optional sorting. Data formats described below. 
Highlight Adaptive Basin(s) | Allows the visualization of basins of accessibility for each peak in the dataset - the peak and all states that can access it through fitness-increasing paths are highlighted
Display Shortest Path(s) | Allows the visualization of the shortest path(s) from some arbitrary state to a selected peak - if multiple paths with equal length are possible, a choice between them is given.

### Data format

Input data may be uploaded in two formats: either as a fitness table, or an adjacency matrix.
1. Fitness table: CSV with columns `species1, species2, Fitness`, where fitness is a numeric value.
```
"Focal","Partner","fitness"
11,11,-0.605299473925827
11,12,-0.574428759712703
11,21,-0.0847747458665471
11,22,-0.0539040316534228
12,11,-0.0710995530232892
12,12,0.435133122477581
12,21,-0.263352403310527
12,22,0.242880272190343
21,11,-0.271682489810885
21,12,-0.338070988431988
21,21,0.0422550966343853
21,22,-0.024133401986718
22,11,-0.307102860757212
22,12,0.0289166217001666
22,21,-0.56703551498691
22,22,-0.231016032529531
```
2. Adjacency matrix: tab-separated matrix with a single header row of format `(ID_1)  (ID_2)  (ID_3)... (ID_N)`, where IDs are identifiers for different species communities. 
```
(000) (001) (010) (011) (100) (101) (110) (111) (000,011) (000,111) (001,010) (001,101) (010,011) (010,111) (101,111)
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
1 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
0 0 1 0 1 0 0 1 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 1 0 0 0 0 0 0 0 1 0 0
0 0 0 0 1 0 0 0 1 0 0 0 0 0 0
1 0 0 0 0 0 0 0 0 0 0 0 1 0 0
1 0 0 1 1 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 1 0 0 1 0 0
0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
```

