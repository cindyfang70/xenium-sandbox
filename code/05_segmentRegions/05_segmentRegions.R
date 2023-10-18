suppressPackageStartupMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(stringr)
    library(BiocSingular)
    library(Matrix)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
})
###########################################################
# Find the connected components within the pruned         #
# triangulation graph, and compute the alpha shape        #
###########################################################

source(here("code","cindy","04_delaunay","delaunay.R"))
sfe <- readRDS(args[[1]])
tri <- readRDS(args[[2]])

# build the adjacency matrix for the triangulation
adj <- delaunayAdjacencyMatrix(sfe=sfe, tri=tri)

# build an igraph graph from the adjacency matrix
g <- graphFromAdjacencyMatrix(adj)

# find the connected components within the graph
comp <- components(adj)
