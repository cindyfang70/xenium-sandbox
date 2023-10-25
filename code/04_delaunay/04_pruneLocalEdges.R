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
    library(alphahull)
    library(sp)
})

sfe <- readRDS("Br2743_Mid_SFE_filt.RDS")
tris <- readRDS("Br2743_Mid_5548-delaunay-lou25.RDS")
source(here("code","04_delaunay","delaunay.R"))

g <- igraph::graph.edgelist(arcs(tri))

paths <- all_simple_paths(g, 1, V(g))

paths.lengths <- lapply(paths, length)
length_two_paths <- paths[which(paths.lengths <=2)]
