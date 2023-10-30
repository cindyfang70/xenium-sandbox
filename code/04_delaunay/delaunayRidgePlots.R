suppressPackageStartupMessages({
    library(Voyager)
    library(patchwork)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(ggplot2)
    library(bluster)
    library(stringr)
    library(BiocSingular)
    library(scater)
    library(rjson)
    library(Matrix)
    library(vroom)
    library(sf)
    library(BiocParallel)
    library(dplyr)
    library(here)
    library(ggforce)
    library(interp)
    library(igraph)
    library(rlist)
    library(gridExtra)
})
###########################################################
# Plot the ridge plots for the whole slide clustering.    #
# Delaunay triangulations.   #
###########################################################
#source(here("code", "cindy", "04_delaunay", "delaunay.R"))
#args <- commandArgs(trailingOnly = TRUE)
source(here("code","04_delaunay","delaunay.R"))
args <- system("ls data/04_delaunay", intern=TRUE)

tri.names <- paste("data", "04_delaunay", args, sep="/")
all.tris <- lapply(tri.names, readRDS)

region.names <- unlist(lapply(strsplit(args, split="-"), "[", 1))

# store the region name for each triangulation in the arcs for easy plotting
for(i in 1:6){
    tris <- all.tris[[i]]
    for (j in 1:length(tris)){
        tri <- tris[[j]]
        tri$arcs$region <- region.names[[i]]
    }
}
