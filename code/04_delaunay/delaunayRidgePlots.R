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
    library(ggridges)
})
###########################################################
# Plot the ridge plots for the whole slide clustering.    #
# Delaunay triangulations.   #
###########################################################
#source(here("code", "cindy", "04_delaunay", "delaunay.R"))
#args <- commandArgs(trailingOnly = TRUE)
source(here("code","04_delaunay","delaunay.R"))
args <- system("ls data/04_delaunay", intern=TRUE)
args <- args[grepl("WholeSlide", args)]

tri.names <- paste("data", "04_delaunay", args, sep="/")
all.tris <- lapply(tri.names, readRDS)

region.names <- unlist(lapply(strsplit(args, split="-"), "[", 1))

# store the region name for each triangulation in the arcs for easy plotting
all.tris.arcs <- list()
for(i in 1:6){
    tris <- all.tris[[i]]
    arcs <- list()
    for (j in 1:length(tris)){
        tri <- tris[[j]]
        tri$arcs$region <- region.names[[i]]
        arcs <- list.append(arcs, tri$arcs)
        
    }
    tris.arcs <- as.data.frame(do.call(rbind, arcs))
    all.tris.arcs <- list.append(all.tris.arcs, tris.arcs)
}
 # need to do it one slide at a time, not both togheter
all.tris.arcs <- do.call(rbind, all.tris.arcs)

tris.5434 <- all.tris.arcs %>%
    filter(grepl("5434", region))%>%
    filter(edge.lengths <=1000)

tris.5548 <- all.tris.arcs %>%
    filter(grepl("5548", region))%>%
    filter(edge.lengths <=1000)

pdf(here("plots", "delaunay", "wholeSlideClusteringRidgePlots.pdf"))
ggplot(tris.5434, aes(x = edge.lengths, y = region, fill =region)) +
    geom_density_ridges() +
    theme_ridges() + 
    facet_wrap(~clust)+
    theme_minimal()+
    theme(legend.position = "none")

ggplot(tris.5548, aes(x = edge.lengths, y = region, fill =region)) +
    geom_density_ridges() +
    theme_ridges() + 
    facet_wrap(~clust)+
    theme_minimal()+
    theme(legend.position = "none")
dev.off()
