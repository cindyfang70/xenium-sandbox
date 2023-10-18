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
})
###########################################################
# Perform Delaunay triangulation on each of the clusters  #
# in the tissue, prune the global long edges, and plot.   #
###########################################################
source(here("code", "cindy", "04_delaunay", "delaunay.R"))
args <- commandArgs(trailingOnly = TRUE)
print(args[[1]])
sfe <- readRDS(args[[1]])
clusterName <- args[[2]] # get the name of the clustering output to do the triangulation on

tris <- list()
for(i in 1:length(unique(colData(sfe)[[clusterName]]))){
    clust <- unique(colData(sfe)[[clusterName]])[[i]]
    print(clust)
    colData(sfe)[clust] <- colData(sfe)[[clusterName]]==clust
    
    # do the triangulation 
    tri <- delaunay(sfe, clust, seed=i)
    summary(tri)
    tris <- rlist::list.append(tris, tri)
}

fname <- paste(sfe$region_id[[1]], "delaunay", clusterName, sep="-")
pdfname <- paste0(fname, ".pdf")

pdf(here("plots", "cindy", "04_delaunay", pdfname))
plotSpatialFeature(sfe, clusterName, colGeometryName="cellSeg")
#print(length(tris))
for (j in 1:length(tris)){
    tri <- tris[[j]]
    print(summary(tri))
    # compute the length of the edges
    tri <- getEdgesSpatialDistance(tri)
    
    longInds <- c()
    for (i in 1:tri$n){
        # compute the global constraint and find the long edges
        gc <- computeGlobalEdgeConstraint(tri, i)
        long <- identifyGlobalLongEdges(tri, i, gc)
        longInds <- c(longInds, long)
    }
    # plot the unpruned triangulation
    p <- plotDelaunay(tri, sfe)
    hist <- plotEdgeLengthHistogram(tri)
    
    # plot the pruned tri
    globalPrunedtri <- tri
    globalPrunedtri$arcs <- globalPrunedtri$arcs[-longInds,]
    pruned.p <- plotDelaunay(globalPrunedtri, sfe)
    pruned.hist <- plotEdgeLengthHistogram(globalPrunedtri, title = "Edge lengths (pruned)")
    
    do.call(grid.arrange, c(list(p, hist, pruned.p, pruned.hist), ncol=2))
    saveRDS(globalPrunedtri, here("processed-data", "04_delaunay", paste0(fname, ".RDS")))
}
    
dev.off()

